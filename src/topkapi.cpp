#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <math.h>
#include <map>
#include <mpi.h>
#include <omp.h>

#include "LossyCountMinSketch.h"
#include "HashFunction.h"

/* Reads data in 64MB chunks */
#define CHUNK_SIZE (64*1024*1024)
#define TRUE 1
#define FALSE 0

/* Allocates a row of topkapi sketch data structure */
void allocate_sketch( LossySketch* sketch,
                      const unsigned range )
{
  int i;
  (*sketch)._b           = range;
  (*sketch).identity     = (char*) malloc(range*MAX_WORD_LENGTH*sizeof(char));
  (*sketch).lossyCount   = (int* ) malloc(range*sizeof(int));
	if ( (*sketch).identity == NULL ||
			 (*sketch).lossyCount == NULL )
	{
		fprintf(stderr, "LossySketch allocation error!\n");
		exit(EXIT_FAILURE);
	}
  /* set counts to -1 to indicate empty counter */
  for (i = 0; i < range; ++i)
    (*sketch).lossyCount[i] = -1;
}

/* Frees a row of topkapi sketch data structure */
void deallocate_sketch( LossySketch* sketch )
{
  free((*sketch).identity);
  free((*sketch).lossyCount);
}

/* Function for merging topkapi sketches between two nodes */
void merge_2( LossySketch* local_sketch,
              LossySketch* sketch_to_merge,
              const int    num_hash_func,
              const int    target_node,
              int          if_sender )
{
  int i;
  /* MPI variables */
  MPI_Status status;
  MPI_Request request_count, request_word;
  
  const int num_elem = local_sketch[0]._b;
  
  /* data transfer of first row */
  if (if_sender)
  {
    MPI_Send(local_sketch[0].lossyCount, num_elem, MPI_INT,
             target_node, 0, MPI_COMM_WORLD);
    MPI_Send(local_sketch[0].identity, (num_elem*MAX_WORD_LENGTH), MPI_CHAR,
             target_node, 0, MPI_COMM_WORLD);

  } else {
    MPI_Recv(sketch_to_merge[0].lossyCount, num_elem, MPI_INT,
             target_node, 0, MPI_COMM_WORLD, NULL);
    MPI_Recv(sketch_to_merge[0].identity, (num_elem*MAX_WORD_LENGTH), MPI_CHAR,
             target_node, 0, MPI_COMM_WORLD, NULL);
  }

  for (i = 1; i < num_hash_func; ++i)
  {
    if (if_sender)
    {
      MPI_Isend(local_sketch[i].lossyCount, num_elem, MPI_INT,
                target_node, 0, MPI_COMM_WORLD, &request_count);
      MPI_Isend(local_sketch[i].identity, (num_elem*MAX_WORD_LENGTH), MPI_CHAR,
                target_node, 0, MPI_COMM_WORLD, &request_word);
      MPI_Wait( &request_count, &status );
      MPI_Wait( &request_word, &status );
    } else {
      MPI_Irecv(sketch_to_merge[i].lossyCount, num_elem, MPI_INT,
                target_node, 0, MPI_COMM_WORLD, &request_count);
      MPI_Irecv(sketch_to_merge[i].identity, (num_elem*MAX_WORD_LENGTH), MPI_CHAR,
                target_node, 0, MPI_COMM_WORLD, &request_word);
      dist_merge_sketch( &local_sketch[i-1], &sketch_to_merge[i-1]);
      MPI_Wait( &request_count, &status );
      MPI_Wait( &request_word, &status );
    }
  }

  /* merging for last row */
  if (!if_sender)
  {
    dist_merge_sketch( &local_sketch[num_hash_func-1], 
                       &sketch_to_merge[num_hash_func-1]);
  }
}

/* Function for performing merging of topkapi sketches over all nodes */
void merge( LossySketch* local_sketch,
            LossySketch* sketch_to_merge,
            const int    num_hash_func,
            const int    num_nodes,
            const int    my_rank )
{
  unsigned i, j;
  int target_node;
  int if_sender;
  unsigned num_reductions = num_nodes;
  int cond_to_participate;
  /* parallel reduction */
  for ( i = 1, j= 0; i < num_nodes; i *= 2, ++j)
  {
    cond_to_participate = (my_rank % i == 0);
    if (num_reductions % 2 != 0)
    {
      cond_to_participate = cond_to_participate && 
                            (my_rank != ((num_reductions - 1) * i));
      num_reductions /= 2;
      ++num_reductions;
    } else {
      num_reductions /= 2;
    }
    if (cond_to_participate)
    {
      if_sender = (my_rank/i) % 2;
      if (if_sender)
        target_node = my_rank - i;
      else
        target_node = my_rank + i;
      merge_2(local_sketch, sketch_to_merge, num_hash_func,
              target_node, if_sender);
    }
  }
}

int main(int argc, char* argv[])
{
  int i, j;
  /* Default Parameter values for Sketch */
  unsigned range = 1024; /* range of buckets for the sketch */
  unsigned log2range = 10; /* log2 value of range */
  unsigned num_hash_func = 4; /* number of hash functions */
  unsigned K = 100; /* Top K */
  unsigned frac_epsilon = 10*K; /* error margin */
  
  
  /* Data for Sketch */
  LossySketch* th_local_sketch; /* array of thread local sketches */
  LossySketch* node_final_sketch; /* final sketch in a node after merging thread local sketches */
  LossySketch* sketch_to_merge; /* required for receiving sketch from other nodes */
  
  /* Random numbers for hash functions */
  char* hash_func_rands;
  int rand_differ;
  srand(time(NULL));

  /* Data structures for sorting the heavy hitters to find topK */
  std::map<int,int> topk_words;
  std::map<int,int>::reverse_iterator rit;

  /* variables for file operation */
  char*  my_file;
  char   default_out_file[] = "topK.out"; /* default output file name */
  char*  output_file = &default_out_file[0]; /* used for storing TopK words */
  char*  str_buff[2]; /* use two buffers alternatively to do parallel read and computation */
  int    str_buff_write_id = 0;
  int    str_buff_read_id = 0;
  size_t read_size;
  size_t compute_size;
  FILE*  fp;
  int    file_end_flag = FALSE;
  int    do_compute_flag = TRUE;

  /* variables for OpenMP */
  int num_threads = omp_get_num_procs(); /* default to number of procs available */

  /* variables for MPI */
  int num_nodes;
  int my_rank;
	int provided_threading;
  int error = 0;
  int reduced_error = 0;
  MPI_Status status;
  double start_time, end_time;
  double total_time = 0.0;

  /* Initialize MPI */
  MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, 
												&provided_threading);
	assert(provided_threading == MPI_THREAD_FUNNELED);
  MPI_Comm_size(MPI_COMM_WORLD, &num_nodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  /* Argument handling */
  if (argc <= 2 && (!strcmp(argv[1], "--help") ||
        !strcmp(argv[1], "-h")))
  {
    if (my_rank == 0)
    {
      fprintf(stdout, "Example usage:\n"
          "\t./topkapi <args> <inputs> <output>\n"
          "<args>:\n\t-b <range value> [--optional] {default 1024}"
          "\n\t-f <number of hash functions> [--optional] {default 4}"
          "\n\t-k <value of K in TopK> [--optional] {default 100}"
          "\n\t-t <number of threads in a node> [--optional] {default 8}"
          "<inputs>:\n\t <input file names> [--required]"
          "<output>:\n\t-o <output file name> [--optional] {default \"topK.out\"}");
    }
    MPI_Finalize();
    exit(EXIT_FAILURE);
  } else {
    for (i = 1; i < argc; ++i)
    {
      if (!strcmp(argv[i], "-b"))
      {
        range = atoi(argv[i+1]);
        log2range = (unsigned) floor(log2(range));
        ++i;
      } else if (!strcmp(argv[i], "-f"))
      {
        num_hash_func = atoi(argv[i+1]);
        ++i;
      } else if (!strcmp(argv[i], "-k"))
      {
        K = atoi(argv[i+1]);
        frac_epsilon = K*10;
        ++i;
      } else if (!strcmp(argv[i], "-t"))
      {
        num_threads = atoi(argv[i+1]);
        ++i;
      } else
        break;
     }
  }
   
  if (argc < (i + num_nodes))
  {
    if (my_rank == 0)
      fprintf(stderr, "Not enough input files!\n"
           "Expecting %d input files, but got "
           " %d input files in arguments.\n", num_nodes,
           (argc - i));
    MPI_Finalize();
    exit(EXIT_FAILURE);
  } else if (argc > (i + num_nodes + 2))
  {
    if (my_rank == 0)
      fprintf(stderr, "More arguments than expected!\n"
           "Please check usage by ./topkapi <-h/--help>\n");
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
   
  my_file = argv[i+my_rank];

  if (argc > (i + num_nodes))
  {
    if (!strcmp(argv[i+num_nodes], "-o"))
    {
      if (argc == (i + num_nodes + 1))
      {
        if (my_rank == 0)
          fprintf(stderr, "Expecting an output file name after argument \"-o\", "
              "none provided!\n");
        MPI_Finalize();
        exit(EXIT_FAILURE);
      } else
        output_file = argv[i+num_nodes+1];
    } else {
      if (my_rank == 0)
        fprintf(stderr, "Unrecognized argument(s) after input file names!\n"
            "Expecting \"-o <output file name>\"\n");
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
  }
  /* end of argument handling */

  /* printing the execution configuration */
  if (my_rank == 0) 
  {
    fprintf (stdout, "*************Parameter Values*************\n"
                     "Range: %d, Number of hash functions: %d\n"
                     "K: %d, Num nodes: %d, Num threads/Node: %d\n"
                     "******************************************\n",
                     range, num_hash_func, K, num_nodes, num_threads);
    fprintf (stdout, "************Output file name************\n"
                     "%s\n"
                     "****************************************\n", 
                     output_file);
  }
 
  /* hash function specific different random number generartion */ 
  hash_func_rands = (char* )malloc(num_hash_func*sizeof(char));
  if (my_rank == 0)
  {
    for (i = 0; i < num_hash_func; ++i)
    {
      rand_differ = FALSE;
      hash_func_rands[i] = (char) (rand() % 47);
      while (!rand_differ)
      {
        rand_differ = TRUE;
        for (j = 0; j < i; ++j)
        {
          if (hash_func_rands[i] == hash_func_rands[j])
          {  
            rand_differ = FALSE;
            break;
          }
        }
        if (!rand_differ)
        {
          hash_func_rands[i] = (char) (rand() % 47);
        }
      }
    }
  }

  MPI_Bcast(hash_func_rands, num_hash_func, MPI_CHAR, 0, MPI_COMM_WORLD);

  /* File operations */
  fp = fopen64(my_file,"r");
  if (fp == NULL)
  {
    fprintf(stderr, "Error opening file!\nNode:%d, Filename:%s\n",
        my_rank, my_file);
    error = 1;
  } 
  MPI_Allreduce(&error, &reduced_error, 1, MPI_INT, MPI_BOR, MPI_COMM_WORLD);
  if (reduced_error)
  {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  str_buff[0] = (char* )malloc(CHUNK_SIZE*sizeof(char));
  str_buff[1] = (char* )malloc(CHUNK_SIZE*sizeof(char));
  if (str_buff[0] == NULL || str_buff[1] == NULL)
  {
    fprintf(stderr, "Error allocating memory for reading files!\nNode:%d, " 
        "Filename:%s\n", my_rank, my_file);
    error = 1;
  }
  MPI_Allreduce(&error, &reduced_error, 1, MPI_INT, MPI_BOR, MPI_COMM_WORLD);
  if (reduced_error)
  {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  read_size = fread(str_buff[str_buff_write_id], 1, CHUNK_SIZE, fp);
  if (read_size != CHUNK_SIZE)
  {
    if (feof(fp))
      file_end_flag = TRUE;
    else
    {
      fprintf(stderr, "Error reading file!\nNode:%d, "
          "Filename:%s\n", my_rank, my_file);
      error = 1;
    }
  }
  MPI_Allreduce(&error, &reduced_error, 1, MPI_INT, MPI_BOR, MPI_COMM_WORLD);
  if (reduced_error)
  {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  compute_size = read_size / num_threads;
  for (i = 0; i < num_threads; ++i)
    str_buff[str_buff_write_id][(i+1)*compute_size - 1] = '\0';
  str_buff_write_id ^= 0x1;


  th_local_sketch = (LossySketch* ) malloc(num_hash_func*num_threads*
                                             sizeof(LossySketch));
  node_final_sketch = th_local_sketch; /* represents the first 
       num_hash_func number of sketches from th_local_sketch array */
  sketch_to_merge = (LossySketch* ) malloc(num_hash_func*
                                             sizeof(LossySketch));
  for (i = 0; i < num_hash_func; ++i)
  {
    allocate_sketch( &sketch_to_merge[i], range);
  }
  
  omp_set_num_threads( num_threads );
  omp_set_dynamic( 0 );
#pragma omp parallel firstprivate(th_local_sketch, range,\
                                  num_hash_func)
  {
    int tid = omp_get_thread_num();
    int th_i;
    /* Allocate and initialize sketch variables */
    for (th_i = 0; th_i < num_hash_func; ++th_i)
    {
      allocate_sketch( &th_local_sketch[tid*num_hash_func+th_i], range);
    }
  }

  while (do_compute_flag)
  {
    if (!file_end_flag)
    {
      read_size = fread(str_buff[str_buff_write_id], 1, CHUNK_SIZE, fp);
      if (read_size != CHUNK_SIZE)
      {
        if (feof(fp))
          file_end_flag = TRUE;
        else
        {
          fprintf(stderr, "Error reading file!\nNode:%d, "
                    "Filename:%s\n", my_rank, my_file);
          exit(EXIT_FAILURE);
        }
      }
      for (i = 0; i < num_threads; ++i)
        str_buff[str_buff_write_id][(i+1)*compute_size - 1] = '\0';
      str_buff_write_id ^= 0x1;
    } else {
      do_compute_flag = FALSE;
    }

    start_time = MPI_Wtime();
    omp_set_num_threads( num_threads );
    omp_set_dynamic( 0 );
#pragma omp parallel firstprivate(range, log2range,\
    num_hash_func, th_local_sketch, str_buff, hash_func_rands) private (i)
    { 
      int tid = omp_get_thread_num();

      char*  word;
      char*  prev_word_ptr;
    
      /* read words from the buffer one by one */
      word = strtok_r (&str_buff[str_buff_read_id][tid*compute_size], " \n", &prev_word_ptr);
      while (word != NULL)
      {
        /* we have a word */
        for (i = 0; i < num_hash_func; ++i)
          update_sketch( &th_local_sketch[tid * num_hash_func + i], hash_func_rands[i], word, log2range );
        word = strtok_r (NULL, " \n", &prev_word_ptr);
      }
    }    
    end_time = MPI_Wtime();
    total_time += (end_time - start_time);
    
    compute_size = read_size / num_threads;
    str_buff_read_id ^= 0x1;
  }
  fclose(fp);

  /* Now merge thread local sketches to produce final sketch for a node */
  start_time = MPI_Wtime();
  for (i = 0; i < num_hash_func; ++i)
    local_merge_sketch(th_local_sketch, num_threads, num_hash_func, i);

  /* Now merge sketches across nodes */
  merge( node_final_sketch, sketch_to_merge, num_hash_func, 
         num_nodes, my_rank );
  end_time = MPI_Wtime();
  total_time += (end_time - start_time);

  /* Print TopK words */  
  if (my_rank == 0)
  {
    int num_heavy_hitter = 0;
    int count;
    char* str;
    int id;
    int* is_heavy_hitter = (int* )malloc(range*sizeof(int));
    int threshold = (int) ((range/K) - 
                           (range/frac_epsilon));
    fp = fopen(output_file,"w");
    if (fp == NULL)
    {
      fprintf(stderr, "Error opening output file!\nFilename:%s\n",
          my_rank, output_file);
      fprintf(stderr, "Printing the TopK words to stdout!\n");
      fp = stdout;
    } 
#pragma omp parallel for schedule(static) firstprivate(threshold,\
    num_hash_func, log2range, is_heavy_hitter, hash_func_rands) private(j, count, str, id)
    for (i = 0; i < range; ++i)
    {
      is_heavy_hitter[i] = FALSE;
      for (j = 0; j < num_hash_func; ++j)
      {
        if ( j == 0)
        {
          str = &node_final_sketch[0].identity[i*MAX_WORD_LENGTH];
          count = node_final_sketch[0].lossyCount[i];
          if (count >= threshold)
          {
            is_heavy_hitter[i] = TRUE;
          }
        } else {
          id = gethash( str, hash_func_rands[j], log2range );
          if (strcmp(&node_final_sketch[j].identity[id*MAX_WORD_LENGTH], str) != 0)
          {
            continue;
          } else if (node_final_sketch[j].lossyCount[id] >= threshold)
          {
            is_heavy_hitter[i] = TRUE;
          }
          if (node_final_sketch[j].lossyCount[id] > count)
            count = node_final_sketch[j].lossyCount[id];
        }
      }
      node_final_sketch[0].lossyCount[i] = count;
    }

    fprintf(fp, "Elapsed time: %fseconds\n", total_time);

    for (i = 0; i < range; ++i)
    {
      if (is_heavy_hitter[i])
      {
        num_heavy_hitter ++;
        topk_words.insert( std::pair<int,int>(node_final_sketch[0].lossyCount[i], i) );
      }
    }

    for (i = 0, rit = topk_words.rbegin(); 
          (i < K) && (rit != topk_words.rend()); 
            ++i, ++rit)
    {
      j = rit->second;
      fprintf(fp, "%s %d\n", &node_final_sketch[0].identity[j*MAX_WORD_LENGTH], 
                rit->first);
    }

    if (fp != stdout) fclose(fp);
    /* free memories */
    free(is_heavy_hitter);
  }
  for (i = 0; i  < num_threads; ++i)
  {
    for (j = 0; j < num_hash_func; ++j)
    {
      if (i == 0)
        deallocate_sketch( &sketch_to_merge[j] );
      deallocate_sketch( &th_local_sketch[i*num_hash_func + j] );
    }
  }
  free(th_local_sketch);
  free(sketch_to_merge);
  free(str_buff[0]);
  free(str_buff[1]);
  free(hash_func_rands);
  MPI_Finalize();
  return 0;
}
