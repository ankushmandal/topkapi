#include <stdlib.h>
#include <string.h>
#include "HashFunction.h"
#include "LossyCountMinSketch.h"

/* This function updates the local sketch 
 * based on input word
 */ 
void update_sketch( LossySketch* _sketch,
                    char hash_func_rand,
                    char* word,
                    unsigned log2range )
{
  if (strlen(word) >= MAX_WORD_LENGTH)
    return;

  unsigned bucket = gethash( word, hash_func_rand, log2range );
  
  int* count_ptr = &((*_sketch).lossyCount[bucket]);
  char* str_ptr = &((*_sketch).identity[bucket*MAX_WORD_LENGTH]);
  if (*count_ptr == -1)
  { /* if the counter is empty */
    strcpy (str_ptr, word);
    *count_ptr = 1;
  }
  else
  { /* if counter is not empty */
    /* we update count based on comparison of strings */
    if (!strcmp(str_ptr, word)) /* strcmp returns 0 if equal */
    { /* if same word */
      (*count_ptr) ++;
    }
    else
    { /* if words are different */
      if (--(*count_ptr) < 0)
      { /* if decreasing the counter makes it's value negative */
        /* replace previous word with new word and set counter */
        strcpy (str_ptr, word);
        *count_ptr = 1;
      }
    }
  }
}

/* This function merges thread local sketches to
 * create the final sketch for a node
 * Note: it is used only when multi-threaded
 * execution happens
 */
void local_merge_sketch( LossySketch*   LCMS,
                         const unsigned num_local_copies,
                         const unsigned num_hash_func,
                         const unsigned hash_func_index )
{
  char* word[num_local_copies];
  int count[num_local_copies];
  unsigned i, j, k, diff_words;
  int max_selected;
  char* current_word;
  int max_count;
  unsigned range = LCMS[0]._b;

  for (i = 0; i < range; ++i)
  {
    word[0] = &(LCMS[hash_func_index].identity[i*MAX_WORD_LENGTH]);
    count[0] = LCMS[hash_func_index].lossyCount[i];
    diff_words = 1;
    for (j = 1; j < num_local_copies; ++j)
    {
      current_word = &(LCMS[j*num_hash_func+hash_func_index].identity[i*MAX_WORD_LENGTH]);
      for ( k = 0; k < diff_words; ++k)
      {
        if (!strcmp(current_word, word[k]) &&
            LCMS[j*num_hash_func+hash_func_index].lossyCount[i] != (-1))
        {
          /* if same word */
          count[k] += LCMS[j*num_hash_func+hash_func_index].lossyCount[i];
          break;
        }
      }
      if (k == diff_words)
      {
        word[diff_words] = current_word;
        count[diff_words] = LCMS[j*num_hash_func+hash_func_index].lossyCount[i];
        diff_words++;
      }
    }
    max_count = -1;
    max_selected = 0;
    k = 0;
    for (j = 0; j < diff_words; ++j)
    {
      if (count[j] != (-1))
      {
        if (max_selected)
        {
          if (count[j] > max_count)
          {
            max_count = (count[j] - max_count);
            k = j;
          } else {
            max_count -= count[j];
          }
        } else {
          max_count = count[j];
          k = j;
          max_selected = 1;
        }
      }
    }
    if (k != 0)
    {
      strcpy(word[0], word[k]);
    }
    LCMS[hash_func_index].lossyCount[i] = max_count;
  }
}

/* This function merges a row of a sketch from other node 
 * to a row of the sketch of current node 
 */
void dist_merge_sketch( LossySketch* dest_sketch,
                        LossySketch* sketch_to_merge )
{
  unsigned i;
  const unsigned range = (*dest_sketch)._b;
  char* dest_word; 
  char* merge_word;
  int* dest_count_ptr; 
  int* merge_count_ptr;
  for (i = 0; i < range; ++i)
  {
    dest_word = &((*dest_sketch).identity[i*MAX_WORD_LENGTH]);
    merge_word = &((*sketch_to_merge).identity[i*MAX_WORD_LENGTH]);
    dest_count_ptr = &((*dest_sketch).lossyCount[i]);
    merge_count_ptr = &((*sketch_to_merge).lossyCount[i]);
    if (*dest_count_ptr != (-1) || *merge_count_ptr != (-1))
    {
      if (!strcmp(dest_word, merge_word))
      {
        /* if words are same */
        *dest_count_ptr += *merge_count_ptr;
      }
      else {
        if (*dest_count_ptr == (-1))
        {
          *dest_count_ptr = *merge_count_ptr;
          strcpy(dest_word, merge_word);
        } else if (*merge_count_ptr != (-1))
        {
          /* if words are different */
          *dest_count_ptr -= *merge_count_ptr;
          if (*dest_count_ptr < 0)
          {
            strcpy(dest_word, merge_word);
            *dest_count_ptr = -(*dest_count_ptr);
          }
        }
      }
    }
  }
}
