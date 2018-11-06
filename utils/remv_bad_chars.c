#include <stdio.h>
#include <stdlib.h>

#define IS_PROPER(c) ((c >= 48 && c <= 57)\
    || (c >= 65 && c <= 90) || (c >= 97 && c <= 122)\
    || (c == 10))

int main(int argc, char* argv[])
{
  char curr_char;
  char* inputfile;
  char* outputfile;

  FILE* in_file;
  FILE* out_file;

  if (argc != 3)
  {
    fprintf(stderr, "Wrong number of arguments!\n");
    fprintf(stdout, "Example Usage:\nremv_bad_char [inputfile_name] [outputfile_name]\n");
    exit(EXIT_FAILURE);
  }
  
  inputfile = argv[1];
	fprintf(stdout, "file to be opened: %s\n", inputfile);
  outputfile = argv[2];  
	fprintf(stdout, "file to be written: %s\n", outputfile);

  in_file = fopen64(inputfile, "r");
  if (in_file == NULL)
  {
    fprintf(stderr, "Can not open input file!\n");
    exit(EXIT_FAILURE);
  }
  out_file = fopen64(outputfile, "w");

	curr_char = fgetc(in_file);
	while (curr_char != EOF)
	{
		if (IS_PROPER(curr_char))
		{
			if (curr_char >= 65 && curr_char <= 90)
				curr_char += 32;
			fputc(curr_char, out_file);
		} 
		else
			fputc(' ', out_file);
		curr_char = fgetc(in_file);
	}

	fclose(in_file);
	fclose(out_file);

  return EXIT_SUCCESS;
}
