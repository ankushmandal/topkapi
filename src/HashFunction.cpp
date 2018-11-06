#include <string.h>
#include "HashFunction.h"
#include "MurmurHash3.h"

/* Details about the API are given in HashFunction.h */
/* Note: this function uses VLA (Variable Length Array)*/
unsigned gethash( const char* word, char hash_rand, unsigned log2range)
{
  unsigned hash;
  uint32_t seed = 0x9747b28c;
  size_t len = strlen( word ) + 1;
  char str_to_hash[len];
  strcpy( str_to_hash, word ) ;
  /* hash_rand is hash function specific random value, 
   * we put it at the end of the string */
  str_to_hash[len-1] = hash_rand;
  MurmurHash3_x86_32( str_to_hash, len*sizeof(char), seed, &hash );
  /* Now, we shift bits to limit the 
   * hash value  within log2range */
  hash = hash << (32 - log2range);
  hash = hash >> (32 - log2range);
  return hash;
}
