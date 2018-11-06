#ifndef _HASH_FUNCTION_H_
#define _HASH_FUNCTION_H_

/* API to get hash value of a string
 * It's a wrapper function for MurmurHash3
 * Arguments:
 * const char* word - the string to be hashed
 * char hash_rand - a variable that is appended at the
 *                  end of 'word' before computing the 
 *                  hashvalue
 * unsigned log2range - log base 2 of range for hashing or number 
 *                      of buckets in hashtable, used for number
 *                      of bits used in final hash value
 * Returns:
 * A hash value between 0 and 2^('log2range')
 */
unsigned gethash( const char* word, char hash_rand, unsigned log2range );

#endif /* _HASH_FUNCTION_H_ */
