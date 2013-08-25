#include <stdlib.h>
#include "packednumbers.h"

PackedNumberArray *NewPackedNumberArray(unsigned int numInts, unsigned int maxInt){
	PackedNumberArray *intArray;
	unsigned long long numBits;
	unsigned int n;
	intArray = (PackedNumberArray *)malloc(sizeof(PackedNumberArray));
	//(intArray->bitsPerWord) = (unsigned char)(sizeof(unsigned long long)*8); // use 64 bit words (8 bytes * 8 bits/byte)
	n = 1; // number of bits needed to store one number
	while( ( (1UL << n) - 1 ) < maxInt ) n++; // n bits per number (stores 2^n numbers, but the last one is (2^n-1))
	(intArray->bitsPerInt) = (unsigned char)n;
	numBits = ( ((unsigned long long)numInts) * ((unsigned long long)n) ); // total number of bits occupied by all the numbers
	if( numBits != 0) numBits--; // if it was a multiple of 64 , it would create an extra unused word
	n = (unsigned int)( (numBits/64ULL) + 1ULL); // number of 64 bit words required to store (numInts) numbers of (bitsPerInt) bits each
	(intArray->numWords) = n;
	(intArray->bitsArray) = (unsigned long long *)calloc(n,sizeof(unsigned long long)); // bit array that will store the numbers
	return intArray;
}

void FreePackedNumberArray(PackedNumberArray *intArray){
	free(intArray->bitsArray);
	free(intArray);
}

void ResetPackedNumberArray(PackedNumberArray *intArray){
	unsigned int n;
	n = (intArray->numWords);
	while( n != 0 ) (intArray->bitsArray)[(--n)] = 0ULL;
}

unsigned int GetPackedNumber(PackedNumberArray *intArray, unsigned int pos){
	unsigned char offset, bitsPerInt;
	unsigned long long temp;
	unsigned int number;
	bitsPerInt = (intArray->bitsPerInt);
	temp = ((unsigned long long)(pos)) * ((unsigned long long)(bitsPerInt)); // number of bits
	offset = (unsigned char)( temp & 63ULL ); // mask by (64-1) to get position of 1st bit inside the word
	pos = (unsigned int)( temp >> 6 ); // divide by 64 to get word position inside bit array
	temp = ( ( 1ULL << bitsPerInt ) - 1ULL ); // mask for the n bits of each number
	number = (unsigned int)( ( (intArray->bitsArray)[pos] >> offset ) & temp );
	offset = ( 64 - offset ); // check if the bits of the number extend to the next word (get number of bits left until the end of this word)
	if( offset < bitsPerInt ) number |= (unsigned int)( ( (intArray->bitsArray)[(++pos)] << offset ) & temp );
	return number;
}

void SetPackedNumber(PackedNumberArray *intArray, unsigned int pos, unsigned int num){
	unsigned char offset, bitsPerInt;
	unsigned long long temp;
	bitsPerInt = (intArray->bitsPerInt);
	temp = ((unsigned long long)(pos)) * ((unsigned long long)(bitsPerInt)); // number of bits
	offset = (unsigned char)( temp & 63ULL ); // mask by (64-1) to get position of 1st bit inside the word
	pos = (unsigned int)( temp >> 6 ); // divide by 64 to get word position inside bit array
	(intArray->bitsArray)[pos] |= ( ((unsigned long long)num) << offset );
	offset = ( 64 - offset ); // check if the bits of the number extend to the next word (get number of bits left until the end of this word)
	if( offset < bitsPerInt ) (intArray->bitsArray)[(++pos)] |= ( ((unsigned long long)num) >> offset );
}

void ReplacePackedNumber(PackedNumberArray *intArray, unsigned int pos, unsigned int num){
	unsigned char offset, bitsPerInt;
	unsigned long long temp;
	bitsPerInt = (intArray->bitsPerInt);
	temp = ((unsigned long long)(pos)) * ((unsigned long long)(bitsPerInt)); // number of bits
	offset = (unsigned char)( temp & 63ULL ); // mask by (64-1) to get position of 1st bit inside the word
	pos = (unsigned int)( temp >> 6 ); // divide by 64 to get word position inside bit array
	temp = ( ((unsigned long long)num) << offset ); // rightmost bits to set
	(intArray->bitsArray)[pos] &= temp; // only keep common 1 bits
	(intArray->bitsArray)[pos] |= temp; // set the missing 1 bits (same as reset and then set)
	offset = ( 64 - offset ); // check if the bits of the number extend to the next word (get number of bits left until the end of this word)
	if( offset < bitsPerInt ){
		pos++; // next word
		temp = ( ((unsigned long long)num) >> offset ); // leftmost bits to set
		(intArray->bitsArray)[pos] &= temp; // only keep common 1 bits
		(intArray->bitsArray)[pos] |= temp; // set the missing 1 bits (same as reset and then set)
	}
}

typedef struct _PackedIncreasingNumberArray {
	unsigned char numTotalBits;
	unsigned char numHighBits;
	unsigned char numLowBits;
	unsigned int highBitsMask;
	unsigned int lowBitsMask;
	unsigned int *highLimits;
	PackedNumberArray *lowBits;
} PackedIncreasingNumberArray;

// TODO: check and fix calculation of optimal number of high/low bits
PackedIncreasingNumberArray *NewPackedIncreasingNumberArray(unsigned int numInts, unsigned int maxInt){
	PackedIncreasingNumberArray *incIntArray;
	unsigned long long numBits;
	unsigned int k;
	unsigned char n;
	incIntArray = (PackedIncreasingNumberArray *)malloc(sizeof(PackedIncreasingNumberArray));
	//k = (unsigned int)(sizeof(unsigned int)*8); // use 32 bits words
	numBits = (unsigned long long)(numInts*32); // total number of bits to store all the words regularly
	n = 0; // get best number of high bits to compact
	while( (((1ULL << n)-1ULL)*32ULL) < numBits ) n++; // for compacting n high bits we need (2^n-1) extra numbers
	if( n != 0 ) n--; // the current value already exceeded the regular size
	(incIntArray->numHighBits) = n;
	n = 0; // number of bits per number
	k = 1; // power of two (2^n)
	while( k && ( (k-1) < maxInt ) ){ // (2^n-1) is the largest number that can be represented by n bits
		n++;
		k <<= 1;
	}
	if( k != 0 ) k--; // mask for all bits
	else k = (~k); // n = 32, whole word (~0UL)
	(incIntArray->numTotalBits) = n; // total bits of each number
	n -= (incIntArray->numHighBits); // number of lower bits
	(incIntArray->numLowBits) = n;
	(incIntArray->lowBitsMask) = ( (1UL << n) - 1UL ); // mask for lower bits
	(incIntArray->highBitsMask) = ( k ^ (incIntArray->lowBitsMask) ); // all bits of number except the low ones
	k = ( ( 1UL << (incIntArray->numHighBits) ) - 1UL ); // number of limits needed for compacting the high bits
	(incIntArray->highLimits) = (unsigned int *)calloc(k,sizeof(unsigned int));
	k = ( ( 1UL << (incIntArray->numLowBits) ) - 1UL ); // highest number representable by the low bits
	(incIntArray->lowBits) = NewPackedNumberArray(numInts,k);
	return incIntArray;
}

void FreePackedIncreasingNumberArray(PackedIncreasingNumberArray *incIntArray){
	FreePackedNumberArray(incIntArray->lowBits);
	free(incIntArray->highLimits);
	free(incIntArray);
}

void ResetPackedIncreasingNumberArray(PackedIncreasingNumberArray *incIntArray){
	unsigned int n;
	n = ( ( 1UL << (incIntArray->numHighBits) ) - 1UL ); // number of limits used
	while( n != 0 ) (incIntArray->highLimits)[(--n)] = 0UL;
	ResetPackedNumberArray(incIntArray->lowBits);
}

unsigned int GetPackedIncreasingNumber(PackedIncreasingNumberArray *incIntArray, unsigned int pos){
	unsigned char n;
	unsigned int number, limit, i, k;
	number = 0UL;
	i = pos; // position in the increasing numbers array
	k = 0; // current position in the limits array
	n = (incIntArray->numHighBits); // current high bit
	while( n && i ){
		limit = (incIntArray->highLimits)[k];
		if( i < limit ){ // bit 0, and advance to limit (2*k+1)
			k <<= 1;
			k++;
		} else { // bit 1, and advance to limit (2*(k+1))
			k++;
			k <<= 1;
			i -= limit; // fix array position relative to the current limit
			number++; // the number has a bit 1 at this position
		}
		n--;
		number <<= 1; // go to next bit of number
	}
	number <<= n; // if we didn't check all high bits, shift the number correctly
	return ( ( number << (incIntArray->numLowBits) ) | GetPackedNumber((incIntArray->lowBits),pos) );
}

void SetPackedIncreasingNumber(PackedIncreasingNumberArray *incIntArray, unsigned int pos, unsigned int num){
	unsigned char n;
	unsigned int mask, k;
	mask = (1UL << ((incIntArray->numTotalBits)-1)); // mask for the highest/leftmost bit of the number
	k = 0; // current position in the limits array
	n = (incIntArray->numHighBits); // current high bits count
	while( n > 0 ){
		if( num & mask ){ // bit 1, and advance to limit (2*(k+1))
			k++;
			k <<= 1;
		} else { // bit 0, and advance to limit (2*k+1)
			((incIntArray->highLimits)[k])++; // increase limit count
			k <<= 1;
			k++;
		}
		n--;
		mask >>= 1; // process next bit at the right
	}
	SetPackedNumber((incIntArray->lowBits),pos,(num & (incIntArray->lowBitsMask)));
}
