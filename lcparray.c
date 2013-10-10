/* ================================================================= *
 *  lcparray.c : Sampled LCP Array                                   *
 *                                                                   *
 *  slaMEM: MUMmer-like tool to retrieve Maximum Exact Matches using *
 *          an FM-Index and a Sampled Longest Common Prefix Array    *
 *                                                                   *
 *  Copyright (c) 2013, Francisco Fernandes <fjdf@kdbio.inesc-id.pt> *
 *  Knowledge Discovery in Bioinformatics group (KDBIO/INESC-ID)     *
 *  All rights reserved                                              *
 *                                                                   *
 *  This file is subject to the terms and conditions defined in the  *
 *  file 'LICENSE', which is part of this source code package.       *
 * ================================================================= */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "lcparray.h"
#include "bwtindex.h"

#define DEBUGLCP 1
//#define BUILDLCP 1 // if we want to build the LCP array here or use the lcparray passed as argument

#ifdef DEBUGLCP
#include <sys/timeb.h>
#endif

/**/
#define BLOCKSIZE 64
#define BLOCKMASK 63
#define BLOCKSHIFT 6
#define BLOCKHALF 32
/**/
/*
#define BLOCKSIZE 128
#define BLOCKMASK 127
#define BLOCKSHIFT 7
#define BLOCKHALF 64
*/

#define BWTBLOCKSIZE 64
#define BWTBLOCKMASK 63
#define BWTBLOCKSHIFT 6

typedef struct _LCPSamplesBlock {
	int baseBwtPos;
	int bigLCPsCount;
	int bigPLPsCount;
	unsigned char sourceLCP[BLOCKSIZE];
	signed char prefixLinkPointer[BLOCKSIZE];
} LCPSamplesBlock;

typedef struct _SampledPosMarks {
	unsigned long long int bits;
	int marksCount;
} SampledPosMarks;

#ifdef DEBUGLCP
// data structure if a representation by intervals was used instead
typedef struct _LCPIntervalTreeBlock {
	unsigned char lcpValue[BLOCKSIZE];
	unsigned char intervalSize[BLOCKSIZE];
	unsigned char parentBwtDistance[BLOCKSIZE];
	int bigValuesCount;
} LCPIntervalTreeBlock;
#endif

typedef struct _IntPair {
	int pos;
	int value;
} IntPair;

static int bwtLength;
static int *fullLCPArray;
static int *fullPLPArray;
static int *fullSuffixArray;
static SampledPosMarks *bwtMarkedPositions;
static int numLCPSamples;
static LCPSamplesBlock *sampledLCPArray;
static LCPSamplesBlock *lastLCPSamplesBlock;
static int numOversizedLCPs, numOversizedPLPs;
static int *extraLCPvalues, *extraPLPvalues;

static unsigned long long int *offsetMasks64bits;
static int *perByteCounts;

#ifdef DEBUGLCP
int testNumCalls = 0;
int testNumFollowedPos = 0;
int testMaxFollowedPos = 0;
#endif

void FreeSampledSuffixArray(){
	free(offsetMasks64bits);
	free(perByteCounts);
	if(fullSuffixArray!=NULL) free(fullSuffixArray);
	if(fullPLPArray!=NULL) free(fullPLPArray);
	if(fullLCPArray!=NULL) free(fullLCPArray);
	free(bwtMarkedPositions);
	free(sampledLCPArray);
	free(extraLCPvalues);
	free(extraPLPvalues);
}

__inline int GetLcpPosFromBwtPos(unsigned int pos){
	SampledPosMarks *bwtBlock;
	unsigned long long int bitsArray;
	bwtBlock = &(bwtMarkedPositions[ (pos >> BWTBLOCKSHIFT) ]);
	bitsArray = (bwtBlock->bits) & offsetMasks64bits[ (pos & BWTBLOCKMASK) ]; // mask = ( 1<<(offset+1) - 1 ) = 0^(64-offset-1)1^(offset+1)
	pos = (bwtBlock->marksCount); // the count is up to but NOT including the 0-th position
	/*
	while( bitsArray ){
		pos++;
		bitsArray &= ( bitsArray - 1ULL );
	}
	*/
	/*
	//bitsArray = ( bitsArray & 0x5555555555555555 ) + ( ( bitsArray >> 1 ) & 0x5555555555555555 ); // 0x5 = 0101b ; 2 bits ( final max count = 2 -> 2 bits )
	bitsArray = ( bitsArray - ( ( bitsArray >> 1 ) & 0x5555555555555555 ) );
	bitsArray = ( bitsArray & 0x3333333333333333 ) + ( ( bitsArray >> 2 ) & 0x3333333333333333 ); // 0x3 = 0011b ; 4 bits ( final max count = 4 -> 3 bits )
	bitsArray = ( ( bitsArray + ( bitsArray >> 4 ) ) & 0x0F0F0F0F0F0F0F0F ); // 0x0F = 00001111b ; 8 bits ( final max count = 8 -> 4 bits )
	//bitsArray = ( bitsArray + ( bitsArray >> 8 ) ); // the final count will be at most 64, which is 7 bits, and since we now have blocks of 8 bits, we can add them without masking because there will not be any overflow
	//bitsArray = ( bitsArray + ( bitsArray >> 16 ) );
	//bitsArray = ( ( bitsArray + ( bitsArray >> 32 ) ) & 0x000000000000007F ); // get the last 7 bits ( final max count = 64 -> 7 bits )
	bitsArray = ( ( bitsArray * 0x0101010101010101 ) >> 56 );
	return ( pos + (int)bitsArray );
	*/
	#ifdef __GNUC__
		return ( pos + __builtin_popcountll( bitsArray ) );
	#else
		pos += perByteCounts[ ( bitsArray & 0x00000000000000FF ) >> 0 ];
		pos += perByteCounts[ ( bitsArray & 0x000000000000FF00 ) >> 8 ];
		pos += perByteCounts[ ( bitsArray & 0x0000000000FF0000 ) >> 16 ];
		pos += perByteCounts[ ( bitsArray & 0x00000000FF000000 ) >> 24 ];
		pos += perByteCounts[ ( bitsArray & 0x000000FF00000000 ) >> 32 ];
		pos += perByteCounts[ ( bitsArray & 0x0000FF0000000000 ) >> 40 ];
		pos += perByteCounts[ ( bitsArray & 0x00FF000000000000 ) >> 48 ];
		pos += perByteCounts[ ( bitsArray & 0xFF00000000000000 ) >> 56 ];
		return pos;
	#endif
}

// Retrieves the LCP value from the specified Sampled LCP Array position
__inline int GetLcpValueFromLcpPos(int pos){
	LCPSamplesBlock *lcpBlock;
	int lcp, extraPos;
	lcpBlock = &(sampledLCPArray[ (pos >> BLOCKSHIFT) ]);
	pos = (pos & BLOCKMASK);
	lcp = (int)(lcpBlock->sourceLCP[pos]);
	if( lcp != UCHAR_MAX ) return lcp; // if not oversized value, return directly
	if( (pos < BLOCKHALF) || (lcpBlock == lastLCPSamplesBlock) ){ // position in the 1st half of the block
		extraPos = (lcpBlock->bigLCPsCount); // the count is up to but NOT including the 0-th position
		pos++; // to include the 0-th position
		while( pos ){
			pos--;
			if( ((lcpBlock->sourceLCP)[pos]) == UCHAR_MAX ) extraPos++;
		}
	} else { // position in the 2nd half of the block
		extraPos = (((LCPSamplesBlock *)(lcpBlock+1))->bigLCPsCount); // get the count in the next lcp block
		pos++; // if there are no large values ahead, the extra position is already extraPos
		while( pos != BLOCKSIZE ){
			if( ((lcpBlock->sourceLCP)[pos]) == UCHAR_MAX ) extraPos--;
			pos++;
		}
	}
	return extraLCPvalues[extraPos];
}

int GetLCP(unsigned int bwtpos){
	int lcpPos = GetLcpPosFromBwtPos(bwtpos);
	if( !(bwtMarkedPositions[(bwtpos>>BWTBLOCKSHIFT)].bits & (1ULL<<(bwtpos&BWTBLOCKMASK))) ) lcpPos++; // if the position is not marked, its LCP is equal to the one of the marked position ahead
	return GetLcpValueFromLcpPos(lcpPos);
}

/*
int GetDepth(int lcpPos){
	int lcp, nextlcp;
	lcp = GetLcpValueFromLcpPos(lcpPos);
	lcpPos++;
	if( lcpPos == numLCPSamples ) return lcp;
	nextlcp = GetLcpValueFromLcpPos(lcpPos);
	if( nextlcp > lcp ) return nextlcp;
	return lcp;
}
*/

int IsTopCorner(int lcpPos){
	return ( (lcpPos != (numLCPSamples-1)) && (GetLcpValueFromLcpPos(lcpPos) < GetLcpValueFromLcpPos(lcpPos+1)) );
}

/*
int IsBottomCorner(int lcpPos){
	return ( (lcpPos == (numLCPSamples-1)) || (GetLcpValueFromLcpPos(lcpPos) > GetLcpValueFromLcpPos(lcpPos+1)) );
}
*/

/*
// Retrieves the position in the full array (BWT) of a position in the sampled array, with the help of a closer known BWT position
int GetBwtPosFromLcpPos(int lcpPos, int bwtPos){
	SampledPosMarks *bwtBlock;
	unsigned long long int bitMask;
	bwtPos = (bwtPos >> BWTBLOCKSHIFT); // block number
	bwtBlock = &(bwtMarkedPositions[bwtPos]);
	bwtPos = (bwtPos << BWTBLOCKSHIFT); // position at the beginning of the block
	if( lcpPos > (bwtBlock->marksCount) ){ // search down/ahead
		while( (bwtPos < bwtLength) && ((bwtBlock->marksCount) < lcpPos) ){ // get the block with the cumulative lcp count closer to our lcp position
			bwtBlock++;
			bwtPos += BWTBLOCKSIZE;
		}
		bwtBlock--;
		bwtPos -= BWTBLOCKSIZE;
	} else { // or search up/behind
		while( (bwtPos > 0) && ((bwtBlock->marksCount) >= lcpPos) ){
			bwtBlock--;
			bwtPos -= BWTBLOCKSIZE;
		}
	}
	lcpPos -= (bwtBlock->marksCount); // we are at the right block, now count this many bits set to 1
	bitMask = 1ULL;
	while( bitMask ){
		if( (bwtBlock->bits) & bitMask ){
			lcpPos--;
			if( lcpPos == 0 ) break;
		}
		bwtPos++;
		bitMask <<= 1;
	}
	return bwtPos;
}
*/

// Retrieves the position in the full array (BWT) of a position in the sampled array
int GetBwtPosFromLcpPos(int lcpPos){
	LCPSamplesBlock *lcpBlock;
	SampledPosMarks *bwtBlock;
	unsigned long long int bitMask;
	int bwtPos;
	lcpBlock = &(sampledLCPArray[ (lcpPos >> BLOCKSHIFT) ]);
	bwtPos = (lcpBlock->baseBwtPos); // BWT pos of the first LCP in this LCP block
	bwtPos = (bwtPos >> BWTBLOCKSHIFT); // BWT block number
	bwtBlock = &(bwtMarkedPositions[bwtPos]); // BWT block of that pos
	bwtPos = (bwtPos << BWTBLOCKSHIFT); // position at the beginning of the block
	bwtBlock++; // check the lcp count of the next block
	bwtPos += BWTBLOCKSIZE;
	if( (bwtPos >= bwtLength) || ((bwtBlock->marksCount) >= lcpPos) ){ // if the BWT position we want is in the initial block
		bwtBlock--;
		bwtPos = (lcpBlock->baseBwtPos); // start searching on the already known position
		lcpPos = (lcpPos & BLOCKMASK); // how many marked positions do we want ahead of that one
		if( lcpPos == 0 ) return bwtPos; // check if this is already the position we want
		bwtPos++; // that position had a sample, but it's not the one we want, so start searching on the next position (it will never be on the next block, othewise we would have (nextBwtBlock->marksCount)<lcpPos)
		bitMask = ( 1ULL << (bwtPos & BWTBLOCKMASK) ); // set the bit mask for our starting position inside the block
	} else { // the LCP sample we want is in a BWT block further ahead
		while( (bwtPos < bwtLength) && ((bwtBlock->marksCount) < lcpPos) ){ // get the block with the cumulative lcp count closer to our lcp position
			bwtBlock++;
			bwtPos += BWTBLOCKSIZE;
		}
		bwtBlock--; // when we exit the loop, we are one block ahead of the one we want
		bwtPos -= BWTBLOCKSIZE;
		lcpPos -= (bwtBlock->marksCount); // we are at the right block, now count this many bits set to 1
		bitMask = 1ULL;
	}
	while( bitMask ){ // get the position inside the block of the n-th set bit, where n is stored in the lcpPos variable
		if( (bwtBlock->bits) & bitMask ){
			lcpPos--;
			if( lcpPos == 0 ) break;
		}
		bwtPos++;
		bitMask <<= 1;
	}
	return bwtPos;
}

/*
void SetPrefixLinkPointer(int sourceBwtPos, int destBwtPos){
	int lcpPos;
	lcpPos = GetLcpPosFromBwtPos(sourceBwtPos);
	( (sampledLCPArray[(lcpPos >> BLOCKSHIFT)]).prefixLinkPointer )[(lcpPos & BLOCKMASK)] = destBwtPos;
}
*/

// Retrieves the prefix link pointer from the specified Sampled LCP Array position
unsigned int GetPrefixLinkFromLcpPos(int pos){
	LCPSamplesBlock *lcpBlock;
	int distance, bwtPos, extraPos;
	bwtPos = GetBwtPosFromLcpPos(pos);
	lcpBlock = &(sampledLCPArray[ (pos >> BLOCKSHIFT) ]);
	pos = (pos & BLOCKMASK);
	distance = (int)(lcpBlock->prefixLinkPointer[pos]);
	if( distance == 0 ){ // get the large value from the extra array
		if( (pos < BLOCKHALF) || (lcpBlock == lastLCPSamplesBlock) ){ // position in the 1st half of the block
			extraPos = (lcpBlock->bigPLPsCount); // the count is up to but NOT including the 0-th position
			pos++; // to include the 0-th position
			while( pos ){
				pos--;
				if( ((lcpBlock->prefixLinkPointer)[pos]) == 0 ) extraPos++;
			}
		} else { // position in the 2nd half of the block
			extraPos = (((LCPSamplesBlock *)(lcpBlock+1))->bigPLPsCount); // get the count in the next lcp block
			pos++; // if there are no large values ahead, the extra position is already extraPos
			while( pos != BLOCKSIZE ){
				if( ((lcpBlock->prefixLinkPointer)[pos]) == 0 ) extraPos--;
				pos++;
			}
		}
		distance = extraPLPvalues[extraPos];
	}
	return (unsigned int)( bwtPos + distance );
}

int GetEnclosingLCPInterval(unsigned int *topptr, unsigned int *bottomptr){
	int lcpPos, destDepth, n;
	unsigned int destTopPtr, destBottomPtr;
	if( (*topptr) != (*bottomptr) ){ // non unitary interval, enlarge the interval using the LCP prefix links
		lcpPos = GetLcpPosFromBwtPos((*topptr));
		destTopPtr = GetPrefixLinkFromLcpPos(lcpPos);
		destDepth = GetLcpValueFromLcpPos(lcpPos); // destination depth = max( LCP(top) , LCP(bottom+1) )
		lcpPos = GetLcpPosFromBwtPos((*bottomptr));
		destBottomPtr = GetPrefixLinkFromLcpPos(lcpPos);
		lcpPos++;
		if( lcpPos != numLCPSamples ) n = GetLcpValueFromLcpPos(lcpPos); // depth of the bottom prefix link
		else n = (-1);
		if( destDepth >= n ){ // follow the link with the highest depth and leave the other pointer as it is
			(*topptr) = destTopPtr;
		}
		if( n >= destDepth ){ // if both depths are equal, both pointers will be changed
			destDepth = n;
			(*bottomptr) = destBottomPtr;
		}
		return destDepth;
	} // else it is an interval with a single position
	destTopPtr = UINT_MAX;
	destBottomPtr = UINT_MAX;
	n = (int)(*topptr);
	if( (bwtMarkedPositions[(n >> BWTBLOCKSHIFT)].bits) & (1ULL << (n & BWTBLOCKMASK)) ){ // if this is a marked LCP position
		lcpPos = GetLcpPosFromBwtPos((*topptr));
		if( IsTopCorner(lcpPos) ){ // top corner
			destTopPtr = (*topptr);
			destDepth = GetLcpValueFromLcpPos(lcpPos+1); // source depth of top corner = LCP(i+1)
		} else { // bottom corner
			destBottomPtr = (*topptr);
			destDepth = GetLcpValueFromLcpPos(lcpPos); // source depth of bottom corner = LCP(i)
		}
	} else { // this bwt pos is not marked by a LCP
		lcpPos = GetLcpPosFromBwtPos((*topptr)); // it will get the closest LCP before/above this bwt pos
		destDepth = GetLcpValueFromLcpPos(lcpPos+1); // source depth is in the next LCP
		if( IsTopCorner(lcpPos) ) destTopPtr = GetBwtPosFromLcpPos(lcpPos); // top corner at the left
		else destBottomPtr = GetPrefixLinkFromLcpPos(lcpPos); // if it's a bottom corner from an interval at the left, set our destination bottom corner
		lcpPos++; // check the next/bellow corner
		if( IsTopCorner(lcpPos) ){ // top corner at the right
			if(destTopPtr==UINT_MAX) destTopPtr = GetPrefixLinkFromLcpPos(lcpPos); // use the prefix link at this top corner at the right to set our top corner
			lcpPos--; // set the LCP pos to the last (possibly only) defined corner
		} else { // bottom corner at the right
			if(destBottomPtr==UINT_MAX) destBottomPtr = GetBwtPosFromLcpPos(lcpPos); // if it's a bottom corner, from an interval at the left, set our destination bottom corner
		}
	}
	if( destTopPtr == UINT_MAX ){ // if we still need to get the top pointer
		lcpPos--; // we are at the bottom pointer, so, go to the left/above
		/*
		while( (n=GetLcpValueFromLcpPos(lcpPos)) > destDepth ) lcpPos--; // get the pos at the left with the LCP lower or equal to ours
		if( n < destDepth ) destTopPtr = GetBwtPosFromLcpPos(lcpPos); // this is the top pos
		else destTopPtr = GetPrefixLinkFromLcpPos(lcpPos); // if it has the same LCP as ours, use its prefix link pointer
		*/
		while( !IsTopCorner(lcpPos) ) lcpPos--; // find the first top corner above
		if( GetLcpValueFromLcpPos(lcpPos) < destDepth ) destTopPtr = GetBwtPosFromLcpPos(lcpPos); // check if the first found top corner is already the one we want
		else { // keep following prefix links until we find a top corner with an LCP value lower than our destination LCP
			destTopPtr = GetPrefixLinkFromLcpPos(lcpPos);
			lcpPos = GetLcpPosFromBwtPos(destTopPtr);
			while( GetLcpValueFromLcpPos(lcpPos) >= destDepth ){
				destTopPtr = GetPrefixLinkFromLcpPos(lcpPos);
				lcpPos = GetLcpPosFromBwtPos(destTopPtr);
			}
		}
	} else if( destBottomPtr == UINT_MAX ){ // if we still need to get the bottom pointer
		lcpPos++; // we are at the top pointer, so, go to the right/bellow
		/*
		if( lcpPos != numLCPSamples ) lcpPos++; // if we are not at the last position, go another position to the right
		n = -1; // set in case we were at the 2 last positions
		while( (lcpPos < numLCPSamples) && ((n=GetLcpValueFromLcpPos(lcpPos)) > destDepth) ) lcpPos++; // get the pos at the right with the LCP lower or equal to ours
		if( (n < destDepth) || (lcpPos == numLCPSamples) ) destBottomPtr = GetBwtPosFromLcpPos((lcpPos-1)); // the previous pos is the top pos
		else { // if it has the same LCP as ours
			if( IsTopCorner(lcpPos) ) destBottomPtr = GetPrefixLinkFromLcpPos(lcpPos-1);// if it's a top corner, use the prefix link of the prev pos, which is a bottom corner
			else destBottomPtr = GetBwtPosFromLcpPos(lcpPos); // if it's a bottom corner, use it
		}
		*/
		while( IsTopCorner(lcpPos) ) lcpPos++; // find the first bottom corner bellow
		if( (lcpPos == numLCPSamples) || (GetLcpValueFromLcpPos(lcpPos+1) < destDepth) ) destBottomPtr = GetBwtPosFromLcpPos(lcpPos); // check if the first found bottom corner is already the one we want
		else {
			destBottomPtr = GetPrefixLinkFromLcpPos(lcpPos); // keep following prefix links until we find a bottom corner that has a position ahead with an LCP value lower than our destination LCP
			lcpPos = GetLcpPosFromBwtPos(destBottomPtr);
			while( (lcpPos != (numLCPSamples-1)) && (GetLcpValueFromLcpPos(lcpPos+1) >= destDepth) ){
				destBottomPtr = GetPrefixLinkFromLcpPos(lcpPos);
				lcpPos = GetLcpPosFromBwtPos(destBottomPtr);
			}
		}
	}
	(*topptr) = destTopPtr;
	(*bottomptr) = destBottomPtr;
	return destDepth;
}

// Get the interval whose depth is equal to the LCP value at lcpPos and that contains that BWT position
int GetEnclosingLCPIntervalFromLCPPos(int lcpPos, unsigned int *topptr, unsigned int *bottomptr){
	int destDepth;
	unsigned int destTopPtr, destBottomPtr;
	destTopPtr = UINT_MAX;
	destBottomPtr = UINT_MAX;
	destDepth = GetLcpValueFromLcpPos(lcpPos);
	if( IsTopCorner(lcpPos) ){ // pos is a top corner
		destTopPtr = GetPrefixLinkFromLcpPos(lcpPos); // the top corner of the interval is found by following the prefix link
		lcpPos++; // we are at the top pointer, so, go to the right/down
		/*
		if( lcpPos != numLCPSamples ) lcpPos++; // if we are not at the last position, go 2 positions to the right/down
		n = -1; // set in case we were at the 2 last positions
		while( (lcpPos < numLCPSamples) && ((n=GetLcpValueFromLcpPos(lcpPos)) > destDepth) ) lcpPos++; // get the pos at the right/down with the LCP lower or equal to ours
		if( (n < destDepth) || (lcpPos == numLCPSamples) ) destBottomPtr = GetBwtPosFromLcpPos((lcpPos-1)); // the previous pos is the top pos
		else { // if it has the same LCP as ours
			if( IsTopCorner(lcpPos) ) destBottomPtr = GetPrefixLinkFromLcpPos(lcpPos-1);// if it's a top corner, use the prefix link of the prev pos, which is a bottom corner
			else destBottomPtr = GetBwtPosFromLcpPos(lcpPos); // if it's a bottom corner, use it
		}
		*/
		while( IsTopCorner(lcpPos) ) lcpPos++; // find the first bottom corner bellow
		if( (lcpPos == numLCPSamples) || (GetLcpValueFromLcpPos(lcpPos+1) < destDepth) ) destBottomPtr = GetBwtPosFromLcpPos(lcpPos); // check if the first found bottom corner is already the one we want
		else {
			destBottomPtr = GetPrefixLinkFromLcpPos(lcpPos); // keep following prefix links until we find a bottom corner that has a position ahead with an LCP value lower than our destination LCP
			lcpPos = GetLcpPosFromBwtPos(destBottomPtr);
			while( (lcpPos != (numLCPSamples-1)) && (GetLcpValueFromLcpPos(lcpPos+1) >= destDepth) ){
				destBottomPtr = GetPrefixLinkFromLcpPos(lcpPos);
				lcpPos = GetLcpPosFromBwtPos(destBottomPtr);
			}
		}
	} else { // pos is a bottom corner
		destBottomPtr = GetBwtPosFromLcpPos(lcpPos); // this is already the bottom corner of the interval
		lcpPos--; // we are at the bottom pointer, so, go to the left/up
		/*
		while( (n=GetLcpValueFromLcpPos(lcpPos)) > destDepth ) lcpPos--; // get the pos at the left with the LCP lower or equal to ours
		if( n < destDepth ) destTopPtr = GetBwtPosFromLcpPos(lcpPos); // this is the top pos
		else destTopPtr = GetPrefixLinkFromLcpPos(lcpPos); // if it has the same LCP as ours, use its prefix link pointer
		*/
		while( !IsTopCorner(lcpPos) ) lcpPos--; // find the first top corner above
		if( GetLcpValueFromLcpPos(lcpPos) < destDepth ) destTopPtr = GetBwtPosFromLcpPos(lcpPos); // check if the first found top corner is already the one we want
		else { // keep following prefix links until we find a top corner with an LCP value lower than our destination LCP
			destTopPtr = GetPrefixLinkFromLcpPos(lcpPos);
			lcpPos = GetLcpPosFromBwtPos(destTopPtr);
			while( GetLcpValueFromLcpPos(lcpPos) >= destDepth ){
				destTopPtr = GetPrefixLinkFromLcpPos(lcpPos);
				lcpPos = GetLcpPosFromBwtPos(destTopPtr);
			}
		}
	}
	(*topptr) = destTopPtr;
	(*bottomptr) = destBottomPtr;
	return destDepth;
}

int GetLCPValuePositionInsideLCPInterval(int lcpvalue, unsigned int topptr, unsigned int bottomptr){
	int startlcppos, endlcppos;
	#ifdef DEBUGLCP
	int testCount = 0;
	testNumCalls++;
	testCount++;
	testNumFollowedPos++;
	#endif
	endlcppos=GetLcpPosFromBwtPos(bottomptr);
	if(GetLcpValueFromLcpPos(endlcppos)==lcpvalue) return endlcppos;
	#ifdef DEBUGLCP
	testCount++;
	testNumFollowedPos++;
	#endif
	startlcppos=GetLcpPosFromBwtPos(topptr);
	startlcppos++;
	if(GetLcpValueFromLcpPos(startlcppos)==lcpvalue) return startlcppos;
	for( (++startlcppos) ; startlcppos<endlcppos ; startlcppos++ ){
		#ifdef DEBUGLCP
		testCount++;
		testNumFollowedPos++;
		#endif
		if(GetLcpValueFromLcpPos(startlcppos)==lcpvalue){
			#ifdef DEBUGLCP
			if(testCount>testMaxFollowedPos) testMaxFollowedPos=testCount;
			#endif
			return startlcppos;
		}
	}
	return (-1);
}

#ifdef DEBUGLCP
int GetTrueDepth(unsigned int topptr, unsigned int bottomptr){
	int depth; // non-singular interval: source depth = min( LCP[top+1] , LCP[bottom] )
	depth=fullLCPArray[bottomptr];
	topptr++;
	if( topptr!=(unsigned int)bwtLength ){
		if( (++bottomptr)==topptr ){
			if( fullLCPArray[topptr]>depth ) depth=fullLCPArray[topptr]; // single position: source depth = max( LCP[pos] , LCP[pos+1] )
		} else {
			if( fullLCPArray[topptr]<depth ) depth=fullLCPArray[topptr];
		}
	}
	return depth;
}
#endif

#ifdef DEBUGLCP
int GetTrueEnclosingLCPInterval(unsigned int *topptr, unsigned int *bottomptr){
	int destdepth; // destination (parent) depth = max( LCP[top] , LCP[bottom+1] )
	destdepth=fullLCPArray[(*topptr)];
	(*bottomptr)++;
	if( (*bottomptr)!=(unsigned int)bwtLength && fullLCPArray[(*bottomptr)]>destdepth ) destdepth=fullLCPArray[(*bottomptr)];
	while( fullLCPArray[(*topptr)]>=destdepth ) (*topptr)--; // find closest pos above with an LCP value lower than this one
	while( (*bottomptr)!=(unsigned int)bwtLength && fullLCPArray[(*bottomptr)]>=destdepth) (*bottomptr)++; // find closest pos bellow with an LCP value lower than this one
	(*bottomptr)--; // go back/up one pos
	return destdepth;
}
#endif

#ifdef DEBUGLCP
int GetTrueEnclosingLCPIntervalWithLcpValue(int lcp, unsigned int *topptr, unsigned int *bottomptr){
	while( (*topptr)!=0 && fullLCPArray[(*topptr)]>=lcp ) (*topptr)--; // find closest pos above with an LCP value lower than this one
	while( (*bottomptr)!=(unsigned int)bwtLength && fullLCPArray[(*bottomptr)]>=lcp) (*bottomptr)++; // find closest pos bellow with an LCP value lower than this one
	(*bottomptr)--; // go back/up one pos
	return lcp;
}
#endif

int CompareIntPair(const void *a, const void * b){
	return ( (((IntPair *)a)->pos) - (((IntPair *)b)->pos) );
}

// TODO: get statistics for number of intervals that do not include all the chars of the alphabet in its BWT range
// TODO: implement SLCP+SV as "lcp-interval-tree":
//  - store info of top corners only
//  - multiple entries if it's a shared top corner
//  - check (parentDistance==0) to distinguish between new corner or next entry of same corner
//  - lcp-value(s), interval size(s), distance(s) to parent interval's pos in BWT (if deep corner, parentDistance=0)
//  - benchmark avg+max results for these 3 fields on large datasets
int BuildSampledLCPArray(char *text, int textsize, unsigned char *lcparray, int verbose){
	int textpos, prevtextpos;
	char *topstring, *bottomstring;
	unsigned char byte;
	unsigned int bwtpos;
	int i, k, lcp, prevlcp, nextlcp, lcppos;
	int numTopCorners, numBottomCorners;
	int maxTopCorners, maxBottomCorners;
	IntPair *topCorners, *bottomCorners;
	IntPair *oversizedCorners;
	unsigned long long int mask;
	SampledPosMarks *bwtBlock;
	LCPSamplesBlock *lcpBlock;
	int maxNumLCPSamples, maxNumOversizedValues;
	long long int sumValues;
	int maxValue, prefixLinkDistance;
	int progressCounter, progressStep;
	#ifdef DEBUGLCP
	unsigned int topptr, bottomptr;
	unsigned int stopptr, sbottomptr;
	struct timeb startTime, endTime;
	double elapsedTime;
	int numLcpIntervals;
	int numBigTopLcps, numBigTopPlps, numBigTopSizes;
	int numSharedTopCorners, avgSharedTopCornersCount, maxSharedTopCornersCount;
	int *sharedTopCornersCount;
	int numIncompleteLcpIntervals;
	unsigned char *charsInsideInterval;
	unsigned char bwtCharMask;
	#endif
	offsetMasks64bits = (unsigned long long int *)malloc(64*sizeof(unsigned long long int));
	mask = 1ULL;
	for(i=0;i<64;i++){
		offsetMasks64bits[i] = mask;
		mask <<= 1;
		mask |= 1ULL;
	}
	perByteCounts = (int *)malloc(256*sizeof(int));
	for(i=0;i<256;i++){
		k = 0;
		byte = (unsigned char)i;
		while( byte ){
			k++;
			byte &= ( byte - (unsigned char)1 );
		}
		perByteCounts[i] = k;
	}
	bwtLength = (textsize+1);
	k = (((bwtLength-1)>>BWTBLOCKSHIFT)+1); // last valid pos, quotient, add one
	bwtMarkedPositions = (SampledPosMarks *)malloc(k*sizeof(SampledPosMarks));
	sampledLCPArray = (LCPSamplesBlock *)malloc(1*sizeof(LCPSamplesBlock));
	extraLCPvalues = (int *)malloc(1*sizeof(int));
	fullSuffixArray = NULL;
	fullPLPArray = NULL;
	fullLCPArray = NULL;
	#ifdef BUILDLCP
	lcparray = NULL; // to fix compiler unused variable warning
	#endif
	/*
	printf("> Building full Suffix Array ... ");
	fflush(stdout);
	fullSuffixArray=(int *)malloc((textsize+1)*sizeof(int));
	n=0;
	bwtpos=0;
	bottomptr=0;
	fullSuffixArray[0]=textsize;
	for(textpos=(textsize-1);textpos>0;textpos--){
		n=(int)FMI_FollowLetter(text[textpos],&bwtpos,&bottomptr);
		if(n!=1) exitMessage("Text matching error");
		fullSuffixArray[bwtpos]=textpos;
	}
	printf("OK\n");
	*/
	if(verbose){ printf("> Building Sampled LCP Array "); fflush(stdout); }
	sumValues=0;
	maxValue=0;
	#ifdef DEBUGLCP
	fullLCPArray=(int *)malloc(bwtLength*sizeof(int));
	fullLCPArray[0]=(-1);
	#endif
	progressStep=(bwtLength/10);
	progressCounter=0;
	k=(-1); // current bwt block id
	bwtBlock=&(bwtMarkedPositions[0]); // first bwt block
	mask=0ULL; // current position mask inside the current bwt block
	lcpBlock=&(sampledLCPArray[0]); // current lcp block
	numLCPSamples=0;
	lastLCPSamplesBlock=0;
	numOversizedLCPs=0;
	maxNumLCPSamples=0;
	maxNumOversizedValues=0;
	prevlcp=(-1); // at the 0-th position there's no position before
	prevtextpos=textsize; // the 0-th BWT position corresponds to the terminator symbol at the last text position
	textpos=0;
	for(bwtpos=1;bwtpos<=(unsigned int)bwtLength;bwtpos++){ // all BWT positions and one fake next-to-last position to fill last position
		if(verbose){
			if(progressCounter==progressStep){ // print progress dots
				putchar('.');
				fflush(stdout);
				progressCounter=0;
			} else progressCounter++;
		}
		if(bwtpos!=(unsigned int)bwtLength){
			//topstring=(char *)(text+fullSuffixArray[(bwtpos-1)]);
			//bottomstring=(char *)(text+fullSuffixArray[bwtpos]);
			#ifdef BUILDLCP
			textpos=(int)FMI_PositionInText(bwtpos);
			topstring=(char *)(text+prevtextpos);
			bottomstring=(char *)(text+textpos);
			lcp=0;
			if(textpos>prevtextpos) prevtextpos=textpos; // use this temporary variable to check when we reach the end of the text, because the '\0' terminator might not be present
			while(((*topstring)==(*bottomstring)) && (prevtextpos!=textsize)){
				topstring++;
				bottomstring++;
				lcp++;
				prevtextpos++;
			}
			#else
			lcp = (int)lcparray[bwtpos];
			if( lcp == (int)UCHAR_MAX ){ // if it is a large (>255) lcp, calculate its value
				prevtextpos=(int)FMI_PositionInText((bwtpos-1));
				textpos=(int)FMI_PositionInText(bwtpos);
				topstring=(char *)( text + prevtextpos + (int)UCHAR_MAX ); // start checking matches 255 positions ahead
				bottomstring=(char *)( text + textpos + (int)UCHAR_MAX );
				if(textpos>prevtextpos) prevtextpos=textpos; // use to check when we reach the end of the text, because the '\0' terminator might not be present
				while(((*topstring)==(*bottomstring)) && (prevtextpos!=textsize)){
					topstring++;
					bottomstring++;
					lcp++;
					prevtextpos++;
				}
			}
			#endif
			#ifdef DEBUGLCP
			fullLCPArray[bwtpos]=lcp; // longest common prefix between positions (bwtpos) and (bwtpos-1)
			#endif
		} else lcp=(-1); // fake next-to-last position
		sumValues+=lcp;
		if(lcp>maxValue) maxValue=lcp;
		if(mask==0ULL){ // advance to new block
			k++;
			if(k!=0) bwtBlock++; // k starts at -1 and we already have the 0-th block loaded
			(bwtBlock->bits)=0ULL; // reset block
			(bwtBlock->marksCount)=(numLCPSamples-1); // the 0-th block is initialized with (-1) here
			mask=1ULL; // in the bit array, this is the position of (bwtpos-1); it will only be incremented at the end of the loop
		}
		if(lcp!=prevlcp){ // set lcp (in the sampled LCP array) of previous position, (bwtpos-1)
			(bwtBlock->bits) |= mask; // mark position (bwtpos-1) as having a sample
			if((numLCPSamples & BLOCKMASK)==0){ // create new lcp block if needed
				if(numLCPSamples==maxNumLCPSamples){
					maxNumLCPSamples += ( (textsize/10) & (~BLOCKMASK) ); // alloc 10% of the text size extra samples each time, and aligned to the block size
					sampledLCPArray = (LCPSamplesBlock *)realloc(sampledLCPArray,((maxNumLCPSamples >> BLOCKSHIFT)+1)*sizeof(LCPSamplesBlock));
				}
				lcpBlock=&(sampledLCPArray[(numLCPSamples >> BLOCKSHIFT)]); // lcpBlock++ would not work because the memory was re-allocated
				(lcpBlock->bigLCPsCount)=(numOversizedLCPs-1);
				(lcpBlock->baseBwtPos)=(bwtpos-1);
			}
			if(prevlcp!=(-1) && prevlcp<UCHAR_MAX) (lcpBlock->sourceLCP)[(numLCPSamples & BLOCKMASK)]=(unsigned char)prevlcp;
			else { // store oversized lcps in a separate array
				(lcpBlock->sourceLCP)[(numLCPSamples & BLOCKMASK)]=UCHAR_MAX;
				if(numOversizedLCPs==maxNumOversizedValues){
					maxNumOversizedValues += (textsize/200); // alloc 0.5% of the text size extra samples each time
					extraLCPvalues=(int *)realloc(extraLCPvalues,(maxNumOversizedValues+1)*sizeof(int));
				}
				extraLCPvalues[numOversizedLCPs]=prevlcp;
				numOversizedLCPs++;
			}
			numLCPSamples++;
		}
		mask <<= 1;
		prevlcp=lcp;
		prevtextpos=textpos;
	}
	sampledLCPArray = (LCPSamplesBlock *)realloc(sampledLCPArray,((numLCPSamples >> BLOCKSHIFT)+1)*sizeof(LCPSamplesBlock));
	extraLCPvalues=(int *)realloc(extraLCPvalues,numOversizedLCPs*sizeof(int)); // shrink the allocated arrays to fit the final number of samples
	lastLCPSamplesBlock = &(sampledLCPArray[(numLCPSamples >> BLOCKSHIFT)]);
	if(verbose){
		printf(" OK\n");
		printf(":: %.2lf%% samples (%d of %d)\n",((double)numLCPSamples/(double)bwtLength)*100.0,numLCPSamples,bwtLength);
		printf(":: %.2lf%% oversized samples (%d of %d)\n",((double)numOversizedLCPs/(double)numLCPSamples)*100.0,numOversizedLCPs,numLCPSamples);
		printf(":: Average LCP value = %d (max=%d)\n",(int)(sumValues/(long long)bwtLength),maxValue);
	}
	#ifdef DEBUGLCP
	if(verbose){ printf("> Testing Sampled LCP Array "); fflush(stdout); }
	ftime(&startTime);
	k=(-1);
	for(bwtpos=0;bwtpos<(unsigned int)bwtLength;bwtpos++){
		if(verbose){
			if(progressCounter==progressStep){ // print progress dots
				putchar('.');
				fflush(stdout);
				progressCounter=0;
			} else progressCounter++;
		}
		if(GetLCP(bwtpos)!=fullLCPArray[bwtpos]) break;
		if(((bwtpos+1)==(unsigned int)bwtLength) || (fullLCPArray[bwtpos]!=fullLCPArray[bwtpos+1])){ // only check marked positions
			k=GetLcpPosFromBwtPos(bwtpos);
			if(GetBwtPosFromLcpPos(k)!=(int)bwtpos) break;
		}
	}
	if((bwtpos==(unsigned int)bwtLength) && (k==(numLCPSamples-1))){
		ftime(&endTime);
		elapsedTime = ( ((endTime.time) + (endTime.millitm)/1000.0) - ((startTime.time) + (startTime.millitm)/1000.0) );
		if(verbose){
			printf(" OK\n");
			printf(":: Done in %.3lf seconds\n",elapsedTime);
		}
	}
	else printf("\n> ERROR: sampledLCP[%u]=%d =!= fullLCP[%u]=%d (bwtPos=%u->lcpPos=%d->bwtPos=%d) \n",bwtpos,GetLCP(bwtpos),bwtpos,fullLCPArray[bwtpos],bwtpos,k,GetBwtPosFromLcpPos(k));
	#endif
	if(verbose){ printf("> Collecting Previous/Next Smaller Values "); fflush(stdout); }
	#ifdef DEBUGLCP
	fullPLPArray=(int *)malloc(numLCPSamples*sizeof(int));
	fullPLPArray[0]=0;
	fullPLPArray[(numLCPSamples-1)]=(bwtLength-1);
	#endif
	progressStep=(numLCPSamples/10);
	progressCounter=0;
	maxNumOversizedValues=(numLCPSamples/100); // allocate 1% of the total number of sampled values
	oversizedCorners=(IntPair *)malloc((maxNumOversizedValues+1)*sizeof(IntPair)); // +1 if textsize is <100 and value above gives 0
	sampledLCPArray[0].prefixLinkPointer[0]=0; // set value for first pos
	//sampledLCPArray[0].bigPLPsCount=(-1);
	oversizedCorners[0].pos=0;
	oversizedCorners[0].value=0;
	numOversizedPLPs=1; // first value already set
	sumValues=0;
	maxValue=0;
	maxTopCorners=1;
	maxBottomCorners=1;
	topCorners=(IntPair *)malloc(maxTopCorners*sizeof(IntPair));
	bottomCorners=(IntPair *)malloc(maxBottomCorners*sizeof(IntPair));
	numTopCorners = 0;
	numBottomCorners = 0;
	#ifdef DEBUGLCP
	numLcpIntervals=0;
	numIncompleteLcpIntervals=0;
	numBigTopLcps=0;
	numBigTopPlps=0;
	numBigTopSizes=0;
	numSharedTopCorners=0;
	avgSharedTopCornersCount=0;
	maxSharedTopCornersCount=0;
	sharedTopCornersCount=(int *)malloc(maxTopCorners*sizeof(int));
	charsInsideInterval=(unsigned char *)malloc(maxTopCorners*sizeof(unsigned char));
	prefixLinkDistance=0;
	#endif
	bwtBlock=&(bwtMarkedPositions[0]);
	mask=1ULL; // mask for offset in current BWT block
	nextlcp=(-1); // this will be the lcp at the 0-th position
	bwtpos=0;
	for( lcppos=0 ; lcppos<numLCPSamples ; lcppos++ ){ // process all LCP samples
		if(verbose){
			if(progressCounter==progressStep){ // print progress dots
				putchar('.');
				fflush(stdout);
				progressCounter=0;
			} else progressCounter++;
		}
		lcp = nextlcp;
		if(lcppos!=(numLCPSamples-1)) nextlcp = GetLcpValueFromLcpPos(lcppos+1);
		else nextlcp=(-1); // simulate next to last LCP
		/*
		if( (lcppos & BLOCKMASK) == 0 ){ // if we are at a new LCP block, save the current extra oversized values count
			if(lcppos!=0) sampledLCPArray[(lcppos >> BLOCKSHIFT)].bigPLPsCount = (numOversizedPLPs-1); // -1 to set the last used position ;  the 0-th position will have (-1) and it was already set before
		}
		*/
		while( ((bwtBlock->bits) & mask)==0 ){ // advance to next marked BWT position
			if(mask==0ULL){
				bwtBlock++;
				mask=1ULL;
				continue;
			}
			#ifdef DEBUGLCP
			bwtCharMask=(unsigned char)FMI_GetCharAtBWTPos(bwtpos);
			if(bwtCharMask=='A') bwtCharMask=0x03;
			else if(bwtCharMask=='C') bwtCharMask=0x0C;
			else if(bwtCharMask=='G') bwtCharMask=0x30;
			else if(bwtCharMask=='T') bwtCharMask=0xC0;
			else bwtCharMask=0x00;
			k=(numTopCorners-1);
			while( k!=(-1) && charsInsideInterval[k]!=0xFF ) charsInsideInterval[k--]|=bwtCharMask;
			#endif
			mask<<=1;
			bwtpos++;
		}
		#ifdef DEBUGLCP
		bwtCharMask=(unsigned char)FMI_GetCharAtBWTPos(bwtpos);
		if(bwtCharMask=='A') bwtCharMask=0x03;
		else if(bwtCharMask=='C') bwtCharMask=0x0C;
		else if(bwtCharMask=='G') bwtCharMask=0x30;
		else if(bwtCharMask=='T') bwtCharMask=0xC0;
		else bwtCharMask=0x00;
		k=(numTopCorners-1);
		while( k!=(-1) && charsInsideInterval[k]!=0xFF ) charsInsideInterval[k--]|=bwtCharMask;
		#endif
		if( nextlcp > lcp ){ // top corner
			lcpBlock = &(sampledLCPArray[(lcppos >> BLOCKSHIFT)]);
			i = (lcppos & BLOCKMASK);
			k = (numTopCorners-1); // the closer upper pos with LCP lower than current LCP is always the last one in the array, because if the LCP was >= it was already removed when going through the bottom corners
			/*
			while( k!=-1 && topCornersLcp[k]>=lcp ) k--; // get closer upper pos with LCP lower than current LCP
			*/
			if( k!=-1 ){ // if we are not at the first position (LCP[0]=-1)
				#ifdef DEBUGLCP
				fullPLPArray[lcppos]=(topCorners[k].pos);
				#endif
				//(lcpBlock->prefixLinkPointer)[i] = topCornersPos[k]; // link current pos with LCP pos above
				prefixLinkDistance = ( topCorners[k].pos - (int)bwtpos ); // =(destination-source)<0
				if( (prefixLinkDistance>(-128)) && (prefixLinkDistance<128) ) (lcpBlock->prefixLinkPointer)[i] = (signed char)prefixLinkDistance; // link current pos with LCP pos above
				else { // if value is <=(-128) or >=(+128) store it in the extra array
					if( numOversizedPLPs == maxNumOversizedValues ){ // realloc array if needed
						maxNumOversizedValues += (numLCPSamples/100); // allocate 1% of the total number of sampled values each time
						oversizedCorners = (IntPair *)realloc(oversizedCorners,(maxNumOversizedValues+1)*sizeof(IntPair));
					}
					(lcpBlock->prefixLinkPointer)[i] = 0; // mark as oversized corner prefix link
					oversizedCorners[numOversizedPLPs].pos = lcppos;
					oversizedCorners[numOversizedPLPs].value = prefixLinkDistance; // add oversized prefix links to separate array
					numOversizedPLPs++;
				}
				if(prefixLinkDistance<0) prefixLinkDistance=(-prefixLinkDistance);
				if(prefixLinkDistance>maxValue) maxValue=prefixLinkDistance;
				sumValues+=prefixLinkDistance;
				//SetPrefixLinkPointer(bwtpos,topCornersPos[k]); // link current pos with lower LCP pos
				//indexLCPInfo[bwtpos].prefixLinkPointer = topCornersPos[k]; // link current pos with lower LCP pos
				//indexLCPInfo[bwtpos].destinationLCP = lcp; // the destination LCP is the highest among the two
			}
			if( numTopCorners == maxTopCorners ){ // realloc arrays if needed
				maxTopCorners++;
				topCorners=(IntPair *)realloc(topCorners,maxTopCorners*sizeof(IntPair));
				#ifdef DEBUGLCP
				sharedTopCornersCount=(int *)realloc(sharedTopCornersCount,maxTopCorners*sizeof(int));
				charsInsideInterval=(unsigned char *)realloc(charsInsideInterval,maxTopCorners*sizeof(unsigned char));
				#endif
			}
			topCorners[numTopCorners].pos = bwtpos; // add this pos to the list of top corners
			topCorners[numTopCorners].value = lcp;
			#ifdef DEBUGLCP
			sharedTopCornersCount[numTopCorners]=0;
			charsInsideInterval[numTopCorners]=bwtCharMask;
			if(prefixLinkDistance>=255) numBigTopPlps++;
			#endif
			numTopCorners++;
		} else if( nextlcp < lcp ){ // bottom corner
			lcp = nextlcp; // now the LCP to consider is the next/bottom one
			k = (numBottomCorners-1);
			while( k!=-1 && (bottomCorners[k].value)>lcp ){ // get all pos behind/above with LCP higher than current LCP
				i = GetLcpPosFromBwtPos((unsigned int)(bottomCorners[k].pos));
				oversizedCorners[numOversizedPLPs].pos = i; // set this here to prevent another call to get the SLCP position; if the PLP is not oversized, it will be replaced by the next oversized one
				#ifdef DEBUGLCP
				fullPLPArray[i]=(int)bwtpos;
				#endif
				lcpBlock = &(sampledLCPArray[(i >> BLOCKSHIFT)]);
				i = (i & BLOCKMASK);
				//(lcpBlock->prefixLinkPointer)[i] = bwtpos; // link LCP pos above with current pos
				prefixLinkDistance = ( (int)bwtpos - (bottomCorners[k].pos) ); // =(destination-source)>0
				if( (prefixLinkDistance>(-128)) && (prefixLinkDistance<128) ) (lcpBlock->prefixLinkPointer)[i] = (signed char)prefixLinkDistance; // link LCP pos above with current pos
				else { // if value is <=(-128) or >=(+128) store it in the extra array
					if( numOversizedPLPs == maxNumOversizedValues ){ // realloc array if needed
						maxNumOversizedValues += (numLCPSamples/100); // allocate 1% of the total number of sampled values each time
						oversizedCorners = (IntPair *)realloc(oversizedCorners,(maxNumOversizedValues+1)*sizeof(IntPair));
					}
					(lcpBlock->prefixLinkPointer)[i] = 0; // mark as oversized corner prefix link
					//oversizedCorners[numOversizedPLPs].pos = GetLcpPosFromBwtPos((unsigned int)(bottomCorners[k].pos));
					oversizedCorners[numOversizedPLPs].value = prefixLinkDistance; // add oversized prefix links to separate array
					numOversizedPLPs++;
				}
				if(prefixLinkDistance<0) prefixLinkDistance=(-prefixLinkDistance);
				if(prefixLinkDistance>maxValue) maxValue=prefixLinkDistance;
				sumValues+=prefixLinkDistance;
				//SetPrefixLinkPointer(bottomCornersPos[k],bwtpos); // link higher LCP pos with current pos
				//indexLCPInfo[(bottomCornersPos[k])].prefixLinkPointer = bwtpos; // link higher LCP pos with current pos
				//indexLCPInfo[(bottomCornersPos[k])].destinationLCP = bottomCornersLcp[k]; // the destination LCP is the highest among the two
				numBottomCorners--; // delete this corner because it has already been linked to another, and no other is going to link to him
				k--;
			}
			if( numBottomCorners == maxBottomCorners ){ // realloc arrays if needed
				maxBottomCorners++;
				bottomCorners=(IntPair *)realloc(bottomCorners,maxBottomCorners*sizeof(IntPair));
			}
			bottomCorners[numBottomCorners].pos = bwtpos; // add this pos to the list of bottom corners
			bottomCorners[numBottomCorners].value = lcp;
			numBottomCorners++;
			i = 0; // number of deleted top corners this time
			k = (numTopCorners-1);
			while( k!=-1 && (topCorners[k].value)>=lcp ){ // get all top corner pos behind/above with LCP higher or equal than current LCP
				numTopCorners--; // delete this top corner because it has already been linked to another, and no other is going to link to him
				#ifdef DEBUGLCP
				sharedTopCornersCount[k]++;
				if(sharedTopCornersCount[k]!=1){
					numSharedTopCorners++;
					avgSharedTopCornersCount+=sharedTopCornersCount[k];
					if(sharedTopCornersCount[k]>maxSharedTopCornersCount) maxSharedTopCornersCount=sharedTopCornersCount[k];
				}
				if(topCorners[k].value>=255) numBigTopLcps++;
				if((bwtpos-topCorners[k].pos)>=255) numBigTopSizes++;
				if(charsInsideInterval[k]!=0xFF) numIncompleteLcpIntervals++;
				#endif
				k--;
				i++;
			}
			#ifdef DEBUGLCP
			numLcpIntervals += i; // if we deleted i top corners, it means i intervals were closed
			if( i==0 || (topCorners[k+1].value)!=lcp ){ // if it's a shared top margin, it's not deleted yet (i=0), but we closed another interval
				numLcpIntervals++; // if we did not end at the same level as the last deleted top margin, we are in the middle of another (shared) top margin, which means another interval
				sharedTopCornersCount[k]++;
				if(topCorners[k].value>=255) numBigTopLcps++;
				if((bwtpos-topCorners[k].pos)>=255) numBigTopSizes++;
			}
			#endif
		}
		/*
		else { // lcp == nextlcp
			//indexLCPInfo[bwtpos].prefixLinkPointer = -1; // undefined values
			//indexLCPInfo[bwtpos].destinationLCP = -1;
		}
		*/
		mask <<= 1;
		bwtpos++;
	}
	if( numTopCorners!=0 || numBottomCorners!=1 ){
		printf("\n> ERROR: Bad corners connection\n");
	}
	free(topCorners);
	free(bottomCorners);
	#ifdef DEBUGLCP
	free(sharedTopCornersCount);
	free(charsInsideInterval);
	#endif
	lcppos = (numLCPSamples-1);  // set value for last pos
	sampledLCPArray[(lcppos >> BLOCKSHIFT)].prefixLinkPointer[(lcppos & BLOCKMASK)] = 0;
	oversizedCorners[numOversizedPLPs].pos = lcppos;
	oversizedCorners[numOversizedPLPs].value = 0; // zero distance from (bwtLength-1);
	numOversizedPLPs++;
	qsort(oversizedCorners,numOversizedPLPs,sizeof(IntPair),CompareIntPair); // sort oversized PLP values by their position in the SLCP array
	extraPLPvalues = (int *)malloc(numOversizedPLPs*sizeof(int)); // final array of oversized values
	lcppos = (-1);
	lcpBlock = &(sampledLCPArray[0]);
	lcpBlock--; // start before the first position so the 0-th block gets updated in the first iteration of the loop
	for( i=0 ; i<numOversizedPLPs ; i++ ){ // fill the extra values array and the current extra values count on each LCP block
		extraPLPvalues[i] = oversizedCorners[i].value;
		k = oversizedCorners[i].pos;
		while( lcppos < k ){
			lcppos++;
			if( (lcppos & BLOCKMASK) == 0 ){ // if we reached the end of this block, go to next block and save the current count of oversized PLPs
				lcpBlock++;
				(lcpBlock->bigPLPsCount) = (i-1); // position in the extra array of the last oversized value before this position in the lcp array
			}
		}
	}
	free(oversizedCorners);
	//SetPrefixLinkPointer(0,0);
	//SetPrefixLinkPointer(textsize,textsize);
	//indexLCPInfo[0].sourceLCP = 0; // set first and last pos values
	//indexLCPInfo[0].destinationLCP = 0;
	//indexLCPInfo[0].prefixLinkPointer = 0;
	//indexLCPInfo[textsize].destinationLCP = 0;
	//indexLCPInfo[textsize].prefixLinkPointer = textsize;
	if(verbose){
		printf(" OK\n");
		printf(":: %.2lf%% oversized values (%d of %d)\n",((double)numOversizedPLPs/(double)numLCPSamples)*100.0,numOversizedPLPs,numLCPSamples);
		printf(":: Average SV distance = %.2lf (max=%d)\n",((double)sumValues/(double)numLCPSamples),maxValue);
		k = (int)( sizeof(*bwtMarkedPositions)*(((bwtLength-1)>>BWTBLOCKSHIFT)+1) + sizeof(*sampledLCPArray)*((numLCPSamples>>BLOCKSHIFT)+1) + sizeof(*extraLCPvalues)*numOversizedLCPs + sizeof(*extraPLPvalues)*numOversizedPLPs );
		printf(":: Total SLCP+SV structure size = %.1lf MB (%d bytes)\n", (double)( (unsigned int)((((bwtLength-1)>>BWTBLOCKSHIFT)+1)*sizeof(SampledPosMarks)) + (unsigned int)((((numLCPSamples-1)>>BLOCKSHIFT)+1)*sizeof(LCPSamplesBlock)) + (unsigned int)((numOversizedLCPs+numOversizedPLPs)*sizeof(int)) ) / (double)1000000U , k);
		#ifdef DEBUGLCP
		k = (int)( sizeof(SampledPosMarks)*(((bwtLength-1)>>BWTBLOCKSHIFT)+1) + sizeof(LCPIntervalTreeBlock)*((numLcpIntervals>>BLOCKSHIFT)+1) + sizeof(*extraLCPvalues)*numOversizedLCPs + sizeof(int)*(numBigTopLcps+numBigTopPlps+numBigTopSizes) );
		printf(":: %d lcp-intervals (tree representation = %.1lf MB)\n",numLcpIntervals, ((double)k)/((double)1000000) );
		printf(":: %.2lf%% with at least one alphabet letter missing\n",((double)numIncompleteLcpIntervals/(double)numLcpIntervals)*100.0);
		printf(":: %d shared top corners (avg count = %d , max count = %d)\n",numSharedTopCorners,(avgSharedTopCornersCount/numSharedTopCorners),maxSharedTopCornersCount);
		#endif
	}
	#ifdef DEBUGLCP
	if(verbose){ printf("> Testing Sampled Smaller Values "); fflush(stdout); }
	progressStep=(numLCPSamples/10);
	progressCounter=0;
	for(lcppos=0;lcppos<numLCPSamples;lcppos++){
		if(verbose){
			if(progressCounter==progressStep){ // print progress dots
				putchar('.');
				fflush(stdout);
				progressCounter=0;
			} else progressCounter++;
		}
		k=GetPrefixLinkFromLcpPos(lcppos);
		if(k!=fullPLPArray[lcppos]) break;
	}
	if(lcppos==numLCPSamples){ if(verbose) printf(" OK\n"); }
	else printf("\n> ERROR: sampledPLP[%d:%d]=%d =!= fullPLP[%d]=%d\n",(lcppos>>BLOCKSHIFT),(lcppos&BLOCKMASK),k,lcppos,fullPLPArray[lcppos]);
	if(verbose){ printf("> Testing Single Parent Intervals "); fflush(stdout); }
	progressStep=(numLCPSamples/10);
	progressCounter=0;
	topptr=0; // prevent compiler warnings
	bottomptr=0;
	for(lcppos=0;lcppos<numLCPSamples;lcppos++){
		if(verbose){
			if(progressCounter==progressStep){ // print progress dots
				putchar('.');
				fflush(stdout);
				progressCounter=0;
			} else progressCounter++;
		}
		bwtpos=GetBwtPosFromLcpPos(lcppos);
		topptr=bwtpos;
		bottomptr=bwtpos;
		k=GetTrueEnclosingLCPInterval(&topptr,&bottomptr);
		stopptr=bwtpos;
		sbottomptr=bwtpos;
		i=GetEnclosingLCPInterval(&stopptr,&sbottomptr);
		if(i!=k || stopptr!=topptr || sbottomptr!=bottomptr) break;
		lcp=GetLcpValueFromLcpPos(lcppos);
		topptr=bwtpos;
		bottomptr=bwtpos;
		k=GetTrueEnclosingLCPIntervalWithLcpValue(lcp,&topptr,&bottomptr);
		stopptr=bwtpos;
		sbottomptr=bwtpos;
		i=GetEnclosingLCPIntervalFromLCPPos(lcppos,&stopptr,&sbottomptr);
		if(i!=k || stopptr!=topptr || sbottomptr!=bottomptr) break;
	}
	if(lcppos==numLCPSamples){
		if(verbose) printf(" OK\n");
	} else printf("\n> ERROR: sampledNPSV[%u]=(%d){%u,%u} =!= fullNPSV[%u]=(%d){%u,%u}\n",bwtpos,i,stopptr,sbottomptr,bwtpos,k,topptr,bottomptr);
	/*
	if(verbose){ printf("> Testing Full Parent Intervals "); fflush(stdout); }
	progressStep=(bwtLength/10);
	progressCounter=0;
	topptr=0;
	bottomptr=0;
	for(bwtpos=0;bwtpos<(unsigned int)bwtLength;bwtpos++){
		if(verbose){
			if(progressCounter==progressStep){ // print progress dots
				putchar('.');
				fflush(stdout);
				progressCounter=0;
			} else progressCounter++;
		}
		stopptr=bwtpos;
		sbottomptr=bwtpos;
		topptr=bwtpos;
		bottomptr=bwtpos;
		i=INT_MAX;
		while(i>0){
			i=GetEnclosingLCPInterval(&stopptr,&sbottomptr);
			k=GetTrueEnclosingLCPInterval(&topptr,&bottomptr);
			if(i!=k || stopptr!=topptr || sbottomptr!=bottomptr) break;
		}
		if(i!=k || stopptr!=topptr || sbottomptr!=bottomptr) break;
	}
	if(bwtpos==(unsigned int)bwtLength){
		if(verbose) printf(" OK\n");
	} else printf("\n> ERROR: sampledNPSV[%u]=(%d){%u,%u} =!= fullNPSV[%u]=(%d){%u,%u}\n",bwtpos,i,stopptr,sbottomptr,bwtpos,k,topptr,bottomptr);
	*/
	#endif
	return numLCPSamples;
}
