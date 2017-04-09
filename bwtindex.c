/* ================================================================= *
 *  bwtindex.c : FM-Index                                            *
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
#include "bwtindex.h"
#include "packednumbers.h"

#define BUILD_LCP 1
//#define UNBOUNDED_LCP 1 // if the lcp values are unbounded (int) or truncated to 255 (unsigned char)
//#define DEBUG_INDEX 1
//#define FILL_INDEX 1

#ifdef DEBUG_INDEX
#include <sys/timeb.h>
#endif

//#define FILEHEADER "IDX0"

// TODO: implement BWT with an wavelet tree
typedef struct _IndexBlock { // 9*32 bits / 32 char per block = 9 bits per char
	unsigned int bwtBits[3]; // 3 bits (x32) for the letters in this order: $NACGT (000 to 101)
	unsigned int letterJumpsSample[5]; // cumulative counts for NACGT (not $) up to but *not* including this block
	unsigned int textPositionSample; // position in the text of the first suffix in this block
} IndexBlock;

#define ALPHABETSIZE 6
//static const char letterChars[ALPHABETSIZE] = { '$' , 'N' , 'A' , 'C' , 'G' , 'T' }; // get letter char from letter id
#define LETTERCHARS "$NACGT"
//static const unsigned int sampleIntervalSize = 32; // = (1<<sampleIntervalShift) = 32
#define SAMPLEINTERVALSIZE 32
//static const unsigned int sampleIntervalShift = 5; // sample interval of 32 positions (2^5=32)
#define SAMPLEINTERVALSHIFT 5
//static const unsigned int sampleIntervalMask = 0x0000001F; // = ((1<<sampleIntervalShift)-1) = (32-1)
#define SAMPLEINTERVALMASK 0x0000001F
//static const unsigned int firstLetterMask = 0x00000001; // lowest bit
#define FIRSTLETTERMASK 0x00000001

// Masks to select only the bit at the offset = (1UL<<offset)
static const unsigned int offsetMasks[32] = {
	0x00000001, // 1st bit
	0x00000002, // 2nd bit
	0x00000004, // 3rd bit
	0x00000008, // 4th bit
	0x00000010, // 5th bit
	0x00000020, // 6th bit
	0x00000040, // 7th bit
	0x00000080, // 8th bit
	0x00000100, // 9th bit
	0x00000200, // 10th bit
	0x00000400, // 11th bit
	0x00000800, // 12th bit
	0x00001000, // 13th bit
	0x00002000, // 14th bit
	0x00004000, // 15th bit
	0x00008000, // 16th bit
	0x00010000, // 17th bit
	0x00020000, // 18th bit
	0x00040000, // 19th bit
	0x00080000, // 20th bit
	0x00100000, // 21st bit
	0x00200000, // 22nd bit
	0x00400000, // 23rd bit
	0x00800000, // 24th bit
	0x01000000, // 25th bit
	0x02000000, // 26th bit
	0x04000000, // 27th bit
	0x08000000, // 28th bit
	0x10000000, // 29th bit
	0x20000000, // 30th bit
	0x40000000, // 31st bit
	0x80000000  // 32nd bit
};
// Masks to select all the bits before (but not at) the offset = ((1UL<<offset)-1UL)
static const unsigned int searchOffsetMasks[33] = {
	0x00000000, // lower 0 letters
	0x00000001, // lower 1 letters
	0x00000003, // lower 2 letters
	0x00000007, // lower 3 letters
	0x0000000F, // lower 4 letters
	0x0000001F, // lower 5 letters
	0x0000003F, // lower 6 letters
	0x0000007F, // lower 7 letters
	0x000000FF, // lower 8 letters
	0x000001FF, // lower 9 letters
	0x000003FF, // lower 10 letters
	0x000007FF, // lower 11 letters
	0x00000FFF, // lower 12 letters
	0x00001FFF, // lower 13 letters
	0x00003FFF, // lower 14 letters
	0x00007FFF, // lower 15 letters
	0x0000FFFF, // lower 16 letters
	0x0001FFFF, // lower 17 letters
	0x0003FFFF, // lower 18 letters
	0x0007FFFF, // lower 19 letters
	0x000FFFFF, // lower 20 letters
	0x001FFFFF, // lower 21 letters
	0x003FFFFF, // lower 22 letters
	0x007FFFFF, // lower 23 letters
	0x00FFFFFF, // lower 24 letters
	0x01FFFFFF, // lower 25 letters
	0x03FFFFFF, // lower 26 letters
	0x07FFFFFF, // lower 27 letters
	0x0FFFFFFF, // lower 28 letters
	0x1FFFFFFF, // lower 29 letters
	0x3FFFFFFF, // lower 30 letters
	0x7FFFFFFF, // lower 31 letters
	0xFFFFFFFF  // lower 32 letters
};
static const unsigned int inverseLetterBitMasks[ALPHABETSIZE][3] = {
	{
		0xFFFFFFFF, // 1st bit mask for '$' (000): ~0...0 = 1...1
		0xFFFFFFFF, // 2nd bit mask for '$' (000): ~0...0 = 1...1
		0xFFFFFFFF, // 3rd bit mask for '$' (000): ~0...0 = 1...1
	},{
		0x00000000, // 1st bit mask for 'N' (001): ~1...1 = 0...0
		0xFFFFFFFF, // 2nd bit mask for 'N' (001): ~0...0 = 1...1
		0xFFFFFFFF  // 3rd bit mask for 'N' (001): ~0...0 = 1...1
	},{
		0xFFFFFFFF, // 1st bit mask for 'A' (010): ~0...0 = 1...1
		0x00000000, // 2nd bit mask for 'A' (010): ~1...1 = 0...0
		0xFFFFFFFF  // 3rd bit mask for 'A' (010): ~0...0 = 1...1
	},{
		0x00000000, // 1st bit mask for 'C' (011): ~1...1 = 0...0
		0x00000000, // 2nd bit mask for 'C' (011): ~1...1 = 0...0
		0xFFFFFFFF  // 3rd bit mask for 'C' (011): ~0...0 = 1...1
	},{
		0xFFFFFFFF, // 1st bit mask for 'G' (100): ~0...0 = 1...1
		0xFFFFFFFF, // 2nd bit mask for 'G' (100): ~0...0 = 1...1
		0x00000000  // 3rd bit mask for 'G' (100): ~1...1 = 0...0
	},{
		0x00000000, // 1st bit mask for 'T' (101): ~1...1 = 0...0
		0xFFFFFFFF, // 2nd bit mask for 'T' (101): ~0...0 = 1...1
		0x00000000  // 3rd bit mask for 'T' (101): ~1...1 = 0...0
	}
};

static IndexBlock *Index = NULL;
static unsigned int bwtSize = 0; // textSize plus counting with the terminator char too
static unsigned int numSamples = 0;
static char *text = NULL;
//static PackedNumberArray *packedText = NULL;
static PackedNumberArray *packedBwt = NULL;
static char *textFilename = NULL;
static unsigned char *letterIds = NULL;

#ifdef BUILD_LCP
	#ifdef UNBOUNDED_LCP
	// use array of ints for storing LCP values
	static int *LCPArray;
	#else
	// use array of chars for storing LCP values (truncate to 255)
	static unsigned char *LCPArray;
	#endif
#endif

#ifdef DEBUG_INDEX
// used to count number of backtracking steps when searching for a text position in function FMI_PositionInText()
static unsigned int numBackSteps;
#endif

// variables and arrays needed to support a text composed of multiple strings
static char **multiStringTexts;
static char *multiStringLastChar;
static unsigned int *multiStringFirstPos;
static unsigned int multiStringPosShift;
static unsigned int *multiStringIdInBlock;

void InitializeIndexArrays(){
	int i;
	letterIds=(unsigned char *)malloc(256*sizeof(unsigned char));
	for(i=0;i<256;i++) letterIds[i]=(unsigned char)1; // ACGT, N and $
	letterIds[(int)'\0']=(unsigned char)0;
	letterIds[(int)'$']=(unsigned char)0;
	letterIds[(int)'N']=(unsigned char)1;
	letterIds[(int)'A']=(unsigned char)2;
	letterIds[(int)'C']=(unsigned char)3;
	letterIds[(int)'G']=(unsigned char)4;
	letterIds[(int)'T']=(unsigned char)5;
	letterIds[(int)'n']=(unsigned char)1;
	letterIds[(int)'a']=(unsigned char)2;
	letterIds[(int)'c']=(unsigned char)3;
	letterIds[(int)'g']=(unsigned char)4;
	letterIds[(int)'t']=(unsigned char)5;
	/*
	unsigned int mask, twoMasks[2];
	offsetMasks = (unsigned int *)malloc(SAMPLEINTERVALSIZE*sizeof(unsigned int));
	searchOffsetMasks = (unsigned int *)malloc((SAMPLEINTERVALSIZE+1)*sizeof(unsigned int));
	searchOffsetMasks[0] = 0x00000000;
	mask = 0x00000001;
	for(i=1;i<=SAMPLEINTERVALSIZE;i++){
		offsetMasks[(i-1)] = mask; // from 0 to 31
		searchOffsetMasks[i] = ( searchOffsetMasks[(i-1)] | mask ); // from 1 to 32
		mask = ( mask << 1 );
	}
	inverseLetterBitMasks = (unsigned int **)malloc(3*sizeof(unsigned int *));
	for(i=0;i<3;i++) inverseLetterBitMasks[i] = (unsigned int *)malloc(ALPHABETSIZE*sizeof(unsigned int));
	twoMasks[0] = 0xFFFFFFFF;
	twoMasks[1] = 0x00000000;
	for(i=0;i<ALPHABETSIZE;i++){
		inverseLetterBitMasks[i][0] = twoMasks[ (i & 0x1) ]; // if 1st bit = 1 , mask = ~0xFFFFFFFF , else mask = ~0x00000000
		inverseLetterBitMasks[i][1] = twoMasks[ ((i & 0x2) >> 1)  ]; // 2nd bit
		inverseLetterBitMasks[i][2] = twoMasks[ ((i & 0x4) >> 2)  ]; // 3rd bit
	}
	*/
}

void FMI_FreeIndex(){
	if(textFilename!=NULL){
		free(textFilename);
		textFilename=NULL;
	}
	if(Index!=NULL){
		free(Index);
		Index=NULL;
	}
	if(letterIds!=NULL){
		free(letterIds);
		letterIds=NULL;
	}
	if(multiStringTexts!=NULL){
		free(multiStringLastChar);
		free(multiStringFirstPos);
		free(multiStringIdInBlock);
		multiStringLastChar=NULL;
		multiStringFirstPos=NULL;
		multiStringIdInBlock=NULL;
		multiStringTexts=NULL;
	}
	/*
	if(offsetMasks!=NULL){
		free(offsetMasks);
		offsetMasks=NULL;
	}
	if(searchOffsetMasks!=NULL){
		free(searchOffsetMasks);
		searchOffsetMasks=NULL;
	}
	if(inverseLetterBitMasks!=NULL){
		for(i=0;i<ALPHABETSIZE;i++) free(inverseLetterBitMasks[i]);
		free(inverseLetterBitMasks);
		inverseLetterBitMasks=NULL;
	}
	*/
}

unsigned int FMI_GetTextSize(){
	return (bwtSize-1); // the bwtSize variable counts the terminator char too
}

unsigned int FMI_GetBWTSize(){
	return bwtSize;
}

char *FMI_GetTextFilename(){
	return textFilename;
}

static __inline void SetCharAtBWTPos( unsigned int bwtpos , unsigned int charid ){
	unsigned int sample = ( bwtpos >> SAMPLEINTERVALSHIFT );
	IndexBlock *block = &(Index[sample]); // get sample block
	unsigned int offset = ( bwtpos & SAMPLEINTERVALMASK );
	unsigned int mask = ( ~ offsetMasks[offset] ); // get all bits except the one at the offset
	unsigned int *letterMasks = (unsigned int *)(inverseLetterBitMasks[charid]);
	(block->bwtBits[0]) &= mask; // reset bits
	(block->bwtBits[1]) &= mask;
	(block->bwtBits[2]) &= mask;
	mask = offsetMasks[offset]; // get only the bit at the offset
	(block->bwtBits[0]) |= ( (~letterMasks[0]) & mask ); // set bits
	(block->bwtBits[1]) |= ( (~letterMasks[1]) & mask );
	(block->bwtBits[2]) |= ( (~letterMasks[2]) & mask );
}

// TODO: check if creating masks on-the-fly is faster than fetching them from array
static __inline unsigned int GetCharIdAtBWTPos( unsigned int bwtpos ){
	unsigned int sample, offset, mask, charid;
	IndexBlock *block;
	sample = ( bwtpos >> SAMPLEINTERVALSHIFT );
	offset = ( bwtpos & SAMPLEINTERVALMASK );
	block = &(Index[sample]); // get sample block
	mask = offsetMasks[offset]; // get only the bit at the offset
	charid = ( ( (block->bwtBits[0]) >> offset ) & FIRSTLETTERMASK ); // get 1st bit
	charid |= ( ( ( (block->bwtBits[1]) >> offset ) & FIRSTLETTERMASK ) << 1 ); // get 2nd bit
	charid |= ( ( ( (block->bwtBits[2]) >> offset ) & FIRSTLETTERMASK ) << 2 ); // get 3rd bit
	return charid;
	/**/
	// TODO: remove this!
	mask = mask;
	/**/
}

__inline char FMI_GetCharAtBWTPos( unsigned int bwtpos ){
	IndexBlock *block;
	unsigned int offset, charid;
	offset = ( bwtpos & SAMPLEINTERVALMASK );
	block = &(Index[( bwtpos >> SAMPLEINTERVALSHIFT )]); // get sample block
	charid = ( ( (block->bwtBits[0]) >> offset ) & FIRSTLETTERMASK ); // get 1st bit
	charid |= ( ( ( (block->bwtBits[1]) >> offset ) & FIRSTLETTERMASK ) << 1 ); // get 2nd bit
	charid |= ( ( ( (block->bwtBits[2]) >> offset ) & FIRSTLETTERMASK ) << 2 ); // get 3rd bit
	return LETTERCHARS[charid]; // get letter char
}

__inline unsigned int BitsSetCount( unsigned int bitArray ){
	/*
	while( bitArray ){
		letterJump++; // if other equal letters exist, increase letter jump
		bitArray &= ( bitArray - 1U ); // clear lower occurrence bit
	}
	*/
	/*
	letterJump += perByteCounts[ ( bitArray & 0x000000FF ) ];
	letterJump += perByteCounts[ ( bitArray & 0x0000FF00 ) >> 8 ];
	letterJump += perByteCounts[ ( bitArray & 0x00FF0000 ) >> 16 ];
	letterJump += perByteCounts[ ( bitArray >> 24 ) ];
	*/
	//bitArray = ( bitArray & 0x55555555 ) + ( ( bitArray >> 1 ) & 0x55555555 ); // sum blocks of 1 bits ; 0x5 = 0101b ; 2 bits ( final max count = 2 -> 2 bits )
	bitArray = bitArray - ( ( bitArray >> 1 ) & 0x55555555 );
	bitArray = ( bitArray & 0x33333333 ) + ( ( bitArray >> 2 ) & 0x33333333 ); // sum blocks of 2 bits ; 0x3 = 0011b ; 4 bits ( final max count = 4 -> 3 bits )
	bitArray = ( ( bitArray + ( bitArray >> 4 ) ) & 0x0F0F0F0F ); // sum blocks of 4 bits ; 0x0F = 00001111b ; 8 bits ( final max count = 8 -> 4 bits )
	//bitArray = ( bitArray + ( bitArray >> 8 ) ); // sum blocks of 8 bits ; the final count will be at most 32, which is 6 bits, and since we now have blocks of 8 bits, we can add them without masking because there will not be any overflow
	//bitArray = ( ( bitArray + ( bitArray >> 16 ) ) & 0x0000003F ); // sum blocks of 16 bits and get the last 6 bits ( final max count = 32 -> 6 bits )
	bitArray = ( ( bitArray * 0x01010101 ) >> 24 );
	return bitArray;
}

// NOTE: if letterId is not at position bwtPos, it considers the jump of the previous occurence behind/above
// NOTE: it assumes we will never try to do a letter jump by the terminator symbol, since there are no jumps stored in the index for it
static __inline unsigned int FMI_LetterJump( unsigned int letterId , unsigned int bwtPos ){
	unsigned int offset, bitArray, letterJump, *letterMasks;
	IndexBlock *block;
	letterMasks = (unsigned int *)(inverseLetterBitMasks[letterId]);
	offset = ( bwtPos & SAMPLEINTERVALMASK );
	block = &(Index[( bwtPos >> SAMPLEINTERVALSHIFT )]);
	bitArray = searchOffsetMasks[(offset+1)]; // all bits bellow and at offset (+1 otherwise it would not include the bit at the offset)
	bitArray &= ( (block->bwtBits[0]) ^ letterMasks[0] ); // keep only positions with the same 1st bit
	bitArray &= ( (block->bwtBits[1]) ^ letterMasks[1] ); // keep only positions with the same 2nd bit
	bitArray &= ( (block->bwtBits[2]) ^ letterMasks[2] ); // keep only positions with the same 3rd bit
	letterJump = (block->letterJumpsSample[(letterId-1)]); // get last letter jump before this block (jumps for '$' are not stored, so it is (letterId-1))
	#if defined(__GNUC__) && defined(__SSE4_2__)
		return ( letterJump + __builtin_popcount( bitArray ) );
	#else
		return ( letterJump + BitsSetCount( bitArray ) );
	#endif
}

// NOTE: returns the size of the BWT interval if a match exists, and 0 otherwise
unsigned int FMI_FollowLetter( char c , unsigned int *topPointer , unsigned int *bottomPointer ){
	/*
	unsigned int charId;
	unsigned int originalTopPointer;
	charId = letterIds[(unsigned char)c];
	originalTopPointer = (*topPointer);
	(*topPointer) = FMI_LetterJump( charId , (*topPointer) );
	if( GetCharIdAtBWTPos( originalTopPointer ) != charId ) (*topPointer)++; // if the letter is not in the topPointer position, its next occurrence is after that
	(*bottomPointer) = FMI_LetterJump( charId , (*bottomPointer) );
	*/
	unsigned int letterId, offset, bitArray, *letterMasks;
	IndexBlock *block;
	letterId = letterIds[(unsigned char)c];
	letterMasks = (unsigned int *)(inverseLetterBitMasks[letterId]);
	offset = ( (*topPointer) & SAMPLEINTERVALMASK );
	block = &(Index[( (*topPointer) >> SAMPLEINTERVALSHIFT )]);
	bitArray = searchOffsetMasks[offset]; // exclusive search mask on top pointer (all bits only bellow offset)
	bitArray &= ( (block->bwtBits[0]) ^ letterMasks[0] );
	bitArray &= ( (block->bwtBits[1]) ^ letterMasks[1] );
	bitArray &= ( (block->bwtBits[2]) ^ letterMasks[2] );
	(*topPointer) = (block->letterJumpsSample[(letterId-1)]);
	#if defined(__GNUC__) && defined(__SSE4_2__)
		(*topPointer) += __builtin_popcount( bitArray );
	#else
		(*topPointer) += BitsSetCount( bitArray );
	#endif
	(*topPointer)++; // if the letter is not in the topPointer position, its next occurrence is after that: LF[top]=count(c,(top-1))+1
	offset = ( (*bottomPointer) & SAMPLEINTERVALMASK );
	block = &(Index[( (*bottomPointer) >> SAMPLEINTERVALSHIFT )]);
	bitArray = searchOffsetMasks[(offset+1)]; // inclusive search mask on bottom pointer (all bits bellow and at offset)
	bitArray &= ( (block->bwtBits[0]) ^ letterMasks[0] );
	bitArray &= ( (block->bwtBits[1]) ^ letterMasks[1] );
	bitArray &= ( (block->bwtBits[2]) ^ letterMasks[2] );
	(*bottomPointer) = (block->letterJumpsSample[(letterId-1)]);
	#if defined(__GNUC__) && defined(__SSE4_2__)
		(*bottomPointer) += __builtin_popcount( bitArray );
	#else
		(*bottomPointer) += BitsSetCount( bitArray );
	#endif
	if( (*topPointer) > (*bottomPointer) ) return 0;
	return ( (*bottomPointer) - (*topPointer) + 1 );
}

unsigned int FMI_PositionInText( unsigned int bwtpos ){
	unsigned int charid, addpos;
	addpos = 0;
	while( bwtpos & SAMPLEINTERVALMASK ){ // move backwards until we land on a position with a sample
		charid = GetCharIdAtBWTPos(bwtpos);
		if( charid == 0 ){ // check if this is the terminator char
			#ifdef DEBUG_INDEX
			numBackSteps = addpos;
			#endif
			return addpos;
		}
		bwtpos = FMI_LetterJump( charid , bwtpos ); // follow the left letter backwards
		addpos++; // one more position away from our original position
	}
	#ifdef DEBUG_INDEX
	numBackSteps = addpos;
	#endif
	return ( (Index[( bwtpos >> SAMPLEINTERVALSHIFT )].textPositionSample) + addpos );
}

// returns the new position in the BWT array after left jumping by the char at the given BWT position
unsigned int FMI_LeftJump( unsigned int bwtpos ){
	unsigned int charid;
	charid = GetCharIdAtBWTPos( bwtpos );
	if( charid == 0 ) return 0U; // terminator symbol jumps to the 0-th position of the BWT
	else return FMI_LetterJump( charid , bwtpos ); // follow the left letter backwards
}

void FMI_GetCharCountsAtBWTInterval(unsigned int topPtr, unsigned int bottomPtr, int *counts){
	unsigned int offset, charid;
	IndexBlock *block;
	counts[0] = 0; // reset counts for chars: A,C,G,T,N
	counts[1] = 0;
	counts[2] = 0;
	counts[3] = 0;
	counts[4] = 0;
	block = &(Index[( topPtr >> SAMPLEINTERVALSHIFT )]);
	offset = ( topPtr & SAMPLEINTERVALMASK );
	while( topPtr <= bottomPtr ){ // process all positions of interval
		if( offset == SAMPLEINTERVALSIZE ){ // go to next sample block if needed
			offset = 0UL;
			block++;
		}
		charid = ( ( (block->bwtBits[0]) >> offset ) & FIRSTLETTERMASK ); // get 1st bit
		charid |= ( ( ( (block->bwtBits[1]) >> offset ) & FIRSTLETTERMASK ) << 1 ); // get 2nd bit
		charid |= ( ( ( (block->bwtBits[2]) >> offset ) & FIRSTLETTERMASK ) << 2 ); // get 3rd bit
		if( charid >= 2 ) counts[(charid - 2)]++; // ACGT
		else counts[4]++; // '$' or 'N'
		offset++;
		topPtr++;
	}
}

void PrintUnsignedNumber( unsigned int number ){
	unsigned int num, denom, quot, rem;
	if( number < 1000 ){
		printf( "%u" , number );
		return;
	}
	denom = 1;
	num = number;
	while( num >= 1000 ){
		num /= 1000;
		denom *= 1000;
	}
	printf( "%u," , num );
	num = ( number - num * denom );
	denom /= 1000;
	while( denom > 1 ){
		quot = ( num / denom );
		rem = ( num - quot * denom );
		printf( "%03u," , quot );
		num = rem;
		denom /= 1000;
	}
	printf( "%03u" , num );
}

/*
void FMI_LoadIndex(char *indexfilename){
	FILE *indexfile;
	size_t readcount;
	fpos_t filepos;
	unsigned int i;
	long long int seqstart, seqend, seqsize;
	char c;
	char fileHeader[5] = FILEHEADER;
	printf("> Loading index from file <%s> ... ",indexfilename);
	fflush(stdout);
	indexfile = fopen(indexfilename,"rb");
	if( indexfile == NULL ){
		printf("\n> ERROR: Cannot open file\n");
		exit(0);
	}
	seqstart=(long long int)ftell(indexfile);
	fseek(indexfile,0L,SEEK_END);
	seqend=(long long int)ftell(indexfile);
	seqsize=(seqend-seqstart);
	rewind(indexfile);
	printf("(");PrintNumber(seqsize);printf(" bytes)");
	fflush(stdout);
	for(i=0;i<4;i++){ // check if header is "IDX0"
		fread( &c , sizeof(char) , (size_t)1 , indexfile );
		if( c != fileHeader[i] ) break;
	}
	if( i != 4 ){
		printf("\n> ERROR: Invalid index file\n");
		exit(0);
	}
	fgetpos(indexfile,&filepos);
	i=0;
	while( c!='\0' && c!=EOF ){ // get text filename size
		fread( &c , sizeof(char) , (size_t)1 , indexfile );
		i++;
	}
	textFilename=(char *)malloc((i)*sizeof(char));
	fsetpos(indexfile,&filepos);
	fread( textFilename , sizeof(char) , (size_t)i , indexfile ); // get text filename chars
	fread( &textSize , sizeof(unsigned int) , (size_t)1 , indexfile ); // counts the terminator char too, so it is actually the BWT size
	fread( &numSamples , sizeof(unsigned int) , (size_t)1 , indexfile );
	#ifdef DEBUG
	printf("\n  [textFilename=\"%s\";textSize=%u;numSamples=%u]",textFilename,textSize,numSamples);
	fflush(stdout);
	#endif
	i = ( ( ( textSize-1 ) >> SAMPLEINTERVALSHIFT ) + 1 );
	if( ((textSize-1) & SAMPLEINTERVALMASK) != 0 ) i++; // extra sample
	if( numSamples != i ){ // check if number of samples is correct based on text size
		printf("\n> ERROR: Invalid index data\n");
		exit(0);
	}
	Index = (IndexBlock *)malloc( numSamples * sizeof(IndexBlock) ); // allocate memory for all index blocks
	if( Index == NULL ){
		printf("\n> ERROR: Not enough memory\n");
		exit(0);
	}
	readcount = fread( Index , sizeof(IndexBlock) , (size_t)(numSamples) , indexfile ); // read all index blocks
	if( readcount != (size_t)numSamples ){
		printf("\n> ERROR: Incomplete index data\n");
		exit(0);
	}
	fclose(indexfile);
	InitializeLetterIdsArray(); // initialize arrays used by index functions
	lastAlignedBwtPos = ((numSamples-1) << SAMPLEINTERVALSHIFT); // set last BWT pos for pattern matching initialization
	printf(" OK\n");
	fflush(stdout);
}

// TODO: ??? use indexBlocksStartInFile, distanceToNextBlockBits
void FMI_SaveIndex(char *indexfilename){
	FILE *indexfile;
	size_t writecount;
	long long int filesize;
	int i;
	printf("> Saving index to file <%s> ... ",indexfilename);
	fflush(stdout);
	indexfile = fopen(indexfilename,"wb");
	if( indexfile == NULL ){
		printf("\n> ERROR: Cannot create file\n");
		exit(0);
	}
	i=0;
	while(textFilename[i]!='\0') i++; // get text filename size
	i++;
	fwrite( (FILEHEADER) , sizeof(char) , (size_t)4 , indexfile );
	fwrite( textFilename , sizeof(char) , (size_t)i , indexfile );
	fwrite( &textSize , sizeof(unsigned int) , (size_t)1 , indexfile );
	fwrite( &numSamples , sizeof(unsigned int) , (size_t)1 , indexfile );
	writecount = fwrite( Index , sizeof(IndexBlock) , (size_t)(numSamples) , indexfile );
	if( writecount != (size_t)(numSamples) ){
		printf("\n> ERROR: Cannot write file\n");
		exit(0);
	}
	filesize=(long long int)ftell(indexfile);
	fclose(indexfile);
	printf("(");PrintNumber(filesize);printf(" bytes) OK\n");
	fflush(stdout);
}
*/


// Function pointer to get char id at corresponding text position
unsigned int (*GetTextCharId)(unsigned int);

// TODO: add code/functions for circular text and for packed binary text
// NOTE: the last pos (pos=textSize) is '\0' which maps to the id of '$' through the letterIds lookup table
unsigned int GetTextCharIdFromPlainText(unsigned int pos){
	return (unsigned int)letterIds[(unsigned char)text[pos]];
}

// Creates a lookup table with entries corresponding to blocks of size floor(log2(smallest_sequence)) and containing the sequence id at the first pos of each block
unsigned int InitializeMultiStringArrays(char **texts, unsigned int *textSizes, unsigned int numTexts){
	unsigned int id, pos, blockSize, blockNum, globalStringSize;
	multiStringTexts = texts;
	multiStringLastChar = (char *)malloc(numTexts*sizeof(char)); // terminator chars for each string
	multiStringFirstPos = (unsigned int *)malloc((numTexts+1)*sizeof(unsigned int)); // start pos in global string, plus one fake position after last one
	blockSize = UINT_MAX; // size of shortest string
	pos = 0; // position in global string
	for (id = 0; id < numTexts; id++){
		multiStringFirstPos[id] = pos;
		multiStringLastChar[id] = 'N';
		if (textSizes[id] < blockSize) blockSize = textSizes[id];
		pos += (textSizes[id] + 1); // string size plus extra terminator char
	}
	multiStringLastChar[(numTexts-1)] = '$'; // terminator char for global string
	multiStringFirstPos[numTexts] = pos; // fake next to last string
	globalStringSize = pos;
	multiStringPosShift = 1;
	while ((1UL << (multiStringPosShift + 1)) < blockSize) multiStringPosShift++; // get highest power of two lower or equal to the minimum size
	blockSize = (1UL << multiStringPosShift); // size of each block
	blockNum = (((globalStringSize - 1) >> multiStringPosShift) + 1); // total number of blocks
	multiStringIdInBlock = (unsigned int *)malloc(blockNum*sizeof(unsigned int)); // string id at the beginning of each block
	id = 0;
	pos = 0;
	blockNum = 0;
	while (pos <= (globalStringSize - 1)){ // fill the string id in all blocks
		if (pos >= multiStringFirstPos[(id + 1)]) id++; // if the 1st pos of this block is at or after the 1st of the next string, set next string as current
		multiStringIdInBlock[blockNum] = id;
		pos += blockSize; // next block
		blockNum++;
	}
	return globalStringSize;
}

// TODO: check implementation with binary search tree of shared and distinct bits of all numbers belonging to the same/different sequences
unsigned int GetTextCharIdFromMultipleStrings(unsigned int pos){
	unsigned int id, lastPos;
	id = multiStringIdInBlock[(pos >> multiStringPosShift)]; // string id at beginning of block containing this pos
	lastPos = (multiStringFirstPos[id + 1] - 1); // last pos of this string
	if (pos >= lastPos){
		if (pos == lastPos) return (unsigned int)letterIds[(unsigned char)multiStringLastChar[id]]; // virtual last char of this string
		id++; // the pos is on the next string
	}
	pos -= multiStringFirstPos[id]; // get pos inside string
	return (unsigned int)letterIds[(unsigned char)(multiStringTexts[id][pos])];
}


void PrintBWT(unsigned int *letterStartPos){
	unsigned int i, n, p;
	printf("%u {", bwtSize);
	for (n = 1; n < ALPHABETSIZE; n++){
		printf(" %c [%02u-%02u] %c", LETTERCHARS[n], letterStartPos[n], (n == (ALPHABETSIZE - 1)) ? (bwtSize - 1) : (letterStartPos[(n + 1)] - 1), (n == (ALPHABETSIZE - 1)) ? '}' : ',');
	}
	printf(" 2^%u=%u %u %#.8X\n", SAMPLEINTERVALSHIFT, SAMPLEINTERVALSIZE, numSamples, SAMPLEINTERVALMASK);
	printf("[ i] (SA) {");
	for (n = 1; n < ALPHABETSIZE; n++){
		printf(" %c%c", LETTERCHARS[n], (n == (ALPHABETSIZE - 1)) ? '}' : ',');
	}
#ifdef BUILD_LCP
	if (LCPArray != NULL) printf(" LCP");
#endif
	printf(" BWT\n");
	for (i = 0; i < bwtSize; i++){ // position in BWT
		p = FMI_PositionInText(i);
		printf("[%02u]%c(%2u) {", i, (i & SAMPLEINTERVALMASK) ? ' ' : '*', p);
		for (n = 1; n < ALPHABETSIZE; n++){
			printf("%02u%c", FMI_LetterJump(n, i), (n == (ALPHABETSIZE - 1)) ? '}' : ',');
		}
#ifdef BUILD_LCP
		if (LCPArray != NULL) printf(" %3d", (int)LCPArray[i]);
#endif
		printf(" %c ", FMI_GetCharAtBWTPos(i));
		n = 0;
		while ((n < (ALPHABETSIZE - 1)) && (i >= letterStartPos[(n + 1)])) n++;
		printf(" %c ", LETTERCHARS[n]);
		if (p != (bwtSize - 1)) p++; // char at the right of the BWT char
		else p = 0;
		//if( text != NULL ) printf("%s", (char *)(text+p) );
		while (p != bwtSize){
			n = GetTextCharId(p);
			printf(" %c ", LETTERCHARS[n]);
			p++;
		}
		printf("\n");
	}
}


typedef struct _LMSPos {
	unsigned int pos;
	#ifdef BUILD_LCP
	int lcp;
	#endif
	int next;
} LMSPos;

static LMSPos *LMSArray;
static int numLMS;

char GetCharType(unsigned int pos){
	char prevType, currentType;
	if(pos==bwtSize) return 'S'; // last position
	if(pos==0 || GetTextCharId(pos-1)<GetTextCharId(pos)) prevType='s';
	else if(GetTextCharId(pos-1)>GetTextCharId(pos)) prevType='l';
	else prevType='?';
	while(pos!=(bwtSize-1) && GetTextCharId(pos)==GetTextCharId(pos+1)) pos++;
	if(pos==(bwtSize-1) || GetTextCharId(pos)>GetTextCharId(pos+1)) currentType='l';
	else currentType='s';
	if(prevType=='?') prevType=currentType;
	if(prevType!=currentType) currentType-=32; // S*-type or L*-type
	return currentType;
}

// NOTE: fills the global LMSArray variable and the input charsCounts and charsFirstPos arrays
void GetLMSs( unsigned int *charsCounts , int *charsFirstPos , char verbose ){
	int arrayMaxSize, arrayGrowSize;
	int *charsLastPos;
	unsigned int i, j, n;
	char type;
	unsigned int progressCounter, progressStep;
	if(verbose){
		printf("> Collecting LMS positions ");
		fflush(stdout);
	}
	progressStep = (bwtSize/10);
	progressCounter = 0;
	charsLastPos = (int *)malloc(ALPHABETSIZE*sizeof(int));
	for( i = 0 ; i < ALPHABETSIZE ; i++ ){ // initialize chars buckets
		charsCounts[i] = 0;
		charsFirstPos[i] = (-1);
		charsLastPos[i] = (-1);
	}
	arrayGrowSize = (bwtSize/20); // how much to expand the LMS array each time we allocate more memory (5%)
	if( arrayGrowSize == 0 ) arrayGrowSize = 20;
	arrayMaxSize = 0;
	LMSArray = NULL;
	n = (bwtSize-1);
	j = GetTextCharId(n); // last char ('$')
	charsCounts[j]++;
	type = 'S'; // type of last char
	numLMS = 0; // current number of LMS (last char is S*-type but it will only be counted next)
	while( n != 0 ){ // process text in reverse from end to start
		if(verbose){
			progressCounter++;
			if(progressCounter==progressStep){ // print progress dots
				printf(".");
				fflush(stdout);
				progressCounter=0;
			}
		}
		n--; // previous position
		i = GetTextCharId(n); // current char
		charsCounts[i]++;
		if( i < j ) type = 'S'; // lower suffix, S-type char
		else if ( i > j ){ // higher suffix, L-type char
			if( type == 'S' ){ // if last char was S-type and this one is L-type, last one was S*-type
				if( numLMS == arrayMaxSize ){
					arrayMaxSize += arrayGrowSize;
					LMSArray = (LMSPos *)realloc(LMSArray,((size_t)arrayMaxSize)*sizeof(LMSPos));
					if( LMSArray == NULL ){
						printf("\n> ERROR: Failed to allocate %lld MB of memory\n",(((long long int)arrayMaxSize)*sizeof(LMSPos))/1000000LL);
						exit(-1);
					}
				}
				LMSArray[numLMS].pos = (n+1); // the previous position (to the right) is an S*-type char
				if( charsFirstPos[j] != (-1) ) LMSArray[ charsLastPos[j] ].next = numLMS; // add to linked list of this char id
				else charsFirstPos[j] = numLMS; // first char with this id
				charsLastPos[j] = numLMS; // current last char with this id
				numLMS++;
			}
			type = 'L';
		} // else (i==j), so use last used type
		j = i; // current char will be previous char on next step
	}
	LMSArray = (LMSPos *)realloc(LMSArray,numLMS*sizeof(LMSPos)); // shrink array to fit exact number of items
	for( i = 0 ; i < ALPHABETSIZE ; i++ ){ // set last position for each char bucket
		if( charsLastPos[i] != (-1) ) LMSArray[ charsLastPos[i] ].next = (-1);
	}
	free(charsLastPos);
	if(verbose){
		printf(" (%d) OK\n",numLMS);
		fflush(stdout);
	}
}

typedef struct _SortDepthState {
	unsigned int depth;
	int numChars;
	int charsFirstPos[ALPHABETSIZE];
	#ifdef DEBUG_INDEX
	int charsIds[ALPHABETSIZE];
	int charsCounts[ALPHABETSIZE];
	#endif
} SortDepthState;

void SortLMSs( int *charsBuckets , char verbose ){
	SortDepthState *sortStates;
	unsigned int textPos, depth;
	int prevSortedPos, prevLowestDepth, sortedCharId;
	int arrayPos, charPos, statePos, maxStatePos;
	int charId, numCharsToSort, nextNumCharsToSort;
	int *firstPos, *lastPos;
	unsigned int progressCounter, progressStep;
	#ifdef DEBUG_INDEX
	int prevDepth;
	int *firstPosIds, *charsCounts;
	int charIdToProcess;
	unsigned int prevTextPos;
	#endif
	if(verbose){
		printf("> Sorting LMS suffixes ");
		fflush(stdout);
	}
	progressStep=(numLMS/10);
	progressCounter=0;
	maxStatePos = 32;
	sortStates = (SortDepthState *)malloc(maxStatePos*sizeof(SortDepthState));
	lastPos = (int *)malloc(ALPHABETSIZE*sizeof(int));
	statePos = 0;
	sortStates[0].depth = 0;
	prevSortedPos = (-1);
	prevLowestDepth = 0;
	arrayPos = 0; // start at the first position of the linked array
	/**/
	depth = sortStates[0].depth;
	firstPos = sortStates[0].charsFirstPos;
	#ifdef DEBUG_INDEX
	firstPosIds = sortStates[0].charsIds;
	charsCounts = sortStates[0].charsCounts;
	charIdToProcess = (-1);
	prevDepth = (-1);
	#endif
	for( charId = 0 ; charId < ALPHABETSIZE ; charId++ ){ // initialize buckets with the already filled LMSs and reset them to save the output
		firstPos[charId] = charsBuckets[charId];
		charsBuckets[charId] = (-1);
		#ifdef DEBUG_INDEX
		firstPosIds[charId] = charId;
		#endif
	}
	goto _skip_char_count; // consider the buckets at depth 0 already filled
	/**/
	while( arrayPos != (-1) ){
		depth = sortStates[statePos].depth;
		firstPos = sortStates[statePos].charsFirstPos;
		#ifdef DEBUG_INDEX
		firstPosIds = sortStates[statePos].charsIds;
		charsCounts = sortStates[statePos].charsCounts;
		#endif
		for( charId = 0 ; charId < ALPHABETSIZE ; charId++ ){ // reset first/last pointers for each char
			firstPos[charId] = (-1);
			#ifdef DEBUG_INDEX
			firstPosIds[charId] = (-1);
			charsCounts[charId] = 0;
			#endif
			lastPos[charId] = (-1);
		}
		while( arrayPos != (-1) ){ // group all positions in this block by their first char (at the current depth)
			textPos = ( LMSArray[arrayPos].pos + depth );
			charId = GetTextCharId(textPos);
			#ifdef DEBUG_INDEX
			charsCounts[charId]++;
			#endif
			if( firstPos[charId] != (-1) ) LMSArray[ lastPos[charId] ].next = arrayPos; // add to already started char bucket
			else {
				firstPos[charId] = arrayPos; // first position in this char bucket
				#ifdef DEBUG_INDEX
				firstPosIds[charId] = charId;
				#endif
			}
			lastPos[charId] = arrayPos; // set current position as the current last one of the char bucket
			arrayPos = LMSArray[arrayPos].next; // next position
			LMSArray[ lastPos[charId] ].next = (-1); // set end of char bucket
		}
		_skip_char_count:
		numCharsToSort = ALPHABETSIZE; // number of char buckets to check on the next loop (if we are at a new depth, check all of them)
		while(1){ // loop for backward steps (we do not need the part above because we have already stored the char buckets at lower depths)
			arrayPos = (-1);
			nextNumCharsToSort = 0;
			#ifdef DEBUG_INDEX
			charIdToProcess = (-1);
			#endif
			for( charId = 0 ; charId < numCharsToSort ; charId++ ){ // check count of each char bucket
				charPos = firstPos[charId];
				if( charPos == (-1) ) continue; // char count = 0
				if( LMSArray[charPos].next == (-1) ){ // char count = 1 , which means this single position is sorted
					if( nextNumCharsToSort == 0 ){ // if there is no non-single bucket to sort before, add to sorted list
						if(verbose){
							progressCounter++;
							if(progressCounter==progressStep){ // print progress dots
								printf(".");
								fflush(stdout);
								progressCounter=0;
							}
						}
						if( prevLowestDepth == 0 ){ // if this is the first sorted position after a passage through depth 0, then set the top-level bucket
							sortedCharId = GetTextCharId( LMSArray[charPos].pos ); // get bucket from first char
							charsBuckets[sortedCharId] = charPos;
							if( prevSortedPos != (-1) ){
								LMSArray[prevSortedPos].next = (-1); // end previous top-level bucket
								prevSortedPos = (-1);
							}
						}
						#ifdef BUILD_LCP
						LMSArray[charPos].lcp = prevLowestDepth; // set LCP for this LMS suffix
						#endif
						prevLowestDepth = depth; // update previously seen minimum depth between two consecutive sorted positions
						if( prevSortedPos != (-1) ) LMSArray[prevSortedPos].next = charPos; // extend sorted list
						prevSortedPos = charPos; // set current last position of sorted list
						firstPos[charId] = (-1); // remove from list
					} else { //  if it's a single position but there's previous non-single buckets, we need to sort those before setting this one
						firstPos[ (nextNumCharsToSort - 1) ] = charPos; // move the unprocessed buckets to the beginning of the array so they will be processed on the next time we come back to this depth
						#ifdef DEBUG_INDEX
						firstPosIds[ (nextNumCharsToSort - 1) ] = firstPosIds[charId];
						#endif
						nextNumCharsToSort++;
					}
					continue;
				} // else, char count >= 2, which means the bucket needs further sorting
				if( nextNumCharsToSort == 0 ){ // if this is the first non-singular bucket, set it to be processed next
					arrayPos = charPos; // set next bucket to sort
					firstPos[charId] = (-1); // set its index as free as it will be already taken care of
					#ifdef DEBUG_INDEX
					charIdToProcess = firstPosIds[charId];
					firstPosIds[charId] = (-1);
					#endif
				} else { // else, there is at least one more non-singular bucket to sort
					firstPos[ (nextNumCharsToSort - 1) ] = charPos; // move it to the beginning of the array to be checked when we later get back to this depth
					#ifdef DEBUG_INDEX
					firstPosIds[ (nextNumCharsToSort - 1) ] = firstPosIds[charId];
					#endif
				}
				nextNumCharsToSort++;
			} // end of loop to check the size of each bucket
			#ifdef DEBUG_INDEX
			for( charId = (nextNumCharsToSort-1) ; charId < ALPHABETSIZE ; charId++ ){
				firstPosIds[charId] = (-1); // reset char ids that do not occur here
			}
			#endif
			sortStates[statePos].numChars = (nextNumCharsToSort - 1); // update the number of unsorted buckets at this depth (the first one will already be finished when we get back to this depth)
			if( nextNumCharsToSort != 0 ){ // if we have buckets left to sort, sort the first one (pointed by the previsouly set arrayPos) at the next depth
				if( nextNumCharsToSort > 1 ){ // if we only have one bucket, we can save the information in the same state, otherwise we have to create a new one
					statePos++;
					if( statePos == maxStatePos ){ // realloc array of states if needed
						maxStatePos += 32;
						sortStates = (SortDepthState *)realloc(sortStates,maxStatePos*sizeof(SortDepthState));
					}
				}
				#ifdef DEBUG_INDEX
				prevDepth = depth;
				#endif
				depth++;
				sortStates[statePos].depth = depth; // save depth for the following sorting step
				break;
			} // else, no more buckets left to sort at this depth, so check previous lower depths
			if( statePos == 0 ){ // check if we reached the end of the array
				arrayPos = (-1);
				break;
			}
			statePos--; // previous state at lower depth
			depth = sortStates[statePos].depth; // restore variables of previous state
			firstPos = sortStates[statePos].charsFirstPos;
			#ifdef DEBUG_INDEX
			firstPosIds = sortStates[statePos].charsIds;
			charsCounts = sortStates[statePos].charsCounts;
			prevDepth = sortStates[(statePos+1)].depth;
			#endif
			numCharsToSort = sortStates[statePos].numChars;
			prevLowestDepth = depth; // update previously seen minimum depth between two sorted positions
		} // end of loop for backward steps
	} // end of loop for all positions of the linked array
	free(lastPos);
	free(sortStates);
	if(verbose){
		printf(" OK\n");
		fflush(stdout);
	}
	#ifdef DEBUG_INDEX
	if(verbose){
		printf("> Checking sort ");
		fflush(stdout);
	}
	statePos = 0;
	prevSortedPos = (-1);
	for( charPos = 0 ; charPos < ALPHABETSIZE ; charPos++ ){
		arrayPos = charsBuckets[charPos];
		while( arrayPos != (-1) ){
			if( prevSortedPos != (-1) ){
				if(verbose){
					progressCounter++;
					if(progressCounter==progressStep){ // print progress dots
						printf(".");
						fflush(stdout);
						progressCounter=0;
					}
				}
				depth = 0;
				while(1){
					prevTextPos = (LMSArray[prevSortedPos].pos + depth);
					textPos = (LMSArray[arrayPos].pos + depth);
					sortedCharId = GetTextCharId(prevTextPos);
					charId = GetTextCharId(textPos);
					if( sortedCharId != charId ) break;
					depth++;
				}
				if( sortedCharId > charId ){
					printf("\n> ERROR: LMS[%d]@text[%d]='%c' > LMS[%d]@text[%d]='%c'\n",statePos,prevTextPos,LETTERCHARS[sortedCharId],(statePos+1),textPos,LETTERCHARS[charId]);
					fflush(stdout);
					charPos = INT_MAX;
					break;
				}
				#ifdef BUILD_LCP
				if( LMSArray[arrayPos].lcp != (int)depth ){
					printf("\n> ERROR: LMS[%d].LCP=%d =!= depth=%d\n",statePos,(LMSArray[arrayPos].lcp),(int)depth);
					fflush(stdout);
					charPos = INT_MAX;
					break;
				}
				#endif
			}
			statePos++;
			prevSortedPos = arrayPos;
			arrayPos = LMSArray[arrayPos].next;
		}
	}
	if( charPos == ALPHABETSIZE ){
		if( statePos != numLMS ) printf("\n> ERROR: #LMS=%d =!= %d\n",numLMS,statePos);
		else printf(" OK\n");
		fflush(stdout);
	}
	#endif
}

void InducedSort( unsigned int *bucketSize , int *bucketStartPos , char verbose ){
	int firstId[ALPHABETSIZE], lastId[ALPHABETSIZE], topSId[ALPHABETSIZE], bottomLId[ALPHABETSIZE];
	int arrayPos, nextArrayPos, charId, leftCharId;
	unsigned int textPos;
	unsigned int bucketPointer[ALPHABETSIZE];
	char processingType, leftType;
	unsigned int progressCounter, progressStep;
	#ifdef BUILD_LCP
	int lcpCharId, lcpValue;
	unsigned int prevTextPos, currentTextPos, lastLSuffixTextPos[ALPHABETSIZE];
	int prevMinLcpValue[ALPHABETSIZE], prevLcpCharLMSPos[ALPHABETSIZE], prevMinLStarLcpValue[ALPHABETSIZE];
	int savedLSBorderLcps[ALPHABETSIZE][ALPHABETSIZE];
	#endif
	if(verbose){
		printf("> Induced Sorting suffixes ");
		fflush(stdout);
	}
	progressStep = (bwtSize/10);
	progressCounter = 0;
	/*
	bucketSize = (unsigned int *)malloc(ALPHABETSIZE*sizeof(unsigned int)); // set pointers to the beginning of the buckets
	for( charId = 0 ; charId < ALPHABETSIZE ; charId++ ) bucketSize[charId] = 0; // reset bucket size
	for( textPos = 0 ; textPos < (bwtSize-1) ; textPos++ ) bucketSize[ GetTextCharId(textPos) ]++; // count number of each alphabet letter in the text (size of each bucket)
	*/
	bucketPointer[0] = 0; // pointers to the location in the suffix array that will be filled up next (for each letter)
	for( charId = 1 ; charId < ALPHABETSIZE ; charId++ ) bucketPointer[charId] = ( bucketPointer[(charId-1)] + bucketSize[(charId-1)] ); // points to the first position in each bucket
	#ifdef BUILD_LCP
	#endif
	for( charId = 0 ; charId < ALPHABETSIZE ; charId++ ){ // initialize pointers that will be used
		firstId[charId] = (-1); // pointers to start and end of array of L-type strings (at the top of buckets)
		lastId[charId] = (-1);
		topSId[charId] = bucketStartPos[charId]; // set this to the start position of S*-type suffixes found after sorting
		bottomLId[charId] = (-1); // pointers to start and end of array of S-type strings (at the bottom of buckets)
		#ifdef BUILD_LCP
		prevMinLcpValue[charId] = (-1); // set to (-1) so the first position in each each bucket is set to 0
		lastLSuffixTextPos[charId] = UINT_MAX;
		#endif
	}
	for( charId = 0 ; charId < ALPHABETSIZE ; charId++ ){ // downwards sweep: induce sort L-type chars from S*-type chars
		processingType = 'L';
		arrayPos = firstId[charId]; // process L-type chars at the top of the bucket
		#ifdef BUILD_LCP
		for( lcpCharId = 0 ; lcpCharId < ALPHABETSIZE ; lcpCharId++ ){ // set lcps to 0 when going to a different char
			if( prevMinLcpValue[lcpCharId] != (-1) ) prevMinLcpValue[lcpCharId] = 0; // it actually should be (-1) when charId=0
			prevMinLStarLcpValue[lcpCharId] = INT_MAX;
		}
		#endif
		L_from_S:
		while( arrayPos != (-1) ){ // newly found L-types will be stored at the bottom of the L-type array (at the top of buckets) pointed by lastId
			textPos = LMSArray[arrayPos].pos;
			if( textPos != 0 ) textPos--; // left position in the text
			else textPos = (bwtSize-1); // terminator char
			leftCharId = GetTextCharId( textPos ); // left char
			if( processingType == 'L' ){ // set the BWT chars when we are processing L-type suffixes
				if(verbose){
					progressCounter++;
					if(progressCounter==progressStep){ // print progress dots
						printf(".");
						fflush(stdout);
						progressCounter=0;
					}
				}
				#ifdef FILL_INDEX
				SetCharAtBWTPos( bucketPointer[charId] , leftCharId ); // fill the BWT array
				if( (bucketPointer[charId] & SAMPLEINTERVALMASK) == 0 ){ // add (textPos+1) sample to BWT Index here if the BWT position is a multiple of the sampling interval
					Index[ (bucketPointer[charId] >> SAMPLEINTERVALSHIFT) ].textPositionSample = LMSArray[arrayPos].pos;
				}
				#else
				SetPackedNumber( packedBwt , bucketPointer[charId] , leftCharId ); // fill the BWT array
				#endif
				#ifdef BUILD_LCP
				lcpValue = LMSArray[arrayPos].lcp;
				#ifdef UNBOUNDED_LCP
				LCPArray[ bucketPointer[charId] ] = lcpValue; // set lcp of L-type chars, since its correct value has already been calculated and stored before
				#else
				LCPArray[ bucketPointer[charId] ] = (lcpValue<UCHAR_MAX)?((unsigned char)lcpValue):(UCHAR_MAX);
				#endif
				if( lcpValue < prevMinLStarLcpValue[charId] ) prevMinLStarLcpValue[charId] = lcpValue; // minimum in interval between L*-type suffixes
				lastLSuffixTextPos[charId] = LMSArray[arrayPos].pos; // save text position of last L-type suffix to compute LCP between that and first S*-type suffix
				#endif
				bucketPointer[charId]++; // fill from the top downwards
				if( leftCharId > charId ) leftType = 'L'; // left char type
				else if( leftCharId < charId ) leftType = 'S';
				else leftType = processingType;
			} else leftType = 'L'; // if we are processing S*-type chars, we already know that the char to the left is L-type
			#ifdef BUILD_LCP
			lcpValue = LMSArray[arrayPos].lcp;
			for( lcpCharId = 0 ; lcpCharId < ALPHABETSIZE ; lcpCharId++ ){ // if there's a smaller lcp value of any char between two occurrences of the same (left) char, always set it to the minimum
				if( lcpValue < prevMinLcpValue[lcpCharId] ) prevMinLcpValue[lcpCharId] = lcpValue; // minimum of both values
			}
			#endif
			if( leftType == 'L' ){ // add to bottom of L-type list at the top of the bucket
				LMSArray[arrayPos].pos = textPos; // go one position backwards (to the left)
				if( firstId[leftCharId] == (-1) ) firstId[leftCharId] = arrayPos; // set first position of L-array (at the top of the bucket)
				else LMSArray[ lastId[leftCharId] ].next = arrayPos; // or connect to the former last position (grow down towards the bottom)
				lastId[leftCharId] = arrayPos; // this is now the current last position
				nextArrayPos = LMSArray[arrayPos].next; // if we are processing the last entry of the bucket and it has the same letter (arrayPos==lastBottomL[leftCharId]), next it will be this last entry again (done in previous line), so save it so it won't get set to (-1) next
				LMSArray[arrayPos].next = (-1); // and has nothing ahead
				#ifdef BUILD_LCP
				LMSArray[arrayPos].lcp = (prevMinLcpValue[leftCharId] + 1); // update lcp of this LMS to be set later at the left jump destination position (but only on L-type, not on L*-type since those already have the correct value set)
				#endif
			} else { // leftType == 'S' , which means current char is L*-type
				nextArrayPos = LMSArray[arrayPos].next; // save next id
				LMSArray[arrayPos].next = bottomLId[charId]; // add to end/bottom of L*-type list at the top of the bucket
				bottomLId[charId] = arrayPos;
				#ifdef BUILD_LCP
				LMSArray[arrayPos].lcp = prevMinLStarLcpValue[charId]; // set minimum between L*-type suffixes because it will be used when inducing S-type suffixes next
				prevMinLStarLcpValue[charId] = INT_MAX;
				#endif
			}
			#ifdef BUILD_LCP
			prevMinLcpValue[leftCharId] = INT_MAX; // since we only want the minimum in the interval [prev_pos_above+1,pos] , set this to max since its lcp value will not be used for next char
			#endif
			arrayPos = nextArrayPos; // next position in array
		}
		if( processingType == 'L' ){
			processingType = 'S';
			arrayPos = topSId[charId]; // process S*-type chars at the bottom of the bucket
			#ifdef BUILD_LCP
			prevTextPos = lastLSuffixTextPos[charId];
			if( (arrayPos != (-1)) && (prevTextPos != UINT_MAX) ){ // compute LCP between last L-type suffix and first S*-type suffix
				textPos = LMSArray[arrayPos].pos;
				lcpValue = 0;
				while( GetTextCharId(prevTextPos) == GetTextCharId(textPos) ){
					lcpValue++;
					textPos++;
					prevTextPos++;
				}
				LMSArray[arrayPos].lcp = lcpValue; // set new lcp value on first S*-type suffix
			}
			for( lcpCharId = 0 ; lcpCharId < ALPHABETSIZE ; lcpCharId++ ){ // save mininum LCP values between the last L-type suffix with this char on the left and the last L-type suffix of all, to be used in the next S-from-L phase
				if( prevMinLcpValue[lcpCharId] != (-1) ) savedLSBorderLcps[charId][lcpCharId] = prevMinLcpValue[lcpCharId]; // copy value
				else savedLSBorderLcps[charId][lcpCharId] = 0; // does not occur in the L-type suffixes left chars
			}
			#endif
			goto L_from_S;
		}
	}
	bucketPointer[0] = bucketSize[0]; // set pointers to the end of the buckets
	for( charId = 1 ; charId < ALPHABETSIZE ; charId++ ) bucketPointer[charId] = ( bucketPointer[(charId-1)] + bucketSize[charId] ); // total number of positions at and before this bucket
	for( charId = 0 ; charId < ALPHABETSIZE ; charId++ ) bucketPointer[charId]--; // points to the last position in each bucket
	for( charId = 0 ; charId < ALPHABETSIZE ; charId++ ){ // reset pointers that will be used
		firstId[charId] = (-1);
		lastId[charId] = (-1);
		topSId[charId] = (-1);
		#ifdef BUILD_LCP
		prevMinLcpValue[charId] = INT_MAX; // initialize with a large value since we will be doing the minimum of the values
		prevLcpCharLMSPos[charId] = (-1); // no LMS positions bellow the last one
		#endif
	}
	for( charId = ALPHABETSIZE ; charId != 0 ; ){ // upwards sweep: induce sort S-type chars from L*-type chars
		charId--; // process letters from bottom to top
		processingType = 'S';
		arrayPos = firstId[charId]; // process S-type chars at the bottom of the bucket
		#ifdef BUILD_LCP
		for( lcpCharId = 0 ; lcpCharId < ALPHABETSIZE ; lcpCharId++ ){ // set lcps to 0 when going to a different char
			if( prevMinLcpValue[lcpCharId] != INT_MAX ) prevMinLcpValue[lcpCharId] = 0;
		}
		#endif
		S_from_L:
		while( arrayPos != (-1) ){
			textPos = LMSArray[arrayPos].pos;
			if( textPos != 0 ) textPos--;
			else textPos = (bwtSize-1); // last position in text (it is the terminal symbol '$', which is an S-type char)
			leftCharId = GetTextCharId( textPos );
			if( processingType == 'S' ){ // set the BWT chars when we are processing S-type suffixes
				if(verbose){
					progressCounter++;
					if(progressCounter==progressStep){ // print progress dots
						printf(".");
						fflush(stdout);
						progressCounter=0;
					}
				}
				#ifdef FILL_INDEX
				SetCharAtBWTPos( bucketPointer[charId] , leftCharId );
				if( (bucketPointer[charId] & SAMPLEINTERVALMASK) == 0 ){ // add (textPos+1) sample to BWT Index here if the BWT position is a multiple of the sampling interval
					Index[ (bucketPointer[charId] >> SAMPLEINTERVALSHIFT) ].textPositionSample = LMSArray[arrayPos].pos;
				}
				#else
				SetPackedNumber( packedBwt , bucketPointer[charId] , leftCharId ); // fill the BWT array
				#endif
				#ifdef BUILD_LCP
				if( (LMSArray[arrayPos].next == (-1)) && ((prevTextPos=lastLSuffixTextPos[charId]) != UINT_MAX) ){ // if this is the last S-suffix and if there are L-suffixes above, explicitely compute the LCP
					currentTextPos = LMSArray[arrayPos].pos; // compute LCP between last S-type suffix and first L-type suffix found before (not necessarily L*-type)
					lcpValue = 0;
					while( GetTextCharId(prevTextPos) == GetTextCharId(currentTextPos) ){
						lcpValue++;
						prevTextPos++;
						currentTextPos++;
					}
					LMSArray[arrayPos].lcp = lcpValue; // fix lcp value of this last S-type suffix
				}
				#ifdef UNBOUNDED_LCP
				LCPArray[ bucketPointer[charId] ] = LMSArray[arrayPos].lcp; // set lcp of S-type chars, since its correct value has already been calculated and stored before
				#else
				lcpValue = LMSArray[arrayPos].lcp;
				LCPArray[ bucketPointer[charId] ] = (lcpValue<UCHAR_MAX)?((unsigned char)lcpValue):(UCHAR_MAX);
				#endif
				#endif
				bucketPointer[charId]--; // fill from the bottom upwards
				if( leftCharId < charId ) leftType = 'S';
				else if( leftCharId > charId ) leftType = 'L';
				else leftType = processingType;
			} else leftType = 'S'; // if we are processing L*-type chars, we already know that the char to the left is S-type
			if( leftType == 'S' ){ // add to top of S-type list at the bottom of the bucket
				LMSArray[arrayPos].pos = textPos;
				if( firstId[leftCharId] == (-1) ) firstId[leftCharId] = arrayPos; // set first position of S-array (at the bottom of the bucket)
				else LMSArray[ lastId[leftCharId] ].next = arrayPos; // or connect to the former last position (grow up towards the top)
				lastId[leftCharId] = arrayPos;
				nextArrayPos = LMSArray[arrayPos].next;
				LMSArray[arrayPos].next = (-1);
				#ifdef BUILD_LCP
				if( prevLcpCharLMSPos[leftCharId] != (-1) ){ // update the lcp of the LMS with the last seen (bellow) occurrence of this same left char with the minimum in the interval until now (exclusive)
					LMSArray[ prevLcpCharLMSPos[leftCharId] ].lcp = (prevMinLcpValue[leftCharId] + 1);
					if( arrayPos == prevLcpCharLMSPos[leftCharId] ){ // if the LMS of the prev (bellow) left char is this one, we already set its value in the LCP array but it is outdated because this step should have been done before
						#ifdef UNBOUNDED_LCP
						LCPArray[ (bucketPointer[charId] + 1) ] = LMSArray[arrayPos].lcp; // fix the value
						#else
						lcpValue = LMSArray[arrayPos].lcp;
						LCPArray[ (bucketPointer[charId] + 1) ] = (lcpValue<UCHAR_MAX)?((unsigned char)lcpValue):(UCHAR_MAX);
						#endif
					}
				}
				#endif
			} else { // leftType == 'L' , which means current char is S*-type
				nextArrayPos = LMSArray[arrayPos].next; // save next id
				LMSArray[arrayPos].next = topSId[charId]; // add to beginning/top of S*-type list at the bottom of the bucket
				topSId[charId] = arrayPos;
			}
			#ifdef BUILD_LCP
			lcpValue = LMSArray[arrayPos].lcp; // this step is also required for L*-type because we already have the min lcps for consecutive L*-type suffixes but not the min lcp for consecutive left char suffixes
			for( lcpCharId = 0 ; lcpCharId < ALPHABETSIZE ; lcpCharId++ ){ // if there's a smaller lcp value of any char between two occurrences of the same char, always set it to the minimum
				if( lcpValue < prevMinLcpValue[lcpCharId] ) prevMinLcpValue[lcpCharId] = lcpValue; // minimum of both values
			}
			prevMinLcpValue[leftCharId] = lcpValue; // since we want the minimum in the interval [next_pos_above+1,pos] , set this to the current lcp to be used by the next char position
			prevLcpCharLMSPos[leftCharId] = arrayPos;
			#endif
			arrayPos = nextArrayPos; // next position in array
		}
		if( processingType == 'S' ){
			processingType = 'L';
			arrayPos = bottomLId[charId]; // process L*-type chars at the top of the bucket
			#ifdef BUILD_LCP
			if( (firstId[charId] != (-1)) && (lastLSuffixTextPos[charId] != UINT_MAX) ){ // if there are S-type suffixes and L-type suffixes above, calculate correct lcp values for the left char pairs that cross the S/L-type border
				for( lcpCharId = 0 ; lcpCharId < charId ; lcpCharId++ ){
					nextArrayPos = prevLcpCharLMSPos[lcpCharId];
					if( nextArrayPos == (-1) ) continue;
					lcpValue = savedLSBorderLcps[charId][lcpCharId]; // minimum lcp after/bellow last occurrence (exclusive) of this left char in L-type suffixes
					if( lcpValue < prevMinLcpValue[lcpCharId] ) prevMinLcpValue[lcpCharId] = lcpValue;
					LMSArray[ nextArrayPos ].lcp = (prevMinLcpValue[lcpCharId] + 1); // update lcp of last S-type suffix with this left char
					prevLcpCharLMSPos[lcpCharId] = (-1); // reset values so the L*-type suffixes can calculate their LCPs based on L*-type suffixes only
					prevMinLcpValue[lcpCharId] = INT_MAX;
				}
			}
			#endif
			goto S_from_L;
		}
		#ifdef BUILD_LCP
		if( (firstId[charId] != (-1)) || (bottomLId[charId] != (-1)) ){ // after processing each char, update the lcps of the last/topmost occorrences of each left char without any pair yet
			for( lcpCharId = 0 ; lcpCharId < charId ; lcpCharId++ ){
				nextArrayPos = prevLcpCharLMSPos[lcpCharId];
				if( nextArrayPos == (-1) ) continue;
				LMSArray[ nextArrayPos ].lcp = 1;
				prevLcpCharLMSPos[lcpCharId] = (-1);
				prevMinLcpValue[lcpCharId] = INT_MAX;
			}
		}
		#endif
	}
	#ifdef BUILD_LCP
	#ifdef UNBOUNDED_LCP
	LCPArray[0] = (-1);
	#else
	LCPArray[0] = UCHAR_MAX;
	#endif
	textPos = 0; // set LCP of all 1st positions of all chars to 0, since the LCP of the 1st occurrence of each char in the BWT was never updated
	for( charId = 0 ; charId < (ALPHABETSIZE-1) ; charId++ ){
		textPos += bucketSize[charId];
		if( textPos < bwtSize ) LCPArray[textPos] = 0;
	}
	#endif
	if(verbose){
		printf(" OK\n");
		fflush(stdout);
	}
}

// TODO: implement BWT as wavelet tree:
//  - select (k) or (n-k) as first[2]={n,k} ; second[2])={k,0} ; result=(first[bit@level]-second[bit@level])
//  - blocks of 64 chars, but SA interval of 32 ( (pos&31==0) , SA_samples[2] @ ((pos>>5)&1) )
//  - sort chars by number of occurrences (or just "$,N" the same, last 2 chars the most frequent ones)
//  - variable length bits per char; process chars bottom up; depth-first while num chars in level is > 2; set chars bit to 0/1 at level
//  - level_size=bwt_size; k=(alphabet_size-1); n=sorted_counts[k]; while(n<(level_size/2)) n+=sorted_counts[--k]; ...
void FMI_BuildIndex(char **inputTexts, unsigned int *inputTextSizes, unsigned int inputNumTexts, unsigned char **lcpArrayPointer, char verbose){
	unsigned int letterId, i, n;
	unsigned int textPos, samplePos;
	unsigned int *letterCounts, *letterStartPos;
	int *letterLMSStartPos;
	IndexBlock *block;
	unsigned int progressCounter, progressStep;
	#ifdef DEBUG_INDEX
	unsigned int bwtPos, prevLetterId, k;
	unsigned int numRuns, sizeRun, longestRun, *runSizesCount;
	struct timeb startTime, endTime;
	double elapsedTime;
	long long unsigned int totalSize;
	#endif
	/*
	if(verbose){
		printf("> Processing reference chars ... ");
		fflush(stdout);
	}
	//textArray = NewPackedNumberArray((size+1),6); // bit array that stores all the chars of the text (plus terminator) in packed bits form
	//if(textArray==NULL){
	//	printf("\n> ERROR: Not enough memory\n");
	//	exit(-1);
	//}
	InitializeLetterIdsArray(); // initialize letter ids array
	for(i=0;i<ALPHABETSIZE;i++) letterCounts[i]=0; // reset counts for each letter ($NACGT)
	textSize=0; // counts all letters plus one sequence terminator symbol
	for(i=0;i<inputTextSize;i++){
		c=text[i];
		letterId = (unsigned int)letterIds[(unsigned char)c]; // get id of this letter
		letterCounts[letterId]++; // increase count of this letter
		textSize++; // one more text char
		if(textSize==UINT_MAX) break; // prevent overflow
		//SetPackedNumber(textArray,i,letterId);
	}
	//SetPackedNumber(textArray,i,0); // set terminator char (id=0) at the end
	if(textSize==0){
		printf("\n> ERROR: No valid characters found\n");
		exit(-1);
	}
	if(textSize==UINT_MAX){
		printf("\n> ERROR: Maximum reference size exceeded\n");
		exit(-1);
	}
	k=letterCounts[1]; // number of invalid chars ('N's)
	if(verbose){
		if(k>0) { printf("(");PrintUnsignedNumber(k);printf(" N's) "); }
		printf("(");PrintUnsignedNumber(textSize);printf(" basepairs) OK\n"); // number of valid chars (not including sequence terminator)
		printf("> Building BWT ");
		fflush(stdout);
	}
	letterCounts[0]++; // terminator char
	textSize++; // count terminator char
	//BWTIS( textArray , textSize , 6 , 0 , verbose );
	//FreePackedNumberArray(textArray); // the packed text array is not needed anymore, only the BWT array
	if(verbose){
		printf(" OK\n");
		fflush(stdout);
	}
	*/
	if (inputNumTexts == 1){ // only one text
		text = inputTexts[0];
		bwtSize = (inputTextSizes[0] + 1); // count terminator char too
		GetTextCharId = GetTextCharIdFromPlainText; // set function pointer
		multiStringTexts = NULL;
	} else { // multiple texts
		text = NULL;
		bwtSize = InitializeMultiStringArrays(inputTexts, inputTextSizes, inputNumTexts);
		GetTextCharId = GetTextCharIdFromMultipleStrings;
	}
	InitializeIndexArrays(); // initialize letter ids array
	letterCounts = (unsigned int *)malloc(ALPHABETSIZE*sizeof(unsigned int));
	letterLMSStartPos = (int *)malloc(ALPHABETSIZE*sizeof(int));

	LMSArray = NULL;
	GetLMSs(letterCounts,letterLMSStartPos,verbose);
	SortLMSs(letterLMSStartPos,verbose);

	#ifdef DEBUG_INDEX
	if(verbose){
		printf(":: ");
		for(i=0;i<ALPHABETSIZE;i++){
			totalSize=(((long long unsigned)letterCounts[i]*100ULL)/(long long unsigned)bwtSize);
			printf("#'%c'=%d(%d%%) ",LETTERCHARS[i],letterCounts[i],(int)totalSize);
		}
		printf("\n");
		// total number of bits in wavelet tree: 3 bits/levels for "$,N,A,C" and 2 bits/levels for "G,T"
		totalSize = 3ULL*(long long unsigned)(letterCounts[0]+letterCounts[1]+letterCounts[2]+letterCounts[3]);
		totalSize += 2ULL*(long long unsigned)(letterCounts[4]+letterCounts[5]);
		// number of bytes in wavelet tree + rank&select samples (interval=64) + SA samples (interval=32)
		printf(":: Size of index with wavelet-tree = %.2lf MB\n",(double)((totalSize/8)+(totalSize/64)*4+(bwtSize/32)*4)/1000000.0);
	}
	#endif
	
	#ifdef FILL_INDEX
	numSamples = ( ( ( bwtSize - 1 ) >> SAMPLEINTERVALSHIFT ) + 1 ); // blocks of 32 chars (the last used pos is (bwtSize-1))
	Index=(IndexBlock *)calloc(numSamples,sizeof(IndexBlock));
	if(Index==NULL){
		printf("> ERROR: Not enough memory to create index\n");
		exit(0);
	}
	packedBwt = NULL;
	#else
	Index = NULL;
	packedBwt = NewPackedNumberArray(bwtSize,ALPHABETSIZE); // bit array that stores all the chars of the BWT in packed bits form
	#endif

	#ifdef BUILD_LCP
	#ifdef UNBOUNDED_LCP
	LCPArray = (int *)malloc(bwtSize*sizeof(int));
	(*lcpArrayPointer) = NULL;
	#else
	LCPArray = (unsigned char *)malloc(bwtSize*sizeof(unsigned char));
	(*lcpArrayPointer) = LCPArray; // output LCP array as pointer in argument
	#endif
	#endif
	InducedSort(letterCounts,letterLMSStartPos,verbose);
	free(LMSArray);
	free(letterLMSStartPos);

	#ifndef FILL_INDEX // allocate index memory now if we did not fill it while building the BWT
	numSamples = ( ( ( bwtSize - 1 ) >> SAMPLEINTERVALSHIFT ) + 1 ); // blocks of 32 chars (if bwtSize was a multiple of 32 it needed +1 additional unused sample, because last used pos is (bwtSize-1))
	Index=(IndexBlock *)calloc(numSamples,sizeof(IndexBlock));
	if(Index==NULL){
		printf("> ERROR: Not enough memory to create index\n");
		exit(0);
	}
	#endif

	if(verbose){
		printf("> Collecting LF samples ");
		fflush(stdout);	
	}
	letterStartPos = (unsigned int *)malloc(ALPHABETSIZE*sizeof(unsigned int));
	letterStartPos[0]=0; // the terminator char is at the top (0-th) position of the BWT (but on the right)
	for(i=1;i<ALPHABETSIZE;i++) letterStartPos[i]=(letterStartPos[(i-1)]+letterCounts[(i-1)]); // where previous letter starts plus number of previous letter occurrences
	for(i=1;i<ALPHABETSIZE;i++) letterCounts[i]=(letterStartPos[i]-1); // initialize all letter jumps with the position before the start of the letter
	progressStep=(bwtSize/10);
	progressCounter=0;
	letterId=0; // just to fix compiler uninitialized warning
	textPos=0;
	samplePos = 0; // start in top position of the BWT and go down
	for( n = 0 ; n < bwtSize ; n++ ){
		if(verbose){
			progressCounter++;
			if(progressCounter==progressStep){ // print progress dots
				printf(".");
				fflush(stdout);
				progressCounter=0;
			}
		}
		#ifdef FILL_INDEX
		letterId = GetCharIdAtBWTPos(n);
		#else
		letterId = GetPackedNumber(packedBwt,n);
		#endif
		if( ( n & SAMPLEINTERVALMASK ) == 0 ){ // if we are over a sample, store here the current letter counts
			block = &(Index[samplePos]);
			#ifndef FILL_INDEX
			(block->bwtBits[0]) = 0U; // reset block
			(block->bwtBits[1]) = 0U;
			(block->bwtBits[2]) = 0U;
			#endif
			for(i=1;i<ALPHABETSIZE;i++) (block->letterJumpsSample[(i-1)]) = letterCounts[i]; // the i-th letter here is the (i-1)-th letter in the index
			samplePos++;
		}
		#ifndef FILL_INDEX
		SetCharAtBWTPos(n,letterId); // copy the current letter from the packed BWT to the BWT in the index
		#endif
		letterCounts[letterId]++;
	}
	#ifndef FILL_INDEX
	FreePackedNumberArray(packedBwt); // the packed BWT array is not needed anymore
	#endif
	if(verbose){
		printf(" OK\n");
		fflush(stdout);
	}
	#ifndef FILL_INDEX
	if(verbose){
		printf("> Collecting SA samples ");
		fflush(stdout);
	}
	textPos=(bwtSize-1); // start with the position of the terminator char in the text
	n=0; // start at the first/topmost BWT position
	progressStep=(bwtSize/10);
	progressCounter=0;
	while(1){
		if(verbose){
			progressCounter++;
			if(progressCounter==progressStep){ // print progress dots
				printf(".");
				fflush(stdout);
				progressCounter=0;
			}
		}
		if( ( n & SAMPLEINTERVALMASK ) == 0 ){ // if we are over a sample, store here the current position of the text
			samplePos = ( n >> SAMPLEINTERVALSHIFT );
			(Index[samplePos].textPositionSample) = textPos;
		}
		if(textPos==0) break;
		i = GetCharIdAtBWTPos(n); // get char at this BWT position (in the left)
		n = FMI_LetterJump(i,n); // follow the letter backwards to go to next position in the BWT
		textPos--;
	}
	if(verbose){
		printf(" OK\n");
		fflush(stdout);
	}
	#endif
	if(verbose){
		printf(":: FM-Index size = %u MB\n",( (unsigned int)((numSamples)*sizeof(IndexBlock))/1000000));
		#ifdef BUILD_LCP
		//printf(":: Short LCP Array size = %u MB\n",( (unsigned int)(bwtSize*sizeof(unsigned char))/1000000));
		#endif
		fflush(stdout);	
	}
	#ifdef DEBUG_INDEX
	if(bwtSize<100) PrintBWT(letterStartPos);
	if(verbose){
		printf("> Checking BWT ");
		#if ( defined(BUILD_LCP) && defined(UNBOUNDED_LCP) )
		printf("and LCP ");
		#endif
		fflush(stdout);
	}
	progressStep=(bwtSize/10);
	progressCounter=0;
	numRuns=1;
	sizeRun=1;
	longestRun=0;
	runSizesCount=(unsigned int *)malloc(1*sizeof(unsigned int));
	runSizesCount[0]=0;
	k = FMI_PositionInText(0); // position in the text of the top BWT position (should be equal to (bwtSize-1))
	for(bwtPos=1;bwtPos<bwtSize;bwtPos++){ // compare current position with position above
		if(verbose){
			progressCounter++;
			if(progressCounter==progressStep){ // print progress dots
				printf(".");
				fflush(stdout);
				progressCounter=0;
			}
		}
		i = FMI_PositionInText(bwtPos); // position in the text of the suffix in this row
		n = 0; // current suffix depth
		while( (prevLetterId=GetTextCharId((k+n))) == (letterId=GetTextCharId((i+n))) ) n++; // keep following suffix chars to the right while their letters are equal
		#if ( defined(BUILD_LCP) && defined(UNBOUNDED_LCP) )
		if( (int)LCPArray[bwtPos] != (int)n ){
			printf("\n> ERROR: LCP[%u]=%d =!= %d\n",bwtPos,(int)LCPArray[bwtPos],(int)n);
			printf("\t[%c] %c|",GetCharType(k),LETTERCHARS[GetCharIdAtBWTPos((bwtPos-1))]);
			for( textPos=k ; textPos<=(k+n) ; textPos++ ) putchar(LETTERCHARS[GetTextCharId(textPos)]);
			putchar('\n');
			printf("\t[%c] %c|",GetCharType(i),LETTERCHARS[GetCharIdAtBWTPos(bwtPos)]);
			for( textPos=i ; textPos<=(i+n) ; textPos++ ) putchar(LETTERCHARS[GetTextCharId(textPos)]);
			putchar('\n');
			fflush(stdout);
			getchar();
		}
		#endif
		if( prevLetterId > letterId ) break; // if the top letter is larger than the bottom letter, it is incorrectly sorted
		prevLetterId = GetCharIdAtBWTPos((bwtPos-1)); // get statistics about run lengths
		letterId = GetCharIdAtBWTPos(bwtPos);
		if( prevLetterId == letterId ) sizeRun++;
		else {
			if( sizeRun > longestRun ){
				runSizesCount = (unsigned int *)realloc(runSizesCount,(sizeRun+1)*sizeof(unsigned int));
				for(n=(longestRun+1);n<=sizeRun;n++) runSizesCount[n] = 0;
				longestRun = sizeRun;
			}
			runSizesCount[sizeRun]++;
			sizeRun = 1;
			numRuns++;
		}
		k = i; // current pos will be prev pos in next step
	}
	if(bwtPos!=bwtSize){
		printf(" FAILED (error at BWT position %u)\n",bwtPos);
		getchar();
		exit(-1);
	}
	if(verbose){
		printf(" OK\n");
		printf(":: Number of char runs = %u (%.2lf%%)\n",numRuns,(((double)numRuns*100.0)/(double)bwtSize));
		printf(":: Average run length = %.2lf (max = %u)\n",((double)bwtSize/(double)numRuns),longestRun);
		n=0;
		sizeRun=0;
		for(k=0;k<=longestRun;k++){
			if(runSizesCount[k]>n){
				n=runSizesCount[k];
				sizeRun=k;
			}
		}
		printf(":: Most frequent run length = %u (%.2lf%%)\n",sizeRun,(((double)n*100.0)/(double)numRuns));
		printf("> Checking LF samples ");
		fflush(stdout);
	}
	free(runSizesCount);
	for(i=0;i<ALPHABETSIZE;i++){ // initialize counters
		if(letterStartPos[i]!=0) letterCounts[i]=(letterStartPos[i]-1);
		else letterCounts[i]=0;
	}
	ftime(&startTime);
	progressCounter=0;
	for(n=0;n<bwtSize;n++){
		if(verbose){
			progressCounter++;
			if(progressCounter==progressStep){ // print progress dots
				printf(".");
				fflush(stdout);
				progressCounter=0;
			}
		}
		if( ( n & SAMPLEINTERVALMASK ) == 0 ){ // if there is a sample at this position, check counts
			samplePos = ( n >> SAMPLEINTERVALSHIFT );
			for(i=1;i<ALPHABETSIZE;i++) if( (Index[samplePos].letterJumpsSample[(i-1)]) != (letterCounts[i]) ) break;
			if(i!=ALPHABETSIZE) break;
		}
		i=GetCharIdAtBWTPos(n); // get letter and update count
		letterCounts[i]++;
		//if(i!=0) if( FMI_LetterJump(i,n) != letterCounts[i] ) break; // check if it is the terminator char because we cannot jump by it
		for(i=1;i<ALPHABETSIZE;i++) if( FMI_LetterJump(i,n) != letterCounts[i] ) break;
	}
	if(n!=bwtSize){
		printf(" FAILED (error at BWT position %u)\n",n);
		getchar();
		exit(-1);
	}
	ftime(&endTime);
	elapsedTime = ( ((endTime.time) + (endTime.millitm)/1000.0) - ((startTime.time) + (startTime.millitm)/1000.0) );
	if(verbose){
		printf(" OK\n");
		printf(":: Done in %.3lf seconds\n",elapsedTime);
		printf("> Checking SA samples ");
		fflush(stdout);
	}
	numRuns=0;
	longestRun=0;
	progressCounter=0;
	textPos=(bwtSize-1); // start at the terminator char
	n=0; // start at the first/topmost BWT position
	while(1){
		if(verbose){
			progressCounter++;
			if(progressCounter==progressStep){ // print progress dots
				printf(".");
				fflush(stdout);
				progressCounter=0;
			}
		}
		numBackSteps = 0;
		if( FMI_PositionInText(n) != textPos ) break;
		if( numBackSteps > longestRun ) longestRun = numBackSteps;
		numRuns += numBackSteps;
		if( ( n & SAMPLEINTERVALMASK ) == 0 ){ // if there is a sample at this BWT position, check text position
			samplePos = ( n >> SAMPLEINTERVALSHIFT );
			if( (Index[samplePos].textPositionSample) != textPos ) break;
		}
		if(textPos==0) break;
		textPos--;
		letterId=GetTextCharId(textPos);
		i=GetCharIdAtBWTPos(n);
		if(letterId!=i) break; // check if it is the same letter at the BWT and at the text
		if(i==0) break; // check if it is the terminator char because it should not be here
		n = FMI_LetterJump(i,n); // follow the letter backwards to go to next position in the BWT
	}
	if( textPos!=0 || FMI_PositionInText(n)!=0 || letterId!=i || i==0 || GetCharIdAtBWTPos(n)!=0 ){
		printf(" FAILED (error at text position %u)\n",textPos);
		getchar();
		exit(-1);
	}
	if(verbose){
		printf(" OK\n");
		printf(":: Average backtracking steps = %.2lf (max = %u)\n",((double)numRuns/(double)bwtSize),longestRun);
		fflush(stdout);
	}
	#endif
	free(letterCounts);
	free(letterStartPos);
}
