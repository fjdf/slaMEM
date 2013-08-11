#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "bwtindex.h"

#define FILEHEADER "IDX0"

#define ALPHABETSIZE 6
#define ALPHABET "$NACGT"

#define BUILD_LCP 1
//#define UNBOUNDED_LCP 1 // if the lcp values are unbounded (int) or truncated to 255 (unsigned char)
//#define DEBUG_INDEX 1
//#define FILL_INDEX 1

typedef struct _IndexBlock {
	unsigned int bwtLowBits;
	unsigned int bwtHighBits;
	unsigned int specialLettersMask;
	//unsigned int bwtBits[3];
	unsigned int letterJumpsSample[5];
	unsigned int textPositionSample;
} IndexBlock;

#ifdef BUILD_LCP
#ifdef UNBOUNDED_LCP
int *LCPArray;
#else
unsigned char *LCPArray;
#endif
#endif

#ifdef DEBUG_INDEX
unsigned int numBackSteps;
#endif

const unsigned int offsetMasks[32] = { // = (1UL<<offset)
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
const unsigned int inverseOffsetMasks[32] = { // = ~(1UL<<offset)
	0xFFFFFFFE, // 1st bit
	0xFFFFFFFD, // 2nd bit
	0xFFFFFFFB, // 3rd bit
	0xFFFFFFF7, // 4th bit
	0xFFFFFFEF, // 5th bit
	0xFFFFFFDF, // 6th bit
	0xFFFFFFBF, // 7th bit
	0xFFFFFF7F, // 8th bit
	0xFFFFFEFF, // 9th bit
	0xFFFFFDFF, // 10th bit
	0xFFFFFBFF, // 11th bit
	0xFFFFF7FF, // 12th bit
	0xFFFFEFFF, // 13th bit
	0xFFFFDFFF, // 14th bit
	0xFFFFBFFF, // 15th bit
	0xFFFF7FFF, // 16th bit
	0xFFFEFFFF, // 17th bit
	0xFFFDFFFF, // 18th bit
	0xFFFBFFFF, // 19th bit
	0xFFF7FFFF, // 20th bit
	0xFFEFFFFF, // 21st bit
	0xFFDFFFFF, // 22nd bit
	0xFFBFFFFF, // 23rd bit
	0xFF7FFFFF, // 24th bit
	0xFEFFFFFF, // 25th bit
	0xFDFFFFFF, // 26th bit
	0xFBFFFFFF, // 27th bit
	0xF7FFFFFF, // 28th bit
	0xEFFFFFFF, // 29th bit
	0xDFFFFFFF, // 30th bit
	0xBFFFFFFF, // 31st bit
	0x7FFFFFFF  // 32nd bit
};
const unsigned int firstLettersMasks[33] = { // all bits before the offset = ((1UL<<offset)-1UL)
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
const unsigned int lastLettersMasks[33] = { // all bits at or after offset = ~((1UL<<offset)-1UL)
	0xFFFFFFFF, // higher 32 letters
	0xFFFFFFFE, // higher 31 letters
	0xFFFFFFFC, // higher 30 letters
	0xFFFFFFF8, // higher 29 letters
	0xFFFFFFF0, // higher 28 letters
	0xFFFFFFE0, // higher 27 letters
	0xFFFFFFC0, // higher 26 letters
	0xFFFFFF80, // higher 25 letters
	0xFFFFFF00, // higher 24 letters
	0xFFFFFE00, // higher 23 letters
	0xFFFFFC00, // higher 22 letters
	0xFFFFF800, // higher 21 letters
	0xFFFFF000, // higher 20 letters
	0xFFFFE000, // higher 19 letters
	0xFFFFC000, // higher 18 letters
	0xFFFF8000, // higher 17 letters
	0xFFFF0000, // higher 16 letters
	0xFFFE0000, // higher 15 letters
	0xFFFC0000, // higher 14 letters
	0xFFF80000, // higher 13 letters
	0xFFF00000, // higher 12 letters
	0xFFE00000, // higher 11 letters
	0xFFC00000, // higher 10 letters
	0xFF800000, // higher 9 letters
	0xFF000000, // higher 8 letters
	0xFE000000, // higher 7 letters
	0xFC000000, // higher 6 letters
	0xF8000000, // higher 5 letters
	0xF0000000, // higher 4 letters
	0xE0000000, // higher 3 letters
	0xC0000000, // higher 2 letters
	0x80000000, // higher 1 letters
	0x00000000  // higher 0 letters
};
const unsigned int searchOffsetMasks[32] = { // all bits before or at the offset, except 1st one = (((1UL<<(offset+1))-1UL)&(~1UL))
	0x00000000, // lower 1 letters
	0x00000002, // lower 2 letters
	0x00000006, // lower 3 letters
	0x0000000E, // lower 4 letters
	0x0000001E, // lower 5 letters
	0x0000003E, // lower 6 letters
	0x0000007E, // lower 7 letters
	0x000000FE, // lower 8 letters
	0x000001FE, // lower 9 letters
	0x000003FE, // lower 10 letters
	0x000007FE, // lower 11 letters
	0x00000FFE, // lower 12 letters
	0x00001FFE, // lower 13 letters
	0x00003FFE, // lower 14 letters
	0x00007FFE, // lower 15 letters
	0x0000FFFE, // lower 16 letters
	0x0001FFFE, // lower 17 letters
	0x0003FFFE, // lower 18 letters
	0x0007FFFE, // lower 19 letters
	0x000FFFFE, // lower 20 letters
	0x001FFFFE, // lower 21 letters
	0x003FFFFE, // lower 22 letters
	0x007FFFFE, // lower 23 letters
	0x00FFFFFE, // lower 24 letters
	0x01FFFFFE, // lower 25 letters
	0x03FFFFFE, // lower 26 letters
	0x07FFFFFE, // lower 27 letters
	0x0FFFFFFE, // lower 28 letters
	0x1FFFFFFE, // lower 29 letters
	0x3FFFFFFE, // lower 30 letters
	0x7FFFFFFE, // lower 31 letters
	0xFFFFFFFE  // lower 32 letters
};
const unsigned int firstLetterMask = 0x00000001; // lowest bit
const unsigned int lastLetterMask  = 0x80000000; // highest bit
const unsigned int secondLetterMask = 0x00000002; // 2nd lowest bit
const unsigned int firstHalfLettersMask = 0x0000FFFF; // lowest 16 bits
const unsigned int lastHalfLettersMask  = 0xFFFF0000; // highest 16 bits
const unsigned int letterLowBitMasks[6] = {
	0x00000000, // unused (for '$' (00))
	0xFFFFFFFF, // unused (for 'N' (01))
	0x00000000, // 1st bit mask for 'A' (00): 0...0
	0xFFFFFFFF, // 1st bit mask for 'C' (01): 1...1
	0x00000000, // 1st bit mask for 'G' (10): 0...0
	0xFFFFFFFF  // 1st bit mask for 'T' (11): 1...1
};
const unsigned int letterHighBitMasks[6] = {
	0x00000000, // unused (for '$' (00))
	0x00000000, // unused (for 'N' (01))
	0x00000000, // 2nd bit mask for 'A' (00): 0...0
	0x00000000, // 2nd bit mask for 'C' (01): 0...0
	0xFFFFFFFF, // 2nd bit mask for 'G' (10): 1...1
	0xFFFFFFFF  // 2nd bit mask for 'T' (11): 1...1
};
const unsigned int inverseLetterLowBitMasks[6] = {
	0xFFFFFFFF, // unused (for '$' (00))
	0x00000000, // unused (for 'N' (01))
	0xFFFFFFFF, // 1st bit mask for 'A' (00): ~0...0 = 1...1
	0x00000000, // 1st bit mask for 'C' (01): ~1...1 = 0...0
	0xFFFFFFFF, // 1st bit mask for 'G' (10): ~0...0 = 1...1
	0x00000000  // 1st bit mask for 'T' (11): ~1...1 = 0...0
};
const unsigned int inverseLetterHighBitMasks[6] = {
	0xFFFFFFFF, // unused (for '$' (00))
	0xFFFFFFFF, // unused (for 'N' (01))
	0xFFFFFFFF, // 2nd bit mask for 'A' (00): ~0...0 = 1...1
	0xFFFFFFFF, // 2nd bit mask for 'C' (01): ~0...0 = 1...1
	0x00000000, // 2nd bit mask for 'G' (10): ~1...1 = 0...0
	0x00000000  // 2nd bit mask for 'T' (11): ~1...1 = 0...0
};
const char letterChars[6] = { '$' , 'N' , 'A' , 'C' , 'G' , 'T' }; // get letter char from letter id
const unsigned int sampleIntervalShift = 5; // sample interval of 32 positions (2^5=32)
const unsigned int sampleIntervalMask = 0x0000001F; // = ((1<<sampleIntervalShift)-1) = (32-1)
const unsigned int sampleIntervalSize = 32; // = (1<<sampleIntervalShift) = 32
const unsigned int sampleIntervalHalfSize = 16; // = (1<<(sampleIntervalShift-1)) = 16

// varibles needed to load/store index
struct _IndexBlock *Index = NULL;
unsigned int textSize = 0; // inside the index functions, textSize always counts with the terminator char
unsigned int numSamples = 0;
unsigned int lastBwtPos = 0; // position in the BWT of the 0-th entry of the last sample (multiple of 32 to speed up search of 1st char of pattern)
char *text = NULL;
struct _PackedNumberArray *packedText = NULL;
struct _PackedNumberArray *packedBwt = NULL;
char *textFilename = NULL;
unsigned char *letterIds = NULL;


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

void SetCharAtBWTPos(unsigned int bwtpos, unsigned int charid){
	unsigned int sample = ( bwtpos >> sampleIntervalShift );
	IndexBlock *block = &(Index[sample]); // get sample block
	unsigned int offset = ( bwtpos & sampleIntervalMask );
	unsigned int mask = inverseOffsetMasks[offset]; // get all bits except the one at the offset
	(block->specialLettersMask) &= mask; // reset special letter mark
	(block->bwtLowBits) &= mask; // reset low bit
	(block->bwtHighBits) &= mask; // reset high bit
	mask = offsetMasks[offset]; // get only the bit at the offset
	if( charid < 2 ) (block->specialLettersMask) |= mask; // set special letter mark if needed (for ids 0 and 1)
	(block->bwtLowBits) |= ( letterLowBitMasks[charid] & mask ); // set low bit
	(block->bwtHighBits) |= ( letterHighBitMasks[charid] & mask ); // set high bit
}

unsigned int GetCharIdAtBWTPos(unsigned int bwtpos){
	unsigned int sample, offset, mask, charid;
	IndexBlock *block;
	sample = ( bwtpos >> sampleIntervalShift );
	offset = ( bwtpos & sampleIntervalMask );
	block = &(Index[sample]); // get sample block
	mask = offsetMasks[offset]; // get only the bit at the offset
	charid = ( ( (block->bwtLowBits) >> offset ) & firstLetterMask ); // get low bit
	charid |= ( ( ( (block->bwtHighBits) >> offset ) & firstLetterMask ) << 1 ); // get high bit
	if( !( (block->specialLettersMask) & mask ) ) charid += 2; // check if special letter bit is set, because regular ids start at id=2
	return charid;
}

char FMI_GetCharAtBWTPos(unsigned int bwtpos){
	IndexBlock *block;
	unsigned int offset, charid;
	offset = ( bwtpos & sampleIntervalMask );
	block = &(Index[( bwtpos >> sampleIntervalShift )]); // get sample block
	charid = ( ( (block->bwtLowBits) >> offset ) & firstLetterMask ); // get low bit
	charid |= ( ( ( (block->bwtHighBits) >> offset ) & firstLetterMask ) << 1 ); // get high bit
	if( (block->specialLettersMask) & offsetMasks[offset] ) return letterChars[charid]; // check if special letter bit is set, because regular ids start at id=2
	return letterChars[(charid + 2)]; // get letter char
}

// TODO: check if when calling function counts of letter 'N' are detected at pos 0 or 4
void FMI_GetCharCountsAtBWTInterval(unsigned int topPtr, unsigned int bottomPtr, int *counts){
	unsigned int offset, charid;
	IndexBlock *block;
	counts[0] = 0; // reset counts for chars: A,C,G,T,N
	counts[1] = 0;
	counts[2] = 0;
	counts[3] = 0;
	counts[4] = 0;
	block = &(Index[( topPtr >> sampleIntervalShift )]);
	offset = ( topPtr & sampleIntervalMask );
	while( topPtr <= bottomPtr ){ // process all positions of interval
		if( offset == sampleIntervalSize ){ // go to next sample block if needed
			offset = 0UL;
			block++;
		}
		charid = ( ( (block->bwtLowBits) >> offset ) & firstLetterMask ); // get low bit
		charid |= ( ( ( (block->bwtHighBits) >> offset ) & firstLetterMask ) << 1 ); // get high bit
		if( (block->specialLettersMask) & offsetMasks[offset] ) counts[4]++; // check if special letter bit is set
		else counts[charid]++;
		offset++;
		topPtr++;
	}
}

/*
unsigned int GetRightCharIdAtBWTPos(unsigned int bwtpos){
	unsigned int charid;
	charid = 5;
	while( (charid != 0) && (bwtpos < letterStartPos[charid]) ) charid--;
	return charid;
}
*/

void InitializeLetterIdsArray(){
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
}

unsigned int FMI_GetTextSize(){
	return (textSize-1); // inside the index functions, the textSize variable counts the terminator char too
}

unsigned int FMI_GetBWTSize(){
	return lastBwtPos; // last valid filled position of the BWT (aligned to a multiple of the sample interval)
}

char *FMI_GetTextFilename(){
	return textFilename;
}

// TODO: find way to remove the decrement in letterJumpsSample[(letterId-1)]
// TODO: move code with the check for special letters to another function, because on normal search we never use jumps by $'s or N's (but on position search, we do by N's)
unsigned int FMI_LetterJump( unsigned int letterId , unsigned int bwtPos ){
	unsigned int offset, bitArray, letterJump;
	IndexBlock *block;
	block = &(Index[( bwtPos >> sampleIntervalShift )]);
	offset = ( bwtPos & sampleIntervalMask );
	bitArray = searchOffsetMasks[offset]; // all bits bellow offset, except 1st one
	if( letterId < 2 ) bitArray &= (block->specialLettersMask); // if we are looking for '$' or 'N', keep only special chars
	else bitArray &= ( ~ (block->specialLettersMask) ); // remove special chars
	bitArray &= ( (block->bwtLowBits) ^ inverseLetterLowBitMasks[letterId] ); // keep only ones with the same low bit
	bitArray &= ( (block->bwtHighBits) ^ inverseLetterHighBitMasks[letterId] ); // keep only ones with the same high bit
	letterJump = (block->letterJumpsSample[(letterId-1)]); // get top letter jump
	while( bitArray ){
		letterJump++; // if other equal letter exist, increase letter jump
		bitArray &= ( bitArray - 1U ); // clear lower occurrence bit
	}
	return letterJump;
}

// TODO: find better way to deal with original toppointer value
// TODO: move code of called functions inside and reuse variables for speed up
// NOTE: returns the size of the BWT interval if a match exists, and 0 otherwise
unsigned int FMI_FollowLetter( char c , unsigned int *topPointer , unsigned int *bottomPointer ){
	unsigned int charid, originaltoppointer;
	originaltoppointer = (*topPointer);
	charid = letterIds[(int)c];
	(*topPointer) = FMI_LetterJump( charid , originaltoppointer );
	if( GetCharIdAtBWTPos(originaltoppointer) != charid ) (*topPointer)++; // if the letter is not in the topPointer position (the initial, not the updated one), its next occurrence is after that
	(*bottomPointer) = FMI_LetterJump( charid , (*bottomPointer) );
	if( (*topPointer) > (*bottomPointer) ) return 0;
	return ( (*bottomPointer) - (*topPointer) + 1 );
	/*
	unsigned int bwtPos, letterId, offset, bitArray, letterJump;
	IndexBlock *block;
	letterId = letterIds[c];
	bwtPos = (*topPointer); // process top pointer
	block = &(Index[( bwtPos >> sampleIntervalShift )]); // get sample block
	bitArray = ( ~ (block->specialLettersMask) ); // remove special chars
	bitArray &= ( (block->bwtLowBits) ^ inverseLetterLowBitMasks[letterId] ); // keep only ones with the same low bit
	bitArray &= ( (block->bwtHighBits) ^ inverseLetterHighBitMasks[letterId] ); // keep only ones with the same high bit
	letterJump = (block->letterJumpsSample[letterId]); // get top letter jump
	offset = ( bwtPos & sampleIntervalMask ); // get offset inside sample
	if( !( bitArray & offsetMasks[offset] ) ) letterJump++; // if the letter is not in the topPointer position, its next occurrence is bellow that
	bitArray &= searchOffsetMasks[offset]; // keep all bits bellow offset, except 1st one
	while( bitArray ){ // while there are occurrences of this letter in the interval
		letterJump++; // if other equal letter exist, increase letter jump
		bitArray &= ( bitArray - 1U ); // clear lower occurrence bit
	}
	(*topPointer) = letterJump;
	bwtPos = (*bottomPointer); // process bottom pointer
	block = &(Index[( bwtPos >> sampleIntervalShift )]);
	offset = ( bwtPos & sampleIntervalMask );
	bitArray = searchOffsetMasks[offset]; // all bits bellow offset, except 1st one
	bitArray &= ( ~ (block->specialLettersMask) ); // remove special chars
	bitArray &= ( (block->bwtLowBits) ^ inverseLetterLowBitMasks[letterId] ); // keep only ones with the same low bit
	bitArray &= ( (block->bwtHighBits) ^ inverseLetterHighBitMasks[letterId] ); // keep only ones with the same high bit
	letterJump = (block->letterJumpsSample[letterId]); // get top letter jump
	while( bitArray ){
		letterJump++; // if other equal letter exist, increase letter jump
		bitArray &= ( bitArray - 1U ); // clear lower occurrence bit
	}
	(*bottomPointer) = letterJump;
	if( (*topPointer) > (*bottomPointer) ) return 0;
	return ( (*bottomPointer) - (*topPointer) + 1 );
	*/
}

// TODO: move code of called functions inside and reuse variables for speed up
unsigned int FMI_PositionInText( unsigned int bwtpos ){
	unsigned int charid, addpos;
	addpos = 0;
	while( bwtpos & sampleIntervalMask ){ // move backwards until we land on a position with a sample
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
	return ( (Index[( bwtpos >> sampleIntervalShift )].textPositionSample) + addpos );
}

// returns the new position in the BWT array after left jumping by the char at the given BWT position
unsigned int FMI_LeftJump( unsigned int bwtpos ){
	unsigned int charid;
	charid = GetCharIdAtBWTPos(bwtpos);
	if( charid == 0 ) return 0U;
	else return FMI_LetterJump( charid , bwtpos ); // follow the left letter backwards
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
	i = ( ( ( textSize-1 ) >> sampleIntervalShift ) + 1 );
	if( ((textSize-1) & sampleIntervalMask) != 0 ) i++; // extra sample
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
	lastBwtPos = ((numSamples-1) << sampleIntervalShift); // set last BWT pos for pattern matching initialization
	printf(" OK\n");
	fflush(stdout);
}

// TODO: use indexBlocksStartInFile, distanceToNextBlockBits
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
}

void PrintBWT(char *text, unsigned int *letterStartPos){
	unsigned int i, n, p;
	printf("%u {", textSize );
	for( n = 1 ; n < ALPHABETSIZE ; n++ ){
		printf(" %c [%02u-%02u] %c", letterChars[n] , letterStartPos[n] , (n==(ALPHABETSIZE-1))?(textSize-1):(letterStartPos[(n+1)]-1) , (n==(ALPHABETSIZE-1))?'}':',' );
	}
	printf(" 2^%u=%u %u %#.8X\n", sampleIntervalShift , sampleIntervalSize , numSamples , sampleIntervalMask );
	printf("[ i] (SA) {");
	for( n = 1 ; n < ALPHABETSIZE ; n++ ){
		printf(" %c%c", letterChars[n] , (n==(ALPHABETSIZE-1))?'}':',' );
	}
	#ifdef BUILD_LCP
	if( LCPArray != NULL ) printf(" LCP");
	#endif
	printf(" BWT\n");
	for( i = 0 ; i < textSize ; i++ ){ // position in BWT
		p = FMI_PositionInText( i );
		printf("[%02u]%c(%2u) {", i , (i & sampleIntervalMask)?' ':'*' , p );
		for( n = 1 ; n < ALPHABETSIZE ; n++ ){
			printf("%02u%c", FMI_LetterJump( n , i ) , (n==(ALPHABETSIZE-1))?'}':',' );
		}
		#ifdef BUILD_LCP
		if( LCPArray != NULL ) printf(" %3d",(int)LCPArray[i]);
		#endif
		printf(" %c ", FMI_GetCharAtBWTPos(i) );
		n = 0;
		while( ( n < (ALPHABETSIZE-1) ) && ( i >= letterStartPos[(n+1)] ) ) n++;
		printf(" %c ", letterChars[n] );
		if( p != (textSize-1) ) p++;
		else p = 0;
		if( text != NULL ) printf("%s", (char *)(text+p) );
		printf("\n");
	}
}


typedef struct _PackedNumberArray {
	unsigned long long *bitsArray;
	unsigned char bitsPerInt;
	unsigned int numWords;
	//unsigned char bitsPerWord; // = 64
} PackedNumberArray;

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

/*
void PrintISState( char *header, PackedNumberArray *text , unsigned int *posLMS , unsigned int *charTypes , unsigned int textSize , int alphabetSize , unsigned int *usedLMS, int level ){
	unsigned int n;
	int i,j,k;
	if(topLMS[0]!=(-1)){ // S-type suffixes
		printf("\n%s\n",header);
		for(i=0;i<16;i++) putchar('-'); putchar('\n');
		k=0;
		for(j=0;j<alphabetSize;j++){
			i=topLMS[j];
			while(i!=(-1)){
				n=posLMS[i];
				printf("[%02d]{%02d} %c %02u %c",k,i,GETCHARTYPE(charTypes,n),n,((usedLMS && GETBIT(usedLMS,i))?'*':' '));
				if(level){
					printf("(%02u)",GetPackedNumber(text,n)); n++;
					if(n==textSize) n=0;
					while(GETCHARTYPE(charTypes,n)!='S'){ printf("(%02u)",GetPackedNumber(text,n)); n++; }
					printf("(%02u)",GetPackedNumber(text,n));
				} else {
					printf("%c",letterChars[GetPackedNumber(text,n)]); n++;
					if(n==textSize) n=0;
					while(GETCHARTYPE(charTypes,n)!='S'){ printf("%c",letterChars[GetPackedNumber(text,n)]); n++; }
					printf("%c",letterChars[GetPackedNumber(text,n)]);
				}
				putchar('\n');
				i=nextId[i];
				k++;
			}
			for(i=0;i<16;i++) putchar('-'); putchar('\n');
		}
	} else { // L-type suffixes
		printf("\n%s\n",header);
		for(i=0;i<16;i++) putchar('-'); putchar('\n');
		k=0;
		for(j=alphabetSize;j!=0;){
			j--;
			i=bottomLML[j];
			while(i!=(-1)){
				n=posLMS[i];
				printf("[%02d]{%02d} %c %02u %c",k,i,GETCHARTYPE(charTypes,n),n,((usedLMS && GETBIT(usedLMS,i))?'*':' '));
				if(level){
					printf("(%02u)",GetPackedNumber(text,n)); n++;
					if(n==textSize) n=0;
					while(GETCHARTYPE(charTypes,n)!='S'){ printf("(%02u)",GetPackedNumber(text,n)); n++; }
					printf("(%02u)",GetPackedNumber(text,n));
				} else {
					printf("%c",letterChars[GetPackedNumber(text,n)]); n++;
					if(n==textSize) n=0;
					while(GETCHARTYPE(charTypes,n)!='S'){ printf("%c",letterChars[GetPackedNumber(text,n)]); n++; }
					printf("%c",letterChars[GetPackedNumber(text,n)]);
				}
				putchar('\n');
				i=nextId[i];
				k++;
			}
			for(i=0;i<16;i++) putchar('-'); putchar('\n');
		}
	}
}
*/

// TODO: add argument and code to get text from file and directly convert it to packed binary
// TODO: change LoadReference function to use filename stored in index
// TODO: on multiple references, add special char to separate them
// TODO: use type double/float for variables used in statistics at the end
// TODO: use inverse of speciallettersmask to prevent extra "~" operation
// TODO: replace original letter*BitsMask arrays by ~inverseLetter*BitsMask
// TODO: move global variables that are not needed for searching, to inside this function
// TODO: initialize mask arrays instead of storing them explicitly
// TODO: in main functions, replace array fetches by defines with the shifts to see if it's faster
// TODO: check speed of FMI_LetterJump with only half 16bits counts
void FMI_BuildIndex(char *text, int size, char verbose){
	text=NULL;
	size=0;
	verbose='\0';
	return;
/*
	char *textPtr, c;
	unsigned int letterId, i, k, n;
	unsigned int bwtPos, samplePos, textPos;
	unsigned int progressCounter, progressStep;
	unsigned int letterCounts[6], letterStartPos[6];
	PackedNumberArray *textArray;
	IndexBlock *block;
	if(verbose){
		printf("> Processing reference chars ... ");
		fflush(stdout);
	}
	textArray = NewPackedNumberArray((size+1),6); // bit array that stores all the chars of the text (plus terminator) in packed bits form
	if(textArray==NULL){
		printf("\n> ERROR: Not enough memory\n");
		exit(-1);
	}
	InitializeLetterIdsArray(); // initialize letter ids array
	for(i=0;i<6;i++) letterCounts[i]=0; // reset counts for each letter (ACGTN)
	textSize=0; // counts all letters plus one sequence terminator symbol
	for(i=0;i<(unsigned int)size;i++){
		c=text[i];
		letterId = (unsigned int)letterIds[(int)c]; // get id of this letter
		letterCounts[letterId]++; // increase count of this letter
		textSize++; // one more text char
		if(textSize==UINT_MAX) break; // prevent overflow
		SetPackedNumber(textArray,i,letterId);
	}
	SetPackedNumber(textArray,i,0); // set terminator char (id=0) at the end
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
	BWTIS( textArray , textSize , 6 , 0 , verbose );
	FreePackedNumberArray(textArray); // the packed text array is not needed anymore, only the BWT array
	if(verbose){
		printf(" OK\n");
		printf("> Allocating memory space for index ... ");
		fflush(stdout);
	}
	numSamples = ( ( ( textSize - 1 ) >> sampleIntervalShift ) + 1 ); // blocks of 32 chars (if textSize was a multiple of 32 it needed +1 additional unused sample, because last used pos is (textSize-1))
	if( ((textSize-1) & sampleIntervalMask) != 0 ) numSamples++; // if the last used position does not fall over a sample (0-th position), add an extra sample after the end (to speed up letter counts on 1st char of pattern search)
	lastBwtPos = ((numSamples-1) << sampleIntervalShift); // so that when getting the sample for this position, if falls in the 0-th pos of the last sample (numSamples-1)
	Index=(IndexBlock *)malloc(numSamples*sizeof(IndexBlock));
	if(Index==NULL){
		printf("\n> ERROR: Not enough memory\n");
		exit(0);
	}
	if(verbose){
		printf("(%u MB) ",( (unsigned int)((numSamples)*sizeof(IndexBlock))/1000000));
		printf("OK\n");
		printf("> Collecting letter jump samples ");
		fflush(stdout);	
	}
	letterStartPos[0]=0; // the terminator char is at the top (0-th) position of the BWT (but on the right)
	for(i=1;i<6;i++) letterStartPos[i]=(letterStartPos[(i-1)]+letterCounts[(i-1)]); // where previous letter starts plus number of previous letter occurrences
	for(i=1;i<6;i++) letterCounts[i] = (letterStartPos[i]-1); // initialize all letter jumps with the position before the start of the letter
	progressStep=(textSize/10);
	progressCounter=0;
	samplePos = 0; // start in top position of the BWT and go down
	for( n = 0 ; n < textSize ; n++ ){
		if(verbose){
			progressCounter++;
			if(progressCounter==progressStep){ // print progress dots
				printf(".");
				fflush(stdout);
				progressCounter=0;
			}
		}
		letterId = GetPackedNumber(bwtArray,n);
		letterCounts[letterId]++;
		if( ( n & sampleIntervalMask ) == 0 ){ // if we are over a sample, store here the current letter counts
			block = &(Index[samplePos]);
			(block->bwtLowBits) = 0U; // reset block
			(block->bwtHighBits) = 0U;
			(block->specialLettersMask) = 0U;
			for(i=1;i<6;i++) (block->letterJumpsSample[(i-1)]) = letterCounts[i]; // the i-th letter here is the (i-1)-th letter in the index
			samplePos++;
		}
		SetCharAtBWTPos(n,letterId); // copy the current letter from the packed BWT to the BWT in the index
	}
	if( ((textSize-1) & sampleIntervalMask) != 0 ){ // if one last extra sample was added at the end, fill it too
		block = &(Index[samplePos]);
		(block->bwtLowBits) = 0U; // reset block
		(block->bwtHighBits) = 0U;
		(block->specialLettersMask) = (~0U); // set special letters block to prevent confusing with real letters
		for(i=1;i<6;i++) (block->letterJumpsSample[(i-1)]) = letterCounts[i];
	}
	FreePackedNumberArray(bwtArray); // the packed BWT array is not needed anymore
	if(verbose){
		printf(" OK\n");
		printf("> Collecting position samples ");
		fflush(stdout);
	}
	textPos=(textSize-1); // start at the terminator char (first text position is 0)
	n=0; // start at the first/topmost BWT position
	progressStep=(textSize/10);
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
		if( ( n & sampleIntervalMask ) == 0 ){ // if we are over a sample, store here the current position of the text
			samplePos = ( n >> sampleIntervalShift );
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
	bwtPos=0; // just so compiler does not complain
	textPtr=NULL;
	#ifdef DEBUG
	if(textSize<100) PrintBWT(text,letterStartPos);
	if(verbose){
		printf("> Checking BWT sort ");
		fflush(stdout);
	}
	progressCounter=0;
	k = FMI_PositionInText(0); // position in the text of the top BWT position (should be equal to (textSize-1))
	for(bwtPos=1;bwtPos<textSize;bwtPos++){ // compare current position with position above
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
		while( text[(k+n)] == text[(i+n)] ) n++; // keep following suffix chars to the right while their letters are equal
		if( letterIds[(int)(text[(k+n)])] > letterIds[(int)(text[(i+n)])] ) break; // if the top letter is larger than the bottom letter, it is incorrectly sorted
		k = i; // current pos will be prev pos in next step
	}
	if(bwtPos!=textSize){
		printf(" FAILED (error at BWT position %u)\n",bwtPos);
		exit(-1);
	}
	if(verbose){
		printf(" OK\n");
		printf("> Checking jump samples ");
		fflush(stdout);
	}
	for(i=0;i<6;i++){ // initialize counters
		if(letterStartPos[i]!=0) letterCounts[i]=(letterStartPos[i]-1);
		else letterCounts[i]=0;
	}
	progressCounter=0;
	for(n=0;n<textSize;n++){
		if(verbose){
			progressCounter++;
			if(progressCounter==progressStep){ // print progress dots
				printf(".");
				fflush(stdout);
				progressCounter=0;
			}
		}
		i=GetCharIdAtBWTPos(n); // get letter and update count
		if(i!=0){ // check if it is the terminator char because we cannot jump by it
			letterCounts[i]++;
			if( FMI_LetterJump(i,n) != letterCounts[i] ) break;
		}
		if( ( n & sampleIntervalMask ) == 0 ){ // if there is a sample at this position, check counts
			samplePos = ( n >> sampleIntervalShift );
			for(i=1;i<6;i++) if( (Index[samplePos].letterJumpsSample[(i-1)]) != (letterCounts[i]) ) break;
			if(i!=6) break;
		}
	}
	if(n!=textSize){
		printf(" FAILED (error at BWT position %u)\n",n);
		exit(-1);
	}
	if(verbose){
		printf(" OK\n");
		printf("> Checking position samples ");
		fflush(stdout);
	}
	progressCounter=0;
	textPos=(textSize-1); // start at the terminator char
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
		if( FMI_PositionInText(n) != textPos ) break;
		if( ( n & sampleIntervalMask ) == 0 ){ // if there is a sample at this BWT position, check text position
			samplePos = ( n >> sampleIntervalShift );
			if( (Index[samplePos].textPositionSample) != textPos ) break;
		}
		if(textPos==0) break;
		textPos--;
		c=text[textPos];
		letterId=letterIds[(int)c];
		i=GetCharIdAtBWTPos(n);
		if(letterId!=i) break; // check if it is the same letter at the BWT and at the text
		if(i==0) break; // check if it is the terminator char because it should not be here
		n = FMI_LetterJump(i,n); // follow the letter backwards to go to next position in the BWT
	}
	if( textPos!=0 || FMI_PositionInText(n)!=0 || letterId!=i || i==0 || GetCharIdAtBWTPos(n)!=0 ){
		printf(" FAILED (error at text position %u)\n",textPos);
		exit(-1);
	}
	if(verbose){
		printf(" OK\n");
		fflush(stdout);
	}
	#endif
	*/
}

// TODO: add code/functions for partial string text, circular text and packed binary text
unsigned int GetTextCharId(unsigned int pos){
	//if(pos==textSize) return (unsigned int)letterIds[(unsigned char)'$'];
	return (unsigned int)letterIds[(unsigned char)text[pos]];
}

typedef struct _LMSPos {
	unsigned int pos;
	#ifdef BUILD_LCP
	int lcp;
	#endif
	int next;
} LMSPos;

LMSPos *LMSArray;
int numLMS;


char GetCharType(unsigned int pos){
	char prevType, currentType;
	if(pos==textSize) return 'S'; // last position
	if(pos==0 || GetTextCharId(pos-1)<GetTextCharId(pos)) prevType='s';
	else if(GetTextCharId(pos-1)>GetTextCharId(pos)) prevType='l';
	else prevType='?';
	while(pos!=(textSize-1) && GetTextCharId(pos)==GetTextCharId(pos+1)) pos++;
	if(pos==(textSize-1) || GetTextCharId(pos)>GetTextCharId(pos+1)) currentType='l';
	else currentType='s';
	if(prevType=='?') prevType=currentType;
	if(prevType!=currentType) currentType-=32; // S*-type or L*-type
	return currentType;
}

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
	progressStep = (textSize/10);
	progressCounter = 0;
	charsLastPos = (int *)malloc(ALPHABETSIZE*sizeof(int));
	for( i = 0 ; i < ALPHABETSIZE ; i++ ){ // initialize chars buckets
		charsCounts[i] = 0;
		charsFirstPos[i] = (-1);
		charsLastPos[i] = (-1);
	}
	arrayGrowSize = (textSize/20); // how much to expand the LMS array each time we allocate more memory (5%)
	if( arrayGrowSize == 0 ) arrayGrowSize = 20;
	arrayMaxSize = 0;
	LMSArray = NULL;
	// TODO: check if n should be textSize or (textSize-1) here
	n = textSize;
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
					LMSArray = (LMSPos *)realloc(LMSArray,arrayMaxSize*sizeof(LMSPos));
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
					printf("\n> ERROR: LMS[%d]@text[%d]='%c' > LMS[%d]@text[%d]='%c'\n",statePos,prevTextPos,ALPHABET[sortedCharId],(statePos+1),textPos,ALPHABET[charId]);
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

// TODO: check if wrap around goes to position (textSize) or (textSize-1)
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
	progressStep = (textSize/10);
	progressCounter = 0;
	/*
	bucketSize = (unsigned int *)malloc(ALPHABETSIZE*sizeof(unsigned int)); // set pointers to the beginning of the buckets
	for( charId = 0 ; charId < ALPHABETSIZE ; charId++ ) bucketSize[charId] = 0; // reset bucket size
	for( textPos = 0 ; textPos < textSize ; textPos++ ) bucketSize[ GetTextCharId(textPos) ]++; // count number of each alphabet letter in the text (size of each bucket)
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
			else textPos = textSize;
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
				if( (bucketPointer[charId] & sampleIntervalMask) == 0 ){ // add (textPos+1) sample to BWT Index here if the BWT position is a multiple of the sampling interval
					Index[ (bucketPointer[charId] >> sampleIntervalShift) ].textPositionSample = LMSArray[arrayPos].pos;
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
			else textPos = textSize; // last position in text (it is the terminal symbol '$', which is an S-type char)
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
				if( (bucketPointer[charId] & sampleIntervalMask) == 0 ){ // add (textPos+1) sample to BWT Index here if the BWT position is a multiple of the sampling interval
					Index[ (bucketPointer[charId] >> sampleIntervalShift) ].textPositionSample = LMSArray[arrayPos].pos;
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
		if( textPos <= textSize ) LCPArray[textPos] = 0;
	}
	#endif
	if(verbose){
		printf(" OK\n");
		fflush(stdout);
	}
}

// TODO: set correct textSize
// TODO: set letterCounts and letterLMSStartPos as global and change names inside functions
// TODO: add stats for average and longest length of runs in the BWT
void FMI_NewBuildIndex(char *inputText, unsigned int inputTextSize, unsigned char **lcpArrayPointer, char verbose){
	unsigned int letterId, i, n;
	unsigned int textPos, samplePos;
	unsigned int *letterCounts, *letterStartPos;
	int *letterLMSStartPos;
	//PackedNumberArray *textArray;
	IndexBlock *block;
	unsigned int progressCounter, progressStep;
	#ifdef DEBUG_INDEX
	unsigned int bwtPos, prevLetterId, k;
	unsigned int numRuns, sizeRun, longestRun;
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
	text = inputText;
	textSize = inputTextSize;
	InitializeLetterIdsArray(); // initialize letter ids array
	letterCounts = (unsigned int *)malloc(ALPHABETSIZE*sizeof(unsigned int));
	letterLMSStartPos = (int *)malloc(ALPHABETSIZE*sizeof(int));

	LMSArray = NULL;
	GetLMSs(letterCounts,letterLMSStartPos,verbose);
	SortLMSs(letterLMSStartPos,verbose);
	
	#ifdef FILL_INDEX
	textSize++; // count terminator char
	numSamples = ( ( ( textSize - 1 ) >> sampleIntervalShift ) + 1 ); // blocks of 32 chars (if textSize was a multiple of 32 it needed +1 additional unused sample, because last used pos is (textSize-1))
	if( ((textSize-1) & sampleIntervalMask) != 0 ) numSamples++; // if the last used position does not fall over a sample (0-th position), add an extra sample after the end (to speed up letter counts on 1st char of pattern search)
	lastBwtPos = ((numSamples-1) << sampleIntervalShift); // so that when getting the sample for this position, if falls in the 0-th pos of the last sample (numSamples-1)
	Index=(IndexBlock *)calloc(numSamples,sizeof(IndexBlock));
	if(Index==NULL){
		printf("> ERROR: Not enough memory to create index\n");
		exit(0);
	}
	textSize--;
	packedBwt = NULL;
	#else
	Index = NULL;
	packedBwt = NewPackedNumberArray((textSize+1),ALPHABETSIZE); // bit array that stores all the chars of the BWT in packed bits form
	#endif

	#ifdef BUILD_LCP
	#ifdef UNBOUNDED_LCP
	LCPArray = (int *)malloc((textSize+1)*sizeof(int));
	(*lcpArrayPointer) = NULL;
	#else
	LCPArray = (unsigned char *)malloc((textSize+1)*sizeof(unsigned char));
	(*lcpArrayPointer) = LCPArray; // output LCP array as pointer in argument
	#endif
	#endif
	InducedSort(letterCounts,letterLMSStartPos,verbose);
	free(LMSArray);
	free(letterLMSStartPos);
	textSize++;

	#ifndef FILL_INDEX // allocate index memory now if we did not fill it while building the BWT
	numSamples = ( ( ( textSize - 1 ) >> sampleIntervalShift ) + 1 ); // blocks of 32 chars (if textSize was a multiple of 32 it needed +1 additional unused sample, because last used pos is (textSize-1))
	if( ((textSize-1) & sampleIntervalMask) != 0 ) numSamples++; // if the last used position does not fall over a sample (0-th position), add an extra sample after the end (to speed up letter counts on 1st char of pattern search)
	lastBwtPos = ((numSamples-1) << sampleIntervalShift); // so that when getting the sample for this position, if falls in the 0-th pos of the last sample (numSamples-1)
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
	progressStep=(textSize/10);
	progressCounter=0;
	letterId=0; // just to fix compiler uninitialized warning
	textPos=0;
	samplePos = 0; // start in top position of the BWT and go down
	for( n = 0 ; n < textSize ; n++ ){
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
		letterCounts[letterId]++;
		if( ( n & sampleIntervalMask ) == 0 ){ // if we are over a sample, store here the current letter counts
			block = &(Index[samplePos]);
			#ifndef FILL_INDEX
			(block->bwtLowBits) = 0U; // reset block
			(block->bwtHighBits) = 0U;
			(block->specialLettersMask) = 0U;
			#endif
			for(i=1;i<ALPHABETSIZE;i++) (block->letterJumpsSample[(i-1)]) = letterCounts[i]; // the i-th letter here is the (i-1)-th letter in the index
			samplePos++;
		}
		#ifndef FILL_INDEX
		SetCharAtBWTPos(n,letterId); // copy the current letter from the packed BWT to the BWT in the index
		#endif
	}
	if( ((textSize-1) & sampleIntervalMask) != 0 ){ // if one last extra sample was added at the end, fill it too
		block = &(Index[samplePos]);
		(block->bwtLowBits) = 0U; // reset block
		(block->bwtHighBits) = 0U;
		(block->specialLettersMask) = (~0U); // set special letters block to prevent confusing with real letters
		for(i=1;i<6;i++) (block->letterJumpsSample[(i-1)]) = letterCounts[i];
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
	textPos=(textSize-1); // start at the terminator char (first text position is 0)
	n=0; // start at the first/topmost BWT position
	progressStep=(textSize/10);
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
		if( ( n & sampleIntervalMask ) == 0 ){ // if we are over a sample, store here the current position of the text
			samplePos = ( n >> sampleIntervalShift );
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
		printf(":: Short LCP Array size = %u MB\n",( (unsigned int)((textSize)*sizeof(unsigned char))/1000000));
		#endif
		fflush(stdout);	
	}
	#ifdef DEBUG_INDEX
	if(textSize<100) PrintBWT(text,letterStartPos);
	if(verbose){
		printf("> Checking BWT ");
		#if ( defined(BUILD_LCP) && defined(UNBOUNDED_LCP) )
		printf("and LCP ");
		#endif
		fflush(stdout);
	}
	progressStep=(textSize/10);
	progressCounter=0;
	numRuns=1;
	sizeRun=1;
	longestRun=0;
	k = FMI_PositionInText(0); // position in the text of the top BWT position (should be equal to (textSize-1))
	for(bwtPos=1;bwtPos<textSize;bwtPos++){ // compare current position with position above
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
			printf("\t[%c] %c|",GetCharType(k),letterChars[GetCharIdAtBWTPos((bwtPos-1))]);
			for( textPos=k ; textPos<=(k+n) ; textPos++ ) putchar(letterChars[GetTextCharId(textPos)]);
			putchar('\n');
			printf("\t[%c] %c|",GetCharType(i),letterChars[GetCharIdAtBWTPos(bwtPos)]);
			for( textPos=i ; textPos<=(i+n) ; textPos++ ) putchar(letterChars[GetTextCharId(textPos)]);
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
			if( sizeRun > longestRun ) longestRun = sizeRun;
			sizeRun = 1;
			numRuns++;
		}
		k = i; // current pos will be prev pos in next step
	}
	if(bwtPos!=textSize){
		printf(" FAILED (error at BWT position %u)\n",bwtPos);
		getchar();
		exit(-1);
	}
	if(verbose){
		printf(" OK\n");
		printf(":: Average run length = %.2lf (max = %u)\n",((double)textSize/(double)numRuns),longestRun);
		printf("> Checking LF samples ");
		fflush(stdout);
	}
	for(i=0;i<ALPHABETSIZE;i++){ // initialize counters
		if(letterStartPos[i]!=0) letterCounts[i]=(letterStartPos[i]-1);
		else letterCounts[i]=0;
	}
	progressCounter=0;
	for(n=0;n<textSize;n++){
		if(verbose){
			progressCounter++;
			if(progressCounter==progressStep){ // print progress dots
				printf(".");
				fflush(stdout);
				progressCounter=0;
			}
		}
		i=GetCharIdAtBWTPos(n); // get letter and update count
		if(i!=0){ // check if it is the terminator char because we cannot jump by it
			letterCounts[i]++;
			if( FMI_LetterJump(i,n) != letterCounts[i] ) break;
		}
		if( ( n & sampleIntervalMask ) == 0 ){ // if there is a sample at this position, check counts
			samplePos = ( n >> sampleIntervalShift );
			for(i=1;i<ALPHABETSIZE;i++) if( (Index[samplePos].letterJumpsSample[(i-1)]) != (letterCounts[i]) ) break;
			if(i!=ALPHABETSIZE) break;
		}
	}
	if(n!=textSize){
		printf(" FAILED (error at BWT position %u)\n",n);
		getchar();
		exit(-1);
	}
	if(verbose){
		printf(" OK\n");
		printf("> Checking SA samples ");
		fflush(stdout);
	}
	numRuns=0;
	longestRun=0;
	progressCounter=0;
	textPos=(textSize-1); // start at the terminator char
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
		if( ( n & sampleIntervalMask ) == 0 ){ // if there is a sample at this BWT position, check text position
			samplePos = ( n >> sampleIntervalShift );
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
		printf(":: Average backtracking steps = %.2lf (max = %u)\n",((double)numRuns/(double)textSize),longestRun);
		fflush(stdout);
	}
	#endif
	free(letterCounts);
	free(letterStartPos);
}
