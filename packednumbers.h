typedef struct _PackedNumberArray {
	unsigned long long *bitsArray;
	unsigned char bitsPerInt;
	unsigned int numWords;
	//unsigned char bitsPerWord; // = 64
} PackedNumberArray;

PackedNumberArray *NewPackedNumberArray(unsigned int numInts, unsigned int maxInt);
void FreePackedNumberArray(PackedNumberArray *intArray);
unsigned int GetPackedNumber(PackedNumberArray *intArray, unsigned int pos);
void SetPackedNumber(PackedNumberArray *intArray, unsigned int pos, unsigned int num);
