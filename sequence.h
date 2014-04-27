typedef struct Sequence {
	unsigned int size;
	char *name;
	char *chars;
	//int rotation;
	int order;
	char strand;
	unsigned char fileid;
	fpos_t sourcefilepos;
	//char *sourcefilename;
	//FILE *sourcefile;
} Sequence;

int numSequences;
Sequence **allSequences;

int LoadSequencesFromFile(char *inputfilename, int loadchars, int mergeseqs, int acgtonly);
void DeleteAllSequences();
void LoadSequenceChars(Sequence *seq);
void FreeSequenceChars(Sequence *seq);
int GetSeqIdFromSeqName(char *seqname);
int GetSeqIdFromMergedSeqsPos(unsigned int *pos);
/*
char CharAt(int pos, int seqid);
char GetNextChar(int seqid);
int GetNextCharCode(int seqid);
*/
void SortSequences(int *seqsizes, int *sortedseqs, int numseqs);
void ReverseComplementSequence(char *text, int textsize);
