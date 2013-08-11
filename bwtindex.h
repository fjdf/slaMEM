unsigned int FMI_PositionInText( unsigned int bwtpos );
unsigned int FMI_FollowLetter( char c , unsigned int *topPointer , unsigned int *bottomPointer );
void FMI_GetCharCountsAtBWTInterval( unsigned int topPtr , unsigned int bottomPtr , int *counts );
unsigned int FMI_LeftJump( unsigned int bwtpos );
char FMI_GetCharAtBWTPos(unsigned int bwtpos);
void FMI_FreeIndex();
void FMI_BuildIndex(char *text, int size, char verbose);
unsigned int FMI_GetTextSize();
unsigned int FMI_GetBWTSize();
char *FMI_GetTextFilename();
void FMI_NewBuildIndex(char *inputText, unsigned int inputTextSize, unsigned char **lcpArrayPointer, char verbose);

