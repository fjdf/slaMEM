FILE *debugfile;

char* AppendToBasename(char *filename, char *extra);
char *NormalizeSeqName(char *name, int mode);
int ParseArgument(int numargs, char** arglist, char *optionchars, int parse);
int ParseNumber(char *numberstring);
void PrintNumber(int number);
void PrintSpace(int spaceval);
void PrintTime(double timeval);
void PrintProgressBar(double percentage, int lineabove);
void exitMessage(char *msg);
