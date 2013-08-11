int BuildSampledLCPArray(char *text, int textsize, unsigned char *lcparray, int verbose);
void FreeSampledSuffixArray();
int GetLCP(unsigned int bwtpos);
int GetEnclosingLCPInterval(unsigned int *topptr, unsigned int *bottomptr);
int GetLCPValuePositionInsideLCPInterval(int lcpvalue, unsigned int topptr, unsigned int bottomptr);
int GetEnclosingLCPIntervalFromLCPPos(int lcppos, unsigned int *topptr, unsigned int *bottomptr);
