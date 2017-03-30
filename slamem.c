/* ================================================================= *
 *  slamem.c : Main application                                      *
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
#include "tools.h"
#include "sequence.h"
#include "bwtindex.h"
#include "lcparray.h"
#include "graphics.h"

#define VERSION "0.8.1"

//#define BENCHMARK 1
#if defined(unix) && defined(BENCHMARK)
#include "unistd.h"
#include "string.h"
#endif

#ifdef _MSC_VER
#define PAUSE_AT_EXIT 1
#endif

#define MATCH_TYPE_CHAR "EAU"

void GetMatches(int numRefs, int numSeqs, int matchType, int minMatchSize, int bothStrands, char *outFilename){
	FILE *matchesOutputFile;
	int i, s, depth, matchSize, numMatches, refId;
	unsigned int j, textsize, refPos;
	long long int sumMatchesSize, totalNumMatches, totalAvgMatchesSize;
	unsigned int topPtr, bottomPtr, prevTopPtr, prevBottomPtr, savedTopPtr, savedBottomPtr, n;
	char c, *text;
	char *refsTexts[1];
	unsigned int refsTextSizes[1];
	unsigned char *lcpArray;
	int progressCounter, progressStep;
	#ifdef DEBUGMEMS
	char *refText;
	int refSize;
	#endif
	#if defined(unix) && defined(BENCHMARK)
	char command[32];
	int commretval;
	#endif
	printf("> Using options: minimum M%cM length = %d\n",MATCH_TYPE_CHAR[matchType],minMatchSize);
	matchesOutputFile=fopen(outFilename,"w");
	if(matchesOutputFile==NULL){
		printf("\n> ERROR: Cannot create output file <%s>\n",outFilename);
		exit(-1);
	}
	printf("> Building index for reference sequence");
	if(numRefs==1) printf(" \"%s\"", (allSequences[0]->name));
	else printf("s");
	printf(" (%u Mbp) ...\n",(allSequences[0]->size)/1000000U);
	fflush(stdout);
	text=(allSequences[0]->chars);
	textsize=(allSequences[0]->size);
	refId=0;
	refsTexts[0]=text;
	refsTextSizes[0]=textsize;
	lcpArray=NULL;
	FMI_BuildIndex(refsTexts,refsTextSizes,1,&lcpArray,1);
	i=BuildSampledLCPArray(text,textsize,lcpArray,minMatchSize,1);
	if(lcpArray!=NULL) free(lcpArray);
	#ifndef DEBUGMEMS
	FreeSequenceChars(allSequences[0]);
	#else
	refText=(allSequences[0]->chars);
	refSize=(allSequences[0]->size);
	#endif
	#if defined(unix) && defined(BENCHMARK)
	sprintf(command,"memusgpid %d &",(int)getpid());
	commretval=system(command);
	#endif
	printf("> Matching query sequences against index ...\n");
	fflush(stdout);
	totalNumMatches=0;
	totalAvgMatchesSize=0;
	for(i=numRefs;i<numSeqs;i++){ // process all queries
		LoadSequenceChars(allSequences[i]);
		text=(allSequences[i]->chars);
		textsize=(allSequences[i]->size);
		progressStep=(textsize/10);
		for(s=0;s<=bothStrands;s++){ // process one or both strands
			if(s==0){ // forward strand
				printf(":: \"%s\" ",(allSequences[i]->name));
				fprintf(matchesOutputFile,">%s\n",(allSequences[i]->name));
			} else { // reverse strand
				ReverseComplementSequence(text,textsize); // convert to reverse strand
				printf(":: \"%s Reverse\" ",(allSequences[i]->name));
				fprintf(matchesOutputFile,">%s Reverse\n",(allSequences[i]->name));
			}
			fflush(stdout);
			progressCounter=0;
			matchSize=0;
			numMatches=0;
			sumMatchesSize=0;
			depth=0;
			topPtr=0;
			bottomPtr=FMI_GetBWTSize();
			prevTopPtr=topPtr;
			prevBottomPtr=bottomPtr;
			for(j=textsize;j!=0;){
				j--;
				if(progressCounter==progressStep){ // print progress dots
					putchar('.');
					fflush(stdout);
					progressCounter=0;
				} else progressCounter++;
				while( (n=FMI_FollowLetter(text[j],&topPtr,&bottomPtr))==0 ){ // when no match exits, follow prefix links to broaden the interval
					topPtr = prevTopPtr; // restore pointer values, because they got lost when no hits exist
					bottomPtr = prevBottomPtr;
					depth = GetEnclosingLCPInterval(&topPtr,&bottomPtr); // get enclosing interval and corresponding destination depth
					if( depth == -1 ) break; // can happen for example when current seq contains 'N's but the indexed reference does not
					prevTopPtr = topPtr; // save pointer values in case the match fails again
					prevBottomPtr = bottomPtr;
				}
				depth++;
				if( depth >= minMatchSize ){
					if(matchType==1 && n!=1) continue; // not a MAM if we are looking for one
					savedTopPtr = topPtr; // save the original interval to restore after finished processing MEMs
					savedBottomPtr = bottomPtr;
					prevTopPtr = (bottomPtr+1); // to process the first interval entirely
					prevBottomPtr = bottomPtr;
					matchSize = depth;
					if( j != 0 ) c = text[j-1]; // next char to be processed (to the left)
					else c = '\0';
					while( matchSize >= minMatchSize ){ // process all parent intervals down to this size limit
						for( n = topPtr ; n != prevTopPtr ; n++ ){ // from topPtr down to prevTopPtr
							if( FMI_GetCharAtBWTPos(n) != c ){
								refPos = FMI_PositionInText(n);
								#ifndef DEBUGMEMS
								if(numRefs!=1){ // multiple refs
									refId = GetSeqIdFromMergedSeqsPos(&refPos); // get ref id and pos inside that ref
									fprintf(matchesOutputFile," %s\t",(allSequences[refId]->name));
								}
								fprintf(matchesOutputFile,"%u\t%d\t%d\n",(refPos+1),(j+1),matchSize);
								#else
								fprintf(matchesOutputFile,"%u\t%d\t%d",(refPos+1),(j+1),matchSize);
								fputc('\t',matchesOutputFile);
								fputc((refPos==0)?('$'):(refText[refPos-1]+32),matchesOutputFile);
								fprintf(matchesOutputFile,"%.*s...%.*s",4,(char *)(refText+refPos),4,(char *)(refText+refPos+matchSize-4));
								fputc(((refPos+matchSize)==refSize)?('$'):(refText[refPos+matchSize]+32),matchesOutputFile);
								fputc('\t',matchesOutputFile);
								fputc((j==0)?('$'):(text[j-1]+32),matchesOutputFile);
								fprintf(matchesOutputFile,"%.*s...%.*s",4,(char *)(text+j),4,(char *)(text+j+matchSize-4));
								fputc(((j+matchSize)==textsize)?('$'):(text[j+matchSize]+32),matchesOutputFile);
								fputc('\n',matchesOutputFile);
								#endif
								numMatches++;
								sumMatchesSize += matchSize;
							}
						}
						for( n = bottomPtr ; n != prevBottomPtr ; n-- ){ // from bottomPtr up to prevBottomPtr
							if( FMI_GetCharAtBWTPos(n) != c ){
								refPos = FMI_PositionInText(n);
								#ifndef DEBUGMEMS
								if(numRefs!=1){ // multiple refs
									refId = GetSeqIdFromMergedSeqsPos(&refPos); // get ref id and pos inside that ref
									fprintf(matchesOutputFile," %s\t",(allSequences[refId]->name));
								}
								fprintf(matchesOutputFile,"%u\t%d\t%d\n",(refPos+1),(j+1),matchSize);
								#else
								fprintf(matchesOutputFile,"%u\t%d\t%d",(refPos+1),(j+1),matchSize);
								fputc('\t',matchesOutputFile);
								fputc((refPos==0)?('$'):(refText[refPos-1]+32),matchesOutputFile);
								fprintf(matchesOutputFile,"%.*s...%.*s",4,(char *)(refText+refPos),4,(char *)(refText+refPos+matchSize-4));
								fputc(((refPos+matchSize)==refSize)?('$'):(refText[refPos+matchSize]+32),matchesOutputFile);
								fputc('\t',matchesOutputFile);
								fputc((j==0)?('$'):(text[j-1]+32),matchesOutputFile);
								fprintf(matchesOutputFile,"%.*s...%.*s",4,(char *)(text+j),4,(char *)(text+j+matchSize-4));
								fputc(((j+matchSize)==textsize)?('$'):(text[j+matchSize]+32),matchesOutputFile);
								fputc('\n',matchesOutputFile);
								#endif
								numMatches++;
								sumMatchesSize += matchSize;
							}
						}
						prevTopPtr = topPtr;
						prevBottomPtr = bottomPtr;
						matchSize = GetEnclosingLCPInterval(&topPtr,&bottomPtr); // get parent interval and its depth
					}
					topPtr = savedTopPtr;
					bottomPtr = savedBottomPtr;
				}
				prevTopPtr=topPtr; // save pointer values in case there's no match on the next char, and they loose their values
				prevBottomPtr=bottomPtr;
			} // end of loop for all chars of seq
			totalNumMatches += numMatches;
			totalAvgMatchesSize += sumMatchesSize;
			matchSize=(int)((numMatches==0)?(0):(sumMatchesSize/(long long)numMatches));
			printf(" (%d M%cMs ; avg size = %d bp)\n",numMatches,MATCH_TYPE_CHAR[matchType],matchSize);
			fflush(stdout);
		} // end of loop for both strands
		FreeSequenceChars(allSequences[i]);
	} // end of loop for all queries
	FMI_FreeIndex();
	FreeSampledSuffixArray();
	if((numSeqs-numRefs)!=1){ // if more than one query, print average stats for all queries
		printf(":: Average %d M%cMs found per query sequence (total = %lld, avg size = %d bp)\n",(int)(totalNumMatches/(numSeqs-numRefs)),MATCH_TYPE_CHAR[matchType],totalNumMatches,(int)(totalAvgMatchesSize/totalNumMatches));
	}
	fflush(stdout);
	printf("> Saving M%cMs to <%s> ... ",MATCH_TYPE_CHAR[matchType],outFilename);
	fclose(matchesOutputFile);
	printf("OK\n");
	fflush(stdout);
}

typedef struct _MEMInfo {
	char refName[65];
	int refPos;
	int queryPos;
	int size;
} MEMInfo;

int MEMInfoSortFunction(const void *a, const void *b){
	char *charsa, *charsb;
	int diff;
	diff = 0;
	charsa = (((MEMInfo *)a)->refName);
	charsb = (((MEMInfo *)b)->refName);
	while((diff=(int)((*charsa)-(*charsb)))==0 && (*charsa)!='\0'){ charsa++; charsb++; }
	if(diff==0){
		diff = ((((MEMInfo *)a)->refPos) - (((MEMInfo *)b)->refPos));
		if(diff==0) diff = ((((MEMInfo *)a)->queryPos) - (((MEMInfo *)b)->queryPos));
	}
	return diff;
}


// NOTE: the spacing of the MUMmer output format is in the form: "  <max_ref_name_size> <9_spaces_number> <9_spaces_number> <9_spaces_number>"
// NOTE: in multi-ref format (4 columns) the name of the ref is considered only up to the 1st space char
void SortMEMsFile(char *memsFilename){
	FILE *memsFile, *sortedMemsFile;
	char c, *sortedMemsFilename, seqname[256];
	int numMems, maxNumMems, numSeqs, refpos, querypos, memsize, formatNumFields, n;
	MEMInfo *memsArray;
	printf("> Sorting MEMs from <%s> ",memsFilename);
	fflush(stdout);
	if((memsFile=fopen(memsFilename, "r"))==NULL){
		printf("\n> ERROR: Cannot read input file\n");
		exit(-1);
	}
	c=fgetc(memsFile);
	if(c!='>'){
		printf("\n> ERROR: Invalid MEMs file\n");
		exit(-1);
	}
	while(c=='>'){
		c=fgetc(memsFile);
		while(c!='\n' && c!=EOF) c=fgetc(memsFile);
		c=fgetc(memsFile);
	}
	if(c==EOF){
		printf("\n> ERROR: No MEMs inside file\n");
		exit(-1);
	}
	formatNumFields=0;
	while(1){
		while(c==' ' || c=='\t') c=fgetc(memsFile);
		if(c!='\n'){
			formatNumFields++;
			while(c!=' ' && c!='\t' && c!='\n' && c!=EOF) c=fgetc(memsFile);
		}
		if(c=='\n' || c==EOF) break;
	}
	if(formatNumFields!=3 && formatNumFields!=4){
		printf("\n> ERROR: Invalid MEMs file format\n");
		exit(-1);
	}
	rewind(memsFile);
	if(formatNumFields==4) printf("(multiple references) ");
	printf("...\n");
	sortedMemsFilename=AppendToBasename(memsFilename,"-sorted.txt");
	if((sortedMemsFile=fopen(sortedMemsFilename,"w"))==NULL){
		printf("> ERROR: Cannot write output file\n");
		exit(-1);
	}
	seqname[0]='\0';
	refpos=-1;
	querypos=-1;
	memsize=-1;
	numSeqs=0;
	numMems=0;
	maxNumMems=0;
	memsArray=NULL;
	while(1){
		c=fgetc(memsFile);
		if(c=='>' || c==EOF){
			if(numSeqs!=0){
				printf("(%d MEMs)\n",numMems);
				fflush(stdout);
				qsort(memsArray,numMems,sizeof(MEMInfo),MEMInfoSortFunction);
				fprintf(sortedMemsFile,">%s\n",seqname);
				while(numMems!=0){
					numMems--;
					refpos=memsArray[numMems].refPos;
					querypos=memsArray[numMems].queryPos;
					memsize=memsArray[numMems].size;
					if(formatNumFields==4) fprintf(sortedMemsFile," %s\t",(memsArray[numMems].refName));
					fprintf(sortedMemsFile,"%d\t%d\t%d\n",refpos,querypos,memsize);
				}
			}
			if(c==EOF) break;
			numMems=0;
			n=fscanf(memsFile," %255[^\n]\n",seqname);
			printf(":: '%s' ... ",seqname);
			fflush(stdout);
			numSeqs++;
			continue;
		} else ungetc(c,memsFile);
		if(numMems==maxNumMems){
			maxNumMems+=1024;
			memsArray=(MEMInfo *)realloc(memsArray,(maxNumMems*sizeof(MEMInfo)));
		}
		if(formatNumFields==4) n=fscanf(memsFile," %64[^\t ]",(memsArray[numMems].refName));
		else memsArray[numMems].refName[0]='\0';
		n=fscanf(memsFile," %d %d %d ",&refpos,&querypos,&memsize);
		if(n!=3){
			printf("\n> ERROR: Invalid format\n");
			getchar();
			exit(-1);
		}
		memsArray[numMems].refPos=refpos;
		memsArray[numMems].queryPos=querypos;
		memsArray[numMems].size=memsize;
		numMems++;
	}
	printf("> Saving sorted MEMs to <%s> ...\n",sortedMemsFilename);
	fflush(stdout);
	fclose(memsFile);
	fclose(sortedMemsFile);
	free(sortedMemsFilename);
	free(memsArray);
	printf("> Done!\n");
	#ifdef PAUSE_AT_EXIT
	getchar();
	#endif
	exit(0);
}

void CreateMemMapImage(char *refFilename, char *queryFilename, char *memsFilename){
	FILE *memsFile;
	char c, seqname[256], **seqsNames, *imageFilename;
	int i, numSeqs, *seqsSizes, numMems, refPos, queryPos, memSize, strand, k;
	numSequences=0;
	numSeqs=LoadSequencesFromFile(refFilename,0,0,0);
	if(numSeqs==0){
		printf("\n> ERROR: Reference sequence not found\n");
		exit(-1);
	}
	if(numSeqs>1){
		printf("\n> ERROR: Support for visualizing multiple reference sequences will be added later but is not available yet\n");
		exit(-1);
	}
	numSeqs=LoadSequencesFromFile(queryFilename,0,0,0);
	if(numSeqs==0){
		printf("\n> ERROR: No query sequences found\n");
		exit(-1);
	}
	printf("> Processing MEMs from <%s> ...\n",memsFilename);
	fflush(stdout);
	if((memsFile=fopen(memsFilename,"r"))==NULL){
		printf("\n> ERROR: Cannot read file\n");
		exit(-1);
	}
	seqsSizes=(int *)malloc(numSequences*sizeof(int));
	seqsNames=(char **)malloc(numSequences*sizeof(char *));
	for(i=0;i<numSequences;i++){
		seqsSizes[i]=(allSequences[i]->size);
		seqsNames[i]=(allSequences[i]->name);
	}
	InitializeRefAlignmentImage(seqsSizes,numSequences);
	seqname[0]='\0';
	refPos=-1;
	queryPos=-1;
	memSize=-1;
	strand=0;
	numSeqs=0;
	numMems=0;
	while(1){
		c=fgetc(memsFile);
		if(c=='>' || c==EOF){
			if(numSeqs!=0){
				printf("(%d MEMs)\n",numMems);
				fflush(stdout);
			}
			if(c==EOF) break;
			numMems=0;
			k=fscanf(memsFile," %255[^\n]\n",seqname);
			printf(":: '%s' ... ",seqname);
			fflush(stdout);
			i=0; // check if seq name ends with string "Reverse"
			while(seqname[i]!='\0') i++;
			if(i>7) i-=7;
			for(k=0;k<8;k++,i++) if(("Reverse"[k])!=seqname[i]) break;
			if(k==8) strand=1; // reverse strand of the previous sequence
			else {
				strand=0;
				numSeqs++;
				if(numSeqs==numSequences){
					printf("\n> ERROR: MEMs file not generated from this query file (too many sequences)\n");
					exit(-1);
				}
			}
			continue;
		} else ungetc(c,memsFile);
		k=fscanf(memsFile," %d %d %d ",&refPos,&queryPos,&memSize);
		if(k!=3){
			printf("\n> ERROR: Invalid MEM format\n");
			exit(-1);
		}
		if( refPos==0 || queryPos==0 || memSize==0 ){
			printf("\n> ERROR: Invalid MEM values\n");
			exit(-1);
		}
		refPos--; // convert from 1-based position to 0-based position
		queryPos--;
		if(strand==0) DrawRefAlignmentBlock(queryPos,refPos,memSize,numSeqs);
		else { // the position is relative to the rev strand, so convert it to fwd strand
			queryPos=((seqsSizes[numSeqs])-(queryPos+memSize)); // pos is on right, add size to go to left, subtract from seq length
			DrawRefAlignmentBlock(queryPos,refPos,(-memSize),numSeqs);
		}
		numMems++;
	}
	fclose(memsFile);
	imageFilename=AppendToBasename(memsFilename,".bmp");
	FinalizeRefAlignmentImage(seqsNames,imageFilename);
	free(imageFilename);
	free(seqsSizes);
	free(seqsNames);
	printf("> Done!\n");
	#ifdef PAUSE_AT_EXIT
	getchar();
	#endif
	exit(0);
}

// Cleans all invalid characters from a FASTA file (leaving only ACGT) and joins multiple sequences if present
void CleanFasta(char *fastafilename){
	FILE *fastafile, *cleanfile;
	char *cleanfilename;
	unsigned int charcount,invalidcharcount,sequencecount,linesize;
	char c;
	printf("> Opening FASTA file <%s> ... ",fastafilename);
	fflush(stdout);
	if((fastafile=fopen(fastafilename,"r"))==NULL){
		printf("\n> ERROR: FASTA file not found\n");
		exit(-1);
	}
	c=fgetc(fastafile);
	if(c!='>'){
		printf("\n> ERROR: Invalid FASTA file\n");
		exit(-1);
	}
	printf("OK\n");
	cleanfilename=AppendToBasename(fastafilename,"-clean.fasta");
	printf("> Creating clean FASTA file <%s> ... ",cleanfilename);
	fflush(stdout);
	if((cleanfile=fopen(cleanfilename,"w"))==NULL){
		printf("\n> ERROR: Can't write clean FASTA file\n");
		exit(-1);
	}
	fprintf(cleanfile,">%s\n",fastafilename); // filename as label
	linesize=0;
	charcount=0;
	invalidcharcount=0;
	sequencecount=0;
	while(c!=EOF){
		if(c=='A' || c=='C' || c=='G' || c=='T'){
			fputc((int)c,cleanfile);
			charcount++;
			linesize++;
			if(linesize==100){ // split sequence by lines of size 100
				fputc('\n',cleanfile);
				linesize=0;
			}
		} else if(c=='a' || c=='c' || c=='g' || c=='t'){
			fputc((int)(c-32),cleanfile);
			charcount++;
			linesize++;
			if(linesize==100){ // split sequence by lines of size 100
				fputc('\n',cleanfile);
				linesize=0;
			}
		} else if(c=='>'){ // new sequence
			while(c!='\n' && c!=EOF) c=fgetc(fastafile); // skip description of sequence
			sequencecount++;
		} else if(c>32 && c<127) invalidcharcount++; // invalid alphanumeric character
		c=fgetc(fastafile);
	}
	if(linesize!=0) fputc('\n',cleanfile); // ending newline
	fclose(fastafile);
	fclose(cleanfile);
	free(cleanfilename);
	printf(" OK\n");
	printf(":: %u total chars",charcount);
	if(invalidcharcount!=0) printf(" (%u non ACGT chars removed)",invalidcharcount);
	if(sequencecount>1) printf(" ; %u sequences merged",sequencecount);
	printf("\n");
	printf("> Done!\n");
	#ifdef PAUSE_AT_EXIT
	getchar();
	#endif
	exit(0);
}

// TODO: add option "-r" to load one single file (or more) and select reference name ("set as reference the sequence containing this string") (create function to move seq to 0-th position of seqs array, load its chars and re-order the remaining seqs)
// TODO: enable option "-v" on normal mode to create image automatically after finding MEMs (directly draw each block inside MEMs finding function) (or if "-v" set, call function to check if any of the input files is a mems file, set new var "memsFile", and pass it to image function)
// TODO: draw gene annotations in image if GFF present (first detect type of all input files, set "memsFile"/"gffFile" variables, add all Fastas to a new array and use it as arg to loadSequences function)
// TODO: if multiple refs exist draw all refs names bellow ref image, one name per row, and vertical lines on each split point extending to corresponding row (box around names or horizontal line above names right extended to the max of name length or split point; grey line when overlapping)
// TODO: remove "baseBwtPos" field from SLCP structure to save memory and benchmark new running times

// TODO: output MUMs (only once in query) and Multi-MEMS (same number in ref and all queries)
int main(int argc, char *argv[]){
	int i, n, numFiles, numSeqsInFirstFile, refFileArgNum, argMatchType, argBothStrands, argNoNs, argMinMemSize;
	char *outFilename;
	printf("[ slaMEM v%s ]\n\n",VERSION);
	if(argc<3){
		printf("Usage:\n");
		printf("\t%s <options> <reference_file> <query_files>\n",argv[0]);
		printf("Options:\n");
		printf("\t-mem\tfind MEMs: any number of occurrences in both ref and query (default)\n");
		printf("\t-mam\tfind MAMs: unique in ref but any number in query\n");
		//printf("\t-mum\tfind MUMs: unique both in ref and query\n");
		printf("\t-l\tminimum match length (default=20)\n");
		printf("\t-o\toutput file name (default=\"*-mems.txt\")\n");
		printf("\t-b\tprocess both strands\n");
		printf("\t-n\tdiscard N's\n");
		printf("Extra:\n");
		printf("\t-v\tgenerate MEMs map image\n");
		//printf("\t-s\tsort MEMs file\n");
		//printf("\t-c\tclean FASTA file\n");
		printf("\n");
		return (-1);
	}
	if( ParseArgument(argc,argv,"S",0) ){ // Sort MEMs
		if(argc!=3){
			printf("Usage: %s -s <mems_file>\n\n",argv[0]);
			return (-1);
		}
		SortMEMsFile(argv[2]);
		return 0;
	}
	if( ParseArgument(argc,argv,"C",0) ){ // Clean FASTA file
		if(argc!=3){
			printf("Usage: %s -c <fasta_file>\n\n",argv[0]);
			return (-1);
		}
		CleanFasta(argv[2]);
		return 0;
	}
	if( ParseArgument(argc,argv,"V",0) ){ // Create MEMs image
		if(argc!=5){
			printf("Usage: %s -v <reference_file> <query_file> <mems_file>\n\n",argv[0]);
			return (-1);
		}
		CreateMemMapImage(argv[2],argv[3],argv[4]);
		return 0;
	}
	argMatchType=0; // MEMs mode
	if( ParseArgument(argc,argv,"MA",0) ) argMatchType=1; // MAMs mode
	argNoNs=ParseArgument(argc,argv,"N",0);
	refFileArgNum=(-1);
	numSeqsInFirstFile=0;
	numFiles=0;
	numSequences=0; // initialize global variable needed by sequence functions
	for(i=1;i<argc;i++){
		if(argv[i][0]=='-'){ // skip arguments for options
			if(argv[i][1]=='L' || argv[i][1]=='l') i++; // skip value of option "-l"
			if(argv[i][1]=='O' || argv[i][1]=='o') i++; // skip value of option "-o"
			continue;
		}
		n=LoadSequencesFromFile(argv[i],((numFiles==0)?1:0),((numFiles==0)?1:0),argNoNs);
		if(n!=0) numFiles++;
		if(numFiles==1){ // reference file
			refFileArgNum=i;
			numSeqsInFirstFile=n;
		}
	}
	if(numFiles==0) exitMessage("No reference or query files provided");
	if(numFiles==1) exitMessage("No query files provided");
	n=(numSequences-numSeqsInFirstFile);
	if(n==0) exitMessage("No valid query sequences found");
	printf("> %d reference%s and %d quer%s successfully loaded\n",numSeqsInFirstFile,((numSeqsInFirstFile==1)?"":"s"),n,((n==1)?"y":"ies"));
	argBothStrands=ParseArgument(argc,argv,"B",0);
	argMinMemSize=ParseArgument(argc,argv,"L",1);
	if(argMinMemSize==(-1)) argMinMemSize=20; // default minimum MEM length is 20
	n=ParseArgument(argc,argv,"O",2);
	if(n==(-1)) outFilename=AppendToBasename(argv[refFileArgNum],"-mems.txt"); // default output base filename is the ref filename
	else outFilename=argv[n];
	GetMatches(numSeqsInFirstFile,numSequences,argMatchType,argMinMemSize,argBothStrands,outFilename);
	if(n==(-1)) free(outFilename);
	DeleteAllSequences();
	printf("> Done!\n");
	#ifdef PAUSE_AT_EXIT
	getchar();
	#endif
	return 0;
}
