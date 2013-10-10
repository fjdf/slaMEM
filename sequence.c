/* ================================================================= *
 *  sequence.c : FASTA files processing                              *
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
#include "sequence.h"

char *charsTable = NULL;

Sequence *AddNewSequence(){
	Sequence *newSeq;
	numSequences++;
	if(numSequences==1) allSequences=(Sequence **)malloc(numSequences*sizeof(Sequence *));
	else allSequences=(Sequence **)realloc(allSequences,numSequences*sizeof(Sequence *));
	newSeq=(Sequence *)calloc(1,sizeof(Sequence));
	allSequences[(numSequences-1)]=newSeq;
	return newSeq;
}

void DeleteAllSequences(){
	Sequence *seq;
	int i;
	for(i=0;i<numSequences;i++){
		seq=allSequences[i];
		if((seq->name)!=NULL) free(seq->name);
		if((seq->chars)!=NULL) free(seq->chars);
		//if((seq->sourcefilename)!=NULL) free(seq->sourcefilename);
		if((seq->sourcefile)!=NULL) fclose(seq->sourcefile);
		free(seq);
	}
	free(allSequences);
	allSequences=NULL;
	if(charsTable!=NULL) free(charsTable);
}

void InitCharsTable(int allowNs){
	int i;
	if(charsTable!=NULL) return;
	charsTable=(char *)calloc(256,sizeof(char)); // default code is 0
	if(allowNs){
		for(i=65;i<=90;i++){ // letters from 'A' to 'Z'
			charsTable[i]='N';
			charsTable[i+32]='N';
		}
	}
	charsTable['A']='A';
	charsTable['C']='C';
	charsTable['G']='G';
	charsTable['T']='T';
	charsTable['a']='A';
	charsTable['c']='C';
	charsTable['g']='G';
	charsTable['t']='T';
	charsTable['>']=(char)EOF;
	charsTable[(unsigned char)EOF]=(char)EOF;
}

// NOTE: returns the number of sequences inside the file
// NOTE: numSequences must be set to 0 before the first invocation of this function
int LoadSequencesFromFile(char *inputfilename, int loadchars, int acgtonly){
	FILE *file;
	char c, *seqchars;
	int k,numseqs,desclen,seqlen,maxseqlen;
	long int filestart,fileend,filesize;
	fpos_t startpos;
	Sequence *seq;
	printf("> Loading sequences from file <%s> ... ",inputfilename);
	if((file=fopen(inputfilename,"r"))==NULL){
		printf("\n> WARNING: Sequence file not found\n");
		return 0;
	}
	filesize=0;
	filestart=ftell(file);
	fseek(file,0L,SEEK_END);
	fileend=ftell(file);
	filesize=(fileend-filestart);
	rewind(file);
	printf("(%ld bytes)\n",filesize);
	c=fgetc(file);
	while(c!=EOF && c!='>') c=fgetc(file);
	if(c==EOF){
		printf("> WARNING: No sequences in file\n");
		return 0;
	}
	InitCharsTable(!acgtonly);
	numseqs=0; // number of sequences inside this file only
	while(1){
		while(c!=EOF && c!='>') c=fgetc(file);
		if(c==EOF) break;
		fgetpos(file,&startpos); // save file position
		printf("# %02d [",(numSequences+1));
		desclen=0;
		while((c=fgetc(file))!=EOF && c!='\n' && c!='\r'){
			if(desclen<40) putchar(c);
			desclen++;
		}
		for(k=desclen;k<40;k++) putchar(' ');
		printf("] ");
		fflush(stdout);
		seqchars=NULL;
		seqlen=0;
		if(loadchars){
			maxseqlen=0;
			while((c=fgetc(file))!='>' && c!=EOF){
				c=charsTable[(unsigned char)c]; // normalize char
				if(c!=0){
					if(seqlen==maxseqlen){
						maxseqlen+=(1<<20); // allocate space in 1MB steps
						seqchars=(char *)realloc(seqchars,maxseqlen*sizeof(char));
					}
					seqchars[seqlen]=c;
					seqlen++;
				}
			}
			if(seqlen!=0){
				seqchars=(char *)realloc(seqchars,(seqlen+1)*sizeof(char)); // shorten allocated space to fit real seq size
				seqchars[seqlen]='\0';
			}
		} else {
			while((c=fgetc(file))!='>' && c!=EOF){
				c=charsTable[(unsigned char)c];
				if(c!=0) seqlen++;
			}
		}
		/*
		while((c=fgetc(file))!=EOF && c!='>'){
			if(c>='a' && c<='z') c-=32;
			if(c>='A' && c<='Z') seqlen++;
		}
		*/
		if(seqlen==0){
			printf("EMPTY\n");
			continue;
		}
		printf("(%d bp) ",seqlen);
		fflush(stdout);
		numseqs++;
		seq=AddNewSequence(); // new sequence ; sets numSequences++
		seq->size=seqlen;
		seq->order=numSequences;
		seq->name=(char *)malloc((desclen+1)*sizeof(char));
		fsetpos(file,&startpos); // restore file position
		k=0;
		while((c=fgetc(file))!=EOF && c!='\n' && c!='\r') (seq->name)[k++]=c;
		(seq->name)[k]='\0';
		if(loadchars){
			(seq->chars)=seqchars;
			if(numseqs==1) loadchars=0; // only load chars for 1st seq inside file
			/*
			(seq->chars)=(char *)malloc((seqlen+1)*sizeof(char));
			k=0;
			while((c=fgetc(file))!=EOF && c!='>'){
				if(c>='a' && c<='z') c-=32;
				if(c>='A' && c<='Z'){
					if(c=='A' || c=='C' || c=='G' || c=='T') (seq->chars)[k++]=c;
					else (seq->chars)[k++]='N';
				}
			}
			(seq->chars)[k]='\0';
			*/
		} else {
			fgetpos(file,&startpos);
			seq->sourcefilepos=startpos;
			/*
			k=0;
			while(inputfilename[k]!='\0') k++;
			seq->sourcefilename=(char *)malloc((k+1)*sizeof(char));
			k=0;
			while((c=inputfilename[k])!='\0') (seq->sourcefilename)[k++]=c;
			(seq->sourcefilename)[k]='\0';
			*/
			seq->sourcefile=fopen(inputfilename,"r");
			fsetpos(seq->sourcefile,&startpos);
			while(c!=EOF && c!='>') c=fgetc(file);
		}
		printf("OK\n");
		fflush(stdout);
	}
	fclose(file);
	return numseqs;
}

void LoadSequenceChars(Sequence *seq){
	FILE *file;
	int i, seqsize;
	char c;
	//if((seq->sourcefilename)==NULL) return;
	//if((file=fopen((seq->sourcefilename),"r"))==NULL) return;
	if((seq->chars)!=NULL || (seq->sourcefile)==NULL) return;
	seqsize=(seq->size);
	(seq->chars)=(char *)malloc((seqsize+1)*sizeof(char));
	file=(seq->sourcefile);
	fsetpos(file,&(seq->sourcefilepos));
	/*
	i=0;
	while((c=fgetc(file))!=EOF && c!='>'){
		if(c>='a' && c<='z') c-=32;
		if(c>='A' && c<='Z' && i<seqsize){
			if(c=='A' || c=='C' || c=='G' || c=='T') (seq->chars)[i++]=c;
			else (seq->chars)[i++]='N';
		}
	}
	*/
	i=0;
	while((c=fgetc(file))!='>' && c!=EOF){
		c=charsTable[(unsigned char)c];
		if(c!=0) (seq->chars)[i++]=c;
	}
	(seq->chars)[i]='\0';
}

void FreeSequenceChars(Sequence *seq){
	if((seq->chars)!=NULL) free(seq->chars);
	seq->chars=NULL;
}

int GetSeqIdFromSeqName(char *seqname){
	char *currentseqname, c, sc;
	int s, k, bestk, bests, i, si;
	bestk=0;
	bests=(-1);
	for(s=0;s<numSequences;s++){
		currentseqname=(allSequences[s]->name);
		k=0;
		//while(seqname[k]!='\0' && seqname[k]==currentseqname[k]) k++;
		i=0;
		si=0;
		c='\0';
		sc='\0';
		while(1){ // get longest match between the two names
			do c=seqname[i++];
			while( c!='\0' && !(c>=48 && c<=57) && !(c>=65 && c<=90) && !(c>=97 && c<=122) ); // alphanumeric chars only
			do sc=currentseqname[si++];
			while( sc!='\0' && !(sc>=48 && sc<=57) && !(sc>=65 && sc<=90) && !(sc>=97 && sc<=122) );
			if(c!=sc) break; // if the query name is a prefix of other incorrect name, the match size is not incremented next
			k++;
			if(c=='\0') break;
		}
		if(k>bestk){
			if(c=='\0' && sc=='\0') return s; // exact name match
			bestk=k; // keep sequence with longest match
			bests=s;
		}
	}
	if(bestk==0) return (-1); // no match
	return bests;
}

char GetNextChar(int seqid){
	char c;
	do {
		c=charsTable[(unsigned char)fgetc(allSequences[seqid]->sourcefile)];
	} while(c==0);
	return c;
	/*
	FILE *file;
	file=(allSequences[seqid]->sourcefile);
	while((c=fgetc(file))!=EOF && c!='>'){
		if(c>='a' && c<='z') c-=32;
		if(c>='A' && c<='Z'){
			if(c=='A' || c=='C' || c=='G' || c=='T') return c;
			else return 'N';
		}
	}
	return (char)EOF;
	*/
}

int GetNextCharCode(int seqid){
	FILE *file;
	char c;
	file=(allSequences[seqid]->sourcefile);
	while((c=fgetc(file))!=EOF && c!='>'){
		if(c>='a' && c<='z') c-=32;
		if(c>='A' && c<='Z'){
			if(c=='A') return 0;
			if(c=='C') return 1;
			if(c=='G') return 2;
			if(c=='T') return 3;
			return 4; // 'N'
		}
		if(c=='-') return 5;
	}
	return (-1);
}

char CharAt(int pos, int seqid){
	Sequence *seq;
	seq=allSequences[seqid];
	if((seq->rotation)!=0){
		pos+=(seq->rotation);
		if(pos>=(seq->size)) pos-=(seq->size);
	}
	return (seq->chars)[pos];
}

// TODO: use quicksort
void SortSequences(int *seqsizes, int *sortedseqs, int numseqs){
	int i, j, k, minsize, minseq;
	for(i=0;i<numseqs;i++) sortedseqs[i]=i;
	for(i=0;i<numseqs;i++){
		k=sortedseqs[i]; // to prevent detection of an already detected minimum
		minseq=i;
		minsize=seqsizes[k];
		for(j=(i+1);j<numseqs;j++){
			k=sortedseqs[j];
			if(seqsizes[k]<minsize){
				minsize=seqsizes[k];
				minseq=j; // seq with min size at pos j of sortedSeqsIds array
			}
		}
		k=sortedseqs[minseq]; // seq with min size at pos k of unsorted seqs array
		sortedseqs[minseq]=sortedseqs[i]; // swap seqs at pos i and pos minseq
		sortedseqs[i]=k;
	}
}

void ReverseComplementSequence(char *text, int textsize){
	int posleft, posright;
	char charleft, charright;
	for(posleft=0,posright=(textsize-1);posleft<=posright;posleft++,posright--){
		charleft=text[posleft];
		charright=text[posright];
		if(charleft=='A') charleft='T';
		else if(charleft=='C') charleft='G';
		else if(charleft=='G') charleft='C';
		else if(charleft=='T') charleft='A';
		if(charright=='A') charright='T';
		else if(charright=='C') charright='G';
		else if(charright=='G') charright='C';
		else if(charright=='T') charright='A';
		text[posleft]=charright;
		text[posright]=charleft;
	}
}
