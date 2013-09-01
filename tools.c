#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

void exitMessage(char *msg){
	printf("> ERROR: %s\n",msg);
	getchar();
	exit(-1);
}

// Checks and parses program arguments sent using the command line parameters
// NOTE: if parse = 0, returns 1 if the argument exists and 0 otherwise
// NOTE: if parse = 1, returns its value if the argument exists and -1 otherwise
// NOTE: if parse = 2, returns its index in the arglist array if the argument exists and -1 otherwise
int ParseArgument(int numargs, char** arglist, char *optionchars, int parse){
	int i;
	char uc[2], lc[2];
	for(i=0;i<2;i++){ // only process the first two chars of the option name
		if(optionchars[i]=='\0'){ // no char
			uc[i] = '\0';
			lc[i] = '\0';
			continue;
		}
		if( ( ((int)optionchars[i]) >= 65 ) && ( ((int)optionchars[i]) <= 90 ) ){ // uppercase
			uc[i] = optionchars[i];
			lc[i] = (char)(optionchars[i]+32);
		} else { // lowercase
			uc[i] = (char)(optionchars[i]-32);
			lc[i] = optionchars[i];
		}
	}
	for( i = 1 ; i < numargs ; i++ ){
		if( arglist[i][0]=='-' &&
				( arglist[i][1]==lc[0] || arglist[i][1]==uc[0] ) &&
				( arglist[i][2]==lc[1] || arglist[i][2]==uc[1] ) ){
			if(parse){
				if(i==(numargs-1)) return (-1); // if it's the last argument, it has nothing in front
				if(parse==1) return atoi(arglist[(i+1)]); // parse argument in front
				else return (i+1);
			}
			return 1; // if we only want to check if it exits, report 1
		}
	}
	if(parse) return (-1); // we wanted a value, but the argument was not found
	return 0; // argument not present
}

// joins the basename of a name of a file with another string
char* AppendToBasename(char *filename, char *extra){
	char *resultfilename;
	int i, n;
	n=(int)strlen(filename);
	for(i=(n-1);i>0;i--){
		if(filename[i]=='.') break;
	}
	if(i==0) i=n;
	n=(int)strlen(extra);
	resultfilename=(char *)calloc((i+n+1),sizeof(char));
	strncpy(resultfilename,filename,i);
	resultfilename[i]='\0';
	strcat(resultfilename,extra);
	return resultfilename;
}

// Returns a new string with all special characters removed or converted
// NOTE: mode=0 removes all non-alphanumeric chars
// NOTE: mode=1 converts all non-alphanumeric chars to underscores
// NOTE: mode=2 keeps all chars except spaces, tabs, newlines and non-printable chars
char *NormalizeSeqName(char *name, int mode){
	char c, *newname;
	int i, n;
	n=0;
	while(name[n]!='\0') n++;
	newname=(char *)malloc((n+1)*sizeof(char));
	i=0;
	n=0;
	while((c=name[n++])!='\0'){
		if( (c>='0' && c<='9') || (c>='A' && c<='Z') || (c>='a' && c<='z') ) newname[i++]=c; // alpha-numeric chars
		else if(mode!=0) { // convert all other chars to underscore
			if(mode==1){ if( i==0 || newname[(i-1)]!='_' ) newname[i++]='_'; } // prevent writing of two followed underscores
			else if( c>='!' && c<='~' ) newname[i++]=c; // non-blank chars
		}
	}
	newname[i]='\0';
	return newname;
}

// prints an integer value in a formated way
void PrintNumber( int number ){
	long long num, denom, quot;
	//div_t divres;
	if( number < 1000 ){
		printf( "%d" , (int)number );
		return;
	}
	denom = 1;
	num = number;
	while( num >= 1000 ){
		//divres = div( num , 1000 );
		//num = divres.quot;
		num /= 1000;
		denom *= 1000;
	}
	printf( "%d," , (int)num );
	num = ( number - num * denom );
	denom /= 1000;
	while( denom > 1 ){
		//divres = div( num , denom );
		//printf( "%03d." , divres.quot );
		//num = divres.rem;
		quot = ( num / denom );
		printf( "%03d," , (int)quot );
		num -= ( quot * denom );
		denom /= 1000;
	}
	printf( "%03d" , (int)num );
}

// prints a space value in a formated way
void PrintSpace( int spaceval ){
	if( (spaceval/1000) == 0 ){ // less than 1KB
		printf( "%d B" , spaceval );
		return;
	}
	spaceval /= 1000;
	if( (spaceval/1000) == 0 ){ // less than 1MB
		printf( "%d KB" , spaceval );
		return;
	}
	spaceval /= 1000;
	if( (spaceval/1000) == 0 ){ // less than 1GB
		printf( "%d MB" , spaceval );
		return;
	}
	printf( "%.1f GB" , ((float)spaceval)/((float)1000) );
}

// prints a time value in a formated way
void PrintTime( double timeval ){
	int time, days, hours, minutes, seconds;
	div_t divres;
	time = (int) timeval;
	if( time == 0 ){
		printf( "0s " );
		return;
	}
	divres = div( time , 60 );
	seconds = divres.rem;
	time = divres.quot;
	divres = div( time , 60 );
	minutes = divres.rem;
	time = divres.quot;
	divres = div( time , 24 );
	hours = divres.rem;
	days = divres.quot;
	if( days != 0 ) printf( "%dD " , days );
	if( hours != 0 ) printf( "%dH " , hours );
	if( minutes != 0 ) printf( "%dm " , minutes );
	if( seconds != 0 ) printf( "%ds " , seconds );
}

// Shows and updates a progress bar (in the current line or in the line above the current one)
void PrintProgressBar(double percentage, int lineabove){
	static int prevbarpos = -1;
	//const char symbols[] = "\x20\xB0\xB1\xB2\xDB\x2D";
	const char symbols[] = "\x2D\x2D\x2D\x2D\x61\x2D";
	const int barsize = 50;
	const int numchars = 5;
	const int leftspaces = 2;
	int barpos, charpos, i;
	if(percentage<0.0) percentage=0.0;
	if(percentage>100.0) percentage=100.0;
	if(percentage==0.0 && prevbarpos==-1){ // initialize progress bar
		for(i=0;i<(leftspaces-1);i++) putchar(' ');
		printf("[");
		for(i=0;i<barsize;i++) putchar(symbols[5]); // fill with gap symbols
		printf("]\n");
		prevbarpos = 0;
		return;
	}
	barpos=(int)floor(percentage*((double)barsize)/100.0); // which pos
	charpos=(int)floor(percentage*(double)(barsize*numchars)/100.0)-(barpos*numchars); // which char (percentage mod numchars)
	if(barpos==barsize){ // even if we are at 100%, the last bar pos is at (barsize-1)
		barpos=(barsize-1);
		charpos=(numchars-1);
	}
	printf("\x1B(0"); // printf("\x1B[11m"); // alternate font
	if(lineabove) printf("\x1B[99D\x1B[K"); // move cursor to beginning of line and delete until the end of the line
	if(lineabove) printf("\x1B[1F"); // previous line (1A = line up)
	i=prevbarpos; // previous filled position in bar
	printf("\x1B[%dC",(leftspaces+i)); // move forward (#G = go to column #)
	while((i++)<barpos) putchar(symbols[4]); // fill incomplete spaces behind if needed
	putchar(symbols[charpos]);
	printf("\x1B[%dD",(leftspaces+i)); // move back (1G = go to column 1)
	//if(lineabove) printf("\x1B[1E"); // next line (1B = line down)
	printf("\x1B(B"); // printf("\x1B[10m"); // back to normal font
	prevbarpos=barpos; // save last bar position
	//if(lineabove) printf("\x1B[99D\x1B[K"); // move cursor to beginning of line and delete until the end of the line
}

// NOTE: returns (-1) if it is not a valid number
int ParseNumber(char *numberstring){
	int i, number, decimals;
	char c;
	number=0;
	decimals=0;
	i=0;
	while((c=numberstring[i])!='\0'){
		if(c>='0' && c<='9'){
			number=(number*10+((int)c-48));
			decimals=(decimals*10);
		} else if(c=='K'){
			number=(number*1000);
			break;
		} else if(c=='M'){
			number=(number*1000000);
			break;
		} else if(c=='G'){
			number=(number*1000000000);
			break;
		} else if(c==',' || c=='.') decimals=1;
		else return (-1);
		i++;
	}
	if(decimals!=0) number=(number/decimals);
	return number;
}
