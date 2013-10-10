/* ================================================================= *
 *  graphics.c : Graphical visualizations                            *
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
#include <math.h>
#include <stdint.h>
#include <limits.h>
#include "bitmap.h"

#define CHAR_WIDTH 5
#define CHAR_HEIGHT 6
#define CHAR_SPACE 2
#define ALPHABET_SIZE 95
const unsigned int alphabet[] = { // least significant bit is the upper left pixel, and so on
	0x00000000, // (space)
	0x08021084, // !
	0x0000014A, // "
	0x15F57D40, // #
	0x1F4717C4, // $
	0x22221111, // %
	0x2C9A88A2, // &
	0x00000084, // '
	0x1821084C, // (
	0x0C842106, // )
	0x2AEFBAA0, // *
	0x084F9080, // +
	0x04400000, // ,
	0x00007C00, // -
	0x08000000, // .
	0x02221110, // /
	0x1D3AD72E, // 0
	0x3E4214C4, // 1
	0x3E17420F, // 2
	0x1F083A1F, // 3
	0x108FA988, // 4
	0x1F083C3F, // 5
	0x1D18BC3E, // 6
	0x0222221F, // 7
	0x1D18BA2E, // 8
	0x1F0F462E, // 9
	0x00400080, // :
	0x04401000, // ;
	0x30609B00, // <
	0x01F07C00, // =
	0x06C83060, // >
	0x0802322E, // ?
	0x05DAF62E, // @
	0x23153944, // A
	0x1F18BE2F, // B
	0x3C10843E, // C
	0x1F18C62F, // D
	0x3E109C3F, // E
	0x02109C3F, // F
	0x3D1C843E, // G
	0x2318FE31, // H
	0x3E42109F, // I
	0x0C94A11E, // J
	0x22928CB9, // K
	0x3E108421, // L
	0x2318D771, // M
	0x239AD671, // N
	0x1D18C62E, // O
	0x0217C62F, // P
	0x2C9AC62E, // Q
	0x2297C62F, // R
	0x1F08383E, // S
	0x0842109F, // T
	0x1D18C631, // U
	0x08A54631, // V
	0x23BAC631, // W
	0x23151151, // X
	0x08422A31, // Y
	0x3E22111F, // Z
	0x1C21084E, // [
	0x20821041, // (backslash)
	0x1C84210E, // ]
	0x00004544, // ^
	0x3E000000, // _
	0x00000082, // `
	0x1CA721C0, // a
	0x1CA53842, // b
	0x1C2109C0, // c
	0x1CA53908, // d
	0x1C2729C0, // e
	0x0842388C, // f
	0x1C8729C0, // g
	0x14A53842, // h
	0x08421004, // i
	0x0C421004, // j
	0x14A32842, // k
	0x08421084, // l
	0x2B5AD560, // m
	0x14A528C0, // n
	0x1CA529C0, // o
	0x042729C0, // p
	0x108729C0, // q
	0x04253840, // r
	0x1C8709C0, // s
	0x0C4211C4, // t
	0x1CA52940, // u
	0x08452940, // v
	0x14AAD6A0, // w
	0x14A22940, // x
	0x1C872940, // y
	0x1C2221C0, // z
	0x18210C4C, // {
	0x08421084, // |
	0x0C846106, // }
	0x009AC800, // ~
	0x3FFFFFFF // invalid
};

uint8_t color_white, color_black, color_grey, color_red, color_green, color_blue;

double pos_per_pixel;
uint8_t **ref_pixels_colors;
int h_margin, num_seqs, seqs_height, *seqs_width, *seqs_start_y, **seqs_pixels_block_size;

void InitColors(){
	color_white=getColorFromPalette(255,255,255);
	color_black=getColorFromPalette(0,0,0);
	//color_grey=getColorFromPalette(128,128,128);
	color_grey=getColorFromPalette(222,222,222);
	color_red=getColorFromPalette(255,0,0);
	color_green=getColorFromPalette(0,255,0);
	color_blue=getColorFromPalette(0,0,255);
}

void DrawChar(char c, int x, int y, uint8_t color){
	unsigned int charpixels;
	unsigned int i, j;
	c-=32; // convert char code to position in alphabet array
	if(c>ALPHABET_SIZE) c=ALPHABET_SIZE; // invalid char
	charpixels=alphabet[(int)c];
	j=CHAR_HEIGHT;
	while(j--){ // all rows
		i=CHAR_WIDTH;
		while(i--){ // all columns of each row
			if(charpixels & 1u) drawPoint(x,y,color);
			charpixels>>=1;
			x++;
		}
		x-=CHAR_WIDTH;
		y++;
	} // draw from left to right, and top to bottom
}

void DrawString(char *s, int maxchars, int x, int y, uint8_t color){
	while((*s)!='\0' && maxchars!=0){
		DrawChar((*s),x,y,color);
		x+=(CHAR_WIDTH+CHAR_SPACE);
		s++;
		maxchars--;
	}
}

int DigitsCount(int number){
	int count;
	count=1;
	while(number>=10){
		number /= 10;
		count++;
	}
	return count;
}

int CharsCount(char *string){
	int count = 0;
	while((*string)){
		count++;
		string++;
	}
	return count;
}

void DrawNumber(int number, int x, int y, uint8_t color){
	int n,d,m;
	n=DigitsCount(number); // print number centered
	n=( n*CHAR_WIDTH + (n-1)*CHAR_SPACE )/2;
	if(x>n) x-=n;
	else x=1;
	n=number;
	d=1;
	while(n>=10){
		n /= 10;
		d *= 10;
	}
	n=number;
	while(d>0){
		m = (n/d); // from 0 to 9
		DrawChar((char)(48+m),x,y,color); // draw char corresponding to this digit
		x+=(CHAR_WIDTH+CHAR_SPACE);
		n -= (m*d);
		d /= 10;
	}
}

void DrawHorizontalLine(int x, int y, int length, uint8_t color){
	while(length--) drawPoint(x++,y,color);
}

void DrawVerticalLine(int x, int y, int length, uint8_t color){
	while(length--) drawPoint(x,y++,color);
}

void DrawEmptyRectangle(int x, int y, int xsize, int ysize, uint8_t color){
	DrawHorizontalLine(x,y,xsize,color);
	DrawHorizontalLine(x,y+ysize-1,xsize,color);
	DrawVerticalLine(x,y,ysize,color);
	DrawVerticalLine(x+xsize-1,y,ysize,color);
}

void DrawFilledRectangle(int x, int y, int xsize, int ysize, uint8_t color){
	while(ysize--) DrawHorizontalLine(x,y++,xsize,color);
}

void DrawArrowTip(int x, int y, int xstep, int ysize, int dir, uint8_t backcolor){
	int xsize;
	xsize=(((ysize-1)/2)*xstep);
	if(dir>=0){ // point to the right
		dir=(+1);
		x=(x-xsize+1);
	} else dir=(-1); // point to the left
	while(xsize>0 && ysize>0){
		DrawHorizontalLine(x,y,xsize,backcolor);
		DrawHorizontalLine(x,(y+ysize-1),xsize,backcolor);
		if(dir==+1) x+=xstep;
		xsize-=xstep;
		y++;
		ysize-=2;
	}
}

void InitializeRefAlignmentImage(int *seqsizes, int numseqs){
	int img_width, img_height;
	int v_margin, header_height, seqs_space, max_seq_size;
	int max_seq_width, num_digits, mark_pos, mark_number;
	int i, n;
	double step;
	v_margin=(2*CHAR_HEIGHT);
	h_margin=(2*CHAR_WIDTH);
	seqs_space=(2*CHAR_HEIGHT);
	seqs_height=(5*CHAR_HEIGHT);
	header_height=(2*CHAR_HEIGHT+2*CHAR_SPACE);
	img_width=1024;
	img_height=( 2*v_margin + header_height + numseqs*seqs_height + (numseqs-1)*seqs_space ); // top/bottom margins, position marks, all seqs, space between seqs
	initializeBitmap(img_width,img_height,0);
	num_seqs=numseqs;
	max_seq_size=seqsizes[0];
	for(i=1;i<num_seqs;i++) if(seqsizes[i]>max_seq_size) max_seq_size=seqsizes[i]; // get size of longest seq
	num_digits=DigitsCount(max_seq_size);
	pos_per_pixel=( (double)(max_seq_size) / (double)(img_width-2*h_margin-(num_digits*(CHAR_WIDTH+CHAR_SPACE)/2)) );
	seqs_width=(int *)malloc(num_seqs*sizeof(int));
	seqs_start_y=(int *)malloc(num_seqs*sizeof(int));
	seqs_pixels_block_size=(int **)malloc(num_seqs*sizeof(int *));
	for(i=0;i<num_seqs;i++){
		seqs_width[i]=(int)ceil(((double)seqsizes[i])/pos_per_pixel);
		seqs_start_y[i]=( v_margin + header_height + i*(seqs_height+seqs_space) );
		if(i!=0) seqs_pixels_block_size[i]=(int *)calloc(seqs_width[i],sizeof(int));
	}
	InitColors();
	for(i=1;i<num_seqs;i++){
		DrawFilledRectangle(h_margin,seqs_start_y[i],seqs_width[i],seqs_height,color_grey);
		DrawEmptyRectangle((h_margin-1),(seqs_start_y[i]-1),(seqs_width[i]+2),(seqs_height+2),color_black);
	}
	ref_pixels_colors=(uint8_t **)malloc(2*sizeof(uint8_t *)); // colors of all pixels of the ref in both strands
	for(i=0;i<2;i++) ref_pixels_colors[i]=(uint8_t *)malloc(seqs_width[0]*sizeof(uint8_t));
	step=(1.0)/((double)(seqs_width[0])); // global progress
	for(i=0;i<seqs_width[0];i++){
		n=(int)floor((128+256)*(i*step)); // forward strand in hot colors
		if(n<128){ // magenta (255,0,128)
			ref_pixels_colors[0][i]=getColorFromPalette(255,0,(128-n));
		} else { // if(n<(128+256)) // red (255,0,0) , orange (255,128,0)
			n-=128;
			ref_pixels_colors[0][i]=getColorFromPalette(255,n,0);
		} // yellow (255,255,0)
		n=(int)floor((2*256+128)*(i*step)); // reverse strand in cold colors
		if(n<256){ // green (0,255,0)
			ref_pixels_colors[1][i]=getColorFromPalette(0,255,n);
		} else if(n<(2*256)){ // cyan (0,255,255)
			n-=256;
			ref_pixels_colors[1][i]=getColorFromPalette(0,(255-n),255);
		} else { // if(n<(2*256+128)) // blue (0,0,255)
			n-=(2*256);
			ref_pixels_colors[1][i]=getColorFromPalette(n,0,255);
		} // purple (128,0,255)
		DrawVerticalLine((h_margin+i),(seqs_start_y[0]),((seqs_height/2)-1),ref_pixels_colors[0][i]);
		DrawVerticalLine((h_margin+i),(seqs_start_y[0]+(seqs_height/2)+1),((seqs_height/2)-1),ref_pixels_colors[1][i]);
	}
	DrawArrowTip((h_margin+seqs_width[0]-1),(seqs_start_y[0]),3,((seqs_height/2)-1),(+1),color_white);
	DrawArrowTip((h_margin),(seqs_start_y[0]+(seqs_height/2)+1),3,((seqs_height/2)-1),(-1),color_white);
	DrawNumber(1,(h_margin),(v_margin),color_black); // print left pos number
	DrawVerticalLine((h_margin),(v_margin+CHAR_HEIGHT+CHAR_SPACE),CHAR_HEIGHT,color_black);
	max_seq_width=(int)ceil(((double)max_seq_size)/pos_per_pixel);
	DrawNumber(max_seq_size,(h_margin+max_seq_width-1),(v_margin),color_black); // print right pos number
	DrawVerticalLine((h_margin+max_seq_width-1),(v_margin+CHAR_HEIGHT+CHAR_SPACE),CHAR_HEIGHT,color_black);
	i=(num_digits*CHAR_WIDTH+(num_digits-1)*CHAR_SPACE); // space occupied by each number
	n=(max_seq_width-i); // space available to print numbers without the numbers already at both ends
	i+=CHAR_WIDTH; // each number needs a space after it
	i=(n/i); // how many numbers can be printed in that space
	n=(max_seq_size/(i+1)); // shortest interval between consecutive marks so that numbers do not overlap
	i=1; // interval between marks
	while(i<n){ // find next larger interval size that is a multiple of 2, 5 or 10
		i=2*i; // x2
		if(i>=n) break;
		i=(i/2)*5; // x5
		if(i>=n) break;
		i=2*i; // x10
	}
	n=(h_margin+max_seq_width-(num_digits*(CHAR_WIDTH+CHAR_SPACE))); // last valid position before overlapping with left number
	mark_number=i;
	while((mark_pos=(int)floor(((double)mark_number)/pos_per_pixel))<n){ // print marks
		DrawNumber(mark_number,(h_margin+mark_pos),(v_margin),color_black);
		DrawVerticalLine((h_margin+mark_pos),(v_margin+CHAR_HEIGHT+CHAR_SPACE),CHAR_HEIGHT,color_black);
		mark_number+=i;
	}
	DrawHorizontalLine((h_margin),(v_margin+CHAR_HEIGHT+CHAR_SPACE+(CHAR_HEIGHT/2)),(max_seq_width),color_black);
}

// NOTE: negative size means the block is in the reverse strand, but the position always refer to the forward strand
void DrawRefAlignmentBlock(int seqpos, int refpos, int size, int seqid){
	int endpos, strand;
	if(size<0){
		strand=1;
		size=(-size);
	} else strand=0;
	endpos=(int)floor(((double)(seqpos+size-1))/pos_per_pixel);
	seqpos=(int)floor(((double)(seqpos))/pos_per_pixel);
	refpos=(int)floor(((double)(refpos))/pos_per_pixel);
	for(;seqpos<=endpos;seqpos++){
		if((seqs_pixels_block_size[seqid][seqpos])<size){
			DrawVerticalLine((h_margin+seqpos),seqs_start_y[seqid],seqs_height,ref_pixels_colors[strand][refpos]);
			seqs_pixels_block_size[seqid][seqpos]=size;
		}
		refpos++;
	}
}

void FinalizeRefAlignmentImage(char **seqnames, char *imagefilename){
	int i, n, m;
	n=CharsCount(seqnames[0]);
	m=(seqs_width[0]-2*CHAR_WIDTH)/(CHAR_WIDTH+CHAR_SPACE); // max chars to print inside seq box
	if(n>m) n=m;
	DrawFilledRectangle((h_margin+(seqs_width[0]/2)-((n*CHAR_WIDTH+(n-1)*CHAR_SPACE)/2)-1),(seqs_start_y[0]+(seqs_height/2)-(CHAR_HEIGHT/2)-1),(n*(CHAR_WIDTH+CHAR_SPACE)-CHAR_SPACE+2),(CHAR_HEIGHT+2),color_white);
	DrawString(seqnames[0],n,(h_margin+(seqs_width[0]/2)-((n*CHAR_WIDTH+(n-1)*CHAR_SPACE)/2)),(seqs_start_y[0]+(seqs_height/2)-(CHAR_HEIGHT/2)),color_black);
	for(i=1;i<num_seqs;i++){
		n=CharsCount(seqnames[i]);
		m=(seqs_width[i]-2*CHAR_WIDTH)/(CHAR_WIDTH+CHAR_SPACE); // max chars to print inside seq box
		if(n>m) n=m;
		DrawFilledRectangle((h_margin+CHAR_WIDTH-1),(seqs_start_y[i]+(seqs_height/2)-(CHAR_HEIGHT/2)-1),(n*(CHAR_WIDTH+CHAR_SPACE)-CHAR_SPACE+2),(CHAR_HEIGHT+2),color_white);
		DrawString(seqnames[i],n,(h_margin+CHAR_WIDTH),(seqs_start_y[i]+(seqs_height/2)-(CHAR_HEIGHT/2)),color_black);
		free(seqs_pixels_block_size[i]);
	}
	free(seqs_pixels_block_size);
	free(seqs_width);
	free(seqs_start_y);
	for(i=0;i<2;i++) free(ref_pixels_colors[i]);
	free(ref_pixels_colors);
	printf("> Saving image to <%s> ... ",imagefilename);
	fflush(stdout);
	if(saveBitmap(imagefilename)==0){
		printf("\n> ERROR: Cannot write file\n");
		exit(0);
	}
	printf("OK\n");
}

/*
void CreateGlobalMap(unsigned int startpos, unsigned int endpos, unsigned int numpos, unsigned int refsize, char *reflabel, char *readsfilename){
		printf("> Creating global reference coverage plot ... ");
	fflush(stdout);
	numfilledpos=0;
	numrows=0;
	for(i=0;i<numpos;i++){
		if(poscount[i]>numrows) numrows=poscount[i]; // find maximum number of rows
		if(poscount[i]!=0) numfilledpos++; // count number of filled positions
	}
	if(numrows>50) numrows=50; // limit maximum number of rows to 50
	rowcount=(unsigned int *)calloc(numrows,sizeof(unsigned int));
	rowheight=4;
	margindim=(2*rowheight); // side and top margins with size 2 times the row thickness
	refheight=(2*rowheight); // height of the reference in pixels
	imgwidth=1024; // fixed image width
	imgheight=( 2*margindim + 2*CHAR_HEIGHT + refheight + 2*numrows*rowheight ); // top margins, position markers, reference, and all rows with a spacing before
	initializeBitmap((int)imgwidth,(int)imgheight,2);
	k=DigitsCount(endpos+1); // make room to print half of a position number on both ends
	refwidth=( imgwidth - 2*( margindim + k*(CHAR_WIDTH+CHAR_SPACE)/2 ) ); // width of the reference in pixels
	posperpixel=(((double)numpos)/((double)refwidth)); // number of ref positions represented by each pixel
	startx=( margindim + (k*(CHAR_WIDTH+CHAR_SPACE))/2 ); // where to start drawing from
	starty=margindim;
	color=BLACK;
	x=startx;
	y=starty;
	DrawNumber((startpos+1),(x-1),y,color); // print left pos number
	DrawVerticalLine((x-1),(y+CHAR_HEIGHT+CHAR_HEIGHT/2),(CHAR_HEIGHT/2),color);
	x=(imgwidth-startx);
	DrawNumber((endpos+1),x,y,color); // print right pos number
	DrawVerticalLine(x,(y+CHAR_HEIGHT+CHAR_HEIGHT/2),(CHAR_HEIGHT/2),color);
	j=(unsigned int)floor(posperpixel*(double)(refwidth-k*(CHAR_WIDTH+CHAR_SPACE))); // last middle position that can have a marked number
	k=( k*CHAR_WIDTH + (k-1)*CHAR_SPACE ); // space occupied by each number
	n=( refwidth - k ); // space available to print numbers without the numbers already at both ends
	k+=CHAR_WIDTH; // each number needs a space after it
	k=(n/k); // how many numbers can be printed in that space
	n=(endpos-startpos+1); // total number of positions
	n=(n/(k+1)); // how many positions can go between each mark (minimum)
	k=1; // how many positions between marks
	while(k<n){
		k=2*k; // x2
		if(k>=n) break;
		k=(k/2)*5; // x5
		if(k>=n) break;
		k=2*k; // x10
	}
	n=(unsigned int)floor(((double)k)/posperpixel); // number of pixels between marks
	x=(startx+n);
	y=starty;
	i=k;
	while(i<j){ // print marks
		DrawNumber(i,x,y,color);
		DrawVerticalLine(x,(y+CHAR_HEIGHT+CHAR_HEIGHT/2),(CHAR_HEIGHT/2),color);
		x+=n;
		i+=k;
	}
}
*/
