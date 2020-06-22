/* msMask.h -- masking for ms ascertainment stuff*/

#ifndef MASK_INC
#define MASK_INC

#define LENGTH	1000000 

#include "stdio.h"

// stringWrap Object Definition 
typedef struct{
  int **bitmap;			// where the masking is stored
  int sampleSize;				// number of alleles represented
} mask;



int biggerlist(int nsam, unsigned nmax, char **list );
char **cmatrix(int nsam, int len);
mask *mask_new(int nsam);	
int readMaskFiles(char *fileName, mask **masks, int nsam);	
int readNextMask(FILE *handle,mask *masks, int nsam);
	

#endif
