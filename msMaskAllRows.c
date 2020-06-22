#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "msMask.h"

#define LINEBUF 1000000

int maxsites = 100000 ;
int main(argc,argv)
	int argc;
char *argv[];
{
	int nsam, i, j, howmany,n  ;
	char **list, **cmatrix(), line[LINEBUF+1]  ;
	char *curLine;
	FILE *fopen(), *pfin, *maskFile ;
	double *posit  ;
	int   segsites, count  , nadv, probflag ;
	char dum[20];
	int  segsub( int nsam, int segsites, char **list ) ;
	mask *aMask;
	
	
/* read in first two lines of output  (parameters and seed) */
	pfin = stdin ;
	fgets( line, LINEBUF, pfin);
	//print line
	printf("%s",line);
	sscanf(line," %s  %d %d", dum,  &nsam, &howmany);
	//now have the info to deal with the maskfiles
	aMask = mask_new(nsam);
	
	//open mask file
	maskFile = fopen(argv[1], "r");
	if (maskFile == NULL){
		fprintf(stderr,"Error opening maskfile! ARRRRR!!!!\n");
		exit(1);
	}

	
	
	fgets( line, LINEBUF, pfin);
	//print line
	printf("%s",line);
	if( argc > 1 ) { 
		nadv = atoi( argv[1] ) ; 
	}

	list = cmatrix(nsam,maxsites+1);
	posit = (double *)malloc( maxsites*sizeof( double ) ) ;

	count=0;
	probflag = 0 ;

	while( howmany-count++ ) {
/* read in a sample */
		do {
			fgets( line, LINEBUF, pfin);
			printf("%s",line);
		}while ( line[0] != '/' );
		fgets( line, LINEBUF, pfin);
		printf("%s",line);
		sscanf(line,"  segsites: %d", &segsites );
	//	fgets( line, LINEBUF, pfin);
		if( segsites >= maxsites){
			maxsites = segsites + 10 ;
			posit = (double *)realloc( posit, maxsites*sizeof( double) ) ;
			biggerlist(nsam,maxsites, list) ;
		}
		//read the next mask
		readNextMask(maskFile,aMask, 1);
		if( segsites > 0) {
			//print positions
			fgets( line, LINEBUF, pfin);
			printf("%s",line);
			//need this stupid dummy ptr to advance buffer
			curLine = line;
			sscanf(curLine, "%*s %lf%n",posit,&n);
			curLine += n;
			for( i=1; i<segsites ; i++){
				sscanf(curLine," %lf%n",posit+i,&n);
				curLine += n;
			}
			for( i=0; i<nsam;i++) fscanf(pfin," %s", list[i] );

			//now swap in the masks
			for(i=0; i<nsam ;i++){
				for(j=0;j<segsites;j++){
					if(aMask->bitmap[0][(int)(posit[j]*LENGTH)] == 1){
						putchar(list[i][j]);
					}
					else {
						putchar('N');
					}
				}
				printf("\n");
			}
		}

	}
	fclose(maskFile);
	return(1);
}


/* allocates space for gametes (character strings) */
char **cmatrix(int nsam,int len){
	int i;
	char **m;

	if( ! ( m = (char **) malloc( (unsigned)( nsam*sizeof( char* )) ) ) )
		perror("alloc error in cmatrix") ;
	for( i=0; i<nsam; i++) {
		if( ! ( m[i] = (char *) malloc( (unsigned) (len*sizeof( char )) )))
			perror("alloc error in cmatric. 2");
	}
	return( m );
}


int biggerlist(int nsam, unsigned nmax, char **list ){
	int i;

	maxsites = nmax  ;
	for( i=0; i<nsam; i++){
		list[i] = (char *)realloc( list[i],maxsites*sizeof(char) ) ;
		if( list[i] == NULL ) perror( "realloc error. bigger");
	}
	return(1);
}                     

//allocs a new mask
mask *mask_new(int nsam){
	mask *tmp; 
	int **tmpBits,i,j;

	if( ! (tmpBits = (int **) malloc((unsigned) (sizeof(int*)))) )
		perror("alloc error in mask_new");
	if( ! ( tmpBits[i] = (int *) malloc( (unsigned) (LENGTH*sizeof( int )) )))
		perror("alloc error in mask_new 2");
	for(j=0;j<LENGTH;j++)tmpBits[0][j] = 1;
	if( !(tmp = malloc(sizeof(mask))))
		perror("alloc error in mask_new 3");
	tmp->bitmap = tmpBits;
	return(tmp);	
}

//reads the mask file and stores them all in the ptr masks
int readMaskFiles(char *fileName,mask **masks, int nsam){
	int i,j,count;
	FILE *infile;
	char line[LINEBUF];
	double start,end;
	
	//open file
	infile = fopen(fileName, "r");
	if (infile == NULL){
		fprintf(stderr,"Error opening maskfile! ARRRRR!!!!\n");
		exit(1);
	}
	//init first mask
	count = 0;
	masks[count] = mask_new(1);
	//start going through the file
	while(fgets( line, LINEBUF,infile)){
		if (line[0] == '/'){
			count++;
			masks[count] = mask_new(1);
		}
		else if(sscanf(line,"%d %lf %lf",&i,&start,&end) == 3){
				//flip some bits
			for(j=start*LENGTH;j<=end*LENGTH;j++)
				masks[count]->bitmap[0][j] = 0;
		}
	

	}
	return(1);
}

//reads the next mask, stores it in ptr
int readNextMask(FILE *handle,mask *aMask, int nsam){
	int i,j;
	char line[LINEBUF];
	double start,end;
	
	//init mask
	for(j=0; j < LENGTH; j++){
		aMask->bitmap[0][j] = 1;
	}
	
	//start going through the file
	fgets( line, LINEBUF,handle);
	while(line[0] != '/'){
		if(sscanf(line,"%d %lf %lf",&i,&start,&end) == 3){
				//flip some bits
			for(j=start*LENGTH;j<=end*LENGTH;j++)
				aMask->bitmap[0][j] = 0;
		}
		fgets( line, LINEBUF,handle);
	}
	return(1);
}

