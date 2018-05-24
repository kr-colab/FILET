/******* maskedStats.c ********
for calculating sample stats from MS output 
after it has been filtered by msMask
********************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "msGeneralStats.h"

#define LINEBUF 1000000
int maxsites = 100000 ;
void usage();

int main(argc,argv)
	int argc;
char *argv[];
{
	int nsam, i,  howmany  ;
	char **list, **cmatrix(), line[LINEBUF+1]  ;
	FILE *fopen(), *pfin ;
	double *posit, *ibsVec, *hetVec;
	int   segsites, count, n1, iss, h, ibsVecLen, hetVecLen, g1Size, g2Size, private1, private2;
	double pi, th,  z, H, tajD, sStarVal;
	double ibsMean, ibsVar, ibsSkew, ibsKurt, ibsMin, ibsMed, ibsMax;
	double hetMean, hetVar, hetSkew, hetKurt, hetMin, hetMed, hetMax;
        double z1, z2, ztot, pi1, pi2, f, snn, dxy_min, dxy_mean, gmin, hetVar1, hetVar2;
	double dd1, dd2, ddRank1, ddRank2, ibsMaxBetween, ibsMeanWithin1, ibsMeanWithin2;
	char dum[20], astr[100] ;
	int starFlag=0;
	int migOption=0;
	int physLen;
//	int bins = 10;
//	double hist[bins];

	if( argc > 2 ) { 
		n1 = atoi( argv[1] ) ;
		physLen = atoi(argv[2]);
		if(argc==4){
			if(argv[3][1] == 'c') migOption=1;
		}
	}
	else{
		usage();
	}

/* read in first two lines of output  (parameters and seed) */
	pfin = stdin ;
	fgets( line, LINEBUF, pfin);
	sscanf(line," %s  %d %d", dum,  &nsam, &howmany);
	fgets( line, LINEBUF, pfin);

	if (n1 <= 0) n1 = nsam;

	list = cmatrix(nsam,maxsites+1);
	posit = (double *)malloc( maxsites*sizeof( double ) ) ;

	count=0;

	//print header line
	printf("pi\tss\tthetaH\ttajd\tH\tHapCount\tZnS\t");
	printf("hetVar\thetSkew\thetKurt\thetMin\thetMed\thetMax\t");
	printf("ibsMean\tibsVar\tibsSkew\tibsKurt\tibsMin\tibsMed\tibsMax\tS*\t");
	printf("Fst\tsnn\tdxy_mean\tdxy_min\tgmin\tzx\tdd1\tdd2\tddRank1\tddRank2\tibsMaxB\tibsMean1\tibsMean2\tprivate1\tprivate2\n");
//	for(i=0;i<bins;i++){
//		printf("\tibs[%d]",i);
//	}
//	printf("\n");
	
	
	while( howmany-count++ ) {

/* read in a sample */
		do {
			fgets( line, LINEBUF, pfin);
		}while ( line[0] != '/' );
		if(line[2] == '*'){
			starFlag = 1;
		}
		else{ 
			starFlag = 0;
		}
		if(migOption==0)starFlag=1;
		
		fscanf(pfin,"  segsites: %d", &segsites );
		if( segsites >= maxsites){
			maxsites = segsites + 10 ;
			posit = (double *)realloc( posit, maxsites*sizeof( double) ) ;
			biggerlist(nsam,maxsites, list) ;
		}
		if( segsites > 0) {
			fscanf(pfin," %s", astr);

			for( i=0; i<segsites ; i++) fscanf(pfin," %lf",posit+i) ;
			for( i=0; i<nsam;i++) fscanf(pfin," %s", list[i] );
		}
		/* analyse sample ( do stuff with segsites and list) */

		pi = nucdivSub(nsam,segsites,0,n1,list);
		iss = segSitesSub(segsites,nsam,0,n1,list);
		th = thetahSub(nsam, segsites,0,n1, list) ;
		h = nHaplotypesSub(segsites,nsam,0,n1,list);
		z = ZnSSub( segsites,  nsam, 0,n1, list);
		H = th-pi;
		tajD = tajd(nsam,iss,pi);
		hetVec = hetVec1Popn(segsites, n1, physLen, &hetVecLen, posit, list);
		//srand ( time(NULL) );//for clustering tiebreaks;
		clusterSeqsFromUnsortedHetVec(hetVec, n1, &g1Size, &g2Size, list);
		statVecMoments(hetVec, hetVecLen, &hetMean, &hetVar, &hetSkew, &hetKurt);
		statVecMinMedMax(hetVec, hetVecLen, &hetMin, &hetMed, &hetMax);
		ibsVec = pairwiseIBSVec1Popn(segsites, n1, &ibsVecLen, posit, list);
		statVecMoments(ibsVec, ibsVecLen, &ibsMean, &ibsVar, &ibsSkew, &ibsKurt);
		statVecMinMedMax(ibsVec, ibsVecLen, &ibsMin, &ibsMed, &ibsMax);
		sStarVal = sStar(nsam, segsites, physLen, posit, list);
		free(hetVec);
		free(ibsVec);

		snn = Snn(segsites,nsam,g1Size,g2Size,list);
                dxy_min = Dxy_min(segsites,nsam,g1Size,g2Size,list);
                dxy_mean = Dxy_mean(segsites,nsam,g1Size,g2Size,list);
                f = fst2Subs(segsites,nsam,0,g1Size,g1Size,g1Size+g2Size,list);
                gmin = dxy_min / dxy_mean;
		z1 = ZnSSub( segsites, nsam, 0, g1Size, list);
		z2 = ZnSSub( segsites, nsam, g1Size, nsam, list);
                ztot = ZnSSub( segsites, nsam, 0, nsam, list);
		pi1 = nucdivSub(nsam,segsites,0,g1Size,list);
                pi2 = nucdivSub(nsam,segsites,g1Size,nsam,list);
                dd1 = dxy_min / pi1;
                dd2 = dxy_min / pi2;
                ddRank1 = pairwiseDistRankAmongSampleRange( segsites, dxy_min, 0, g1Size, &hetVar1, list);
                ddRank2 = pairwiseDistRankAmongSampleRange( segsites, dxy_min, g1Size, g2Size, &hetVar2, list);
		privateSegSitesInTwoPopns(segsites, nsam, g1Size, &private1, &private2, list);
		if(starFlag==1){
                        ibsMaxBetween = pairwiseIBSMax2Popn(segsites, nsam, g1Size, posit, list);
                        ibsMeanWithin1 = pairwiseIBSMeanWithin(segsites,0, g1Size, posit, list);
                        ibsMeanWithin2 = pairwiseIBSMeanWithin(segsites, g1Size, nsam, posit, list);
			printf("%lf\t%d\t%lf\t%f\t%f\t%d\t%f\t",pi, iss, th, tajD, H, h, z) ;
			printf("%g\t%f\t%f\t%f\t%f\t%f\t", hetVar, hetSkew, hetKurt, hetMin, hetMed, hetMax);
			printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t", ibsMean, ibsVar, ibsSkew, ibsKurt, ibsMin, ibsMed, ibsMax, sStarVal);
			printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t",f, snn, dxy_mean, dxy_min, gmin, (z1+z2)/2.0/ztot, dd1, dd2, ddRank1, ddRank2);
                        printf("%f\t%f\t%f\t%d\t%d\n", ibsMaxBetween, ibsMeanWithin1, ibsMeanWithin2, private1, private2);
		}
	}
	free(posit);
	cmatrix_free(list,nsam);
	return(0);
}

void usage(){
	printf("onePopnStats_forGhostIntroML n1 physLen\n");
	printf("Returns analysis of Hudson style output assuming two subpops; first is of size n1 and second is ignored. ");
	printf("Physical sequence length (physLen) is required for Plagnol and Wall's S* statistic\n");
	printf("options:\n\t-c <condition on migration>\n");
	exit(1);
}

