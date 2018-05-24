//pgStats for a bedFile rather than windows
// calculates fst and related

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "miscCode/stringWrap.h"
#include "miscCode/sequenceMatrix.h"
#include "pgSummaryStats.h"
#include "miscCode/bedFile.h"


void usage();

int main(int argc, char *argv[]){
	struct sequenceMatrix *fullSeqMat, *aSeqMat, *ancSeqMat, *bSeqMat;
	int i, bedElNumber, h, ss, private1=0, private2=0;
	long int numSites;
	struct bedEl data[80000];
	double pi, theta_h, z, H, tajD, sStarVal, *hetVec, *ibsVec;
	double ibsMean, ibsVar, ibsSkew, ibsKurt, ibsMin, ibsMed, ibsMax;
	double hetMean, hetVar, hetSkew, hetKurt, hetMin, hetMed, hetMax;
	double hetVar1, hetVar2;
	double z1, z2, pi1, pi2, f, snn, dxy_min, dxy_mean, gmin;
	double dd1, dd2, ddRank1, ddRank2, ibsMaxBetween, ibsMeanWithin1, ibsMeanWithin2;
	int hetVecLen, ibsVecLen;
	float maxFracMissingData;
	
	if(argc < 5){
		usage();
		exit(1);
	}

	//open fastaFile and bedFile
	fullSeqMat = sequenceMatrix_importFasta(argv[1]);
	aSeqMat = malloc(sizeof(sequenceMatrix));
	bSeqMat = malloc(sizeof(sequenceMatrix));
	ancSeqMat = sequenceMatrix_importFasta(argv[2]);
	bedElNumber = bedFileImport3(argv[3], data);
	maxFracMissingData = atof(argv[4]);
	
	sequenceMatrix_NOutSitesWithTooMuchMissingData(fullSeqMat,maxFracMissingData);

	//print header
	printf("chrom\tchromStart\tchromEnd\tnumSites\tpi\tss\tthetaH\ttajd\tH\tHapCount\tZnS\t");
        printf("hetVar\thetSkew\thetKurt\thetMin\thetMed\thetMax\t");
        printf("ibsMean\tibsVar\tibsSkew\tibsKurt\tibsMin\tibsMed\tibsMax\tS*\t");
	printf("Fst\tsnn\tdxy_mean\tdxy_min\tgmin\tzx\tdd1\tdd2\tddRank1\tddRank2\tibsMaxB\tibsMean1\tibsMean2\tprivate1\tprivate2\n");
	

	//loop through beds; the adjustments to end are to honor the zero indexed half open bed convention
	for(i=0;i<bedElNumber;i++){
		numSites = numColumnsNotNedOutFromTo(fullSeqMat, data[i].chromStart, data[i].chromEnd+1);
		printf("%s\t%ld\t%ld\t%ld\t",data[i].chrom,data[i].chromStart,data[i].chromEnd, numSites);

		hetVec = hetVec1PopnFromTo(fullSeqMat, &hetVecLen, data[i].chromStart, data[i].chromEnd+1);
		clusterSeqsFromUnsortedHetVec(hetVec, fullSeqMat, aSeqMat, bSeqMat);
		statVecMoments(hetVec, hetVecLen, &hetMean, &hetVar, &hetSkew, &hetKurt);
		statVecMinMedMax(hetVec, hetVecLen, &hetMin, &hetMed, &hetMax);
		ibsVec = pairwiseIBSVec1PopnFromTo(fullSeqMat, &ibsVecLen, data[i].chromStart, data[i].chromEnd+1);
		statVecMoments(ibsVec, ibsVecLen, &ibsMean, &ibsVar, &ibsSkew, &ibsKurt);
		statVecMinMedMax(ibsVec, ibsVecLen, &ibsMin, &ibsMed, &ibsMax);
		sStarVal = sStarFromTo(fullSeqMat, data[i].chromStart, data[i].chromEnd+1);
		//ingroup 1 stats
		pi = nucdivFromTo(fullSeqMat, data[i].chromStart, data[i].chromEnd+1);
		ss = segSiteCountFromTo(fullSeqMat, data[i].chromStart, data[i].chromEnd+1);
		theta_h = thetaHFromTo(fullSeqMat,ancSeqMat,data[i].chromStart, data[i].chromEnd+1);
		h = nHaplotypes(fullSeqMat,data[i].chromStart, data[i].chromEnd+1);
		z = ZnSFromTo(fullSeqMat,data[i].chromStart, data[i].chromEnd+1);
		tajD = tajd(fullSeqMat->sampleSize,ss,pi);
		H = theta_h-pi;
		printf("%f\t%d\t%f\t%f\t%f\t%d\t%f\t",pi,ss,theta_h,tajD,H,h,z);
		printf("%f\t%f\t%f\t%f\t%f\t%f\t", hetVar, hetSkew, hetKurt, hetMin, hetMed, hetMax);
		printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t", ibsMean, ibsVar, ibsSkew, ibsKurt, ibsMin, ibsMed, ibsMax, sStarVal);
		
		dxy_min = Dxy_minFromTo(aSeqMat, bSeqMat, data[i].chromStart, data[i].chromEnd+1);
		ddRank1 = pairwiseDistRankAmongSampleRange(aSeqMat, dxy_min, &hetVar1, data[i].chromStart, data[i].chromEnd+1);
		ddRank2 = pairwiseDistRankAmongSampleRange(bSeqMat, dxy_min, &hetVar2, data[i].chromStart, data[i].chromEnd+1);
		f = fstFromTo(aSeqMat,bSeqMat, fullSeqMat, data[i].chromStart, data[i].chromEnd+1);
		snn = SnnFromTo(aSeqMat, bSeqMat, data[i].chromStart, data[i].chromEnd+1);
		dxy_mean = DxyFromTo(aSeqMat, bSeqMat, data[i].chromStart, data[i].chromEnd+1);
		gmin = dxy_min / dxy_mean; 
		pi1 = nucdivFromTo(aSeqMat, data[i].chromStart, data[i].chromEnd+1);
		pi2 = nucdivFromTo(bSeqMat, data[i].chromStart, data[i].chromEnd+1);
		dd1 = dxy_min / pi1;
		dd2 = dxy_min / pi2;
		
		z1 = ZnSFromTo(aSeqMat,data[i].chromStart, data[i].chromEnd+1);
		z2 = ZnSFromTo(bSeqMat,data[i].chromStart, data[i].chromEnd+1);
		printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t",f,snn,dxy_mean,dxy_min,gmin,(z1+z2)/2.0/z,dd1,dd2,ddRank1,ddRank2) ;

		ibsMaxBetween = pairwiseIBSMax2PopnFromTo(aSeqMat, bSeqMat, data[i].chromStart, data[i].chromEnd+1);
		ibsMeanWithin1 = pairwiseIBSMeanWithinFromTo( aSeqMat, data[i].chromStart, data[i].chromEnd+1);
		ibsMeanWithin2 = pairwiseIBSMeanWithinFromTo( bSeqMat, data[i].chromStart, data[i].chromEnd+1);
		privateSegSitesInTwoPopnsFromTo(aSeqMat, bSeqMat, &private1, &private2, data[i].chromStart, data[i].chromEnd+1);
		printf("%f\t%f\t%f\t%d\t%d",ibsMaxBetween,ibsMeanWithin1,ibsMeanWithin2,private1,private2);
		
		
		printf("\n");
		free(hetVec);
		free(ibsVec);
	}
	sequenceMatrix_free(aSeqMat);
	sequenceMatrix_free(bSeqMat);
	free(fullSeqMat);//the parts of this guy were already freed in the above two lines -- yes, this is ugly
	return(0);
}	

void usage(){
	printf("pgStatsBedOnePop_forGhostML ingroupFastaFile ancestorFastaFile bedFile maxFractionMissingData\n");
}
