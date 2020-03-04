import sys
from miscFuncs import *

vcfFileName1, vcfFileName2, repeatGffFileName, refFaFileName, qualThresh, sampSizePerc, outMaskedFaFileName = sys.argv[1:]
qualThresh = int(qualThresh)
sampSizePerc = float(sampSizePerc)

def writeFa(seqData, outFaFileName, upper=False):
    with open(outFaFileName, "w") as outFaFile:
        for arm in sorted(seqData):
            outStr = ">" + arm
            for i in xrange(len(seqData[arm])):
                if i % 60 == 0:
                    outStr += "\n"
                outStr += seqData[arm][i]
            outStr += "\n\n"
            outFaFile.write(outStr)

genome = readFaAsLists(refFaFileName, upper=True)

with open(repeatGffFileName) as repeatGffFile:
    for line in repeatGffFile:
        if not line.startswith("#"):
            line = line.strip().split("\t")
            arm, prog, basis, s, e, score, strand, blah, subject = line
            arm = "chr" + arm
            s, e = int(s), int(e)
            for pos in range(s, e+1):
                genome[arm][pos-1] = 'N'

altAlleles = {}
qualFiltered = {}
nFiltered = {}
indelFiltered = {}
missingFiltered = {}
totalSnps = {}
prevArm = False
for vcfFileName in [vcfFileName1, vcfFileName2]:
    with open(vcfFileName) as vcfFile:
        for line in vcfFile:
            if line.startswith("#CHROM\tPOS\t"):
                isInbred = [False]*9
                line = line.strip().split("\t")
                for i in range(9, len(line)):

                    #TODO: the line below is meant to flag individuals that are
                    # inbred and should thus not have heterozygous base calls.
                    # In the sim-sech study this applied to all simulans genomes
                    # which began either with the prefix MD or NS (while all of
                    # the sechellia genome identifiers began with SECH.
                    if line[i][:2] in ("MD", "NS"):
                        isInbred.append(True)
                    else:
                        isInbred.append(False)
            elif not line.startswith("#"):
                line = line.strip().split("\t")
                arm, pos, varId, ref, alt, qual, varFilter, info, fmt = line[:9]
                pos = int(pos)
                if arm == prevArm and pos > prevPos+1:
                    for missingPos in range(prevPos+1, pos):
                        genome[arm][missingPos-1] = 'N'
                        missingFiltered[(arm, pos)] = 1
                #if arm != "chr2L" or pos > 10000:
                #    sys.exit()
                if varFilter == "LowQual":
                    genome[arm][pos-1] = 'N'
                    qualFiltered[(arm, pos)] = 1
                elif not alt in list("ACGT."):
                    genome[arm][pos-1] = 'N'
                    indelFiltered[(arm, pos)] = 1
                elif not ref in list("ACGT."):
                    if genome[arm][pos-1] == 'N':
                        nFiltered[(arm, pos)] = 1
                    else:
                        indelFiltered[(arm, pos)] = 1
                        genome[arm][pos-1] = 'N'
                else:
                    genoFormat = fmt.split(":")
                    genoIndex = None
                    gqIndex = None
                    for j in range(len(genoFormat)):
                        if genoFormat[j] in ["GQ","RGQ"]:
                            gqIndex = j
                        elif genoFormat[j] == "GT":
                            genoIndex = j
                    goodSamps, totalSamps = 0, 0
                    for i in range(9, len(line)):
                        gtLs = line[i].split(":")
                        if gtLs[genoIndex] == "./." or (len(gtLs) == 1 and gtLs[genoIndex] == "0/0") or (isInbred[i] and gtLs[genoIndex] == "0/1"):
                            gq = 0
                        else:
                            try:
                                gq = int(gtLs[gqIndex])
                            except Exception:
                                print line
                                print gtLs
                                print genoIndex, gqIndex
                                raise Exception
                        if gq >= qualThresh:
                            goodSamps += 1
                        totalSamps += 1
                    if alt != ".":
                        if not (arm, pos) in altAlleles:
                            altAlleles[(arm, pos)] = {}
                        altAlleles[(arm, pos)][alt] = 1
                        totalSnps[(arm, pos)] = 1
                    if goodSamps / float(totalSamps) < sampSizePerc:
                        genome[arm][pos-1] = 'N'
                        if alt != ".":
                            qualFiltered[(arm, pos)] = 1
                prevArm = arm
                prevPos = pos

triallelicFiltered = 0
for arm, pos in altAlleles:
    if len(altAlleles[(arm, pos)]) > 1:
        triallelicFiltered += 1
        genome[arm][pos-1] = 'N'

writeFa(genome, outMaskedFaFileName, upper=True)
print "masked for qual", len(qualFiltered)
print "already N in ref", len(nFiltered)
print "masked for indel", len(indelFiltered)
print "masked for missing", len(missingFiltered)
print "masked for triallelic", triallelicFiltered
print "total snp count:", len(totalSnps)
