import sys
from miscFuncs import *

vcfFileName, maskedRefFileName, fastaFileName1, fastaFileName2 = sys.argv[1:]

arm = vcfFileName.split("/")[-1].split(".phased.vcf")[0]

#TODO: this pulls out all sechellia individuals included in the analysis for
# the FILET paper. These genomes had identifiers ranging from SECH_17 to
# SECH_23
def isSechToKeep(name):
    if "_" in name:
        name = int(name.split("_")[1])
        return name >= 17 and name <= 23
    else:
        return False

#TODO: this pulls out all simulans individuals included in the analysis for
# the FILET paper. These genomes had identifiers beginning with "MD" or "NS".
def isSimToKeep(name):
    return name.startswith("MD") or name.startswith("NS")

def alleleNumToBase(alleleNum, ref, alt):
    assert alleleNum in ["0","1"]
    if alleleNum == "0":
        return ref
    else:
        return alt

def getInbredHap(genos, ref, alt):
    allele1, allele2 = genos.split("|")
    if allele1 == allele2:
        return alleleNumToBase(allele1, ref, alt)
    else:
        return 'N'

def getDiploidHaps(genos, ref, alt):
    allele1, allele2 = genos.split("|")
    return alleleNumToBase(allele1, ref, alt), alleleNumToBase(allele2, ref, alt)

def writeFastaFileWithRefAndSnpGenos(names, snpGenos, arm, maskedData, fastaFileName):
    with open(fastaFileName, "w") as outFile:
        for name in names:
            sys.stderr.write("writing %s for %s\n" %(name, arm))
            outS = ">%s\n" %(name)
            for i in range(len(maskedData[arm])):
                pos = i+1
                if snpGenos.has_key((arm, pos)) and maskedData[arm][i] != 'N':
                    outS += snpGenos[(arm, pos)][name]
                else:
                    outS += maskedData[arm][i]
                if pos % 60 == 0:
                    outS += '\n'
            outFile.write(outS + '\n')
            if pos % 60 != 0:
                outFile.write('\n')

def readSnpHapsFromPhasedSimSechVcf(vcfFileName):
    simIndices, sechIndices = [], []
    simNames, sechNames = [], []
    snpGenosSim, snpGenosSech = {}, {}
    with open(vcfFileName) as vcfFile:
        for line in vcfFile:
            if line.startswith("#CHROM"):
                line = line.strip().split()
                for i in range(9, len(line)):
                    #TODO: the line below checks to see if we are examining a
                    # sechellia genome in our data set. See function definition
                    # at the top of the file and modify as necessary to pull
                    # out individuals in population 2 of your data set
                    if isSechToKeep(line[i]):
                        sechIndices.append(i)
                        sechNames.append(line[i]+".1")
                        sechNames.append(line[i]+".2")

                    #TODO: the line below checks to see if we are examining a
                    # simulans genome in our data set. See function definition
                    # at the top of the file and modify as necessary to pull
                    # out individuals in population 1 of your data set.
                    elif isSimToKeep(line[i]):
                        simIndices.append(i)
                        simNames.append(line[i])
                header = line
            elif not line.startswith("#"):
                line = line.strip().split()
                c, pos, varId, ref, alt = line[:5]
                assert ref in ['G','T','C','A'] and alt in ['G','T', 'C', 'A']
                pos = int(pos)
                snpGenosSim[(c, pos)] = {}
                snpGenosSech[(c, pos)] = {}
                for i in simIndices:
                    # TODO: here we assume our simulans data (i.e. pop 1) are inbred
                    # you may wish to replace with getDiploidHaps if both of your
                    # populations consist of outbred diploid individuals.
                    snpGenosSim[(c, pos)][header[i]] = getInbredHap(line[i], ref, alt)
                for i in sechIndices:
                    snpGenosSech[(c, pos)][header[i]+".1"], snpGenosSech[(c, pos)][header[i]+".2"] = getDiploidHaps(line[i], ref, alt)
    return simNames, snpGenosSim, sechNames, snpGenosSech

headersSim, snpGenosSim, headersSech, snpGenosSech = readSnpHapsFromPhasedSimSechVcf(vcfFileName)

maskedData = readFa(maskedRefFileName, upper=True)
writeFastaFileWithRefAndSnpGenos(headersSim, snpGenosSim, arm, maskedData, fastaFileName1)
writeFastaFileWithRefAndSnpGenos(headersSech, snpGenosSech, arm, maskedData, fastaFileName2)
