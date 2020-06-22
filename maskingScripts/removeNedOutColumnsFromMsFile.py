#!/usr/bin/env python
import sys, gzip

msFile = sys.argv[1]

def isSnp(samples, i):
    alleleH = {}
    for j in range(len(samples)):
        if samples[j][i] != "N":
            alleleH[samples[j][i]] = 1
    return len(alleleH.keys()) > 1

def processSimulation(samples, positions):
    newPositions = []
    newSamples = [""]*len(samples)
    if positions:
        for i in range(len(samples[0])):
            if isSnp(samples, i):
                newPositions.append(positions[i])
                for j in range(len(samples)):
                    newSamples[j] += samples[j][i]
        if len(newPositions) > 0:
            newPositionsStr = "positions: " + " ".join([str(newPositions[x]) for x in range(len(newPositions))])
        else:
            newPositionsStr = ""
        newSamplesStr = "\n".join([str(newSamples[x]) for x in range(len(newSamples))])
        print "\n//\nsegsites: %s\n%s\n%s" %(len(newPositions), newPositionsStr, newSamplesStr)
    else:
        print "\n//\nsegsites: 0\n"

if msFile == "stdin":
    isFile = False
    msStream = sys.stdin
else:
    isFile = True
    if msFile.endswith(".gz"):
        msStream = gzip.open(msFile)
    else:
        msStream = open(msFile)

header = msStream.readline()
program,numSamples,numSims = header.strip().split()[:3]
numSamples,numSims = int(numSamples),int(numSims)

processedSims = 0
#advance to first simulation
line = msStream.readline().strip()
print header.strip()
while not line.strip().startswith("//"):
    if line.strip() != "":
        print line.strip()
    line = msStream.readline()
while line:
    if not line.strip().startswith("//"):
        sys.exit("Malformed ms-style output file: read '%s' instead of '//'. AAAARRRRGGHHH!!!!!\n" %(line.strip()))
    segsitesBlah,segsites = msStream.readline().strip().split()
    segsites = int(segsites)
    if segsitesBlah != "segsites:":
        sys.exit("Malformed ms-style output file. AAAARRRRGGHHH!!!!!\n")

    positionsLine = msStream.readline().strip()
    if not positionsLine.startswith("positions:"):
        if segsites == 0:
            processSimulation([], [])
        else:
            sys.exit("Malformed ms-style output file. AAAARRRRGGHHH!!!!!\n")
    else:
        positionsLine = positionsLine.split()
        positions = [float(x) for x in positionsLine[1:]]

        samples = []
        for i in range(numSamples):
            sampleLine = msStream.readline().strip()
            if len(sampleLine) != segsites:
                sys.exit("Malformed ms-style output file %s segsites but %s columns in line: %s; line %s of %s samples AAAARRRRGGHHH!!!!!\n" %(segsites,len(sampleLine),sampleLine,i,numSamples))
            samples.append(sampleLine)
        if len(samples) != numSamples:
            raise Exception
        processSimulation(samples,positions)
    processedSims += 1
    line = msStream.readline()
    #advance to the next non-empty line or EOF
    while line and line.strip() == "":
        line = msStream.readline()
if processedSims != numSims:
    sys.exit("Malformed ms-style output file: %s of %s sims processed. AAAARRRRGGHHH!!!!!\n" %(processedSims,numSims))

if isFile:
    msStream.close()
