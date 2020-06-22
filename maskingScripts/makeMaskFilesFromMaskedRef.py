import sys,os
import miscFuncs

maskedRefFaFileName, maskFileDir, winSize = sys.argv[1:]
winSize = int(winSize)
os.system("mkdir -p %s" %(maskFileDir))

maskedGenome = miscFuncs.readFa(maskedRefFaFileName)

def getMaskedRuns(seq, winSize):
    inRun = False
    runs = []
    for i in range(len(seq)):
        if inRun:
            if seq[i] == 'N':
                end = i
            else:
                runs.append((start/float(winSize), end/float(winSize)))
                inRun = False
        else:
            if seq[i] == 'N':
                start = i
                end = i
                inRun = True
            else:
                pass
    if inRun:
        runs.append((start/float(winSize), end/float(winSize)))
    return runs

def writeMaskedRegionsToFile(seq, maskFileName, winSize):
    assert len(seq) == winSize
    with open(maskFileName, "w") as maskFile:
        maskedRuns = getMaskedRuns(seq, winSize)
        for s, e in maskedRuns:
            maskFile.write("0 %f %f\n" %(s, e))
        maskFile.write("//\n")

for arm in maskedGenome:
    prevEnd = 0
    while prevEnd < len(maskedGenome[arm]):
        currStart = prevEnd+1
        currEnd = currStart + winSize-1
        if currEnd <= len(maskedGenome[arm]):
            maskFileName = maskFileDir + "/%s.%d.%d.maskedRegions" %(arm, currStart, currEnd)
            writeMaskedRegionsToFile(maskedGenome[arm][currStart-1:currEnd], maskFileName, winSize)
        prevStart = currStart
        prevEnd = currEnd
