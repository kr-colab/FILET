def readFa(faFileName, upper=False):
    seqData = {}
    with open(faFileName) as faFile:
        reading = False
        for line in faFile:
            if line.startswith(">"):
                if reading:
                    if upper:
                        seqData[currChr] = seq.upper()
                    else:
                        seqData[currChr] = seq
                else:
                    reading = True
                currChr = line[1:].strip()
                seq = ""
            else:
                seq += line.strip()
    if upper:
        seqData[currChr] = seq.upper()
    else:
        seqData[currChr] = seq
    return seqData

def readFaAsLists(faFileName, upper=False):
    seqData = {}
    with open(faFileName) as faFile:
        reading = False
        for line in faFile:
            if line.startswith(">"):
                if reading:
                    if upper:
                        seqData[currChr] = list(seq.upper())
                    else:
                        seqData[currChr] = list(seq)
                else:
                    reading = True
                currChr = line[1:].strip()
                seq = ""
            else:
                seq += line.strip()
    if upper:
        seqData[currChr] = list(seq.upper())
    else:
        seqData[currChr] = list(seq)
    return seqData
