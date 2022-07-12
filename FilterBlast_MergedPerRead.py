from __future__ import division
import sys, os


def mergeLineOut(line1, line2):
    led1 = line1.rstrip('\n').split('\t')
    led2 = line2.rstrip('\n').split('\t')
    newLine = []
    for i in range(len(led1)):
        if led1[i] == led2[i]:
            newLine.append(led1[i])
        else:
            newLine.append(led1[i] + '&' + led2[i])
    return '\t'.join(newLine) + '\n'

def countLines(blast):
    nrLines = 0
    with open(blast) as b:
        for l in b:
            nrLines += 1
    return nrLines


def FilterBlast(BlastOut, out, TotalNrLines):
    result = {}
    with open(BlastOut) as f:
        qid = ''
        id = 0
        cov = 0
        rec = 0
        lineOut = ''
        nrLine = 0
        for l in f:
            nrLine += 1
            led = l.rstrip('\n').split('\t')
            qidLed = led[0]
            idLed = float(led[1])
            covLed = int(led[3])/int(led[2])
            lineage = led[-1]


            BestData = result.get(idLed, [100, 0, 0, 0])
            BestEval = BestData[0]
            Bestid = BestData[1]
            Bestcov = BestData[2]
            Bestrec = BestData[3]




            if qidLed != qid:
                if rec > 0:
                    out.write(lineOut)
                else:
                    rec += 1

                qid = qidLed
                id = idLed
                cov = covLed
                lineOut = l



            elif nrLine == TotalNrLines:

                    if idLed > id:
                        id = idLed
                        cov = covLed
                        lineOut = l
                        out.write(lineOut)
                    elif idLed == id:
                        if cov < covLed:
                            cov = covLed
                            lineOut = l
                            out.write(lineOut)
                        elif cov == covLed:
                            lineOut = mergeLineOut(lineOut, l)
                            out.write(lineOut)
                    else:
                        out.write(lineOut)

            else:
                if idLed > id:
                    id = idLed
                    cov = covLed
                    lineOut = l
                elif idLed == id:
                    if cov < covLed:
                        cov = covLed
                        lineOut = l
                    elif cov == covLed:
                        lineOut = mergeLineOut(lineOut, l)


inDir = sys.argv[1]
outDir = sys.argv[2]


for f in os.listdir(inDir):
    if f.endswith(".blastout"):
        print("processing " + f)
        out = open(outDir + f.replace(".blastout", ".Best.blastout"), "w")
        nrLines = countLines(inDir + f)
        FilterBlast(inDir + f, out, nrLines)
        out.close()
