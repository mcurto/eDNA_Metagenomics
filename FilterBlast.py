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


def FilterBlast(BlastOut, out, minId , minEval, Overlap, TotalNrLines):
    result = {}
    with open(BlastOut) as f:
        qid = ''
        eval = 100
        id = 0
        cov = 0
        rec = 0
        lineOut = ''
        nrLine = 0
        for l in f:
            nrLine += 1
            led = l.rstrip('\n').split('\t')
            qidLed = led[0]
            idLed = float(led[3])
            covLed = int(led[5])/int(led[4])
            evalLed = float(led[6])
            lineage = led[-1]

            if idLed >= minId and evalLed <= minEval and covLed >= Overlap:
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
                    eval = evalLed
                    id = idLed
                    cov = covLed
                    lineOut = l


                elif nrLine == TotalNrLines:

                    if evalLed < eval:
                        eval = evalLed
                        id = idLed
                        cov = covLed
                        lineOut = l
                        out.write(lineOut)
                    elif evalLed == eval:
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
                    if evalLed < eval:
                        eval = evalLed
                        id = idLed
                        cov = covLed
                        lineOut = l
                    elif evalLed == eval:
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
            elif nrLine == TotalNrLines:
                out.write(lineOut)

inDir = sys.argv[1]
outDir = sys.argv[2]
eval = sys.argv[3]
id = sys.argv[4]
Overlap = sys.argv[5]

for f in os.listdir(inDir):
    if f.endswith(".blastout"):
        print("processing " + f)
        out = open(outDir + f.replace(".blastout", ".filt.ID" + id + "Eval" + str(eval).split(".")[-1] + "Ov" + Overlap + ".blastout"), "w")
        nrLines = countLines(inDir + f)
        FilterBlast(inDir + f, out, float(id), float(eval), float(Overlap), nrLines)
        out.close()
