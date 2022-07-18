from __future__ import division
import sys, os

# Merge results from two equally good matches
def mergeLineOut(line1, line2):
    # make a list of the match record
    led1 = line1.rstrip('\n').split('\t')
    led2 = line2.rstrip('\n').split('\t')
    newLine = []
    #iterate simultaneously through the elements in both records and compare them. If they show the same value, keep it; otherwise merge them togetehr seperated with "&"
    for i in range(len(led1)):
        if led1[i] == led2[i]:
            newLine.append(led1[i])
        else:
            newLine.append(led1[i] + '&' + led2[i])
    return '\t'.join(newLine) + '\n'

#Get the number of lines in the blast output file
def countLines(blast):
    nrLines = 0
    with open(blast) as b:
        for l in b:
            nrLines += 1
    return nrLines

#Read the blast results per read and find the best match
def FilterBlast(BlastOut, out, minId , minEval, Overlap, TotalNrLines):
    result = {}
    with open(BlastOut) as f:
        #Set initial reference blast parameters to wors values as possible
        qid = ''
        eval = 100
        id = 0
        cov = 0
        rec = 0
        lineOut = ''
        nrLine = 0
        for l in f:
            #getting blast parameters from current match
            nrLine += 1
            led = l.rstrip('\n').split('\t')
            qidLed = led[0]
            idLed = float(led[3])
            covLed = int(led[5])/int(led[4])
            evalLed = float(led[6])
            lineage = led[-1]

            #Check if match passes the threashold values
            if idLed >= minId and evalLed <= minEval and covLed >= Overlap:

                #Check if the match corresponds to a new record. If that is the case, it sets the reference parameter to the ones found for the match and saves the results from the previous match
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

                # In case it is a suboptimal match and it corresponds to the last line of the file, it compares the parameters to the best match so far
                elif nrLine == TotalNrLines:

                    #First the evalue and if this value is lower than the reference consider the match from the last line as the best one and save it
                    if evalLed < eval:
                        eval = evalLed
                        id = idLed
                        cov = covLed
                        lineOut = l
                        out.write(lineOut)
                    #if evalue is the same it compares the identity value
                    elif evalLed == eval:
                        if idLed > id:
                            id = idLed
                            cov = covLed
                            lineOut = l
                            out.write(lineOut)
                        # if id value is the same it compares the query coverage
                        elif idLed == id:
                            if cov < covLed:
                                cov = covLed
                                lineOut = l
                                out.write(lineOut)
                            # if query coverage is the same merge current match with best match
                            elif cov == covLed:
                                lineOut = mergeLineOut(lineOut, l)
                                out.write(lineOut)
                    #if the parameters from the match of the last line are not better than the reference save the best match so far
                    else:
                        out.write(lineOut)

                #if the match does not correspond to the last line use the same conditions to decide if it is the match
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

#Import information from positional variables
inDir = sys.argv[1]
outDir = sys.argv[2]
eval = sys.argv[3]
id = sys.argv[4]
Overlap = sys.argv[5]

#Iterate through multiple files and save filtered results to the output directory
for f in os.listdir(inDir):
    if f.endswith(".blastout"):
        print("processing " + f)
        out = open(outDir + f.replace(".blastout", ".filt.ID" + id + "Eval" + str(eval).split(".")[-1] + "Ov" + Overlap + ".blastout"), "w")
        nrLines = countLines(inDir + f)
        FilterBlast(inDir + f, out, float(id), float(eval), float(Overlap), nrLines)
        out.close()
