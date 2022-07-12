library(dada2); packageVersion("dada2")

#processing MiFish sequences
path <- "CutadaptOut_NoQC/"  
list.files(path)
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.ad.pr.fastq.gz and SAMPLENAME_R2_001.ad.pr.fastq.gz
fnFs <- sort(list.files(path, pattern="_R1_001.ad.pr.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.ad.pr.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
# Plots
#R1
plotQualityProfile(fnFs[1:2])
#R2 
plotQualityProfile(fnRs[1:2])
#Filtering
# Place filtered files in filtered/ subdirectory  
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(2,2), truncLen=c(150,120), rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
# Plot filt
#R1
plotQualityProfile(filtFs[1:2])
#R2
plotQualityProfile(filtRs[1:2])
#Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
#plot errors    
plotErrors(errF, nominalQ=TRUE)
#Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]
#Merge Pairs
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#saving ASF table

write.table(seqtab.nochim, file = "Dada_CutadaptOut_NoQC_AVS.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)

# saving track information
write.table(track, file = "Dada_CutadaptOut_NoQC_Counts.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)

