#Adapted from the following resources:
#https://benjjneb.github.io/dada2/ITS_workflow.html
#https://benjjneb.github.io/dada2/tutorial.html
#https://benjjneb.github.io/dada2/bigdata_paired.html

#Adapting for Novaseq data:
#https://github.com/benjjneb/dada2/issues/1307
#https://github.com/ErnakovichLab/dada2_ernakovichlab/tree/split_for_premise


library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(Biostrings); packageVersion("Biostrings")
library(magrittr)
library(dplyr)

# Import available CPU number from SLURM (determined from jobscript)
numCPUs <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))

# Import raw data.
path <- "/work/bcp30/LZ_scord/miseq/PE250"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

######################################
#### Remove primers with Cutadapt ####
######################################

# Define primers.
FWD <- "AACMGGATTAGATACCCKG"
REV <- "ACGTCATCCCCACCTTCC"

# Ensure primers are correct and correct orientation.
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna),
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

# Pre-filter sequences by removing ambiguous bases (Ns) to make mapping of short primer sequences less difficult.
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = numCPUs) # on windows, set multithread = FALSE

# Count number of times the primers appear in F and R reads.
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

# Remove primers with cutadapt.
cutadapt <- "cutadapt" # cutadapt conda environment must be activated, otherwise CHANGE ME to the cutadapt path on your machine.
system2(cutadapt, args = "--version") # Run shell commands from R, will print version if correctly configured.

# Create output file names for the cutadapt-ed files.
path.cut <- file.path(path, "cutadapt-0.2_min200")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC, "--minimum-length 200")
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC, "--minimum-length 200")
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-e", 0.2, # Max error rate of overlapping region
                             "--match-read-wildcards", # Allow IUPAC wildcards in reads (default: False)
                             "--discard-untrimmed", # Discard reads that do not contain primers
                             "-j", numCPUs, # On Windows set "numCPUs" to "0" to enable multithreading on all avaibable cpus
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

# Uncomment line below to save object.
#saveRDS(sample.names, file="sample.names.rds")

# Visualize quality profiles of forward reads.
pdf(file="cutFs_plot.pdf")
plotQualityProfile(cutFs[1:2])
dev.off()

# Visualize quality profiles of reverse reads.
pdf(file="cutRs_plot.pdf")
plotQualityProfile(cutRs[1:2])
dev.off()


# Assign filenames.
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))


# Filter.
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, truncLen=c(220,220), maxN = 0, maxEE = c(2, 2),
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = numCPUs)  # on windows, set multithread = FALSE
head(out)

# Uncomment lines below to save objects.
#saveRDS(filtFs, file="filtFs.rds")
#saveRDS(filtRs, file="filtRs.rds")
#saveRDS(out, file="out.rds")

#################################################
#### Learn error rates with default settings ####
#################################################

errF <- learnErrors(filtFs, multithread = numCPUs)
errR <- learnErrors(filtRs, multithread = numCPUs)

# Visualize estimated error rates.
errF_plot <- plotErrors(errF, nominalQ = TRUE)
errR_plot <- plotErrors(errR, nominalQ = TRUE)

pdf(file="errF_plot.pdf")
errF_plot
dev.off()

pdf(file="errR_plot.pdf")
errR_plot
dev.off()


#########################################
#### Relearn errors for NovaSeq data ####
#########################################

## Note: This section can be skipped for MiSeq data because the standard error model is sufficient.

##Relearning errors for NovaSeq data: "Trial 1"
##Because error is underestimated with NovaSeq data (binned q-scores), re-learn error rates with altered weight and span in loess and with forced monotonicity. From JacobRPrice: https://github.com/benjjneb/dada2/issues/1307

loessErrfun_mod1 <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # Gulliem Salazar's solution
        # https://github.com/benjjneb/dada2/issues/938
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

# check what this looks like
errF_1 <- learnErrors(
  filtFs,
  multithread = numCPUs,
  errorEstimationFunction = loessErrfun_mod1,
  verbose = TRUE
)

errR_1 <- learnErrors(
  filtRs,
  multithread = numCPUs,
  errorEstimationFunction = loessErrfun_mod1,
  verbose = TRUE
)

# Uncomment lines below to save objects.
#saveRDS(errF_1, file="errF_1.rds")
#saveRDS(errR_1, file="errR_1.rds")

##Visualize newly estimated error rates: Trial 1 (alter span and weight in loess, enforce montonicity)
errF_1plot <- plotErrors(errF_1, nominalQ = TRUE)
errR_1plot <-plotErrors(errR_1, nominalQ = TRUE)

pdf(file="errF_1plot.pdf")
errF_1plot
dev.off()

pdf(file="errR_1plot.pdf")
errR_1plot
dev.off()

##Relearning errors for NovaSeq data: "Trial 4"
##Because error is underestimated with NovaSeq data (binned q-scores), re-learn error rates with altered weight, span AND DEGREE in loess and with forced monotonicity. From JacobRPrice: https://github.com/benjjneb/dada2/issues/1307

#loessErrfun_mod4 <- function(trans) {
#  qq <- as.numeric(colnames(trans))
#  est <- matrix(0, nrow=0, ncol=length(qq))
#  for(nti in c("A","C","G","T")) {
#    for(ntj in c("A","C","G","T")) {
#      if(nti != ntj) {
#        errs <- trans[paste0(nti,"2",ntj),]
#        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
#        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
#        rlogp[is.infinite(rlogp)] <- NA
#        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
#        
#        # original
#        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
#        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
#        # #        mod.lo <- loess(rlogp ~ q, df)
#        
#        # jonalim's solution
#        # https://github.com/benjjneb/dada2/issues/938
#        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),degree = 1, span = 0.95)
#        
#        pred <- predict(mod.lo, qq)
#        maxrli <- max(which(!is.na(pred)))
#        minrli <- min(which(!is.na(pred)))
#        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
#        pred[seq_along(pred)<minrli] <- pred[[minrli]]
#        est <- rbind(est, 10^pred)
#      } # if(nti != ntj)
#    } # for(ntj in c("A","C","G","T"))
#  } # for(nti in c("A","C","G","T"))
#  
#  # HACKY
#  MAX_ERROR_RATE <- 0.25
#  MIN_ERROR_RATE <- 1e-7
#  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
#  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
#  
#  # enforce monotonicity
#  # https://github.com/benjjneb/dada2/issues/791
#  estorig <- est
#  est <- est %>%
#    data.frame() %>%
#    mutate_all(funs(case_when(. < X40 ~ X40,
#                              . >= X40 ~ .))) %>% as.matrix()
#  rownames(est) <- rownames(estorig)
#  colnames(est) <- colnames(estorig)
#  
#  # Expand the err matrix with the self-transition probs
#  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
#               est[4,], 1-colSums(est[4:6,]), est[5:6,],
#               est[7:8,], 1-colSums(est[7:9,]), est[9,],
#               est[10:12,], 1-colSums(est[10:12,]))
#  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
#  colnames(err) <- colnames(trans)
#  # Return
#  return(err)
#}
#
# check what this looks like
#errF_4 <- learnErrors(
#  filtFs,
#  multithread = numCPUs,
#  errorEstimationFunction = loessErrfun_mod4,
#  verbose = TRUE
#)
#
#errR_4 <- learnErrors(
#  filtRs,
#  multithread = numCPUs,
#  errorEstimationFunction = loessErrfun_mod4,
#  verbose = TRUE
#)
#
###Visualize newly estimated error rates: Trial 4 (alter loess (span, weight, and degree) and enforce monotonicity)
#errF_4plot <- plotErrors(errF_4, nominalQ = TRUE)
#errR_4plot <-plotErrors(errR_4, nominalQ = TRUE)
#
#pdf(file="errF_4plot.pdf")
#errF_4plot
#dev.off()
#
#pdf(file="errR_plot4.pdf")
#errR_4plot
#dev.off()

###################################################################
#### Dereplication, inference, and merging of paired-end reads ####
###################################################################

# Note: Resume here for MiSeq data. Adapted from: https://github.com/ErnakovichLab/dada2_ernakovichlab/tree/split_for_premise

# Dereplicate forward reads
derepFs <- derepFastq(filtFs)
names(derepFs) <- sample.names

# Uncomment line below to save object.
#saveRDS(derepFs, file="derepFs.rds")

# Dereplicate reverse reads
derepRs <- derepFastq(filtRs)
names(derepRs) <- sample.names

# Uncomment line below to save object.
#saveRDS(derepRs, file="derepRs.rds")

#########################################################################
#### Make sequence frequency table using standard error model (errF) ####
#########################################################################

# Infer sequences for forward reads
dadaFs <- dada(derepFs, err = errF, multithread = numCPUs)
names(dadaFs) <- sample.names

# Uncomment line below to save object.
#saveRDS(dadaFs, file="dadaFs.rds")

# Infer sequences for reverse reads
dadaRs <- dada(derepRs, err = errR, multithread = numCPUs)
names(dadaRs) <- sample.names

# Uncomment line below to save object.
#saveRDS(dadaRs, file="dadaRs.rds")

# Merge reads together
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)

#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Save table as an R data object file
saveRDS(seqtab, file="seqtab.rds")

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=numCPUs, verbose=TRUE)

#Inspect distribution of sequence lengths:
table(nchar(getSequences(seqtab.nochim)))

# Print percentage of our seqences that were not chimeric.
100*sum(seqtab.nochim)/sum(seqtab)

# Save table as an r data object file
saveRDS(seqtab.nochim, file="seqtab.nochim.rds")

##Track reads through pipeline:
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,
                                                                       getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged",
                     "nonchim")
rownames(track) <- sample.names
head(track)

head(seqtab.nochim)

write.csv(track, file="track.csv")
write.csv(seqtab.nochim, file="seqtab.nochim.csv")

# Export for to QIIME2 for further analysis.
write.table(t(seqtab.nochim), "q2-seqtab-nochim.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
uniquesToFasta(seqtab.nochim, fout='q2-rep-seqs.fna', ids=colnames(seqtab.nochim))

# Uncomment line below to save environment.
#save.image(file='env.RData')

#####################################################################################
#### Make sequence frequency table using modified error model (errF_1 or errF_4) ####
#####################################################################################

# Note: This section can be skipped for MiSeq data where relearning errors was not necessary. For proceeding with NovaSeq data, use either errF_1 or errF_4, which ever error model fit the data better.

# Infer sequences for forward reads
dadaFs_1 <- dada(derepFs, err = errF_1, multithread = numCPUs)
names(dadaFs_1) <- sample.names

# Uncomment line below to save object.
#saveRDS(dadaFs_1, file="dadaFs_1.rds")

# Infer sequences for reverse reads
dadaRs_1 <- dada(derepRs, err = errR_1, multithread = numCPUs)
names(dadaRs_1) <- sample.names

# Uncomment line below to save object.
#saveRDS(dadaRs_1, file="dadaRs_1.rds")

# Merge reads together
mergers_1 <- mergePairs(dadaFs_1, derepFs, dadaRs_1, derepRs)

#Construct sequence table
seqtab_1 <- makeSequenceTable(mergers_1)
dim(seqtab_1)

# Save table as an R data object file
saveRDS(seqtab_1, file="seqtab_1.rds")

# Remove chimeras
seqtab_1.nochim <- removeBimeraDenovo(seqtab_1, method="consensus", multithread=numCPUs, verbose=TRUE)

#Inspect distribution of sequence lengths:
table(nchar(getSequences(seqtab_1.nochim)))

# Print percentage of our seqences that were not chimeric.
100*sum(seqtab_1.nochim)/sum(seqtab_1)

# Save table as an r data object file
saveRDS(seqtab_1.nochim, file="seqtab_1.nochim.rds")

##Track reads through pipeline:
getN_1 <- function(x) sum(getUniques(x))
track_1 <- cbind(out, sapply(dadaFs_1, getN_1), sapply(dadaRs_1, getN_1), sapply(mergers_1,
                                                                                 getN_1), rowSums(seqtab_1.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs_1, getN_1) with getN_1(dadaFs_1)
colnames(track_1) <- c("input", "filtered", "denoisedF", "denoisedR", "merged",
                       "nonchim")
rownames(track_1) <- sample.names
head(track_1)

head(seqtab_1.nochim)

write.csv(track_1, file="track_1.csv")
write.csv(seqtab_1.nochim, file="seqtab_1.nochim.csv")

# Export for to QIIME2 for further analysis.
write.table(t(seqtab_1.nochim), "q2-seqtab_1-nochim.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
uniquesToFasta(seqtab_1.nochim, fout='q2-rep-seqs_1.fna', ids=colnames(seqtab_1.nochim))

# Uncomment line below to save environment.
save.image(file='env.RData')
