
suppressPackageStartupMessages(library(GenomeInfoDb))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(Rsamtools))
library(yaml)
library(readr)

config <- read_yaml("config.yaml")

# params from workflow
ref_fasta <- config$ref$sequence
outfile <- "grouped_contigs.tsv"

# make GenomicRanges from the genome
ref_gr <- as(seqinfo(FaFile(ref_fasta)), "GRanges")
all_chroms <- seqlevels(ref_gr)
std_chroms <- standardChromosomes(ref_gr)
unplaced_contigs <- all_chroms[which(!all_chroms %in% std_chroms)]

df <- data.frame(name=standardChromosomes(ref_gr), contigs=standardChromosomes(ref_gr))

if (length(unplaced_contigs) > 0){
    df <- rbind(df, data.frame(name="unplaced_contigs", contigs=paste(unplaced_contigs, collapse=",")))
}

df %>% write_tsv(., outfile)
