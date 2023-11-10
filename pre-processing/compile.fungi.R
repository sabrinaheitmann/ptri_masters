#Compile fungal OTUs into single phyloseq object and collapse to 99% similarity
#Also filter non-target seqs and add taxonomy

#### 10/12/21 - dada2 is not available for this version of R - need to upload new traits_meta for phyloseq object
 
library(tidyverse)
library(magrittr)
library(foreach)
library(dada2)
library(Biostrings)
library(phyloseq)

#make scratch directory
dir.create("output/scratch")

#identify paths to fungal data
in.paths <- list.files("output/dada","seqTab.rds",recursive = T, full.names = T) %>% 
  grep("fungi",.,value = T)

#merge sequence tables
seqTabs <- foreach(tab.in=in.paths) %do% {readRDS(tab.in)}
seqTab <- mergeSequenceTables(tables=seqTabs)

#Create data frame for summary of processing
compile.summary <- data.frame(SampID=rownames(seqTab),
                              denoised=rowSums(seqTab))

#Extract ITS sequences and write to scratch directory
seqs <- getSequences(seqTab) %>% dada2::rc() %>% DNAStringSet
OTU.names.tmp <- as.character(seqs) %>% openssl:::md5()
names(seqs) <- OTU.names.tmp
colnames(seqTab) <- OTU.names.tmp
writeXStringSet(seqs,file="output/scratch/raw.fa",width=600)

# Trim conserved regions with ITSx
# ITSx parameters
ITSx.flags <- paste("-i output/scratch/raw.fa",
                    "-t 'all'",
                    "--preserve T",
                    "--complement F",
                    "--summary T",
                    "--cpu 12",
                    "--save_regions 'ITS2'",
                    "-only_full T",
                    "-o output/scratch/ITSx",
                    "-E 1e-2")
system2("ITSx", args = ITSx.flags)

# Remove OTUs from seqTab that are not in ITSx output (and longer than 75 bp)
seqs.ITS2 <- readDNAStringSet("output/scratch/ITSx.ITS2.fasta") %>% .[.@ranges@width > 75]
seqTab %<>% .[,names(seqs.ITS2)] 
#Update summary file
compile.summary$ITSx <- seqTab %>% rowSums()

#remove host contamination with Bowtie2 - uses a Ptrichocarpa v3 genome assembly pre-processed with bowtie2-build
writeXStringSet(seqs.ITS2,file="output/scratch/ITS2.fasta",width=600)
system("~/miniconda3/bin/bowtie2 --threads 12 -x data/Ptri_genome/CopyOfPtri.v.3.db -f output/scratch/ITS2.fasta --un output/scratch/noHost.fa --al output/scratch/host.fa --sensitive-local")
seqs.noHost <- readDNAStringSet("output/scratch/noHost.fa")
seqTab %<>% .[,names(seqs.noHost)] 
compile.summary$noHost <- seqTab %>% rowSums()

#Collapse any remaining identical sequences
colnames(seqTab) <- as.character(seqs.noHost)
seqTab %<>% collapseNoMismatch()

#convert to phyloseq object
seqs.noHost.collapse <- getSequences(seqTab) %>% DNAStringSet
OTU.names <- as.character(seqs.noHost.collapse) %>% openssl:::md5()
names(seqs.noHost.collapse) <- OTU.names
colnames(seqTab) <- OTU.names
otuTab <- seqTab %>% as.data.frame() %>% otu_table(taxa_are_rows = F)
phy <- phyloseq(otuTab, refseq(seqs.noHost.collapse))

#collapse sequences with >= 99% similarity
cluster <- function(phy.in,method="single",dissimilarity=0.01){
  require(DECIPHER)
  clusters <- DistanceMatrix(refseq(phy.in), includeTerminalGaps = T, processors=NULL) %>%
    IdClusters(method=method, cutoff=dissimilarity, processors=NULL) 
  clusters[,1] %<>% as.character()
  tax_table(phy.in) <- clusters %>% as.matrix %>% tax_table
  phy.in %<>% speedyseq::tax_glom("cluster")
  tax_table(phy.in) <- NULL
  return(phy.in)
}  
phy %<>% cluster  

#Remove non-fungal sequences by searching against the UNITE database - all eukaryotes including singletons (manually removed entries with no Kingdom level assignment).
writeXStringSet(refseq(phy), file="output/scratch/noHost.fasta",width=600)
vsearch.flags <- paste("--usearch_global output/scratch/noHost.fasta",
                       "--db data/taxonomy_db/sh_general_release_dynamic_s_all_04.02.2020.fasta",
                       "--id 0.50",
                       "--userout output/scratch/UNITEmatches.txt",
                       "--userfields query+target+id",
                       "--notmatched output/scratch/noMatch.fasta",
                       "--maxhits 5",
                       "--maxaccepts 500 --maxrejects 0"
)
system2("vsearch", args = vsearch.flags)
# Read in search results and calculate the proportion of top matches that are fungal
UNITEmatches <- read.table("output/scratch/UNITEmatches.txt",as.is=T) %>%
  group_by(V1) %>% 
  dplyr::summarise(isFungi= ifelse(sum(!grepl("k__Fungi",V2))==0,1,
                               sum(grepl("k__Fungi",V2))/sum(!grepl("k__Fungi",V2))))
# Make list of OTUs to remove from the data
fungal <- UNITEmatches$V1[UNITEmatches$isFungi>0.5] %>% 
  c(readDNAStringSet("output/scratch/noMatch.fasta") %>% names)

# Prune phyloseq object
phy %>% prune_taxa(fungal,.)

#Predict taxonomy with UNITE fungal database
taxa <- assignTaxonomy(refseq(phy),"data/taxonomy_db/sh_general_release_dynamic_04.02.2020.fasta", multithread = T)
rownames(taxa) %<>% names

#Add taxonomy to phyloseq object
tax_table(phy) <- taxa %>% as.matrix %>% tax_table

#Add meta data
#read in sample data
meta <- read.csv("data/gardenData/trait_full_meta.csv",as.is=T)

#convert meta to phyloseq sample data
meta %<>% mutate(SampID=substring(Sample_ID,10)) %>%
  column_to_rownames("SampID") %>%
  sample_data

#add meta data to phyloseq object
phy %<>% merge_phyloseq(meta)

#write phyloseq object to output
saveRDS(phy,"output/compiled/phy.rds")
  
#delete temp files
unlink("output/scratch",recursive=T)


