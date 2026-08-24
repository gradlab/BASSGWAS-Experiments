library(tidyverse)
library(ape)
args <- commandArgs(TRUE)
#args <- c("test", "test/annotated_hits.csv", "test/sequences.csv", "test")

#args <- c("../ciplow/a1/RUN_2/", "../../data/lab_cip_ctr/sequences.csv", "references.txt", "anno.csv")
ifile <- args[1]
ddir <- args[2]
odir <- args[3]
reffile <-args[4]

data_in <- vroom::vroom(ifile, show_col_types = FALSE, progress = FALSE, altrep = FALSE)
unitigs_all <- vroom::vroom(paste0(ddir, "/unitigs.rtab"), show_col_types = FALSE, progress = FALSE, altrep = FALSE)
tree <- read.tree(paste0(ddir, "/tree.tre"))

tfile <- paste0(odir, "/tree.tre")
ufile <- paste0(odir, "/unitigs.rtab")
cfile <- paste0(odir,"/phylogeny_similarity.tsv")
dfile <- paste0(odir, "/data.tsv")
pfile <- paste0(odir, "/kmer_patterns.txt")
ofile <- paste0(odir, "/hit_kmers.txt")
sigfile <- paste0(odir, "/significant_kmers.txt")
afile <- paste0(odir, "/annotated_kmers.txt")
afile2 <- paste0(odir, "/annotated_kmers_chisqfilt.txt")
sfile <- paste0(odir, "/gene_hits.txt")
sfile2 <- paste0(odir, "/gene_hits_chisqfilt.txt")

ids <- data_in[[1]]
unitigs_all<-unitigs_all %>% 
    select(all_of(c("Unitig_sequence",ids)))
rs <- apply(unitigs_all[,-1], 1, sum)
to_keep <- which((rs>0) & (rs<nrow(unitigs_all))) 

unitigs_all[to_keep,] %>%
    write_tsv(ufile)

unitigs_all<- NULL
tree <- keep.tip(tree, ids)

write.tree(tree,tfile)
system(paste0("python scripts/phylogeny_distance.py --lmm ", tfile, " > ", cfile))
data_in[[2]] <- as.numeric(data_in[[2]])
write_tsv(data_in, dfile)

system(paste0("pyseer --lmm --phenotypes ", dfile, 
    " --pres ", ufile,
    " --similarity ", cfile, 
    " --output-patterns ", pfile,
    " --min-af 0.01",
    " --max-af 0.99",
    " --cpu 8 > ", ofile))
countout <- system(paste0("python scripts/count_patterns.py ", pfile),intern=T)
thresh <- str_extract(countout[2], "\\d+([.]{0,1}\\d+){0,1}E[-]\\d+")
system(paste("bash -c", shQuote(paste0("cat <(head -1 ",ofile,") <(awk '$4<",thresh," {print $0}' ",ofile,") > ",sigfile))))

sighits <- read_delim(sigfile) 
annotations<-read_delim(paste0(ddir,"/annotated.txt"),col_names=c("variant","ann"))
sighits %>% left_join(annotations) %>%
        write_tsv(afile, col_names=NA)
sighits %>% left_join(annotations) %>% 
	filter(is.na(notes) | !str_detect(notes,"bad-chisq")) %>% 
	write_tsv(afile2, col_names=NA)

#system(paste0("annotate_hits_pyseer ",sigfile," ",reffile," ",afile))
system(paste0("python scripts/summarise_annotations.py ",afile," > ",sfile))
system(paste0("python scripts/summarise_annotations.py ",afile2," > ",sfile2))

#system("rm remaining_kmers.fa")
#system("rm remaining_kmers.txt")
