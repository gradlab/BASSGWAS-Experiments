library(tidyverse)
args <- commandArgs(TRUE)
#args <- c("../ciplow/a1/RUN_2/", "../../data/lab_cip_ctr/sequences.csv", "references.txt", "anno.csv")
sfile <- args[1]
reffile <- args[2]
ofile <- args[3]
dir <- args[4]

seqs <- read_csv(sfile)

ln1 <- c("")
lns <- sapply(seqs[["sequence"]], paste0)
lns <- unique(unname(lns))
#for (i in 1:nrow(seqs_f))
#{
#    lns <- c(lns, paste0(seqs_f[i,2]))
#}

writeLines(c(ln1,lns), paste0(dir,"/hits.txt"))
system(paste0("annotate_hits_pyseer ", dir, "/hits.txt ", reffile, " ", dir, "/annotated.txt"))
anno <- read_delim(paste0(dir,"/annotated.txt"),col_names=c("sequence","ann"),delim="\t")
anno <- anno %>% mutate(ann = str_extract_all(anno$ann, "[^;,]+:\\d+-\\d+[^,]+")) %>% unnest_longer(ann)
anno <- anno %>% mutate(contig = str_extract(ann, "[^:]+(?=:)"), range = str_extract(ann, "\\d+-\\d+"), ann = str_extract_all(ann,"(?<=((;)(\")?))[\\w\\.]*(?=(\")?($|[;,]))")) %>%
	unnest_wider(ann, names_sep ="_")
anno %>% write_csv(ofile)

