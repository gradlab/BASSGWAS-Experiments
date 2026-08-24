library(tidyverse)
library(ape)

tigs <- read_delim("unitigs.rtab")
tree <- read.tree("tree.tre")
svd_threshold = 0.99
af_threshold = 0

tig_seqs <- tigs$Unitig_sequence
stopifnot(all(tree$tip.label %in% colnames(tigs)))
stopifnot(all(colnames(tigs)[-1] %in% tree$tip.label))

tigs <- tigs %>% relocate(all_of(c("Unitig_sequence",tree$tip.label)))
stopifnot(all(colnames(tigs)[-1] == tree$tip.label))

patts <- t(as.matrix(tigs[,-1]))

#Convert to minor allele frequency
n <- nrow(patts)
signs <- rep(1L, ncol(patts))
for (i in 1:ncol(patts))
{
    v <- patts[, i]
    sv <- sum(v)
    af <- sv/n
    if(af < af_threshold || af > (1-af_threshold))
    {
        patts[, i] <- 0L
    }
    else if (sv > n/2)
    {
	signs[i] <- -1L
        patts[, i] <- 1L-v
    }
}

pattern_ids <- apply(patts, 2, paste0, collapse="")
patts_minor <- as_tibble(t(patts)) %>% 
    mutate(pattern = pattern_ids, 
        sequence = tig_seqs,
	direction = signs)

patts_minor <- patts_minor %>% 
    filter(!if_all(!all_of(c("pattern", "sequence","direction")), ~ .x == 0)) %>%
    mutate(pattern = as.integer(factor(pattern))) 

seqs <- patts_minor %>% 
    select(c("pattern", "sequence","direction"))
patts_minor <- patts_minor %>% 
    select(!all_of(c("sequence","direction"))) %>% 
    distinct(.keep_all=T)

patts_minor <- patts_minor %>% 
    mutate(pattern_new = 1:nrow(patts_minor))
seqs <- seqs %>% 
    left_join(patts_minor %>% select(c("pattern", "pattern_new")), by="pattern") 

sequences_out <- seqs %>% select(c("pattern_new", "sequence","direction")) %>% 
    rename(pattern_id=pattern_new, sequence=sequence)

patterns_out <- patts_minor %>% 
    select(!all_of(c("pattern")))%>%#,"pattern_new"))) %>%
    rename(pattern_id=pattern_new) %>% 
    relocate("pattern_id") %>%
    as.matrix() %>% t()

colnames(patterns_out) <- paste0(patterns_out[1, ])
patterns_out <- patterns_out[-1, ]
patterns_ctr <- patterns_out

for (i in 1:ncol(patterns_ctr))
{
	patterns_ctr[,i] <- patterns_ctr[,i] - mean(patterns_ctr[,i])
}

popVCV <- vcv(tree, corr=T)

SV <- svd(popVCV)
ssum <- sum(SV$d)
truncpt <- which(cumsum(SV$d)/ssum > svd_threshold)[1]

U = SV$u[,1:truncpt]
S = SV$d[1:truncpt]

write.csv(patterns_out, "patterns.csv")
write.csv(patterns_ctr, "patterns_centred.csv")
write.csv(sequences_out, "sequences.csv",row.names = FALSE)
write.csv(U, "U.csv",row.names = FALSE)
write.csv(S, "S.csv",row.names = FALSE)
