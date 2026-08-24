library(tidyverse)
args <- commandArgs(TRUE)
#args <- c("test", "test/annotated_hits.csv", "test/sequences.csv", "test")

#args <- c("../ciplow/a1/RUN_2/", "../../data/lab_cip_ctr/sequences.csv", "references.txt", "anno.csv")
dir <- args[1]
afile <- args[2]
sfile <- args[3]
ofile <- args[4]
Sys.setenv(VROOM_CONNECTION_SIZE = 1500072)

pips <- vroom::vroom(paste0(dir, "/pips_draws.csv"), show_col_types = FALSE, progress = FALSE, altrep = FALSE)
pars <- vroom::vroom(paste0(dir, "/cpars_draws.csv"), show_col_types = FALSE, progress = FALSE, altrep = FALSE)
pars <- pars %>% mutate(Z=Z/sum(Z),iteration=row_number())
nr <- nrow(pars)
pipsM <- as.matrix(pips[,3:ncol(pips)])
pips <- NULL
gc()

pipsE <- rep(0, ncol(pipsM))
names(pipsE) <- colnames(pipsM)

for (i in 1:ncol(pipsM))
{
    pipsE[i] <- sum(pipsM[,i]*pars$Z)
}

pipsM <- NULL
gc()
seqs <- vroom::vroom(sfile, show_col_types = FALSE, progress = FALSE, altrep = FALSE)

pipsS <- data.frame(pattern_id = as.numeric(str_extract(names(pipsE), "\\d+")), pip = pipsE)
#seqs <- seqs %>% right_join(pipsS,by="pattern_id") #seqs %>% filter(pattern_id %in% str_extract(nmret, "\\d+"))
anno <- vroom::vroom(afile,show_col_types = FALSE) %>% 
    pivot_longer(c(ann_1,ann_2,ann_3),values_to="locus") %>% 
    select(!c(name)) %>% 
    distinct() %>%
    mutate(locus = ifelse(is.na(locus), "unannotated", locus)) %>%
    mutate(locus2 = 
        lapply(locus, 
            function(x) if(x%in%c("mexB","mexA","acrR","mtrE")) {
                list("mtrCDE", x)} else 
                {list(x)})
    ) %>%
    unnest_longer(locus2) %>% 
    mutate(locus=locus2) %>%
    select(!locus2) %>%
    distinct() %>%
    mutate(down=as.integer(str_extract(range, "^\\d+")), 
        up=as.integer(str_extract(range, "\\d+$")))
seqs <- seqs %>% left_join(anno,by="sequence", relationship = "many-to-many") %>%
    mutate(locus = ifelse(is.na(contig), "unmapped", locus))
anno <- NULL
seqs_with_anno <- seqs %>% 
    ungroup() %>%
        select(pattern_id, contig, locus) %>%
        distinct()
#seqs <- NULL
gc()

annopips <- seqs_with_anno %>% 
    right_join(pipsS,by="pattern_id") 
annopips %>%
    write_csv(paste0(ofile,".patterns.annotated.csv"))

annopips <- NULL
seqs_with_anno <- NULL
gc()

seqs4anno <- seqs %>% 
    ungroup() %>%
    select(pattern_id, locus) %>%
    distinct() %>%
    filter(!(locus %in% c("unannotated", "unmapped")))

seqs4pos <- seqs %>% 
    ungroup() %>%
    filter(!is.na(range)) %>%
    filter(str_detect(contig, "WHO_N"))%>%
    select(pattern_id, contig, down, up) %>%
    distinct()

gammas <- vroom::vroom(paste0(dir, "/gammas_draws.csv"), show_col_types = FALSE, progress = FALSE, altrep = FALSE)
gammas <- gammas[,3:ncol(gammas)] > 0 
stopifnot(nrow(gammas)==nr)

annopipcombs <- seqs4anno %>% 
	    group_by(locus) %>%
   	    summarise(pincl = sum(pars$Z*apply(gammas[,pattern_id,drop=F],1,any))) %>%
        ungroup()
annopipcombs %>% 
    write_csv(paste0(ofile,".annotation.summaries.csv"))
annopipcombs <- NULL
gc()

bwds <- c(500L,1000L,10000L)
for (bwd in bwds)
{
    seqs4pos.2 <- seqs4pos %>% 
        mutate(down2=trunc(down/(bwd/2))-1,up2=trunc(up/(bwd/2))) %>%
        mutate(down2=ifelse(down2 >= 0, down2, 0)) %>%
        mutate(sp=lapply(1:n(), function(i) down2[i]:up2[i])) %>%
        unnest_longer(sp) %>%
        group_by(contig, pattern_id) %>% 
        distinct(sp, .keep_all=T) %>%
        mutate(sp=sp*bwd/2) %>%
	    ungroup() %>%
        select(contig,pattern_id,sp)

    pospipsums <- seqs4pos.2 %>% 
        right_join(pipsS,by="pattern_id") %>%
	    group_by(sp, contig) %>%
   	    summarise(pip_sum=sum(pip))
    
    pospipcombs <- seqs4pos.2 %>% 
	    group_by(sp, contig) %>%
   	    summarise(pincl = sum(pars$Z*apply(gammas[,pattern_id,drop=F],1,any))) %>%
        ungroup()
    
    pospipcombs <- pospipcombs %>% 
        full_join(pospipsums, by=c("sp","contig")) 
    pospipcombs %>%
        write_csv(paste0(ofile,".bwd",bwd,".position.summaries.csv"))
    gc()
}

