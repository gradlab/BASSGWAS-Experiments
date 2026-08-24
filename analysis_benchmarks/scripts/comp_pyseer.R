library(tidyverse)

get_files <- function(patt,dir)
{
    pfiles <- list.files(path=dir, pattern = patt)
    dfs <- list()
    for (f in pfiles)
    {
        s <- str_extract(f, "[AR]\\d+_\\d+")
        rep_idx <- str_extract(s, "[AR]\\d+")
        run_idx <- str_extract(s,paste0("(?<=_)\\d+"))
        dat <- vroom::vroom(paste0(dir,"/",f), show_col_types = FALSE, progress = FALSE, altrep = FALSE)
        if(nrow(dat) > 1)
        {
            dat <- dat %>% mutate(run = run_idx, 
                    replicate = rep_idx)
                    
            dfs[[length(dfs)+1]] <- dat
        }
    }
    do.call(rbind,dfs)
}

true_hits <- "(mexA)|(mexB)|(acrR)|(mtrE)"
genehits <- get_files("gene_hits.*\\.txt","azi8_pyseer_genehits") %>%
    mutate(gene=str_extract(gene, "[^\"]+"))

gs1 <- genehits %>% 
    mutate(gene = ifelse(str_detect(gene, true_hits), "mtrCDE", gene)) %>%
    select(gene,replicate,run) %>%
    distinct(gene,replicate,run) %>%
    mutate(gene = ifelse(gene=="mtrCDE", "mtrCDE", "Other")) %>%
    group_by(run,gene,replicate) %>% 
    summarise(ng=n()) %>%
    ungroup()%>%
    mutate(is_hit=gene=="mtrCDE") %>% 
    group_by(run,replicate) %>% 
    summarise(n_pos=sum(is_hit*ng),n_neg=sum((!is_hit)*ng))%>%
    ungroup() %>%
    complete(run,replicate,fill=list(n_pos=0,n_neg=0))%>%
    group_by(run,replicate) %>% 
    summarise(sens=n_pos, spec=n_neg/(n_pos+n_neg)) %>%
    mutate(spec=ifelse(is.nan(spec), NA,spec))%>%
    mutate(type=str_extract(replicate,"[AR]"))%>%
    group_by(type,run) %>%
    summarise(TPR=mean(sens),FPR=mean(spec,na.rm=T))



true_hits <- "(gyrA)|(parC)"
genehits <- get_files("gene_hits.*\\.txt","cip8_pyseer_genehits") %>%
    mutate(gene=str_extract(gene, "[^\"]+"))

gs2 <- genehits %>% 
    select(gene,replicate,run) %>%
    distinct(gene,replicate,run) %>%
    mutate(gene = ifelse(str_detect(gene,true_hits), gene, "Other")) %>%
    group_by(run,gene,replicate) %>% 
    summarise(ng=n()) %>%
    ungroup()%>%
    mutate(is_hit=gene!="Other") %>% 
    group_by(run,replicate) %>% 
    summarise(n_pos=sum(is_hit*ng),n_neg=sum((!is_hit)*ng))%>%
    ungroup() %>%
    complete(run,replicate,fill=list(n_pos=0,n_neg=0))%>%
    group_by(run,replicate) %>% 
    summarise(sens=n_pos/2, spec=n_neg/(n_pos+n_neg)) %>%
    mutate(spec=ifelse(is.nan(spec), NA,spec))%>%
    mutate(type=str_extract(replicate,"[AR]"))%>%
    group_by(type,run) %>%
    summarise(TPR=mean(sens),FPR=mean(spec,na.rm=T))