library(tidyverse)
library(ggplot2)
library(rtracklayer)
library(viridis)
library(RColorBrewer)
library(patchwork)
#dir
source("../plotting/plotting_funcs.R")
adir <- "../plotting/annotations"

freq <- "10"
fdir <- paste0("mapping/tet",freq,"_full_pips")

afiles <- list.files(path=adir, pattern = paste0("\\",".gff","$"))
bwd <- 1000

gdfs <- list()
for (af in afiles)
{
    gdfs[[length(gdfs)+1]] <- readGFF(paste0(adir,"/",af)) %>% as_tibble()
}

cns <- colnames(gdfs[[1]])
cns <- cns[!(cns %in% c("colour","inference","note","EC_number"))]

gff_df <- do.call(rbind,lapply(gdfs, select, all_of(cns))) %>% 
    unnest_longer(product)

#genes <- "(tetM)|(gyrA)"
genes <- "(tetM)|(pezT)|(WHO_N_pConjugative_00008)|(WHO_N_pConjugative_00010)|(WHO_N_pConjugative_00011)"
gene_names_ann <- c("tetM")
ids <- "(WHO_N_pConjugative_00008)|(WHO_N_pConjugative_00010)|(WHO_N_pConjugative_00011)"


regions_whoN <- gff_df %>%
    filter(str_detect(seqid, "WHO_N")) %>%
    mutate(hit = str_detect(gene, genes) | str_detect(ID, genes)) %>%
    mutate(hit = ifelse(is.na(hit), FALSE, hit)) %>%
    mutate(gene=str_extract(gene,"[^\"]+"))%>%
    mutate(gene=ifelse(is.na(gene)&str_detect(ID,ids), ID, gene)) %>%
    mutate(gene=ifelse(is.na(gene), "Off-Target", gene)) %>%
    rename(seqid="contig") %>%
    select(contig, hit, type, product, gene, start, end, ID) 

hits_whoN <- regions_whoN %>% 
    mutate(lo=(trunc(start/(bwd/2))-1L), hi=(trunc(end/(bwd/2)))) %>%
    mutate(lo=ifelse(lo >= 0, lo, 0)) %>%
    mutate(sp=lapply(1:n(), function(i) (lo[i]):(hi[i]))) %>%
    unnest_longer(sp) %>%
    mutate(sp=as.integer(sp*bwd/2))

mapping_pos_full <- vroom::vroom(paste0(fdir,"/","pips.bwd",bwd,".position.summaries.csv"), 
    show_col_types = FALSE,progress = FALSE, altrep = FALSE) %>% 
    mutate(sp=as.integer(sp))

mapping_ann_full <- vroom::vroom(paste0(fdir,"/","pips.annotation.summaries.csv"), 
    show_col_types = FALSE,progress = FALSE, altrep = FALSE)%>% 
    rename(locus="gene") 

patterns_full <- vroom::vroom(paste0(fdir,"/","pips.patterns.annotated.csv"), 
    show_col_types = FALSE,progress = FALSE, altrep = FALSE)%>% 
    rename(locus="gene") 

canonical_hits <- mapping_pos_full %>% 
    left_join(hits_whoN, by =c("contig","sp")) %>% 
        mutate(hit=ifelse(is.na(hit), FALSE, hit), 
        gene=ifelse(is.na(gene), "Off-Target", gene)) %>%
        group_by(contig, sp) %>%
        mutate(has_hit = any(hit,na.rm=T)) %>%
        group_by(hit, gene, contig) %>%
        mutate(pincl_full = ifelse(hit, max(pincl), NA)) %>%
        select(contig, sp, gene, hit, has_hit, pincl_full) %>%
        distinct(gene, sp, contig, .keep_all = T) 

mapping_ann_full <- mapping_ann_full %>%
    rename(pincl="PIP") %>%
    arrange(-PIP) %>%
    mutate(gene=str_replace(gene, "WHO_N_pConjugative","pConj")) %>%
    mutate(gene=str_replace(gene, "_0+","_"))

hit_hilight <- regions_whoN %>% 
    filter(gene %in% gene_names_ann)%>%
    filter(hit) %>%
    select(gene, ID, start,end, contig) %>% 
    distinct

hit_hilight<-hit_hilight %>%
    group_by(contig,ID, gene)%>%
    summarise(pos = trunc(min(start)+(max(end)-min(start))/2))

mapping_pos_full <- mapping_pos_full %>% rename(pincl="PIP")

pfull_ann <- plot_ann(mapping_ann_full, minp=.25, hits=gene_names_ann)
ggsave(paste0("plots/tet",freq,"_ann_full.pdf"),pfull_ann,width=180/3, height=228/4, units="mm",device=cairo_pdf)

pfull_map <- plot_pos_WHON(mapping_pos_full, 
    hit_names=hit_hilight$gene, 
    hit_pos=hit_hilight$pos, 
    hit_contig=hit_hilight$contig %>% as.character)
ggsave(paste0("plots/tet",freq,"_map_full.pdf"),pfull_map,width=2*180/3, height=228/4, units="mm",device=cairo_pdf)
ggsave(paste0("plots/tet",freq,"_comb_full.pdf"),((pfull_ann+theme(axis.title.y = element_text(size=8),axis.title.x = element_text(size=8)))|pfull_map+
    theme(axis.text.y = element_blank(),axis.title.y = element_blank(),axis.title.x = element_text(size=8)))+plot_layout(widths = c(3, 5)),width=180, height=228/6, units="mm",device=cairo_pdf)

pyseer_hits <- read_delim("mapping/tet10_full_pyseer/annotated_kmers.txt",col_names=NA) %>%
     select(c(X1,X4,X9)) 
colnames(pyseer_hits) <- c("kmer", "pval", "ann")
pyseer_hits<-pyseer_hits %>% 
    mutate(mapping=str_extract_all(ann, "WHO_N[_a-zA-Z]*:\\d+\\-\\d+;")) %>% 
    unnest_longer(mapping) %>% 
    select(pval,mapping) %>% 
    mutate(contig=str_extract(mapping, "[^:]+(?=:)")) %>%
    mutate(hi=str_extract(mapping, "(?<=\\-)\\d+"), lo=str_extract(mapping, "\\d+(?=\\-)")) %>%
    mutate(hi=as.integer(hi),lo=as.integer(lo)) %>% 
    mutate(mid = (hi-lo)/2+lo) %>%
    mutate(logp = -log10(pval))%>%
    mutate(sp=mid) %>%
    select(sp, logp, contig)

p<-plot_pyseer_WHON(pyseer_hits, hit_names=hit_hilight$gene, 
    hit_pos=hit_hilight$pos, 
    hit_contig=hit_hilight$contig %>% as.character)
ggsave("plots/tet10_map_pyseer.pdf",p,width=180, height=228/8, units="mm",device=cairo_pdf)

get_files <- function(patt,dir)
{
    pfiles <- list.files(path=dir, pattern = paste0("\\",patt,"$"))
    dfs <- list()
    for (f in pfiles)
    {
        rep_idx <- str_extract(f, "[AR]\\d+")
        run_idx <- str_extract(f,paste0("(?<=_)\\d+(?=\\",patt,")"))
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


    dir <- paste0("mapping/tet",freq,"_map_all_2")

    posdf <- get_files(paste0(".bwd",bwd,".position.summaries.csv"),dir)
    anndf <- get_files(paste0(".annotation.summaries.csv"),dir) %>% 
        rename(locus="gene")

    posplt_df <- posdf %>% 
        mutate(sp=as.integer(sp)) %>%
        filter(str_detect(contig,"WHO_N")) %>%
        left_join(canonical_hits,relationship = "many-to-many") %>%
        mutate(type=ifelse(str_detect(replicate,"A"), "Adaptive", "Random")) %>%
        mutate(gene=ifelse(hit, gene, "Off-Target"))

    post_plt2 <- posplt_df %>% 
        group_by(gene, replicate, run, type) %>%
        summarise(max_p_incl = max(pincl * (hit | !has_hit)),
            test=all(hit) | !any(hit)) %>%
        ungroup()
    stopifnot(all(post_plt2$test)) 

    plt <- post_plt2 %>% 
        filter((gene%in%gene_names_ann) | !str_detect(gene, genes))%>%
        rename(max_p_incl = "pip")%>%
        plot_heatmaps_pos
    ggsave(paste0("plots/tet",freq,"_map_plot.pdf"),plt, width=2*180/3, height=228/4, units="mm", device=cairo_pdf)

    adf <- anndf %>% 
        filter((gene%in%gene_names_ann) | !str_detect(gene, genes))%>%
        mutate(gene=ifelse(str_detect(gene,genes),  gene, "Off-Target")) %>%
        group_by(gene,replicate,run)%>%
        summarise(pip=max(pincl)) %>%
        ungroup() 

    adf %>% write_csv(paste0("summary_data/tet",freq, "_bench.csv"))

    plt_ann <- adf %>%
        plot_heatmaps_ann
    ggsave(paste0("plots/tet",freq,"_ann_plot.pdf"), plt_ann, width=180/2, height=228/4, units="mm", device=cairo_pdf)
    
    plt_pow <- adf %>% plot_power
    ggsave(paste0("plots/tet",freq,"_pow_plot.pdf"),plt_pow, width=180/2, height=228/3, units="mm", device=cairo_pdf)

