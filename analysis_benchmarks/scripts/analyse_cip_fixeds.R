library(tidyverse)
library(ggplot2)
library(rtracklayer)
library(viridis)
library(RColorBrewer)
library(patchwork)
#dir
source("../plotting/plotting_funcs.R")
adir <- "../plotting/annotations"
fdir <- "mapping/cipfs_full"

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

genes <- "(parC)|(gyrA)"
gene_names_ann <- c("gyrA", "parC")

regions_whoN <- gff_df %>%
    filter(str_detect(seqid, "WHO_N")) %>%
    mutate(hit = str_detect(gene, genes)) %>%
    mutate(hit = ifelse(is.na(hit), FALSE, hit)) %>%
    mutate(gene=str_extract(gene,"[^\"]+"))%>%
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
    arrange(-PIP)

hit_hilight <- regions_whoN %>% 
    filter(hit) %>%
    select(gene, ID, start,end, contig) %>% 
    distinct

hit_hilight<-hit_hilight %>%
    group_by(contig,ID, gene)%>%
    summarise(pos = trunc(min(start)+(max(end)-min(start))/2))

mapping_pos_full <- mapping_pos_full %>% rename(pincl="PIP")

pfull_ann <- plot_ann(mapping_ann_full, minp=.25, hits=gene_names_ann)
ggsave("plots/cipfs_ann_full.pdf",pfull_ann,width=180/3, height=228/4, units="mm",device=cairo_pdf)

pfull_map <- plot_pos_WHON(mapping_pos_full, 
    hit_names=hit_hilight$gene, 
    hit_pos=hit_hilight$pos, 
    hit_contig=hit_hilight$contig %>% as.character)

ggsave("plots/cipfs_map_full.pdf",pfull_map,width=2*180/3, height=228/4, units="mm",device=cairo_pdf)
ggsave("plots/cipfs_comb_full.pdf",((pfull_ann+theme(axis.title.y = element_text(size=8),axis.title.x = element_text(size=8)))|pfull_map+
    theme(axis.text.y = element_blank(),axis.title.y = element_blank(),axis.title.x = element_text(size=8)))+plot_layout(widths = c(3, 5)),width=180, height=228/6, units="mm",device=cairo_pdf)

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

dir <- paste0("mapping/cip_fixedsig")

posdf <- get_files(paste0(".bwd",bwd,".position.summaries.csv"),dir)
anndf <- get_files(paste0(".annotation.summaries.csv"),dir) %>% 
    rename(locus="gene")

    adf <- anndf %>% 
        mutate(gene=ifelse(str_detect(gene,genes),  gene, "Off-Target")) %>%
        group_by(gene,replicate,run)%>%
        summarise(pip=max(pincl)) %>%
        ungroup() 
    adf %>% write_csv(paste0("summary_data/cipfs_bench.csv"))


    plt_ann <- adf %>%
        plot_heatmaps_ann
    ggsave(paste0("plots/cipfs_ann_plot.pdf"),plt_ann,width=3*180/4, height=228/6, units="mm", device=cairo_pdf)
    
    plt_pow <- adf %>% plot_power
    ggsave(paste0("plots/cipfs_pow_plot.pdf"),plt_pow, width=180/2, height=228/6, units="mm", device=cairo_pdf)
