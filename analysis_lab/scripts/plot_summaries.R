library(tidyverse)
library(ggplot2)
library(rtracklayer)
library(viridis)
library(RColorBrewer)
library(gridExtra)
library(patchwork)

adir <- "../plotting/annotations"
dir <- "mapping/"#"cip8_map_all_2"

hits <- ""
ddir <- ""
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

regions_whoN <- gff_df %>%
    filter(str_detect(seqid, "WHO_N")) %>%
    mutate(gene=str_extract(gene,"[^\"]+"))%>%
    rename(seqid="contig") %>%
    select(contig, type, product, gene, start, end, ID) 

get_files <- function(patt,dir)
{
    pfiles <- list.files(path=dir, pattern = paste0("\\",patt,"$"))
    dfs <- list()
    for (f in pfiles)
    {
        run_idx <- as.integer(str_extract(f,paste0("(?<=_)\\d+(?=\\",patt,")")))
        dat <- vroom::vroom(paste0(dir,"/",f), show_col_types = FALSE, progress = FALSE, altrep = FALSE)
        if(nrow(dat) > 1)
        {
            dat <- dat %>% mutate(run = run_idx)
                    
            dfs[[length(dfs)+1]] <- dat
        }
    }
    do.call(rbind,dfs)
}


posdf <- get_files(paste0(".bwd",bwd,".position.summaries.csv"),dir)
anndf <- get_files(paste0(".annotation.summaries.csv"),dir)
pattdf <- get_files(paste0(".patterns.annotated.csv"),dir)

hit_genes <- c("parC")
anndf <- anndf %>%
    rename(locus="gene")%>%
    rename(pincl="PIP") %>%
    arrange(-PIP)

hit_hilight <- regions_whoN %>% 
    filter(gene=="parC") %>%
    select(gene, ID, start,end, contig) %>% 
    distinct %>%
    group_by(contig,ID, gene)%>%
    summarise(pos = trunc(min(start)+(max(end)-min(start))/2))

posdf<-posdf %>% rename(pincl="PIP")

source("../plotting/plotting_funcs.R")
plts_ann <- list()
for(i in 1:8) {
    plts_ann[[i]] <- plot_ann(anndf %>% filter(run==i), minp=.1, hits=hit_genes)
}

plts_pos <- list()
for(i in 1:8) {
    plts_pos[[i]] <- plot_pos_WHON(posdf %>% filter(run==i), 
    hit_names=hit_hilight$gene, 
    hit_pos=hit_hilight$pos, 
    hit_contig=hit_hilight$contig %>% as.character)
}


ggsave("plots/annotations_experiments.pdf",  
    marrangeGrob(grobs = plts_ann, nrow=4, ncol=2, top=NULL), width=180, height=228, units="mm",device=cairo_pdf)
ggsave("plots/position_experiments.pdf",  
    marrangeGrob(grobs = plts_pos, nrow=8, ncol=1, top=NULL), width=180, height=228, units="mm",device=cairo_pdf)
p <- wrap_plots(plts_pos)+ plot_layout(guides="collect",axes = "collect",ncol = 1)
ggsave("plots/position_experiments_comb.pdf",  
    p, width=180, height=208, units="mm",device=cairo_pdf)


plts_pos_7 <- plot_pos_WHON(posdf %>% filter(run==8,contig=="WHO_N"), 
    hit_names=hit_hilight$gene, 
    hit_pos=hit_hilight$pos, 
    hit_contig=hit_hilight$contig %>% as.character)

ggsave("plots/batch_7.pdf", 
    plts_pos_7 + 
    theme(legend.position = "none", 
            axis.title.y = element_text(size = 6),
            axis.title.x = element_text(size = 6)),
        width=58, height=32, units="mm", device=cairo_pdf)

annotation_hmap <- anndf %>%
    group_by(gene) %>%
    mutate(toinclude=any(PIP>.1)) %>%
    ungroup() %>%
    filter(toinclude) %>%
    mutate(run=factor(run-1,levels=1:max(run)-1)) %>%
    mutate(gene=factor(gene,levels=unique(gene))) %>%
    ggplot(aes(x=gene,y=run,fill=PIP))+
    geom_tile() +
    geom_text(aes(label=scales::label_number(accuracy=.01)(round(PIP,digits=2))),size=4,fontface = "bold", color="gray15")+
    scale_fill_distiller(palette="RdBu",limits=c(0,1),guide = guide_colourbar(title.position = "left"))+
    labs(x="Locus", y="Batch", fill="Posterior Incl. Probability")+
    theme_minimal()+
    theme(aspect.ratio=1,
        axis.text.x = element_text(angle = 60, size = 10, hjust=1,vjust = 1),
        axis.text.y = element_text(size=10),
        legend.title = element_text(angle = 90, size = 12, hjust=.5,vjust = 1))

ggsave("plots/heatmap_experiments.pdf", annotation_hmap, width=180, height=160, units="mm",device=cairo_pdf)