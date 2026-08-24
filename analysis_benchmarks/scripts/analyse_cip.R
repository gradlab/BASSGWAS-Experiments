library(tidyverse)
library(ggplot2)
library(rtracklayer)
library(viridis)
library(RColorBrewer)
library(patchwork)
#dir
source("../plotting/plotting_funcs.R")
adir <- "../plotting/annotations"
fdir <- "mapping/cip_full_pips"

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
ggsave("plots/cip_ann_full.pdf",pfull_ann,width=180/3, height=228/4, units="mm",device=cairo_pdf)

pfull_map <- plot_pos_WHON(mapping_pos_full, 
    hit_names=hit_hilight$gene, 
    hit_pos=hit_hilight$pos, 
    hit_contig=hit_hilight$contig %>% as.character)

ggsave("plots/cip_map_full.pdf",pfull_map,width=2*180/3, height=228/4, units="mm",device=cairo_pdf)
ggsave("plots/cip_comb_full.pdf",((pfull_ann+theme(axis.title.y = element_text(size=8),axis.title.x = element_text(size=8)))|pfull_map+
    theme(axis.text.y = element_blank(),axis.title.y = element_blank(),axis.title.x = element_text(size=8)))+plot_layout(widths = c(3, 5)),width=180, height=228/6, units="mm",device=cairo_pdf)

pyseer_hits <- read_delim("mapping/cip_full_pyseer/annotated_kmers.txt",col_names=NA) %>%
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
ggsave("plots/cip_map_pyseer.pdf",p,width=180, height=228/8, units="mm",device=cairo_pdf)

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

batches <- c(4,8,16)
ann_plts <-list()
for (bsz in batches)
{
    dir <- paste0("mapping/cip",bsz,"_map_all_2")

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
        rename(max_p_incl = "pip")%>%
        plot_heatmaps_pos
    ggsave(paste0("plots/cip",bsz,"_map_plot.pdf"),plt,width=3*180/4, height=228/4, units="mm", device=cairo_pdf)

    adf <- anndf %>% 
        mutate(gene=ifelse(str_detect(gene,genes),  gene, "Off-Target")) %>%
        group_by(gene,replicate,run)%>%
        summarise(pip=max(pincl)) %>%
        ungroup() 
    adf %>% write_csv(paste0("summary_data/cip",bsz, "_bench.csv"))


    plt_ann <- adf %>%
        plot_heatmaps_ann
    ggsave(paste0("plots/cip",bsz,"_ann_plot.pdf"),plt_ann,width=3*180/4, height=228/4, units="mm", device=cairo_pdf)
    
    ann_plts[[length(ann_plts)+1]] <- plt_ann

    plt_pow <- adf %>% plot_power
    ggsave(paste0("plots/cip",bsz,"_pow_plot.pdf"),plt_pow, width=180/2, height=228/3, units="mm", device=cairo_pdf)
}

dfs_ann <- list()
dfs_pos <- list()
plt_4_16 <- ((ann_plts[[1]]) / (ann_plts[[3]])) + plot_layout(guides="collect")
ggsave("plots/cip_4_16_ann_plot.pdf",plt_4_16, width=130, height=140, units="mm", device=cairo_pdf)

for (bsz in batches)
{
    dir <- paste0("mapping/cip",bsz,"_map_all_2")

    dfs_pos[[length(dfs_pos)+1]] <- get_files(paste0(".bwd",bwd,".position.summaries.csv"),dir) %>%
        mutate(batch_sz=bsz)
    dfs_ann[[length(dfs_ann)+1]]  <- get_files(paste0(".annotation.summaries.csv"),dir) %>% 
        rename(locus="gene") %>%
        mutate(batch_sz=bsz)
}

ann_df_all <- do.call(rbind, dfs_ann) 
ann_sum <-ann_df_all %>% 
    filter(str_detect(gene,genes)) %>%
        group_by(gene,replicate,run,batch_sz)%>%
        summarise(pip=max(pincl)) %>%
        ungroup() %>%
        mutate(type=ifelse(str_detect(replicate,"A"), "Adaptive", "Random")) %>%
        mutate(PIP = lapply(1:n(), function(x) c(.5,.7,.9)))%>%
        unnest_longer(PIP) %>%
        group_by(gene, run, PIP,type,batch_sz) %>%
        summarise(num_replicates = 100*sum(pip >= PIP)/n(),.keep_all=T) %>%
        mutate(n_samples=as.integer(0+as.integer(run))*batch_sz) 
p <- ann_sum %>%
    filter(type=="Adaptive")%>%
    filter(n_samples %in% as.integer(seq(from=4,to=124,by=4))) %>%
    mutate(PIP=factor(PIP, levels=c(.5,.7,.9), labels=paste0("PIP \u2265 ",c(.5,.7,.9))))%>%
    mutate(batch_sz= factor(batch_sz)) %>%
    ggplot(aes(x=n_samples, y=num_replicates, color=batch_sz,shape=batch_sz)) +
        geom_line(linewidth=.5) +
        geom_point(alpha=.8)+
        scale_color_brewer(palette="Dark2")+
        scale_shape_manual(values=c(5,2,4))+
        scale_y_continuous(breaks = 100*seq(0,1,length.out=11), 
            labels=scales::percent_format(scale = 1),
            expand = c(0,2.5))+
        scale_x_continuous(breaks=seq(from=0,to=116,by=8))+
        facet_grid(cols=vars(gene),rows=vars(PIP)) +
        labs(x="Additional Strains",y="Power", color="Batch Size", shape="Batch Size")+
        theme_bw()+
        theme(
            aspect.ratio=1,
            text=element_text(size=8),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, size=6),
            legend.title = element_text(angle = 0, size = 8, hjust=.5, vjust = 1))
ggsave("plots/comp_batch.pdf",p,width=180, height=180, units="mm", device=cairo_pdf)
        
