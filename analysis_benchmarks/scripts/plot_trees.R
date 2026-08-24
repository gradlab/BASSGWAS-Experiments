library(tidyverse)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(viridis)
library(RColorBrewer)
library(patchwork)
library(ape)

plot_tree <- function(tr, dat)
{
    dat2<-dat[,2]
    dat2 <- dat2 %>% as.data.frame()
    rownames(dat2) <- dat[[1]]
    p <- ggtree(tr %>% makeNodeLabel,linewidth=0.1,layout="circular")
    p <- gheatmap(p, dat2, offset=rel(8), width=0.25, colnames=F, legend_title="MIC \u2265 Cutoff",color=NA)  +
        scale_fill_brewer(palette = "Dark2",name="MIC \u2265 Cutoff")
    p
}

tr_azi <- read.tree("../data/azi_mic_2/tree.tre")
dat_azi <- read_delim("../data/azi_mic_2/data_bin.csv")
p1 <- plot_tree(tr_azi,dat_azi) +  ggtitle('Azithromycin (MIC \u2265 2), AF=0.1')# %>% 
#ggsave(filename = "plots/tree_azi.pdf", width=180/4, height=228/8, units="mm",device=cairo_pdf)

tr_cip <- read.tree("../data/cip/tree.tre")
dat_cip <- read_delim("../data/cip/data_bin.csv")
p2 <- plot_tree(tr_cip,dat_cip) +ggtitle('Ciprofloxacin (MIC \u2265 1), AF=0.1')

tr_tet10 <- read.tree("../data/tet10/tree.tre")
dat_tet10 <- read_delim("../data/tet10/data_bin.csv")
p3 <- plot_tree(tr_tet10,dat_tet10)+  ggtitle('Tetracycline (MIC \u2265 8), AF=0.1')

tr_tet2.5 <- read.tree("../data/tet2.5/tree.tre")
dat_tet2.5 <- read_delim("../data/tet2.5/data_bin.csv")
p4 <- plot_tree(tr_tet2.5,dat_tet2.5) + ggtitle('Tetracycline (MIC \u2265 8), AF=0.025')

p_comb <- ((p3 | p4) / (p2|p1))+ plot_layout(guides="collect")
ggsave("plots/trees_bench.pdf",p_comb,width=180, height=160, units="mm",device=cairo_pdf)
