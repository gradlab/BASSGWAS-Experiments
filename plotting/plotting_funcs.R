plot_ann <- function(ann_df, minp, hits=c())
{
    ann_df <- ann_df %>% 
        mutate(hit=(gene%in%hits))%>%
        filter(PIP >= minp) %>%
        arrange(-PIP, gene) %>%
        mutate(gene = factor(gene, levels=gene, ordered=T))

        ggplot(ann_df, aes(x=gene,y=PIP))+#,fill=hit))+
            geom_col() +
            theme_minimal() +
            labs(x="Locus", y="Posterior Incl. Probability")+
            #scale_fill_manual(values = 
            #    c("#276FBF", "gray35"),
            #    breaks=c(T,F),
            #)+
            scale_y_continuous(breaks = seq(from=.1, to=1, by=.1),limits=c(0,1.02), expand=c(0,0))+
            theme(legend.position = "none",
                axis.text.x = element_text(size = 6, hjust=0.5, vjust = 1),
                axis.text.y =  element_text(size = 6),
                axis.title.x = element_text(size = 8),
                axis.title.y =  element_text(size = 8,hjust=.5,vjust=.5))
}

#pconj only
plot_pos_WHON <- function(pos_df, hit_names=c(), hit_pos=c(), hit_contig=c())
{
    WHO_N_len <- 2172826
    WHO_N_pConjugative_len <- 42004
    WHO_N_pBlaTEM_len <- 7449  
    WHO_N_pCryptic_len <- 4207

    offs <- c(WHO_N=0, 
        WHO_N_pConjugative=WHO_N_len)#, 
        #WHO_N_pBlaTEM=WHO_N_len+WHO_N_pConjugative_len,
        #WHO_N_pCryptic=WHO_N_len+WHO_N_pConjugative_len+WHO_N_pBlaTEM_len)
    
    xlims <- c(0, WHO_N_len+WHO_N_pConjugative_len)#+WHO_N_pBlaTEM_len+WHO_N_pCryptic_len)
    xlabs <- c("WHO_N", "WHO_N_pConj")#,"WHO_N_pBlaTEM","WHO_N_pCryptic")
    xbreaks <- c(offs[1] + WHO_N_len/2, 
        offs[2] + WHO_N_pConjugative_len/2, 
        offs[3] + WHO_N_pBlaTEM_len/2,
        offs[4] + WHO_N_pCryptic_len/2)

    contigs_inord <- c("WHO_N", "WHO_N_pConjugative")# "WHO_N_pBlaTEM", "WHO_N_pConjugative","WHO_N_pCryptic")

    annotations <- data.frame(hit_names=hit_names, sp=hit_pos, contig=hit_contig) %>%
        mutate(sp=sp+offs[contig])

    pos_df %>% 
        filter(!is.na(sp)) %>%
        filter(contig%in%c("WHO_N","WHO_N_pConjugative")) %>%
        mutate(sp=sp+offs[contig]) %>%
        mutate(contig=factor(contig,levels=contigs_inord)) %>%
        ggplot(aes(x=sp,y=PIP,color=contig,fill=contig)) +
        geom_col(alpha=.5) +
        geom_vline(data=annotations, aes(xintercept=sp), 
            color="darkred", linetype="dashed", alpha=.5, linewidth=.5) +
        geom_text(data=annotations, aes(x=sp, y=.9, label=hit_names),
            hjust=1, nudge_x=-xlims[2]*0.005, vjust=0, size=2.0, angle=90, color="darkred") +
        scale_y_continuous(
            breaks = seq(from=.1, to=1, by=.1),
            expand = c(0,0), 
            limits = c(0,1.02)
        ) +
        scale_x_continuous(expand=c(0,0),limits=c(xlims[1]-10, xlims[2]+10),breaks=seq(from=0,to=2172826, by =.5*1e6),labels=seq(from=0,to=2172826, by =.5*1e6)/1e3)+
        scale_color_brewer(palette="Dark2")+
        scale_fill_brewer(palette="Dark2")+
        labs(x="Postion (kb)", y="Posterior Incl. Probability",fill="Chromosome",color="Chromosome")+
        theme_minimal() +
        theme(
            legend.title = element_text(size=8),
            legend.text = element_text(size=6),
            panel.grid.major.x = element_blank(),
            #panel.grid.minor.x = element_blank(),
            axis.text.y =  element_text(size = 6),
            axis.text.x = element_text(size=6),
            axis.title.x=element_text(size=8),
            axis.title.y=element_text(size=8)
        )
}


plot_pyseer_WHON <- function(pos_df, hit_names=c(), hit_pos=c(), hit_contig=c())
{
    WHO_N_len <- 2172826
    WHO_N_pConjugative_len <- 42004
    WHO_N_pBlaTEM_len <- 7449  
    WHO_N_pCryptic_len <- 4207

    offs <- c(WHO_N=0, 
        WHO_N_pConjugative=WHO_N_len)#, 
        #WHO_N_pBlaTEM=WHO_N_len+WHO_N_pConjugative_len,
        #WHO_N_pCryptic=WHO_N_len+WHO_N_pConjugative_len+WHO_N_pBlaTEM_len)
    
    xlims <- c(0, WHO_N_len+WHO_N_pConjugative_len)#+WHO_N_pBlaTEM_len+WHO_N_pCryptic_len)
    xlabs <- c("WHO_N", "WHO_N_pConj")#,"WHO_N_pBlaTEM","WHO_N_pCryptic")
    xbreaks <- c(offs[1] + WHO_N_len, 
        offs[2] + WHO_N_pConjugative_len)#, 
        #offs[3] + WHO_N_pBlaTEM_len,
        #offs[4] + WHO_N_pCryptic_len)

    contigs_inord <- c("WHO_N", "WHO_N_pConjugative")

    annotations <- data.frame(hit_names=hit_names, sp=hit_pos, contig=hit_contig) %>%
        mutate(sp=sp+offs[contig])

    pos_df %>% 
        filter(!is.na(sp)) %>%
        filter(contig%in%c("WHO_N","WHO_N_pConjugative")) %>%
        mutate(sp=sp+offs[contig]) %>%
        mutate(contig=factor(contig,levels=contigs_inord,labels=contigs_inord)) %>%
        ggplot(aes(x=sp,y=logp,color=contig)) +
        geom_text(data=annotations, aes(x=sp, y=.9*max(pos_df$logp), label=hit_names),
            hjust=1, nudge_x=-xlims[2]*0.005, vjust=0, size=2.0, angle=90, color="darkred") +
        geom_point(alpha=.5,show.legend = TRUE) +
        geom_vline(data=annotations, aes(xintercept=sp), 
            color="darkred", linetype="dashed", alpha=.5,linewidth=.5) +
        scale_y_continuous(
            expand=c(0,0),
            limits=c(0,1.02*max(pos_df$logp))
        ) +
        scale_x_continuous(expand=c(0,0),limits=c(xlims[1]-10, xlims[2]+10),,breaks=seq(from=0,to=xbreaks[length(xbreaks)], by =.5*1e6),labels=seq(from=0,to=xbreaks[length(xbreaks)], by =.5*1e6)/1e3)+
        scale_color_brewer(palette="Dark2",drop = FALSE,breaks=contigs_inord)+
        labs(x="Position (kb)", y="-log p",fill="Chromosome",color="Chromosome")+
        theme_minimal() +
        theme(
            legend.title = element_text(size=8),
            legend.text = element_text(size=6),
            axis.line.x = element_line(color="gray75"),
            axis.text.y = element_text(size = 6),
            axis.text.x = element_text(size=6),
            axis.title.x = element_text(size=8),
            axis.title.y = element_text(size=8)
        )
}

.plot_heatmaps <- function(df)
{
    max_b <- max(as.integer(df$run))
    df <- .prepare_df(df)

    df %>% ggplot(aes(x=run,y=PIP,fill=num_replicates)) + 
        geom_tile()+
        scale_fill_distiller(palette="RdBu",
            labels=scales::percent_format(scale = 1),
            guide = guide_colourbar(title.position = "left")) +
        scale_x_continuous(breaks=seq(from=0, to=max_b, by = (1+max_b)/10),
            labels=paste0(0+seq(from=0, to=max_b, by=(1+max_b)/10)),
            expand = c(0,0)) +
        scale_y_continuous(breaks = seq(.1,.9,by=.2), 
            labels = paste0("\u2265",seq(.1,.9,by=.2)),
            expand = c(0,0))+
        facet_grid(cols=vars(gene),rows=vars(type))+
        theme_bw()+
        theme(
            aspect.ratio=1,
            text=element_text(size=8),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, size=6),
            legend.title = element_text(angle = 90, size = 8, hjust=.5, vjust = 1))
}

plot_heatmaps_ann <- function(df)
{
    df %>% .plot_heatmaps + 
        labs(fill="% Replicates", y="Posterior Inclusion Probability", x="Batch")
}

plot_heatmaps_pos <- function(df)
{
    df %>% .plot_heatmaps + 
        labs(fill="% Replicates", y="Max Posterior Inclusion Probability", x="Batch")
}


.prepare_df <- function(df,levs=seq(.1,.9,length.out=9))
{
    gls <- unique(df$gene)
    gls <- c(gls[!(gls=="Off-Target")],"Off-Target")
    df %>% 
        mutate(gene=factor(gene, levels=gls))%>%
        mutate(type=ifelse(str_detect(replicate,"A"), "Adaptive", "Random")) %>%
        mutate(PIP = lapply(1:n(), function(x) levs))%>%
        unnest_longer(PIP) %>%
        group_by(gene, run, PIP,type) %>%
        summarise(num_replicates = 100*sum(pip >= PIP)/n(),.keep_all=T) %>%
        mutate(run=as.integer(run))
}

plot_power <- function(df,alphas =c(.5,.7,.9)) {
    max_b <- max(as.integer(df$run))
    .prepare_df(df,levs=alphas) %>% 
        filter(gene!="Off-Target") %>%
        filter(PIP%in%alphas)%>%
        mutate(PIP=factor(PIP, levels=alphas, labels=paste0("\u2265 ",alphas)))%>%
        ggplot(aes(x=run,y=num_replicates,color=PIP)) + 
        geom_line()+
        scale_color_brewer(palette="Dark2")+
        scale_x_continuous(breaks=seq(from=0, to=max_b, by = (1+max_b)/10),
            labels=paste0(0+seq(from=0, to=max_b, by=(1+max_b)/10)),
            expand = c(0,0)) +
        scale_y_continuous(breaks = 100*seq(0,1,length.out=11), 
            labels=scales::percent_format(scale = 1),
            expand = c(0,2.5))+
        labs(x="Batch", y="Power", color="PIP")+
        facet_grid(cols=vars(gene),rows=vars(type))+
        theme_bw()+
        theme(
            aspect.ratio=1,
            text=element_text(size=8),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, size=6),
            legend.title = element_text(angle = 0, size = 8, hjust=.5, vjust = 1))
}