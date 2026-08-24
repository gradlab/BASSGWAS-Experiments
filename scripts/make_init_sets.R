library(tidyverse)
args <- commandArgs(trailingOnly=TRUE) 

ifile <- args[1]
ninit <- as.integer(args[2])
nset <- as.integer(args[3])
seed <- as.integer(args[4])
odir <- args[5]

set.seed(seed)
dat <- read_csv(ifile)

foo <- function(x)
{
    idx1 <- which(dat[[2]])
    idx0 <- which(!dat[[2]])
    s <- dat[c(sample(idx1,1), sample(idx0,1)),][[1]]
    iset <- sample(dat[!(dat[[1]] %in% s),][[1]], ninit-2)
    c(s,iset)
}

sample_sets <- lapply(1:nset, foo)

for (i in seq_along(sample_sets))
{
    s <- sample_sets[[i]]
    d <- dat[dat[[1]] %in% s,]
    dir.create(paste0(odir, "/R",i))
    dir.create(paste0(odir, "/R",i,"/RUN_0"))
    dir.create(paste0(odir, "/A",i))
    dir.create(paste0(odir, "/A",i,"/RUN_0"))
    write_csv(d, paste0(odir, "/A",i,"/RUN_0/data_bin",".csv"))
    write_csv(d, paste0(odir, "/R",i,"/RUN_0/data_bin",".csv"))
}
warnings()
