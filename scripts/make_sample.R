library(tidyverse)

data_set_sz <- 800
cutoffs <- list(azi = 2, cip = 1, tet = 8)
freq <- list(azi = .1, cip = .1, tet10 = .1, tet2.5 = 0.025)
counts <- lapply(freq, function(x) x*data_set_sz)

meta <- read_delim("data/meta.tsv") %>%
    filter(reference == "Reimche2023") %>% 
    select(wgs_id, azithromycin, ciprofloxacin, tetracycline) %>%
    rename(azi_mic=azithromycin, cip_mic=ciprofloxacin, tet_mic=tetracycline) %>%
    mutate(azi_mic=as.numeric(azi_mic),
        cip_mic=as.numeric(cip_mic),
        tet_mic=as.numeric(tet_mic)) %>%
    mutate(azi = azi_mic >= cutoffs["azi"],
        cip = cip_mic >= cutoffs["cip"],
        tet = tet_mic >= cutoffs["tet"])

downsample <- function(x,y,a,b)
{
    pl <- sample(which(y), a, replace=F)
    mi <- sample(which(!y), b-a, replace=F)
    idx <- sort(c(pl,mi))
    return(x[idx])
}

data_azi <- downsample(meta$wgs_id, meta$azi, counts[["azi"]], data_set_sz)
data_cip <- downsample(meta$wgs_id, meta$cip, counts[["cip"]], data_set_sz)
data_tet10 <- downsample(meta$wgs_id, meta$tet, counts[["tet10"]], data_set_sz)
data_tet2.5 <- downsample(meta$wgs_id, meta$tet, counts[["tet2.5"]], data_set_sz)

write_delim(meta, "out/meta.tsv")
writeLines(data_azi, "out/azi/isolates.txt")
write_csv(meta %>% filter(wgs_id %in% data_azi) %>%select(wgs_id, azi), "out/azi/data_bin.csv")
writeLines(data_cip, "out/cip/isolates.txt")
write_csv(meta %>% filter(wgs_id %in% data_cip) %>%select(wgs_id, cip), "out/cip/data_bin.csv")

writeLines(data_tet10, "out/tet10/isolates.txt")
write_csv(meta %>% filter(wgs_id %in% data_tet10) %>%select(wgs_id, tet), "out/tet10/data_bin.csv")
writeLines(data_tet2.5, "out/tet10/isolates.txt")
write_csv(meta %>% filter(wgs_id %in% data_tet2.5) %>%select(wgs_id, tet), "out/tet2.5/data_bin.csv")

