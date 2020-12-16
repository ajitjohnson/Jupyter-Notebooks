# Metagene2 plot

library(metagene2)

# load BAM files
setwd("D:/Metagene plots")

# import bam files
bam_files <- file.path(getwd(), list.files(pattern = ".bam$", recursive = TRUE))
names(bam_files) = c("H3K27ac", "IKZF1", "ZFP91", "ctrl")

# import regions
bed_files <- file.path(getwd(), list.files(pattern = ".bed$", recursive = TRUE))
#bed_files <- bed_files[1]

mg <- metagene2$new(regions = bed_files, 
                    bam_files = bam_files,
                    assay='chipseq',
                    paired_end=T,
                    bin_count=200)

mg$produce_metagene(title = "MetaGene Plot")
mg$plot()


