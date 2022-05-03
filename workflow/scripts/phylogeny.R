#!/usr/bin/env Rscript
# -*- mode: R -*-

suppressPackageStartupMessages(
    {
        ## require(tidyverse)
        require(ape)
        require(phyclust)
        ## require(ggfortify)
        require(ggtree)
    })

args <- commandArgs(trailingOnly = TRUE)
input=args[1]
output=args[2]


y <- read.delim(input)

# correct factor order
y$group.name <- factor(as.character(y$group.name), levels=c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10", "S2-55s")

colfunc <- colorRampPalette(c("red", "yellow"))
phage.colors = c(colfunc(10), rainbow(8)[3:7])
phage.colors[11] <- "#B6FF00"


ggplot(y, aes(x=path.length, color=group.name)) + geom_density() + scale_color_manual("passage",values=c(phage.colors, "black"))
ggsave(paste(output, "path.length.density.pdf", sep="."), height=6, width=10)
ggplot(y, aes(x=path.length, fill=group.name)) + geom_histogram(binwidth=50) + scale_fill_manual("passage",values=c(phage.colors, "black"))
ggsave(paste(output, "path.length.hist.pdf", sep="."), height=6, width=10)
