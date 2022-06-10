#!/usr/bin/env Rscript
# -*- mode: R -*-

suppressPackageStartupMessages(
    {
        require(tidyverse)
        require(ape)
        require(phyclust)
        require(ggfortify)
        require(ggtree)
    })

args <- commandArgs(trailingOnly = TRUE)

fastani_out=args[1]
codes_txt=args[2]
output=args[3]
title=args[4]

pdf(file=paste(dirname(output), 'Rplots.pdf', sep = '/'))

# COLORS
codes <- read.delim(codes_txt, header = FALSE)
names(codes) <- c('run', 'label', 'id')
colfunc <- colorRampPalette(c("red", "yellow"))
not_P <- sort(levels(factor(codes$label[-grep('^P', codes$label)])))
is_P <- paste0('P',1:10)
col_pallet <- c( rainbow(length(not_P)+ 2)[3:8], colfunc(length(is_P)))
names(col_pallet) <- c(not_P, is_P)


y <- read.delim(fastani_out, header = FALSE, sep = '\t')

names(y) <- c('seq.A', 'seq.B', 'ANI')

y %>% pivot_wider(names_from = seq.B, values_from=ANI ) %>% replace(is.na(.), 0) -> y.dist

y.tree <- nj(as.dist(y.dist[, !names(y.dist) %in% c("seq.A")]))


group.info <- split(y.tree$tip.label, gsub("#.+", "", y.tree$tip.label))
y.tree <- groupOTU(y.tree, group.info)

ggtree(y.tree, branch.length = 'none', layout = 'rectangular') + geom_tippoint(aes(color=group), size=1) + scale_color_manual('Passages' , values=col_pallet) +  ggtitle(title)

ggsave(output, height=60, width=10, limitsize = FALSE)

ggtree(y.tree, branch.length = 'none', layout = 'circular') + geom_tippoint(aes(color=group), size=1) + scale_color_manual('Passages' , values=col_pallet) + ggtitle(title)

ggsave(gsub('rectangular', 'circular', output), height=15, width=15)
