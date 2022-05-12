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

# args <- commandArgs(trailingOnly = TRUE)
# input=args[1]
# output=args[2]

# input <- '/home/hugo/projects/phage-evo-paper/results/single/pggb/distance_matrix_removed_ecoli_removed_phages.tsv'
input <- '/home/hugo/projects/phage-evo-paper/phages.random.500.each.passage.out.gz'
output <- 'plots'

colfunc <- colorRampPalette(c("red", "yellow"))
# phage.groups <- c('Ec_DE3', 'Ec_LE', 'S1-55L', 'S2-36s', 'S2-55s', paste0('P',1:10))
phage.groups <- c('S1-55L', 'S2-36s', 'S2-55s', paste0('P',1:10))
phage.colors <- c(rainbow(8)[5:8], colfunc(10))
names(phage.colors) <- phage.groups

y <- read.delim(input, header = FALSE, sep = ' ')

names(y) <- c('seq.A', 'seq.B', 'ANI')

y %>% pivot_wider(names_from = seq.B, values_from=ANI ) %>% replace(is.na(.), 0) -> y.dist
y.tree <- nj(as.dist(y.dist[, !names(y.dist) %in% c("seq.A")]))

# write.tree(y.tree, 'newick.dnd')
# y.tree <- read.tree('/home/hugo/projects/phage-evo-paper/newick.dnd')

group.info <- split(y.tree$tip.label, gsub("#.+", "", y.tree$tip.label))
y.tree <- groupOTU(y.tree, group.info)
p <- ggtree(y.tree, branch.length = 'none', layout = 'daylight') + geom_tippoint(aes(color=group), size=1) + scale_color_manual('Passages' , values=phage.colors)
# p <- ggtree(y.tree, branch.length = 'none') + geom_tippoint(aes(color=group), size=0.5) + scale_color_manual('Passages' , values=phage.colors)

ggsave(paste(output, "ggtree.passage.pdf", sep="."), height=40, width=9)
