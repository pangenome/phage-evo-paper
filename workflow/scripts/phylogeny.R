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


input <- '/home/hugo/projects/phage-evo-paper/results/single/pggb/distance_matrix_removed_ecoli_removed_phages.tsv'
output <- 'plots'

colfunc <- colorRampPalette(c("red", "yellow"))
# phage.colors=c(colfunc(10), rainbow(8)[3:7])

passages = paste0("P",1:10)
phage.colors = structure(colfunc(length(phage_and_passages)), names = passages)
print(phage.colors)


x <- read.delim(input)
x %>% mutate(path.a=path.a, path.b=path.b, jaccard=jaccard, path.a.length=NULL, path.b.length=NULL, intersection=NULL, euclidean=NULL) %>% pivot_wider(names_from=path.b, values_from=jaccard) %>% replace(is.na(.), 0) -> x.dist
x.tree <- nj(as.dist(x.dist[, !names(x.dist) %in% c("path.a")]))
p <- ggtree(x.tree, options(ignore.negative.edge=TRUE),  branch.length='none') # + xlim(0, 0.025 ) + geom_treescale()

teste <- data.frame(x.tree["tip.label"])
teste$prefix <- str_split_fixed(teste$tip.label, "#", 2)[,1]
teste$color <- vapply(teste$prefix, (function(x) phage.colors[x]), FUN.VALUE = 'character', )
teste[["color"]][is.na(teste[["color"]])] <- 'black'
