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

# colfunc <- colorRampPalette(c("red", "yellow"))
# phage.colors=c(colfunc(10), rainbow(8)[3:7])

y <- read.delim(input)
y %>% mutate(path.a=path.a, path.b=path.b, jaccard=jaccard, path.a.length=NULL, path.b.length=NULL, intersection=NULL, euclidean=NULL) %>% pivot_wider(names_from=path.b, values_from=jaccard) %>% replace(is.na(.), 0) -> y.dist
y.tree <- nj(as.dist(y.dist[, !names(y.dist) %in% c("path.a")]))

p <- ggtree(y.tree, layout = 'rectangular', branch.length = 'none') %<+% data_frame(
  labels = y.tree$tip.label,
  passages = str_split_fixed(y.tree$tip.label, "#", 2)[,1]
) + aes(color=passages)  + geom_tree()

# Eriks code
# ggtree(y.tree) %<+% data.frame(
#   node=1:nrow(y.tree$edge),
#   group.name=factor()
# ) + aes(color=group.name)
#   + geom_tree()
#   + scale_color_manual("passage",values=c(phage.colors, "black"))

# p <- ggtree(y.tree) %<+% data.frame(
#   node=1:nrow(y.tree$edge),
#   group.name=factor(
#     c(as.character(y$path.a),rep("internal",nrow(y.tree$edge) - nrow(y))),
#     levels=c(levels(y$path.a),"internal") 
#   )
# )
#   + aes(color=group.name)
#   + geom_tree()
#   + scale_color_manual("passage",values=c(phage.colors, "black"))

# ggsave(paste(output, "ggtree.passage.pdf", sep="."), height=40, width=9)