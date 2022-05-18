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

# Read user input table
args <- commandArgs(trailingOnly = TRUE)
input=args[1]
output=args[2]

# Get an color function object that aceppts an interger 
colfunc <- colorRampPalette(c("red", "yellow"))

# Generate labels for the colors 
phage.groups <- c('E_coli_LE', 'E_coli_bl21_DE3_polished', 'S1-55L', 'S2-36s', 'S2-55s', paste0('P',1:10))

# Get a list of colors 
phage.colors <- c(rainbow(8)[3:7], colfunc(10))

# name the colors with the labels 
names(phage.colors) <- phage.groups

# Read input label 
y <- read.delim(input)

y %>% mutate(path.a=path.a, path.b=path.b, jaccard=jaccard, path.a.length=NULL, path.b.length=NULL, intersection=NULL, euclidean=NULL) %>% pivot_wider(names_from=path.b, values_from=jaccard) %>% replace(is.na(.), 0) -> y.dist
y.tree <- nj(as.dist(y.dist[, !names(y.dist) %in% c("path.a")]))

group.info <- split(y.tree$tip.label, gsub("#.+", "", y.tree$tip.label))
y.tree <- groupOTU(y.tree, group.info)

ggtree(y.tree, branch.length = 'none') + geom_tippoint(aes(color=group), size=1) +
  scale_color_manual('Passages' , values=phage.colors)

ggsave(output, height=40, width=9)

ggtree(y.tree, branch.length = 'none', layout = 'daylight') + geom_tippoint(aes(color=group), size=1) +
  scale_color_manual('Passages' , values=phage.colors)

ggsave(gsub('rectangular', 'daylight', output), height=10, width=20)
