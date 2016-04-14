# Enrichment plot
Plots GSEA-like enrichment plot using ggplot2 and gridExtra

Usage:
```R
head(scores) # named vector
# 0610007P14Rik 0610009B22Rik 0610010B08Rik 0610010K14Rik 0610012G03Rik 0610037L13Rik 
#   -0.09936889   -0.09320537   -0.18495076    0.01866837   -0.37919820   -0.03912840 
head(gene_set) # character vector
# "Strn4"   "Atp13a1" "Tom1"    "Nmt2"    "Nmt1"    "Pkd2" 
png("example.png", width=8, height=5, units="in", res=300)
plot_enrichment(scores, gene_set, name="Example enrichment")
dev.off()
```

![alt text](https://github.com/konsolerr/utils/blob/master/enrichment_plot/example.png "Example of enrichment plot")
