library(ggplot2)
library(gridExtra)
library(parallel)

plot_enrichment <- function(scores, gene_set, power = 1, name="GSEA-like Enrichment", plot=T, nperm=10000, nproc=1) {
  scores <- scores[order(-as.numeric(scores))]
  selected <- which(names(scores) %in% gene_set)
  
  m <- length(selected)
  n <- length(scores)
  
  estimate_score <- function(selected, power, only_score=F) {
    scores_adjusted <- abs(scores[selected])^power
    nr <- sum(scores_adjusted)
    
    score_cum_sum <- cumsum(scores_adjusted) / nr
    
    tops <- score_cum_sum - (selected - seq_along(selected)) / (n - m)
    bottoms <- tops - scores_adjusted / nr
    
    if (only_score) {
      if (max(tops) > -min(bottoms))
        return(max(tops))
      else
        return(min(bottoms))
    } else {
      return(list(tops=tops, bottoms=bottoms))
    }
  }
  
  
  values <- mclapply(1:nperm, function(k) {
    rand_sample <- sort(sample(seq_along(scores), size=length(selected)))
    estimate_score(rand_sample, power, only_score=T)
  }, mc.cores=nproc)
  values = unlist(values)
  
  
  res <- estimate_score(selected, power)
  bottoms = res$bottoms
  tops = res$tops
  
  xs <- as.vector(rbind(selected - 1, selected))
  ys <- as.vector(rbind(bottoms, tops))
  to_plot <- data.frame(x=c(0, xs, n + 1), y=c(0, ys, 0))
  
  diff <- (max(tops) - min(bottoms)) / 16
  bmin <- min(min(bottoms) - 2*diff, min(values))
  tmax <- max(max(tops), max(values))
  
  
  # bar code at the bottom
  to_plot_lines <- data.frame(
    x=selected,
    y=rep(min(bottoms)-2 * diff, length(selected)),
    xend=selected,
    yend=rep(min(bottoms)-diff, length(selected))
  )
  
  

  yl <- c(bmin, tmax)
  
  g <- ggplot(to_plot, aes(x=x, y=y)) +
    geom_point(color="green", size=0.1) +
    geom_hline(yintercept=max(tops), colour="red", linetype="dashed") + 
    geom_hline(yintercept=min(bottoms), colour="red", linetype="dashed") +
    geom_hline(yintercept=0, colour="black") +
    geom_line(color="green") + theme_bw() +
    geom_segment(data=to_plot_lines, mapping=aes(x=x, y=y, xend=xend, yend=yend), size=0.1) +
    theme(axis.text.x=element_blank(), axis.title.x=element_blank()) + labs(title=name, y="enrichment score") +
    scale_y_continuous(expand=c(0, 0), limits=yl)

  distribution <- ggplot(mapping=aes(values)) + geom_histogram(binwidth=0.01) + scale_x_continuous(expand=c(0, 0), limits=yl) + coord_flip() +
    theme_bw() +
    theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank()) + 
      labs(title="Distribution")
  
  
  grid.arrange(g, distribution, ncol=2, nrow=1, widths=c(5, 3), heights=c(1))
  
}
