
rm(list=ls())
library(ggplot2)
library(cowplot)
library(ggridges)
require(reshape2)
library(colorspace)
library(dplyr)
library(ggbeeswarm)
library(mclust)
library(DataExplorer)
library(data.table)

#draws PCA Figure 5 and S5

d_104 = read.table("data/104/Final_ABC_save_parameters.txt",header=TRUE)
#nice labels
n<-c(expression('log c'[tran]),expression('log c'[com]),expression('log c'[circ]),expression('log c'[rep~pos]),
     expression('log c'[rep~neg]),expression('log c'[pack]),expression('log com'[max]),expression('log rep'[max]),expression('log c'["3A"]),expression('c'[stay]))


d_104<-d_104[1:10]
colnames(d_104)<-n


normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}


#change viz enviroment a little
theme_dviz_vgrid <- function(font_size = 14, font_family = dviz_font_family, line_size = .5,
                             rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14,
                             colour = "grey90") {
  half_line <- font_size / 2

  cowplot::theme_minimal_vgrid(font_size = font_size, font_family = font_family, line_size = line_size,
                               rel_small = rel_small, rel_tiny = rel_tiny, rel_large = rel_large,
                               colour = colour)  %+replace%
    theme(
      plot.margin = margin(half_line/2, 1.5, half_line/2, 1.5),
      complete = TRUE
    )
}


#rehash of what data explore does, just draws it nicer
my_plot_prcomp <- function(data, variance_cap = 0.8, maxcat = 50L, prcomp_args = list("scale." = TRUE), title = NULL, nrow = 3L, ncol = 3L, parallel = FALSE) {
  ## Declare variable first to pass R CMD check
  pc <- pct <- cum_pct <- Feature <- variable <- value <- NULL
  ## Check if input is data.table
  data <- data.table(data)
  ## Dummify data
  dt <- suppressWarnings(split_columns(dummify(data, maxcat = maxcat))$continuous)
  prcomp_args_list <- list("x" = dt, "retx" = FALSE)
  ## Analyze principle components
  pca <- tryCatch(
    do.call("prcomp", c(prcomp_args_list, prcomp_args)),
    error = function(e) {
      message(e$message)
      if (grepl("missing", e$message)) stop("\nConsider passing `na.omit(data)` as input.")
    }
  )

  ## Calcualte principle components standard deviation
  var_exp <- pca$sdev ^ 2
  pc_var <- data.table(
    "pc" = paste0("PC ",seq_along(pca$sdev)),
    "var" = var_exp,
    "pct" = var_exp / sum(var_exp),
    "cum_pct" = cumsum(var_exp) / sum(var_exp)
  )
  min_cum_pct <- min(pc_var$cum_pct)
  pc_var2 <- pc_var[cum_pct <= max(variance_cap, min_cum_pct)]
  ## Create explained variance plot
  varexp_plot <- ggplot(pc_var2, aes(x = reorder(pc, pct), y = pct)) +
    geom_col(fill = "#999999", alpha = 0.9) +
    geom_text(aes(label = paste("cumulative: ", round(cum_pct*100,2),"%",sep="")), color = "white", hjust = 1.1) +
    scale_y_continuous(limits = c(0, .22),
                       expand = c(0, 0),
                       breaks = c(0, 0.05, 0.10, 0.15,0.20),
                       labels = c("0%","5%", "10%", "15%","20%"),
                       name = "variance explained") +
    scale_x_discrete(expand = c(0, 0.5), name="principle components") +
    coord_flip(clip = "off") +
    theme(
      #axis.ticks.length = grid::unit(0, "pt"),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.x = element_blank(),
      panel.grid.major.x = element_line(colour="gray"),
      axis.ticks.x=element_line(colour="gray")
    )

  return(varexp_plot)

}

my_plot_prcomp_features <- function(data, variance_cap = 0.8, maxcat = 50L, prcomp_args = list("scale." = TRUE), title = NULL, nrow = 3L, ncol = 3L, parallel = FALSE,PC=1) {
  ## Declare variable first to pass R CMD check
  pc <- pct <- cum_pct <- Feature <- variable <- value <- NULL
  ## Check if input is data.table
  data <- data.table(data)
  ## Dummify data
  dt <- suppressWarnings(split_columns(dummify(data, maxcat = maxcat))$continuous)
  prcomp_args_list <- list("x" = dt, "retx" = FALSE)
  ## Analyze principle components
  pca <- tryCatch(
    do.call("prcomp", c(prcomp_args_list, prcomp_args)),
    error = function(e) {
      message(e$message)
      if (grepl("missing", e$message)) stop("\nConsider passing `na.omit(data)` as input.")
    }
  )

  ## Calcualte principle components standard deviation
  var_exp <- pca$sdev ^ 2
  pc_var <- data.table(
    "pc" = paste0("PC ",seq_along(pca$sdev)),
    "var" = var_exp,
    "pct" = var_exp / sum(var_exp),
    "cum_pct" = cumsum(var_exp) / sum(var_exp)
  )
  min_cum_pct <- min(pc_var$cum_pct)
  pc_var2 <- pc_var[cum_pct <= max(variance_cap, min_cum_pct)]
  rotation_dt <- data.table(
    "Feature" = rownames(pca$rotation),
    data.table(pca$rotation)[, seq.int(nrow(pc_var2)), with = FALSE]
  )
  melt_rotation_dt <- melt.data.table(rotation_dt, id.vars = "Feature")
  feature_names <- rotation_dt[["Feature"]]
  print(melt_rotation_dt)
  ## Calculate number of pages
  ## Create list of ggplot objects
  n<-c(expression('log c'[tran]),expression('log c'[com]),expression('log c'[circ]),expression('log c'[rep~pos]),
       expression('log c'[rep~neg]),expression('log c'[pack]),expression('log com'[max]),expression('log rep'[max]),expression('log c'["3A"]),expression('c'[stay]))
  df<-melt_rotation_dt[variable==paste("PC",PC,sep="")]
  print(df)
  d<-as.data.frame(df)
  print(d)
  d$lab<-as.list(n)
  print(d)

  d<- d[order(d$value),]
  d$Feature <- factor(d$Feature, levels = d$Feature)
  d$lab <- factor(d$lab, levels = d$lab)

   p1<-ggplot(d, aes(x = Feature, y = value))+
    geom_col(fill = "#999999", alpha = 0.9) +
        coord_flip() +
        scale_y_continuous(name = "relative importance",limits =c(-.8,.8), labels=function(l) parse(text=l))+
        scale_x_discrete(name = "feature", labels=parse(text=levels(d$lab)))+
        ggtitle(paste("PC ",PC,sep=""))+
  theme(
    #axis.ticks.length = grid::unit(0, "pt"),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_blank(),
    panel.grid.major.x = element_line(colour="gray"),
    axis.ticks.x=element_line(colour="gray")
  )





  return(p1)
}

norm<-d_104

#normal pca
p<-my_plot_prcomp(norm)

#features
p1<-my_plot_prcomp_features(norm,PC=1)+theme(axis.title.x = element_blank())
p2<-my_plot_prcomp_features(norm,PC=2)+theme(axis.title.y = element_blank())
p3<-my_plot_prcomp_features(norm,PC=3)+theme(axis.title.y = element_blank(),axis.title.x = element_blank())
p4<-my_plot_prcomp_features(norm,PC=4)+theme(axis.title.x = element_blank())
p5<-my_plot_prcomp_features(norm,PC=5)+theme(axis.title.y = element_blank())
p6<-my_plot_prcomp_features(norm,PC=6)+theme(axis.title.y = element_blank(),axis.title.x = element_blank())


pp<- plot_grid( p1,
                p2,
                p3,
                p4,
                p5,
                p6,
                align = 'vh',
                labels = c("B","C","D","E","F","G"),
                hjust = -1,
                nrow = 2)
pp

ppp<-plot_grid(p,labels=c("A"),align='vh',nrow=1,hjust=-1)

pppp<-plot_grid(ppp,pp,nrow=2,align='vh',hjust=-1,rel_heights = c(.5,1))
pppp

ggplot2::ggsave("pca.pdf",plot=pppp,width = 12, height = 12,units = "in")

