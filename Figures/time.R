#draws timing of events Figure 4

rm(list = ls())
library(ggplot2)
library(cowplot)
library(ggridges)
require(reshape2)

#read in data
filenames <-
  list.files(path = "C:/Users/User/Documents/Polio_data_104/events",
             pattern = "events*",
             full.names = TRUE)

RNA_P <- NULL
RNA_N <- NULL
GFP <- NULL

for (i in filenames) {
  name <- paste(i)
  temp <- read.csv(name, sep = ",")
  #all times need to be divided by 2, to get to half hour mark
  #get the min value, because we want to know when the event first happened
  RNA_P <- c(RNA_P, min(which(temp$RNA_P != 0)) / 2)
  RNA_N <- c(RNA_N, min(which(temp$RNA_N != 0)) / 2)
  GFP <- c(GFP, min(which(temp$GFP != 0)) / 2)
}

RNA_P <- RNA_P[which(is.finite(RNA_P))]
RNA_N <- RNA_N[which(is.finite(RNA_N))]
protein <- GFP[which(is.finite(GFP))]

#everything gray because i just did this for the no drug case
cbPalette <- c("#999999", "#999999", "#999999", "#999999")

all_alpha = 0.2
all_size = 1

RNA_P <- na.omit(as.data.frame(RNA_P))
g2 <-
  ggplot(RNA_P, aes(RNA_P, colour = cbPalette[2], fill = cbPalette[2])) +
  geom_density(alpha = all_alpha, size = all_size) + scale_fill_manual(values =
                                                                         cbPalette[2]) +
  scale_colour_manual(values = cbPalette[2]) + scale_y_continuous(expand = c(0.001, 0)) +
  scale_x_continuous(expand = c(0.001, 0), limits = c(0, 10)) + xlab("hours until first positive sense RNA observed") +
  geom_vline(size = .5,
             xintercept = c(median(RNA_P$RNA_P)),
             linetype = "dashed") + theme(legend.position = "none")

RNA_N <- na.omit(as.data.frame(RNA_N))
g3 <-
  ggplot(RNA_N, aes(RNA_N, colour = cbPalette[3], fill = cbPalette[3])) +
  geom_density(alpha = all_alpha, size = all_size) + scale_fill_manual(values =
                                                                         cbPalette[3]) +
  scale_colour_manual(values = cbPalette[3]) + scale_y_continuous(expand = c(0.001, 0)) +
  scale_x_continuous(expand = c(0.001, 0), limits = c(0, 10)) + xlab("hours until first negative sense RNA observed") +
  geom_vline(size = .5,
             xintercept = c(median(RNA_N$RNA_N)),
             linetype = "dashed") + theme(legend.position = "none")


protein <- na.omit(as.data.frame(protein))
g4 <-
  ggplot(protein, aes(protein, colour = cbPalette[4], fill = cbPalette[4])) +
  geom_density(alpha = all_alpha,
               size = all_size,
               adjust = 2) + scale_fill_manual(values = cbPalette[4]) +
  scale_colour_manual(values = cbPalette[4]) + scale_y_continuous(expand = c(0.001, 0)) +
  scale_x_continuous(expand = c(0.001, 0), limits = c(0, 10)) + xlab("hours until first protein observed") +
  geom_vline(size = .5,
             xintercept = c(median(protein$protein)),
             linetype = "dashed") + theme(legend.position = "none")

pp <- plot_grid(
  g4,
  g2,
  g3,
  align = 'vh',
  hjust = -1,
  nrow = 3,
  labels = "AUTO"
)
pp <- pp + theme(plot.margin = unit(c(.2, .2, .2, .2), "cm"))

ggplot2::ggsave(
  "timing2.pdf",
  plot = pp,
  width = 12,
  height = 12,
  units = "in"
)