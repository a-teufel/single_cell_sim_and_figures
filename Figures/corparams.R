#draws correlations between parameters for a given data set (Fig 5,S3, S6, & S9)

rm(list = ls())
library(ggplot2)
library(cowplot)
library(ggridges)
require(reshape2)
library(colorspace)
library(dplyr)
hcl_palettes(plot = TRUE)

df <- NULL

d_104 = read.table(
  "C:/Users/User/Documents/Polio_data_104/estimation_iterations/104_Up2/Final_ABC_save_parameters_labs_new3.txt",
  header = TRUE
)


d_104 <- cbind(d_104)

#pretty labels
n <-
  c(
    expression('log c'[tran]),
    expression('log c'[com]),
    expression('log c'[circ]),
    expression('log c'[rep ~ pos]),
    expression('log c'[rep ~ neg]),
    expression('log c'[pack]),
    expression('log com'[max]),
    expression('log rep'[max]),
    expression('log c'["3A"]),
    expression('c'[stay])
  )

d_104 <- d_104[1:10]
colnames(d_104) <- n

#correlations
cm <- cor(d_104)

#cluster
df_wide <- as.data.frame(cm)
df_long <- stack(df_wide)
names(df_long) <- c("cor", "var1")
df_long <-
  cbind(df_long, var2 = rep(rownames(cm), length(rownames(cm))))
clust <- hclust(1 - as.dist(cm), method = "ward.D2")

levels <- clust$labels[clust$order]
df_long$var1 <- factor(df_long$var1, levels = levels)
df_long$var2 <- factor(df_long$var2, levels = levels)
filter(df_long, as.integer(df_long$var1) < as.integer(df_long$var2))
df_long <-
  subset(df_long, as.integer(df_long$var1) < as.integer(df_long$var2))

#plot
p <-
  ggplot(df_long, aes(
    df_long$var1,
    df_long$var2,
    fill = cor,
    size = abs(cor)
  )) +
  geom_point(shape = 21,
             stroke = 1,
             color = "black") +
  scale_x_discrete(
    position = "top",
    name = NULL,
    expand = c(0, 0.6),
    labels = function(l)
      parse(text = l)
  ) +
  scale_y_discrete(
    position = "left",
    name = NULL,
    expand = c(0, 0.6),
    labels = function(l)
      parse(text = l)
  ) +
  scale_size_area(max_size = 25,
                  limits = c(0, 0.6),
                  guide = "none") +
  scale_fill_gradient2(
    low = "#D53E4F",
    mid = "#FEE08B",
    high = "#3288BD",
    name = "correlation",
    limits = c(-.55, .55),
    breaks = c(-.5, 0, .5),
    guide = guide_colorbar(
      direction = "vertical",
      label.position = "right",
      title.position = "left",
      label.theme = element_text(angle = 90, hjust =
                                   .5),
      title.theme = element_text(angle = 90),
      barheight = grid::unit(140, "pt"),
      barwidth  = grid::unit(17.5, "pt"),
      ticks.linewidth = 1
    )
  ) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.ticks.length = grid::unit(3, "pt"),
    legend.position = c(.97, .0),
    legend.justification = c(1, 0),
    legend.title.align = 0.5,
    axis.text.x.top = element_text(
      angle = 90,
      vjust = .5,
      hjust = 0
    ),
    axis.text.y = element_text(vjust = .5)
  )

ggplot2::ggsave(
  "C:/Users/User/Documents/Polio_data_104/cor104_labs_small2.pdf",
  plot = p,
  width = 7,
  height = 7,
  units = "in"
)
