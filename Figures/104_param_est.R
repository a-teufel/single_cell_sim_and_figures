#draws estimates of parameter values from final rounds of abc for the no drug case. Figure 3

rm(list = ls())
library(ggplot2)
library(cowplot)
library(ggridges)
require(reshape2)


data <- NULL
names <- NULL

df <- NULL

#read in data
d_104 = read.table(
  "C:/Users/User/Documents/Polio_data_104/estimation_iterations/Final_ABC_save_parameters_labs_new3.txt",
  header = TRUE
)


lower <- apply(d_104[1:10], 2, min)
upper <- apply(d_104[1:10], 2, max)


D_104 <- as.data.frame(d_104)
D_104$Experiment <- "no drug"

#palettes for drug and no drug cases
darkPalette <-
  c(
    "#999999",
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7"
  )
cbPalette <-
  c(
    "#999999",
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7"
  )

lighten <- function(color, factor = 1.2) {
  col <- col2rgb(color)
  col <- col / factor
  col <- rgb(t(col), maxColorValue = 255)
  col
}


cbPalette <- c(cbPalette[1], lighten(cbPalette[1]))


all_alpha = 0.1
all_size = 1.5
all_adjust = 1

leg_pos <- theme(legend.position = "top")
non_pos <- theme(legend.position = 'none')


total <- rbind(D_104)
p1 <-
  ggplot(total, aes(c.trans., colour = Experiment, fill = Experiment)) +
  geom_density(alpha = all_alpha,
               size = all_size,
               adjust = all_adjust) + scale_fill_manual(values = cbPalette) + scale_colour_manual(values =
                                                                                                    cbPalette) + scale_y_continuous(expand = c(0.001, 0)) + scale_x_continuous(expand = c(0.001, 0),
                                                                                                                                                                               limits = c(lower[1], upper[1])) + xlab(expression('log c'[tran]))
p2 <-
  ggplot(total, aes(c.com., colour = Experiment, fill = Experiment)) + geom_density(alpha =
                                                                                      all_alpha,
                                                                                    size = all_size,
                                                                                    adjust = all_adjust) + scale_fill_manual(values = cbPalette) + scale_colour_manual(values =
                                                                                                                                                                         cbPalette) + scale_y_continuous(expand = c(0.001, 0)) + scale_x_continuous(expand = c(0.001, 0),
                                                                                                                                                                                                                                                    limits = c(lower[2], upper[2])) + xlab(expression('log c'[com]))

#annotate(geom = 'text', label = c(paste("  p =",p[2])), x = -Inf, y = Inf, hjust = -4, vjust = 1,parse=FALSE)
p3 <-
  ggplot(total, aes(c.circ., colour = Experiment, fill = Experiment)) + geom_density(alpha =
                                                                                       all_alpha,
                                                                                     size = all_size,
                                                                                     adjust = all_adjust) + scale_fill_manual(values = cbPalette) + scale_colour_manual(values =
                                                                                                                                                                          cbPalette) + scale_y_continuous(expand = c(0.001, 0)) + scale_x_continuous(expand = c(0.001, 0),
                                                                                                                                                                                                                                                     limits = c(lower[3], upper[3])) + xlab(expression('log c'[circ]))
#annotate(geom = 'text', label = c(paste("  p =",p[3])), x = -Inf, y = Inf, hjust = -4, vjust = 1,parse=FALSE)
p4 <-
  ggplot(total, aes(c.rep.., colour = Experiment, fill = Experiment)) + geom_density(alpha =
                                                                                       all_alpha,
                                                                                     size = all_size,
                                                                                     adjust = all_adjust) + scale_fill_manual(values = cbPalette) + scale_colour_manual(values =
                                                                                                                                                                          cbPalette) + scale_y_continuous(expand = c(0.001, 0)) + scale_x_continuous(expand = c(0.001, 0),
                                                                                                                                                                                                                                                     limits = c(lower[4], upper[4])) + xlab(expression('log c'[rep ~ pos]))
#annotate(geom = 'text', label = c(paste("  p =",p[4])), x = -Inf, y = Inf, hjust = -4, vjust = 1,parse=FALSE)
p5 <-
  ggplot(total, aes(c.rep...1, colour = Experiment, fill = Experiment)) +
  geom_density(alpha = all_alpha,
               size = all_size,
               adjust = all_adjust) + scale_fill_manual(values = cbPalette) + scale_colour_manual(values =
                                                                                                    cbPalette) + scale_y_continuous(expand = c(0.001, 0)) + scale_x_continuous(expand = c(0.001, 0),
                                                                                                                                                                               limits = c(lower[5], upper[5])) + xlab(expression('log c'[rep ~ neg]))
#annotate(geom = 'text', label = c(paste("  p =",p[5])), x = -Inf, y = Inf, hjust = 0, vjust = 1,parse=FALSE)
p6 <-
  ggplot(total, aes(c.pack., colour = Experiment, fill = Experiment)) + geom_density(alpha =
                                                                                       all_alpha,
                                                                                     size = all_size,
                                                                                     adjust = all_adjust) + scale_fill_manual(values = cbPalette) + scale_colour_manual(values =
                                                                                                                                                                          cbPalette) + scale_y_continuous(expand = c(0.001, 0)) + scale_x_continuous(expand = c(0.001, 0),
                                                                                                                                                                                                                                                     limits = c(lower[6], upper[6])) + xlab(expression('log c'[pack]))
#annotate(geom = 'text', label = c(paste("  p =",p[6])), x = -Inf, y = Inf, hjust = 0, vjust = 1,parse=FALSE)

p7 <-
  ggplot(total, aes(com.max., colour = Experiment, fill = Experiment)) + geom_density(alpha =
                                                                                        all_alpha,
                                                                                      size = all_size,
                                                                                      adjust = all_adjust) + scale_fill_manual(values = cbPalette) + scale_colour_manual(values =
                                                                                                                                                                           cbPalette) + scale_y_continuous(expand = c(0.001, 0)) + scale_x_continuous(expand = c(0.001, 0),
                                                                                                                                                                                                                                                      limits = c(lower[7], upper[7]))  + xlab(expression('log com'[max]))
#annotate(geom = 'text', label = c(paste("  p =",p[7])), x = -Inf, y = Inf, hjust = 0, vjust = 1,parse=FALSE)
p8 <-
  ggplot(total, aes(rep.max., colour = Experiment, fill = Experiment)) + geom_density(alpha =
                                                                                        all_alpha,
                                                                                      size = all_size,
                                                                                      adjust = all_adjust) + scale_fill_manual(values = cbPalette) + scale_colour_manual(values =
                                                                                                                                                                           cbPalette) + scale_y_continuous(expand = c(0.001, 0)) + scale_x_continuous(expand = c(0.001, 0),
                                                                                                                                                                                                                                                      limits = c(lower[8], upper[8]))  + xlab(expression('log rep'[max]))
#annotate(geom = 'text', label = c(paste("  p =",p[8])), x = -Inf, y = Inf, hjust = 0, vjust = 1,parse=FALSE)

p9 <-
  ggplot(total, aes(c.3A., colour = Experiment, fill = Experiment)) + geom_density(alpha =
                                                                                     all_alpha,
                                                                                   size = all_size,
                                                                                   adjust = all_adjust) + scale_fill_manual(values = cbPalette) + scale_colour_manual(values =
                                                                                                                                                                        cbPalette) + scale_y_continuous(expand = c(0.001, 0)) + scale_x_continuous(expand = c(0.001, 0),
                                                                                                                                                                                                                                                   limits = c(lower[9], upper[9])) + xlab(expression('log c'["3A"]))
#annotate(geom = 'text', label = c(paste("  p =",p[9])), x = -Inf, y = Inf, hjust = 0, vjust = 1,parse=FALSE)
p10 <-
  ggplot(total, aes(c.stay., colour = Experiment, fill = Experiment)) + geom_density(alpha =
                                                                                       all_alpha,
                                                                                     size = all_size,
                                                                                     adjust = all_adjust) + scale_fill_manual(values = cbPalette) + scale_colour_manual(values =
                                                                                                                                                                          cbPalette) + scale_y_continuous(expand = c(0.001, 0)) + scale_x_continuous(expand = c(0.001, 0),
                                                                                                                                                                                                                                                     limits = c(lower[10], upper[10]))  + xlab(expression('c'[stay]))
#annotate(geom = 'text', label = c(paste("  p =",p[10])), x = -Inf, y = Inf, hjust = -4.5, vjust = 1,parse=FALSE)
p1 <- p1 + coord_cartesian(clip = "off") + non_pos
p2 <- p2 + coord_cartesian(clip = "off") + non_pos
p3 <- p3 + coord_cartesian(clip = "off") + non_pos
p4 <- p4 + coord_cartesian(clip = "off") + non_pos


p5 <- p5 + coord_cartesian(clip = "off") + non_pos
p6 <- p6 + coord_cartesian(clip = "off") + non_pos
p7 <- p7 + coord_cartesian(clip = "off") + non_pos
p8 <- p8 + coord_cartesian(clip = "off") + non_pos
p9 <- p9 + coord_cartesian(clip = "off") + non_pos
p10 <- p10 + coord_cartesian(clip = "off") + non_pos



pp <- plot_grid(
  p1,
  p2,
  p3,
  p4,
  p5,
  p6,
  p7,
  p8,
  p9,
  p10,
  labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"),
  align = 'vh',
  hjust = -1,
  nrow = 5
)


ggplot2::ggsave(
  "C:/Users/User/Documents/Polio_data_104/just104_range.pdf",
  plot = pp,
  width = 12,
  height = 12,
  units = "in"
)