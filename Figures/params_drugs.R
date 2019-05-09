#draws comparison of no drug estimated parameters to each of the drug cases (Fig. S1,S2,S4,s5,s7,s8) to get the variance parameters just change the labels in the ggplot graph part

rm(list = ls())
library(ggplot2)
library(cowplot)
library(ggridges)
require(reshape2)


data <- NULL
names <- NULL

df <- NULL


d_104 = read.table(
  "C:/Users/User/Documents/Polio_data_104/estimation_iterations/104_Up2/Final_ABC_save_parameters_labs_new3.txt",
  header = TRUE
)
d_105 = read.table(
  "C:/Users/User/Documents/Polio_data_104/estimation_iterations/105_Up2/Final_ABC_save_parameters_labs_new3.txt",
  header = TRUE
)
d_106 = read.table(
  "C:/Users/User/Documents/Polio_data_104/estimation_iterations/106_Up2/Final_ABC_save_parameters_labs_new2.txt",
  header = TRUE
)
d_107 = read.table(
  "C:/Users/User/Documents/Polio_data_104/estimation_iterations/107_Up2/Final_ABC_save_parameters_labs_new3.txt",
  header = TRUE
)


D_104 <- as.data.frame(d_104)
D_104$Experiment <- "no drug"

D_105 <- as.data.frame(d_105)
D_105$Experiment <- "rupintrivir"

D_106 <- as.data.frame(d_106)
D_106$Experiment <- "2'-C-meA"

D_107 <- as.data.frame(d_107)
D_107$Experiment <- "ganetespib"

#need to change the p-value calulations depending on what data sets are being compared
p_vals <- NULL
for (i in 1:10) {
  test <- ks.test(d_104[, i], d_107[, i])
  p_vals <- c(p_vals, test$p.value)
}

p <- format(p.adjust(p_vals, "bonferroni"), digits = 2)

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

#change color for each data set
cbPalette <- c(cbPalette[4], cbPalette[1])
all_alpha = 0.1
all_size = 1.5

leg_pos <- theme(legend.position = "top")
non_pos <- theme(legend.position = "none")


total <- rbind(D_104, D_107)
p1 <-
  ggplot(total, aes(c.trans., colour = Experiment, fill = Experiment)) +
  geom_density(alpha = all_alpha, size = all_size) + scale_fill_manual(values =
                                                                         cbPalette) + scale_colour_manual(values = cbPalette) + scale_y_continuous(expand = c(0.001, 0)) +
  scale_x_continuous(expand = c(0.001, 0)) + xlab(expression('log c'[tran])) +
  annotate(
    geom = 'text',
    label = c(paste("  p =", p[1])),
    x = -Inf,
    y = Inf,
    hjust = 0,
    vjust = 1,
    parse = FALSE
  )

p2 <-
  ggplot(total, aes(c.com., colour = Experiment, fill = Experiment)) + geom_density(alpha =
                                                                                      all_alpha, size = all_size) + scale_fill_manual(values = cbPalette) + scale_colour_manual(values =
                                                                                                                                                                                  cbPalette) + scale_y_continuous(expand = c(0.001, 0)) + scale_x_continuous(expand = c(0.001, 0)) +
  xlab(expression('log c'[com])) +
  annotate(
    geom = 'text',
    label = c(paste("  p =", p[2])),
    x = -Inf,
    y = Inf,
    hjust = -4,
    vjust = 1,
    parse = FALSE
  )
p3 <-
  ggplot(total, aes(c.circ., colour = Experiment, fill = Experiment)) + geom_density(alpha =
                                                                                       all_alpha, size = all_size) + scale_fill_manual(values = cbPalette) + scale_colour_manual(values =
                                                                                                                                                                                   cbPalette) + scale_y_continuous(expand = c(0.001, 0)) + scale_x_continuous(expand = c(0.001, 0)) +
  xlab(expression('log c'[circ])) +
  annotate(
    geom = 'text',
    label = c(paste("  p =", p[3])),
    x = -Inf,
    y = Inf,
    hjust = -4,
    vjust = 1,
    parse = FALSE
  )
p4 <-
  ggplot(total, aes(c.rep.., colour = Experiment, fill = Experiment)) + geom_density(alpha =
                                                                                       all_alpha, size = all_size) + scale_fill_manual(values = cbPalette) + scale_colour_manual(values =
                                                                                                                                                                                   cbPalette) + scale_y_continuous(expand = c(0.001, 0)) + scale_x_continuous(expand = c(0.001, 0)) +
  xlab(expression('log c'[rep ~ pos])) +
  annotate(
    geom = 'text',
    label = c(paste("  p =", p[4])),
    x = -Inf,
    y = Inf,
    hjust = -4,
    vjust = 1,
    parse = FALSE
  )
p5 <-
  ggplot(total, aes(c.rep...1, colour = Experiment, fill = Experiment)) +
  geom_density(alpha = all_alpha, size = all_size) + scale_fill_manual(values =
                                                                         cbPalette) + scale_colour_manual(values = cbPalette) + scale_y_continuous(expand = c(0.001, 0)) +
  scale_x_continuous(expand = c(0.001, 0)) + xlab(expression('log c'[rep ~
                                                                       neg])) +
  annotate(
    geom = 'text',
    label = c(paste("  p =", p[5])),
    x = -Inf,
    y = Inf,
    hjust = 0,
    vjust = 1,
    parse = FALSE
  )
p6 <-
  ggplot(total, aes(c.pack., colour = Experiment, fill = Experiment)) + geom_density(alpha =
                                                                                       all_alpha, size = all_size) + scale_fill_manual(values = cbPalette) + scale_colour_manual(values =
                                                                                                                                                                                   cbPalette) + scale_y_continuous(expand = c(0.001, 0)) + scale_x_continuous(expand = c(0.001, 0)) +
  xlab(expression('log c'[pack])) +
  annotate(
    geom = 'text',
    label = c(paste("  p =", p[6])),
    x = -Inf,
    y = Inf,
    hjust = 0,
    vjust = 1,
    parse = FALSE
  )

p7 <-
  ggplot(total, aes(com.max., colour = Experiment, fill = Experiment)) + geom_density(alpha =
                                                                                        all_alpha, size = all_size) + scale_fill_manual(values = cbPalette) + scale_colour_manual(values =
                                                                                                                                                                                    cbPalette) + scale_y_continuous(expand = c(0.001, 0)) + scale_x_continuous(expand = c(0.001, 0))  +
  xlab(expression('log com'[max])) +
  annotate(
    geom = 'text',
    label = c(paste("  p =", p[7])),
    x = -Inf,
    y = Inf,
    hjust = 0,
    vjust = 1,
    parse = FALSE
  )
p8 <-
  ggplot(total, aes(rep.max., colour = Experiment, fill = Experiment)) + geom_density(alpha =
                                                                                        all_alpha, size = all_size) + scale_fill_manual(values = cbPalette) + scale_colour_manual(values =
                                                                                                                                                                                    cbPalette) + scale_y_continuous(expand = c(0.001, 0)) + scale_x_continuous(expand = c(0.001, 0))  +
  xlab(expression('log rep'[max])) +
  annotate(
    geom = 'text',
    label = c(paste("  p =", p[8])),
    x = -Inf,
    y = Inf,
    hjust = 0,
    vjust = 1,
    parse = FALSE
  )

p9 <-
  ggplot(total, aes(c.3A., colour = Experiment, fill = Experiment)) + geom_density(alpha =
                                                                                     all_alpha, size = all_size) + scale_fill_manual(values = cbPalette) + scale_colour_manual(values =
                                                                                                                                                                                 cbPalette) + scale_y_continuous(expand = c(0.001, 0)) + scale_x_continuous(expand = c(0.001, 0)) +
  xlab(expression('log c'["3A"])) +
  annotate(
    geom = 'text',
    label = c(paste("  p =", p[9])),
    x = -Inf,
    y = Inf,
    hjust = 0,
    vjust = 1,
    parse = FALSE
  )
p10 <-
  ggplot(total, aes(c.stay., colour = Experiment, fill = Experiment)) + geom_density(alpha =
                                                                                       all_alpha, size = all_size) + scale_fill_manual(values = cbPalette) + scale_colour_manual(values =
                                                                                                                                                                                   cbPalette) + scale_y_continuous(expand = c(0.001, 0)) + scale_x_continuous(expand = c(0.001, 0))  +
  xlab(expression('c'[stay])) +
  annotate(
    geom = 'text',
    label = c(paste("  p =", p[10])),
    x = -Inf,
    y = Inf,
    hjust = -4.5,
    vjust = 1,
    parse = FALSE
  )
legend <- get_legend(p1 + theme(legend.position = "left"))
p1 <- p1 + non_pos
p2 <- p2 + non_pos
p3 <- p3 + non_pos
p4 <- p4 + non_pos


p5 <- p5 + non_pos
p6 <- p6 + non_pos
p7 <- p7 + non_pos
p8 <- p8 + non_pos
p9 <- p9 + non_pos
p10 <- p10 + non_pos



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

p <- plot_grid(pp, legend, rel_heights = c(3, .2), nrow = 2)


ggplot2::ggsave(
  "C:/Users/User/Documents/Polio_data_104/104vs107.pdf",
  plot = p,
  width = 12,
  height = 12,
  units = "in"
)