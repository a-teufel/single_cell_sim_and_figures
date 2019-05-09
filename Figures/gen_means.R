#draws generation means (Fig. S10)
rm(list = ls())
library(ggplot2)
library(cowplot)
library(ggridges)
require(reshape2)
library(plyr)

name <-
  paste(
    "C:/Users/User/Documents/Polio_data_104/estimation_iterations/105_Up2/Est_gen_ABC.txt"
  )
temp <- scan(name)
gen <- rbind(temp)
#the error code was zero, so just take these out
g <- gen[gen > 0]


name <-
  paste(
    "C:/Users/User/Documents/Polio_data_104/estimation_iterations/104_Up2/Est_gen_ABC.txt"
  )
temp <- scan(name)
gen1 <- rbind(temp)
g1 <- gen1[gen1 > 0]


name <-
  paste(
    "C:/Users/User/Documents/Polio_data_104/estimation_iterations/106_Up2/Est_gen_ABC.txt"
  )
temp <- scan(name)
gen2 <- rbind(temp)
g2 <- gen2[gen2 > 0]



name <-
  paste(
    "C:/Users/User/Documents/Polio_data_104/estimation_iterations/107_Up2/Est_gen_ABC.txt"
  )
temp <- scan(name)
gen3 <- rbind(temp)
g3 <- gen3[gen3 > 0]


D_104 <- data.frame(generations = g1)
D_105 <- data.frame(generations = g)
D_106 <- data.frame(generations = g2)
D_107 <- data.frame(generations = g3)

D_104 <- as.data.frame(D_104)
D_104$experiment <- "no drug  "

D_105 <- as.data.frame(D_105)
D_105$experiment <- "rupintrivir  "

D_106 <- as.data.frame(D_106)
D_106$experiment <- "2'-C-meA  "

D_107 <- as.data.frame(D_107)
D_107$experiment <- "ganetespib  "

#set up colors
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
darkPalette <-
  c("#4c4c4c",
    "#734f00",
    "#12608d",
    "#004f39",
    "#8d840b",
    "#003859",
    "#6a2f00")

all_alpha = 0.1
all_size = .75
line_size = 1

#do all combinations and add in p-values
combined <- rbind.fill(D_104, D_105)

p <-
  format(ks.test(D_104$generations, D_105$generations)$p.value,
         digits = 2)

g1 <-
  ggplot(combined,
         aes(generations, colour = experiment, fill = experiment)) + geom_density(alpha =
                                                                                    all_alpha, size = all_size) + scale_fill_manual(values = cbPalette) + scale_colour_manual(values =
                                                                                                                                                                                cbPalette) + scale_y_continuous(expand = c(0.001, 0)) + scale_x_continuous(expand = c(0.001, 0)) +
  xlim(c(0, 10)) +
  geom_vline(
    size = line_size,
    xintercept = c(median(D_104$generations), median(D_105$generations)),
    linetype = "dashed",
    color = c(cbPalette[1], cbPalette[2])
  ) +
  annotate(
    geom = 'text',
    label = c(paste("  p =", p)),
    x = -Inf,
    y = Inf,
    hjust = 0,
    vjust = 1,
    parse = FALSE
  )


combined <- rbind.fill(D_104, D_106)

p <-
  format(ks.test(D_104$generations, D_106$generations)$p.value,
         digits = 2)

g2 <-
  ggplot(combined,
         aes(generations, colour = experiment, fill = experiment)) + geom_density(alpha =
                                                                                    all_alpha, size = all_size) + scale_fill_manual(values = c(cbPalette[3], cbPalette[1])) + scale_colour_manual(values =
                                                                                                                                                                                                    c(cbPalette[3], cbPalette[1])) + scale_y_continuous(expand = c(0.001, 0)) +
  scale_x_continuous(expand = c(0.001, 0)) + xlim(c(0, 10)) +
  geom_vline(
    size = line_size,
    xintercept = c(median(D_106$generations), median(D_104$generations)),
    linetype = "dashed",
    color = c(cbPalette[3], cbPalette[1])
  ) +
  annotate(
    geom = 'text',
    label = c(paste("  p =", p)),
    x = -Inf,
    y = Inf,
    hjust = 0,
    vjust = 1,
    parse = FALSE
  )


combined <- rbind.fill(D_104, D_107)

p <-
  format(ks.test(D_104$generations, D_107$generations)$p.value,
         digits = 2)

g3 <-
  ggplot(combined,
         aes(generations, colour = experiment, fill = experiment)) + geom_density(alpha =
                                                                                    all_alpha, size = all_size) + scale_fill_manual(values = c(cbPalette[4], cbPalette[1])) + scale_colour_manual(values =
                                                                                                                                                                                                    c(cbPalette[4], cbPalette[1])) + scale_y_continuous(expand = c(0.001, 0)) +
  scale_x_continuous(expand = c(0.001, 0)) + xlim(c(0, 10)) +
  geom_vline(
    size = line_size,
    xintercept = c(median(D_107$generations), median(D_104$generations)),
    linetype = "dashed",
    color = c(cbPalette[4], cbPalette[1])
  ) +
  annotate(
    geom = 'text',
    label = c(paste("  p =", p)),
    x = -Inf,
    y = Inf,
    hjust = 0,
    vjust = 1,
    parse = FALSE
  )


combined <- rbind.fill(D_105, D_106)

p <-
  format(ks.test(D_105$generations, D_106$generations)$p.value,
         digits = 2)

g4 <-
  ggplot(combined,
         aes(generations, colour = experiment, fill = experiment)) + geom_density(alpha =
                                                                                    all_alpha, size = all_size) + scale_fill_manual(values = c(cbPalette[3], cbPalette[2])) + scale_colour_manual(values =
                                                                                                                                                                                                    c(cbPalette[3], cbPalette[2])) + scale_y_continuous(expand = c(0.001, 0)) +
  scale_x_continuous(expand = c(0.001, 0)) + xlim(c(0, 10)) +
  geom_vline(
    size = line_size,
    xintercept = c(median(D_106$generations), median(D_105$generations)),
    linetype = "dashed",
    color = c(cbPalette[3], cbPalette[2])
  ) +
  annotate(
    geom = 'text',
    label = c(paste("  p =", p)),
    x = -Inf,
    y = Inf,
    hjust = 0,
    vjust = 1,
    parse = FALSE
  )


combined <- rbind.fill(D_105, D_107)

p <-
  format(ks.test(D_105$generations, D_107$generations)$p.value,
         digits = 2)

g5 <-
  ggplot(combined,
         aes(generations, colour = experiment, fill = experiment)) + geom_density(alpha =
                                                                                    all_alpha, size = all_size) + scale_fill_manual(values = c(cbPalette[4], cbPalette[2])) + scale_colour_manual(values =
                                                                                                                                                                                                    c(cbPalette[4], cbPalette[2])) + scale_y_continuous(expand = c(0.001, 0)) +
  scale_x_continuous(expand = c(0.001, 0)) + xlim(c(0, 10)) +
  geom_vline(
    size = line_size,
    xintercept = c(median(D_105$generations), median(D_107$generations)),
    linetype = "dashed",
    color = c(cbPalette[2], cbPalette[4])
  ) +
  annotate(
    geom = 'text',
    label = c(paste("  p =", p)),
    x = -Inf,
    y = Inf,
    hjust = 0,
    vjust = 1,
    parse = FALSE
  )

combined <- rbind.fill(D_106, D_107)

p <-
  format(ks.test(D_106$generations, D_107$generations)$p.value,
         digits = 2)

g6 <-
  ggplot(combined,
         aes(generations, colour = experiment, fill = experiment)) + geom_density(alpha =
                                                                                    all_alpha, size = all_size) + scale_fill_manual(values = c(cbPalette[3], cbPalette[4])) + scale_colour_manual(values =
                                                                                                                                                                                                    c(cbPalette[3], cbPalette[4])) + scale_y_continuous(expand = c(0.001, 0)) +
  scale_x_continuous(expand = c(0.001, 0)) + xlim(c(0, 10)) +
  geom_vline(
    size = line_size,
    xintercept = c(median(D_106$generations), median(D_107$generations)),
    linetype = "dashed",
    color = c(cbPalette[3], cbPalette[4])
  ) +
  annotate(
    geom = 'text',
    label = c(paste("  p =", p)),
    x = -Inf,
    y = Inf,
    hjust = 0,
    vjust = 1,
    parse = FALSE
  )


pp <-
  plot_grid(
    g1 + theme(legend.position = "bottom"),
    g2 + theme(legend.position = "bottom"),
    g3 + theme(legend.position = "bottom"),
    g4 + theme(legend.position = "bottom"),
    g5 + theme(legend.position = "bottom"),
    g6 + theme(legend.position = "bottom"),
    labels = c("A", "B", "C", "D", "E", "F"),
    align = 'vh',
    hjust = -1,
    nrow = 3
  )

ggplot2::ggsave(
  "C:/Users/User/Documents/Polio_data_104/generations_means.pdf",
  plot = pp,
  width = 12,
  height = 12,
  units = "in"
)
