#draws a comparison of the known distributions of slope, max, mid, and lysis to the estimated distributions for each of the experiments (Fig. S1)

rm(list = ls())
library(ggplot2)
library(cowplot)
library(ggridges)
require(reshape2)


slope_est <- NULL
mid_est <- NULL
max_est <- NULL
lysis_est <- NULL


#gets estimated final distributions for slope, max, mid, lysis time
get_data <- function(exp) {
  slope_est <- NULL
  mid_est <- NULL
  max_est <- NULL
  lysis_est <- NULL
  score_est <- NULL

  name <-
    paste(
      "data/",exp,"/Final_ABC_B11.txt",
      sep = ""
    )
  temp <- scan(name)
  data <- data.frame(d = as.vector(temp))
  slope_est <- rbind(slope_est, data)


  name <-
    paste(
      "data/",exp,"/Final_ABC_M11.txt",
      sep = ""
    )
  temp <- scan(name)
  data <- data.frame(d = as.vector(temp))
  mid_est <- rbind(mid_est, data)

  name <-
    paste(
      "data/",exp,"/Final_ABC_Ka.txt",
      sep = ""
    )
  temp <- scan(name)
  data <- data.frame(d = as.vector(temp))
  max_est <- rbind(max_est, data)

  name <-
    paste(
      "data/",exp,"/Final_ABC_lysis.txt",
      sep = ""
    )
  temp <- scan(name)
  data <- data.frame(d = as.vector(temp))
  lysis_est <- rbind(lysis_est, data)

  d <-
    data.frame(
      slope = slope_est,
      mid = mid_est,
      max = max_est,
      lysis = lysis_est
    )
  colnames(d) <- c("Slope", "Midpoint", "Maximum", "Lysis")
  return(d)
}

#get experimental data for drug of intrest
get_csv_data <- function(drug) {
  RDA_file <-
    read.csv("data/all.csv", header = TRUE)

  RDA_data <- RDA_file

  rows_where_5_fit_best <- NULL
  rows_where_10_fit_best <- NULL
  rows_where_line_fit_best <- NULL

  #get just the lines with data
  for (i in 1:length(RDA_data[, 1])) {
    if (RDA_data[i, ]$decision_bio == "infection&lysis" &
        RDA_data[i, ]$dataSet_bio == drug) {
      rows_where_10_fit_best <- c(rows_where_10_fit_best, i)
    }
    if (RDA_data[i, ]$decision_bio == "infection" &
        RDA_data[i, ]$dataSet_bio == drug) {
      rows_where_5_fit_best <- c(rows_where_5_fit_best, i)
    }
  }


  RDA_data_just_5_mod <- RDA_data[rows_where_5_fit_best, ]

  RDA_data_just_10_mod <- RDA_data[rows_where_10_fit_best, ]

  max <-
    c(RDA_data_just_5_mod$COMB_maximum_y,
      RDA_data_just_10_mod$COMB_maximum_y)
  slope <-
    c(RDA_data_just_5_mod$COMB_slope1,
      RDA_data_just_10_mod$COMB_slope1)
  midpoint <-
    c(RDA_data_just_5_mod$COMB_midPoint1_x,
      RDA_data_just_10_mod$COMB_midPoint1_x)
  lysis <-
    c(rep(24, nrow(RDA_data_just_5_mod)),
      RDA_data_just_10_mod$COMB_startDeclinePoint_x)

  print(max)
  print(slope)
  print(midpoint)
  print(lysis)
  d <- data.frame(slope, midpoint, max, lysis)

  colnames(d) <- c("Slope", "Midpoint", "Maximum", "Lysis")
  return(d)


}

#draws 2 distributions together
draw_exp <- function(exp, obs, col, name, Tit) {
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


  data_k <- data.frame(d = as.vector(obs), name = "observed")
  slope_est <- data.frame(d = as.vector(exp), name = "estimated")
  df <- rbind(data_k, slope_est)

  set_scale <- 3
  d <- data.frame(df)

  colnames(d) <- c("x", "y")


  d$y <- factor(d$y,
                levels = c("estimated", "observed"))

  all_alpha = 0.2
  all_size = 1

  number_ticks <- function(n) {
    function(limits)
      pretty(limits, n)
  }

  p1 <-
    ggplot(d, aes(x = x, colour = y, fill = y)) + geom_density(alpha = all_alpha,
                                                               size = all_size,
                                                               adjust = .65) + scale_fill_manual(values = c(darkPalette[col], (cbPalette[col]))) + scale_colour_manual(values =
                                                                                                                                                                         c(darkPalette[col], (cbPalette[col]))) + scale_y_continuous(expand = c(0.01, 0)) +
    scale_x_continuous(expand = c(0.01, 0),
                       breaks = scales::pretty_breaks(n = 5)) + xlab(name)

  return(p1)

}


#draws all of an experiment together
all_exp <- function(exp, obs, color, label, letters) {
  g1 <- draw_exp(exp$Slope, obs$Slope, color, "slope")
  g2 <- draw_exp(exp$Maximum, obs$Maximum, color, "maximum")
  g3 <- draw_exp(exp$Midpoint, obs$Midpoint, color, "midpoint")
  g4 <- draw_exp(exp$Lysis, obs$Lysis, color, "lysis")

  leg_pos <-
    theme(legend.position = "right", legend.title = element_blank())
  non_pos <- theme(legend.position = "none")
  g0 <- g1

  g1 <- g1 + non_pos
  g2 <- g2 + non_pos
  g3 <- g3 + non_pos
  g4 <- g4 + non_pos
  g5 <- g0 + leg_pos


  legend <- get_legend(g5)

  pp <- plot_grid(
    g1,
    g2,
    g3,
    g4,
    align = 'vh',
    hjust = -1,
    nrow = 1,
    labels = letters
  )

  p <- plot_grid(pp, legend, rel_widths = c(4, .5), nrow = 1)

  title <- ggdraw() + draw_label(label, fontface = 'bold')
  plot_grid(title, p, ncol = 1, rel_heights = c(0.2, 1)) # rel_heights values control title margins

  No_drug <- plot_grid(title, p, ncol = 1, rel_heights = c(0.2, 1))

  return(No_drug)

}

exp_104 <- get_data(104)
exp_105 <- get_data(105)
exp_106 <- get_data(106)
exp_107 <- get_data(107)

obs_104 <- get_csv_data("no drug")
obs_105 <- get_csv_data("10 nM rupintrivir")
obs_106 <- get_csv_data("50 µM 2'-C-meA")
obs_107 <- get_csv_data("2 μM Ganetespib")


letters1 <- c("A", "B", "C", "D")
letters2 <- c("E", "F", "G", "H")
letters3 <- c("I", "J", "K", "L")
letters4 <- c("M", "N", "O", "P")

no_drug <- all_exp(exp_104, obs_104, 1, "no drug", letters1)
drug1 <- all_exp(exp_105, obs_105, 2, "rupintrivir", letters2)
drug2 <- all_exp(exp_106, obs_106, 3, "2'-C-meA", letters3)
drug3 <- all_exp(exp_107, obs_107, 4, "Ganetespib", letters4)

pp <- plot_grid(
  no_drug,
  drug1,
  drug2,
  drug3,
  align = 'vh',
  hjust = -1,
  nrow = 4
)

ggplot2::ggsave(
  "all_con.pdf",
  plot = pp,
  width = 12,
  height = 9,
  units = "in"
)
