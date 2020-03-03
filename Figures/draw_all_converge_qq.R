#draw QQ plot comparing each of the fitted distributions (Fig. S2)

rm(list = ls())
library(ggplot2)
library(cowplot)
library(ggridges)
require(reshape2)


slope_est <- NULL
mid_est <- NULL
max_est <- NULL
lysis_est <- NULL
score_est <- NULL

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

get_csv_data <- function(drug) {
  RDA_file <-
    read.csv("data/all.csv", header = TRUE)

  RDA_data <- RDA_file


  #record the rows where mod 10 fit best
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

  d <- data.frame(slope, midpoint, max, lysis)

  colnames(d) <- c("Slope", "Midpoint", "Maximum", "Lysis")
  return(d)


}

#do QQ
refine <- function(x, y) {
  sx <- sort(x)
  sy <- sort(y)
  lenx <- length(sx)
  leny <- length(sy)
  if (leny < lenx)
    sx <- approx(1L:lenx, sx, n = leny)$y
  if (leny > lenx)
    sy <- approx(1L:leny, sy, n = lenx)$y
  data <- data.frame(obs = sx, exp = sy)
  return(data)
}

#read estimted data
exp_104 <- get_data(104)
exp_105 <- get_data(105)
exp_106 <- get_data(106)
exp_107 <- get_data(107)

#read experimental data
obs_104 <- get_csv_data("no drug")
obs_105 <- get_csv_data("10 nM rupintrivir")
obs_106 <- get_csv_data("50 µM 2'-C-meA")
obs_107 <- get_csv_data("2 μM Ganetespib")


make_plot <- function(exp, obs, color, label) {
  d <- refine(exp, obs)
  reg <- lm(obs ~ exp, data = d)

  m <- (ks.test(obs, exp))

  #do Ks test
  lm_eqn <- function(df) {
    m <- (ks.test(df$obs, df$exp))
    m$p
    eq <- substitute("  p =" ~ r2,
                     list(r2 = format(m$p, digits = 3)))
    as.character(as.expression(eq))

  }

  xl <- paste("est.\n", label)
  yl <- paste("obs.\n", label)

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

  g <-
    ggplot(data = d, aes(x = obs, y = exp)) + geom_point(colour = darkPalette[color]) + xlab(xl) + ylab(yl) +
    xlim(c(min(d), max(d))) + ylim(c(min(d), max(d))) +
    geom_abline(
      intercept = 0,
      slope = 1,
      size = 1,
      col = "black"
    ) +
    annotate(
      geom = 'text',
      label = lm_eqn(d),
      x = -Inf,
      y = Inf,
      hjust = 0,
      vjust = 1,
      parse = TRUE
    )
  return(g)

}



#draw all experiment together
all_exp <- function(exp, obs, color, label, letters) {
  g1 <- make_plot(exp$Slope, obs$Slope, color, "slope")
  g2 <- make_plot(exp$Maximum, obs$Maximum, color, "maximum")
  g3 <- make_plot(exp$Midpoint, obs$Midpoint, color, "midpoint")
  g4 <- make_plot(exp$Lysis, obs$Lysis, color, "lysis")


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

  p <- plot_grid(pp, nrow = 1)

  title <- ggdraw() + draw_label(label, fontface = 'bold')
  plot_grid(title, p, ncol = 1, rel_heights = c(0.2, 1)) # rel_heights values control title margins

  No_drug <- plot_grid(title, p, ncol = 1, rel_heights = c(0.2, 1))

  return(No_drug)

}

letters1 <- c("A", "B", "C", "D")
letters2 <- c("E", "F", "G", "H")
letters3 <- c("I", "J", "K", "L")
letters4 <- c("M", "N", "O", "P")
no_drug <- all_exp(exp_104, obs_104, 1, "no drug", letters1)
drug1 <- all_exp(exp_105, obs_105, 2, "rupintrivir", letters2)
drug2 <- all_exp(exp_106, obs_106, 3, "2'-C-meA", letters3)
drug3 <- all_exp(exp_107, obs_107, 4, "ganetespib", letters4)

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
  "all_qq_ks.pdf",
  plot = pp,
  width = 12,
  height = 12,
  units = "in"
)
