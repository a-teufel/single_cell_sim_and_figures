rm(list = ls())
library(ggplot2)
library(cowplot)
library(ggridges)
require(reshape2)
library(ggpubr)

#makes comparions of drug estimates to no drug estimates shown in figure 6 and S4
data <- NULL
names <- NULL

df <- NULL


d_104 = read.table(
  "data/104/Final_ABC_save_parameters.txt",
  header = TRUE
)
d_105 = read.table(
  "data/105/Final_ABC_save_parameters.txt",
  header = TRUE
)
d_106 = read.table(
  "data/106/Final_ABC_save_parameters.txt",
  header = TRUE
)
d_107 = read.table(
  "data/107/Final_ABC_save_parameters.txt",
  header = TRUE
)

all<-rbind(d_104,d_105,d_106,d_107)

colmins<-apply(all,2,min)
colmax<-apply(all,2,max)

n<-c(expression('log c'[tran]),expression('log c'[com]),expression('log c'[circ]),expression('log c'[rep~pos]),
     expression('log c'[rep~neg]),expression('log c'[pack]),expression('log com'[max]),expression('log rep'[max]),expression('log c'["3A"]),expression('c'[stay]))


d_104<-d_104[1:10]
#colnames(d_104)<-n

d_105<-d_105[1:10]
#colnames(d_105)<-n

d_106<-d_106[1:10]
#colnames(d_106)<-n


d_107<-d_107[1:10]
#colnames(d_107)<-n


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
  test <- ks.test(d_104[, i], d_105[, i])
  p_vals <- c(p_vals, test$p.value)
}

p <- format(p.adjust(p_vals, "bonferroni"), digits = 2)

all<-rbind(D_104,D_105,D_106,D_107)
names<-c("c.trans.","c.com.","c.circ.","c.rep..","c.rep...1","c.pack.","com.max.","rep.max.","c.3A.","c.stay.","Experiment")
colnames(all)<-names

all$Experiment<-factor(all$Experiment,levels=c("no drug","rupintrivir","2'-C-meA","ganetespib"))

m<-melt(all)
t<-m[m$variable == names[1],]

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


onepan<-function(i,m){

t<-m[m$variable == names[i],]

base<-t[t$Experiment == 'no drug',]

e105<-t[t$Experiment == 'rupintrivir',]

e106<-t[t$Experiment == "2'-C-meA",]

e107<-t[t$Experiment == "ganetespib",]

p_vals<-NULL
test <- ks.test(base$value, e105$value)
p_vals <- c(p_vals, test$p.value)
test <- ks.test(base$value, e106$value)
p_vals <- c(p_vals, test$p.value)
test <- ks.test(base$value, e107$value)
p_vals <- c(p_vals, test$p.value)

my_comparisons<-list(c('no drug', 'rupintrivir'),c('no drug', "2'-C-meA"),c('no drug', 'ganetespib'))


n<-c(expression('log c'[tran]),expression('log c'[com]),expression('log c'[circ]),expression('log c'[rep~pos]),
     expression('log c'[rep~neg]),expression('log c'[pack]),expression('log com'[max]),expression('log rep'[max]),expression('log c'["3A"]),expression('c'[stay]))

#n<-c(expression('var. log c'[tran]),expression('var. log c'[com]),expression('var. log c'[circ]),expression('var. log c'[rep~pos]),
#     expression('var. log c'[rep~neg]),expression('var. log c'[pack]),expression('var. log com'[max]),expression('var. log rep'[max]),expression('var. log c'["3A"]),expression('var. c'[stay]))


p <- ggplot(t, aes(Experiment,value,fill=Experiment))+ geom_violin() +scale_fill_manual(values = cbPalette) + scale_colour_manual(values = cbPalette)+
  stat_compare_means(aes(label = ..p.signif..),method = "wilcox.test", ref.group = "no drug",vjust =-.02,textsize = 5)+xlab(n[i])+ylab("estimate")+coord_cartesian(clip = 'off')


  #scale_fill_manual(values = cbPalette) + scale_colour_manual(values = cbPalette)+


#+stat_compare_means(comparisons = p_vals, label = "p.signif")
return(p)
}


m<-melt(all)


t<-m[m$variable == names[i],]

base<-t[t$Experiment == 'no drug',]

e105<-t[t$Experiment == 'rupintrivir',]

e106<-t[t$Experiment == "2'-C-meA",]

e107<-t[t$Experiment == "ganetespib",]

p_vals<-NULL
test <- ks.test(base$value, e105$value)
p_vals <- c(p_vals, test$p.value)
test <- ks.test(base$value, e106$value)
p_vals <- c(p_vals, test$p.value)
test <- ks.test(base$value, e107$value)
p_vals <- c(p_vals, test$p.value)

my_comparisons<-list(c('no drug', 'rupintrivir'),c('no drug', "2'-C-meA"),c('no drug', 'ganetespib'))
t<-m[m$variable == names[1],]


p <- ggplot(t, aes(Experiment,value,fill=Experiment))+ geom_violin() +scale_fill_manual(values = cbPalette) + scale_colour_manual(values = cbPalette)+
  stat_compare_means(aes(label = ..p.signif..),
                     method = "wilcox.test", ref.group = "no drug")

  geom_signif(comparisons = my_comparisons)

p

p1<-onepan(1,m)
p2<-onepan(2,m)
p3<-onepan(3,m)
p4<-onepan(4,m)
p5<-onepan(5,m)
p6<-onepan(6,m)
p7<-onepan(7,m)
p8<-onepan(8,m)
p9<-onepan(9,m)
p10<-onepan(10,m)


leg_pos <- theme(legend.position = "left")
non_pos <- theme(legend.position = "none")

sam<-p1

legend <- get_legend(sam + theme(legend.position = "top",legend.justification="top" ))
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



ppp <- plot_grid(pp, legend, rel_heights = c(3, .2),nrow = 2)


save_plot("violincmp.pdf",plot=ppp,base_height=12,units = "in")
