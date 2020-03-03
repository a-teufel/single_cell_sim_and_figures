#makes figure of inital input data from experiments, The bottom part of figure 2 in the paper

rm(list=ls())
library(ggplot2)
library(cowplot)
library(ggridges)
require(reshape2)

data<-NULL
names<-NULL
slope_est<-NULL
mid_est<-NULL
max_est<-NULL
lysis_est<-NULL
score_est<-NULL

#read in data and get data just for the drug case we want
get_csv_data<-function(drug){

  RDA_file<-read.csv("data/all.csv",header=TRUE)

  RDA_data<-RDA_file

  #record the rows where mod 10 fit best
  rows_where_5_fit_best<-NULL
  rows_where_10_fit_best<-NULL
  rows_where_line_fit_best<-NULL

  #get just the lines with data
  for( i in 1:length(RDA_data[,1])){

    if(RDA_data[i,]$decision_bio == "infection&lysis" & RDA_data[i,]$dataSet_bio== drug  ){
      rows_where_10_fit_best<-c(rows_where_10_fit_best,i)
    }
    if(RDA_data[i,]$decision_bio == "infection" & RDA_data[i,]$dataSet_bio== drug  ){
      rows_where_5_fit_best<-c(rows_where_5_fit_best,i)
    }
  }


  RDA_data_just_5_mod<-RDA_data[rows_where_5_fit_best,]

  RDA_data_just_10_mod<-RDA_data[rows_where_10_fit_best,]

  max<-c(RDA_data_just_5_mod$COMB_maximum_y,RDA_data_just_10_mod$COMB_maximum_y)
  slope<-c(RDA_data_just_5_mod$COMB_slope1,RDA_data_just_10_mod$COMB_slope1)
  midpoint<-c(RDA_data_just_5_mod$COMB_midPoint1_x,RDA_data_just_10_mod$COMB_midPoint1_x)
  lysis<-c(rep(24,nrow(RDA_data_just_5_mod)), RDA_data_just_10_mod$COMB_startDeclinePoint_x)

  d<-data.frame(slope,midpoint,max,lysis)

  colnames(d)<-c("Slope","Midpoint","Maximum","Lysis")
  return(d)
}



all_alpha=0.2
all_size=1
#inside fill colors (each drug has a diffrent color, so if you want to draw the input data for the drug cases you need to change the indexs in the scale_fill_manual calls)
cbPalette <- c("#999999","#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#line fill colors
darkPalette<-cbPalette#c("#4c4c4c","#734f00","#12608d","#004f39","#8d840b","#003859","#6a2f00")


#Slope
obs_104<-get_csv_data("no drug")
slope<-obs_104$Slope
slope<-as.data.frame(slope)
d<-melt(slope)
name<-"slope"
col=1
p1<-ggplot(d, aes(x=value,fill=variable,color = variable)) +geom_density(alpha=all_alpha,size=all_size,adjust=.65)+scale_fill_manual(values=c(darkPalette[col],(cbPalette[col])))+ scale_colour_manual(values=c(darkPalette[col],(cbPalette[col]))) +scale_y_continuous(expand = c(0.001, 0))+scale_x_continuous(expand = c(0.001, 0),breaks = scales::pretty_breaks(n = 5)) +xlab(name)


#Midpoint
Midpoint<-obs_104$Midpoint
Midpoint<-as.data.frame(Midpoint)
d<-melt(Midpoint)
name<-"midpoint"
p2<-ggplot(d, aes(x=value,fill=variable,color=variable)) +geom_density(alpha=all_alpha,size=all_size,adjust=.65)+scale_fill_manual(values=c(darkPalette[col],(cbPalette[col])))+ scale_colour_manual(values=c(darkPalette[col],(cbPalette[col]))) +scale_y_continuous(expand = c(0.001, 0))+scale_x_continuous(expand = c(0.001, 0),breaks = scales::pretty_breaks(n = 5)) +xlab(name)


#Max
Maximum<-obs_104$Maximum
Maximum<-as.data.frame(Maximum)
d<-melt(Maximum)
name<-"maximum"
p3<-ggplot(d, aes(x=value,fill=variable,color = variable)) +geom_density(alpha=all_alpha,size=all_size,adjust=.65)+scale_fill_manual(values=c(darkPalette[col],(cbPalette[col])))+ scale_colour_manual(values=c(darkPalette[col],(cbPalette[col]))) +scale_y_continuous(expand = c(0.001, 0))+scale_x_continuous(expand = c(0.001, 0),breaks = scales::pretty_breaks(n = 5)) +xlab(name)

#Lysis
Lysis<-obs_104$Lysis
Lysis<-as.data.frame(Lysis)
d<-melt(Lysis)
name<-"lysis"
p4<-ggplot(d, aes(x=value,fill=variable,color = variable)) +geom_density(alpha=all_alpha,size=all_size,adjust=.65)+scale_fill_manual(values=c(darkPalette[col],(cbPalette[col])))+ scale_colour_manual(values=c(darkPalette[col],(cbPalette[col]))) +scale_y_continuous(expand = c(0.001, 0))+scale_x_continuous(expand = c(0.001, 0),breaks = scales::pretty_breaks(n = 5)) +xlab(name)


#remove from other plots
non_pos<- theme(legend.position = "none")
g1<-p1+non_pos
g2<-p2+non_pos
g3<-p3+non_pos
g4<-p4+non_pos

pp<- plot_grid( g1,g3,g2,g4,
                align = 'vh',
                hjust = -1,
                nrow = 2,labels=c("A","B","C","D"))

#look at plot to make sure it looks right
pp

ggplot2::ggsave("just_input.pdf",plot=pp,width = 12, height = 9,units = "in")
