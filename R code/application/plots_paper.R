load("~/college_application.Rdata")
load("~/community_college_application.Rdata")
load("~/high_school_application.Rdata")
df <- read.csv("~/AUM_edu_pm25_X_county.csv")

#libraries
library(ggplot2)
library(dplyr)
library(forcats)
library(hrbrthemes)
library(viridis)
library(radarchart)
library(fmsb)

#clean dataset
high_school_application<-application
df <- subset(df, select = -c(X,Unnamed..0))
df <- na.omit(df)
#standardize values that are not percentages
df$med_hhinc1990 <- scale(df$med_hhinc1990)
df$popdensity2000 <- scale(df$popdensity2000)
df$mean_winter_prcp <- scale(df$mean_winter_prcp)
df$mean_winter_tmin <- scale(df$mean_winter_tmin)
df$mean_summer_prcp <- scale(df$mean_summer_prcp)
df$mean_summer_tmax <- scale(df$mean_summer_tmax)
#create the treatment
df$pm25_1982_binary <- ifelse(df$pm25_1982 >= median(df$pm25_1982),1,0)
# create confounders
X <- cbind(1,df$frac_coll_plus2000, df$med_hhinc1990, df$popdensity2000, df$poor_share1990,
           df$share_black2000, df$share_white2000, df$share_hisp2000, df$share_asian2000,
           df$emp2000, df$census_region, df$mean_winter_prcp, df$mean_winter_tmin,
           df$mean_summer_prcp, df$mean_summer_tmax)

#############################################################################

M_hs <- rowMeans(high_school_application$P1_POST)-
  rowMeans(high_school_application$P0_POST)
ATE_hs<- rowMeans(high_school_application$Y1_POST)-
  rowMeans(high_school_application$Y0_POST)
matrix_hs <- as.data.frame(cbind(M_hs, ATE_hs, X))

dissociative_hs <- matrix_hs[which(M_hs<=0.01&M_hs>=-0.01),]
associative_plus_hs <- matrix_hs[which(M_hs > 0.01),]
associative_minus_hs <- matrix_hs[which(M_hs <= -0.01),]

strata_hs<-rep("EDE", 3009)
strata_hs[which(M_hs > 0.01)]<-"EAE+"
strata_hs[which(M_hs <= -0.01)]<-"EAE-"

M_cc <- rowMeans(community_college_application$P1_POST)-
  rowMeans(community_college_application$P0_POST)
ATE_cc<- rowMeans(community_college_application$Y1_POST)-
  rowMeans(community_college_application$Y0_POST)
matrix_cc <- as.data.frame(cbind(M_cc, ATE_cc, X))

dissociative_cc <- matrix_cc[which(M_cc<=0.01 & M_cc>=-0.01),]
associative_plus_cc <- matrix_cc[which(M_cc > 0.01),]
associative_minus_cc <- matrix_cc[which(M_cc <= -0.01),]

strata_cc<-rep("EDE", 3009)
strata_cc[which(M_cc > 0.01)]<-"EAE+"
strata_cc[which(M_cc <= -0.01)]<-"EAE-"

M_c <- rowMeans(college_application$P1_POST)-
  rowMeans(college_application$P0_POST)
ATE_c <- rowMeans(college_application$Y1_POST)-
  rowMeans(college_application$Y0_POST)
matrix_c <- as.data.frame(cbind(M_c, ATE_c, X))

dissociative_c <- matrix_c[which(M_c<=0.01 & M_c>=-0.01),]
associative_plus_c <- matrix_c[which(M_c > 0.01),]
associative_minus_c <- matrix_c[which(M_c <= -0.01),]

strata_c<-rep("EDE", 3009)
strata_c[which(M_c > 0.01)]<-"EAE+"
strata_c[which(M_c <= -0.01)]<-"EAE-"

#############################################################################
effects<- c("EAE+","EAE-","EDE")

post_pcf_hs<-sapply(effects, function(i) colMeans(high_school_application$Y1_POST[which(strata_hs==i),]-
                      high_school_application$Y0_POST[which(strata_hs==i),]))
post_pcf_c<-sapply(effects, function(i) colMeans(college_application$Y1_POST[which(strata_c==i),]-
                                                   college_application$Y0_POST[which(strata_c==i),]))
post_pcf_cc<-sapply(effects, function(i) colMeans(community_college_application$Y1_POST[which(strata_cc==i),]-
                                                    community_college_application$Y0_POST[which(strata_cc==i),]))


dataset_posterior<- as.data.frame(cbind(values=c(post_pcf_hs, post_pcf_c, post_pcf_cc),
                                      stratum=rep(rep(effects, each=5000), 3),
                                      school=c(rep("high school",15000),
                                               rep("community college",15000),
                                               rep("college",15000))))
dataset_posterior$values<-as.numeric(dataset_posterior$values)
dataset_posterior$stratum=ordered(dataset_posterior$stratum, levels=c("EAE-","EDE","EAE+"))
dataset_posterior$school=ordered(dataset_posterior$school, levels=c("high school","community college","college"))


cbPalette <- c("#D90224","#EBB600","#00944B")

ggplot(dataset_posterior, aes(x=school, y=values, fill=stratum)) + 
  scale_fill_manual(values=cbPalette, name=" ")+   ####REMOVE FOR standard color
  #scale_color_manual(values=cbPalette, name="")+   ####REMOVE FOR standard color
  geom_boxplot(lwd=0.3,fatten = 1.5, outlier.size = 0.3)+
  geom_hline(yintercept = 0, col="#80dfff", size=0.4) +
  theme(panel.background = element_rect(fill='white'),
        plot.background = element_rect(fill ="white"),
        #panel.grid.minor = element_line(color = "grey"),
        axis.title = element_text(size=14),
        legend.position = "right",
        legend.text=element_text(size=10),
        legend.title=element_text(size=10),
        plot.title = element_text(hjust = 0.2),
        title =element_text(size=18),
        legend.background = element_rect(fill='transparent'),
        panel.grid.major = element_line(color = "grey",size = 0.1))+
  ylab("Principal causal effects") +
  xlab("") 


# post_princCE

summary(post_pcf_hs)
summary(post_pcf_cc)
summary(post_pcf_c)
apply(post_pcf_hs,2,quantile,prob=c(0.025,0.975))
apply(post_pcf_cc,2,quantile,prob=c(0.025,0.975))
apply(post_pcf_c,2,quantile,prob=c(0.025,0.975))


##############################################################
# ---    spiderplot ----
##############################################################


spiderplots_schools<-function(X, causal_effects, title_edu, value_max){
  
  X_matrix <- cbind(eae_n = apply(X[which(causal_effects=='EAE-'),-1],2, mean),
                    ade = apply(X[which(causal_effects=='EDE'),-1],2, mean),
                    eae_p = apply(X[which(causal_effects=='EAE+'),-1],2, mean))
  
  data_spider=as.data.frame(rbind(max=rep(1,dim(X_matrix)[1]),min=rep(-1,dim(X_matrix)[1]),
                                  t(X_matrix),
                                  mean=apply(X[,-1],2,mean)))
  colnames(data_spider)=c('frac_coll_plus2000', 'med_hhinc1990', 'popdensity2000', 'poor_share1990',
                          'share_black2000', 'share_white2000', 'share_hisp2000', 'share_asian2000',
                          'emp2000', 'census_region', 'mean_winter_prcp', 'mean_winter_tmin',
                          'mean_summer_prcp', 'mean_summer_tmax')
  rownames(data_spider)=c("max","min","EAE -","EDE","EAE +","mean pop.")
  
  data_spider<-data_spider[,-10]
  
  data_spider[3,]<-(data_spider[3,]-data_spider[6,])*100
  data_spider[4,]<-(data_spider[4,]-data_spider[6,])*100
  data_spider[5,]<-(data_spider[5,]-data_spider[6,])*100
  data_spider[6,]<- 0
  data_spider[1,]<- (-value_max)
  data_spider[2,]<- value_max
  
  seq_perc<-c(paste0('-',value_max,'%'),paste0('-',value_max/2,'%'),
              ' ',paste0(value_max,'%'),paste0(value_max/2,'%'))
  
  # a spider plot for each group:
  par(mfrow=c(1,1))
  for(cl in 1:3){
    pdf(file=paste0("spiderplot_",title_edu,cl,"perc.pdf"))
    data_cl= data_spider[c(1,2,6,cl+2),]
    radarchart( data_cl  , axistype=1 , 
                pcol=c("grey",cbPalette[cl]), plwd=4 ,plty=1,
                cglcol="grey", cglty=1, 
                axislabcol="black", caxislabels=seq_perc, calcex=0.7,
                cglwd=0.8, vlcex=0.8 
    )
    legend(x=0.97, y=-1.05, legend=rownames(data_spider[c(cl+2,6),]), bty = "n", pch=20 , 
           col=c(cbPalette[cl],"grey") ,  cex=0.8, pt.cex=3)
    title(main=title_edu)
    dev.off()
  }
}

spiderplots_schools(X, causal_effects=strata_hs, title_edu='High School', value_max=30)
spiderplots_schools(X, causal_effects=strata_cc, title_edu='Comunity College', value_max=50)
spiderplots_schools(X, causal_effects=strata_hs, title_edu='College', value_max=30)

