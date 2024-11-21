load("~/Downloads/college_application.Rdata")
load("~/Downloads/community_college_application.Rdata")
load("~/Downloads/high_school_application.Rdata")
df <- read.csv("~/Documents/GitHub/StrataBayes/R code/data/AUM_edu_pm25_X_county.csv")

high_school_application<-application
#clean dataset
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

dataset_boxplot<- as.data.frame(cbind(values=c(ATE_hs, ATE_cc, ATE_c),
                        stratum=c(strata_hs,strata_cc,strata_c),
                        school=c(rep("high school",3009),
                                 rep("community college",3009),
                                 rep("private college",3009))))
dataset_boxplot$values<-as.numeric(dataset_boxplot$values)
dataset_boxplot$stratum=ordered(dataset_boxplot$stratum, levels=c("EAE-","EDE","EAE+"))

cbPalette <- c("#D90224","#EBB600","#00944B")

ggplot(dataset_boxplot, aes(x=school, y=values, fill=stratum)) + 
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
