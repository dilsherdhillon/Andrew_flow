## Import the required libraries 
library(tidyverse)
library(corrplot)
## read in data 
dta <- read_csv("~/Box/Dilsher/Andrew/Flow project/Data Files/CD4_Master_MFI.csv")
dta<-dta%>%
  mutate(Sample=X1)%>%
  select(-c(X1))
## DATA WRANGLING ##
#sapply(strsplit(dta$Sample,"_"), "[", 1) ## can use this instead of gsub 
## parse out the subject names, create a new column that indicates whether the sample is nil or not ## 
dta<-
  dta%>%
  mutate(sample_group=gsub( "_.*$", "",Sample))%>%
  mutate(sample_group_2=gsub(" .*$", "",sample_group))%>%
  mutate(nil=ifelse(grepl("nil",Sample,ignore.case = TRUE),"Yes","No"))%>%
  modify_at(c(1:9),as.numeric)  %>%
{.}

## Some samples have been run twice - why were they run twice and which values should I take?
## 01082019, met with andrew to discuss this and he suggested I send him the sample numbers 
## and he will look into them. 
##
dta%>%
  mutate(sample_group=gsub( "_.*$", "",Sample))%>%
  mutate(sample_group_2=gsub(" .*$", "",sample_group))%>%
  count(sample_group_2)%>%
  filter(n>5)%>%
  inner_join(.,dta)%>%
  write_csv(.,path="~/Box/Dilsher/Andrew/Flow project/Data Files/Samples_Run_MoreThanOnce.csv")
### Sent this file to andrew - he gave me the samples that need to be excluded ###
### Re-read the data and process new data ### January 9th 2019 



dta_v2 <- read_csv("~/Box/Dilsher/Andrew/Flow project/Data Files/CD4_Master_MFI_v2.csv")
dta_v2<-dta_v2%>%
  mutate(Sample=X1)%>%
  select(-c(X1))
## DATA WRANGLING ##
#sapply(strsplit(dta$Sample,"_"), "[", 1) ## can use this instead of gsub 
## parse out the subject names, create a new column that indicates whether the sample is nil or not ## 
dta_v2<-
  dta_v2%>%
  mutate(sample_group=gsub( "_.*$", "",Sample))%>%
  mutate(sample_group_2=gsub(" .*$", "",sample_group))%>%
  mutate(nil=ifelse(grepl("nil",Sample,ignore.case = TRUE),"Yes","No"))%>%
  modify_at(c(1:9),as.numeric)  %>%

    group_by(sample_group)%>%
    arrange(desc(nil),.by_group=TRUE)%>%
    mutate_at(.funs = funs(fc=./.[1]),.vars = (1:9))%>%
  ungroup()


## Now, we create a new df, which only has the 
## 1. correlations and sample info
## 2. Remove the "nil" samples 
## 3. Remove the sample group columns as well 
fc_dta<-dta_v2%>%
  select(-c(1:9,sample_group,sample_group_2))%>%
  filter(nil=="No")%>%
  select(-c(nil))%>%
  mutate(stimulus=case_when(grepl("BCG",Sample,ignore.case = TRUE)~"BCG",
                            grepl("SEB",Sample,ignore.case = TRUE)~"SEB",
                                      grepl("HIV",Sample,ignore.case = TRUE)~"HIV",
                                                grepl("E,2f,C",Sample,ignore.case = TRUE)~"E2FC",
                                                          FALSE~"NA"))
## create individual dfs for each stimulus 
bcg<-fc_dta%>%
  filter(stimulus=="BCG")%>%
  select(-c(stimulus,Sample))
seb<-fc_dta%>%
  filter(stimulus=="SEB")%>%
  select(-c(stimulus,Sample))
hiv<-fc_dta%>%
  filter(stimulus=="HIV")%>%
  select(-c(stimulus,Sample))
e2fc<-fc_dta%>%
  filter(stimulus=="E2FC")%>%
  select(-c(stimulus,Sample))
## Individual dfs created ##################



###############################################################################################
####################### #######################PLOTTING ####################### #######################

## create correlation dfs ##
c<-cor(bcg,use="pairwise.complete.obs",method = "spearman")
#corrplot(c, type = "upper",order="hclust",hclust.method="ward")

## p-value 
p<-cor.mtest(bcg)[[1]]
p_adj<-p.adjust(p,method = "BH") ### BH correction multiple testing 
p_adj_mat<-matrix(p_adj,nrow = 9,ncol=9) ## convert into matrix 
rownames(p_adj_mat)<-rownames(c) ## 
colnames(p_adj_mat)<-colnames(c)
#corrplot(p_adj_mat,type="upper")

## clustering
dd<-as.dist((1-abs(c)))
#plot(hclust(dd,method = "average"))

par(mfrow=c(2,2))
corrplot(c, type = "upper")
corrplot(p_adj_mat,type="upper",method="number",cl.lim = c(0, 1))
plot(hclust(dd,method = "average"))
graphics.off()
rm(c,p,p_adj,p_adj_mat,dd)
###############################################################################################




###############################################################################################
###################### play with colors ###################### ###################### 
col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
                           "cyan", "#007FFF", "blue", "#00007F"))
corrplot(p_adj_mat,type="upper",cl.lim = c(0, 1),col=col1(50))


## add p-values to plot? 
corrplot(c, type = "upper",order="hclust",p.mat = p_adj_mat, sig.level = 0.05,insig = "p-value")

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(c, method="color", col=col(200),  
         diag=FALSE, # tl.pos="d", 
         type="upper",
         title="BCG Stimulus", 
         addCoef.col = "black", # Add coefficient of correlation
         # Combine with significance
         p.mat = p_adj_mat, sig.level = 0.05, insig = "blank" 
         # hide correlation coefficient on the principal diagonal
)


###############################################################################################








###############################################################################################
#### write a function ###
dd_correlation_plots<-function(df){
  ## The data frame df should be a numeric df with the variables listed and the rows are simply the rownames ##
  
  c<-cor(df,use="pairwise.complete.obs",method = "spearman")
  p<-cor.mtest(df)[[1]]
  p_adj<-p.adjust(p,method = "BH") ### BH correction multiple testing 
  p_adj_mat<-matrix(p_adj,nrow = 9,ncol=9) ## convert into matrix 
  rownames(p_adj_mat)<-rownames(c) ## 
  colnames(p_adj_mat)<-colnames(c)
  ## clustering
  dd<-as.dist((1-abs(c)))
  
  #dir<-getwd()
  #dir.create(getwd)
  col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
                             "cyan", "#007FFF", "blue", "#00007F"))
  
  pdf(file.path(getwd(),paste0(deparse(substitute(df)),".pdf")))
  par(mfrow=c(2,2))
  corrplot(c, type = "upper",method="number",col=col1(20),title="Spearman Correlation",mar=c(0,0,1,0),number.cex = .7)
  corrplot(p_adj_mat,type="upper",method="number",cl.lim = c(0, 1),col=col1(20),title="BH Adjusted p-values",mar=c(0,0,1,0),number.cex = .7)
  plot(hclust(dd,method = "ward.D"),main="Heirarchial Clustering")
  title(paste0(toupper(deparse(substitute(df)))," Stimulus Correlation"), line = -21, outer = TRUE)
  graphics.off()
  rm(c,p,p_adj,p_adj_mat,dd,col1)

}
###############################################################################################

##### Use the df objects created earlier #####
dd_correlation_plots(e2fc)
dd_correlation_plots(hiv)
dd_correlation_plots(seb)
dd_correlation_plots(bcg)
#############################################