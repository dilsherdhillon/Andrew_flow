library(tidyverse)


############################################### META DATA ##############################################
## IMPORT META DATA 
meta <- read_csv("~/Box/Dilsher/Andrew/Flow project/Data Files/20181209_Meta_Data_with_contactScore.csv")
meta<-meta%>%
  mutate(tb_status=`TB Status`,helminth=`_Composite_Helminth__PCR_and_serology_`,schisto=`_201805_schisto_ELISA`,
         ascaris=`_201805_ascaris_ELISA`,Sample=ID_stim)%>%
  select(Sample,tb_status,HIV,SEX,helminth,schisto,ascaris)
## Tracking data (experiment dates)
tracking_file <- read_csv("~/Box/Dilsher/Andrew/Flow project/Data Files/tracking_file.csv", 
                          +     col_types = cols(`CD14 MFI Exported` = col_skip(), 
                                                 +         `CD4 MFI Exported` = col_skip(), 
                                                 +         `CD56 MFI Exported` = col_skip(), 
                                                 +         `CD8 MFI Exported` = col_skip(), 
                                                 +         ID = col_skip(), `NKT Cells` = col_skip(), 
                                                 +         `NL MFI Exported` = col_skip(), NOTES = col_skip(), 
                                                 +         `Patient visit ID` = col_skip(), 
                                                 +         X12 = col_skip(), X13 = col_skip(), 
                                                 +         X14 = col_skip()))

### DATA ###
#### THIS IS CD4 DATA ######
dta_v2 <-
  read_csv("~/Box/Dilsher/Andrew/Flow project/Data Files/CD4_Master_MFI_v2.csv")
dta_v2 <- dta_v2 %>%
  mutate(Sample = X1) %>%
  select(-c(X1))
## DATA WRANGLING ## sapply(strsplit(dta$Sample,"_"), "[", 1) ## can use this
#instead of gsub # parse out the subject names, create a new column that
#indicates whether the sample is nil or not ##
dta_v2 <-
  dta_v2 %>%
  mutate(sample_group = gsub("_.*$", "", Sample)) %>%
  mutate(sample_group_2 = gsub(" .*$", "", sample_group)) %>%
  mutate(nil = ifelse(grepl("nil", Sample, ignore.case = TRUE), "Yes", "No")) %>%
  modify_at(c(1:9), as.numeric)  %>%
  group_by(sample_group) %>%
  arrange(desc(nil), .by_group = TRUE) %>%
  mutate_at(.funs = funs(fc = . / .[1]), .vars = (1:9)) %>%
  ungroup()
## Now, we create a new df, which only has the 
## 1. correlations and sample info
## 2. Remove the "nil" samples 
## 3. Remove the sample group columns as well 
## 4. Calculates fold change of stimulus using the nil 
fc_dta <- dta_v2 %>%
  select(-c(1:9, sample_group, sample_group_2)) %>%
  filter(nil == "No") %>%
  select(-c(nil)) %>%
  mutate(stimulus = case_when(
    grepl("BCG", Sample, ignore.case = TRUE) ~ "BCG",
    grepl("SEB", Sample, ignore.case = TRUE) ~ "SEB",
    grepl("HIV", Sample, ignore.case = TRUE) ~ "HIV",
    grepl("E,2f,C", Sample, ignore.case = TRUE) ~
      "E2FC",
    FALSE ~ "NA"
  ))



## select bcg from fc_data including the sample ID ##
bcg<-fc_dta%>%
  filter(stimulus=="BCG")%>%
  select(-c(stimulus))%>%
  {.}

## HEATMAP +HEIRARCHIAL CLUSTERING ## 
library("gplots")
library("heatmap.plus")
library("RColorBrewer")
## prep the matrix 
bcg_mat<-as.matrix(t(bcg))
colnames(bcg_mat)<-bcg_mat[1,]
names_bcg<-rownames(bcg_mat)
names_bcg<-names_bcg[c(-1,-11)]
rownames(bcg_mat)<-NULL
bcg_mat<-bcg_mat[c(-1,-11),]
rownames(bcg_mat)<-names_bcg
bcg_mat<-apply(bcg_mat, 2, as.numeric)
rownames(bcg_mat)<-names_bcg

## prep the annotations 

labels<-meta%>%
  right_join(.,bcg,by="Sample")%>%
  left_join(.,tracking_file)%>%
  distinct(Sample,.keep_all = TRUE)%>%
  mutate_at(vars(c(2:7)), funs(ifelse(.=="#N/A",NA,.)))
labels%>%
  #count(helminth)
  #count(schisto)
  #count(HIV)
  #count(ascaris)
  #count(SEX)
  #count(Experiment_Date)

colors<-labels%>%
  mutate(colors_helminth=case_when(helminth=="Negative" ~'#FFC0CB',
            helminth=="Positive" ~'#CC0000',
            is.na(helminth)~'#808080'),
         colors_schisto=case_when(schisto=="Negative" ~'#FFC0CB',
                                  schisto=="Positive" ~'#CC0000',
                                  is.na(schisto)~'#808080'),
         colors_ascaris=case_when(ascaris=="Negative" ~'#FFC0CB',
                                  ascaris=="Positive" ~'#CC0000',
                                  is.na(ascaris)~'#808080'),
         colors_hiv=case_when(HIV=="Negative" ~'#FFC0CB',
                              HIV=="Positive" ~'#CC0000',
                              is.na(HIV)~'#808080'),
         colors_gender=case_when(SEX=="Male" ~'#FFC0CB',
                                 SEX=="Female" ~'#CC0000',
                                 is.na(SEX)~'#808080'))%>%
  select(colors_helminth,colors_schisto,colors_ascaris,colors_hiv,colors_gender)
colors<-as.matrix(colors)

pdf(file="/Users/dilsherdhillon/Box/Dilsher/Andrew/Flow project/Andrew_Flow/Outputs/BCG_CD4_Heatmap_v2_10012018_averageMethod.pdf")
par(cex.main=0.8,mar=c(1,1,1,1),cex.axis=0.5)
heatmap.plus(bcg_mat, scale="row",col=blues9,cexRow=1,cexCol=0.2, margins = c(20,13), main="BCG Stimulus for CD4",
             ColSideColors=colors,hclust=function(x)hclust(x,method="complete"))
legend(0.009,0.98,legend=c("Negative","Positive","NA"),fill=c('#FFC0CB','#CC0000','#808080'),cex=0.5)
legend(0.90,0.98,legend=c("Male","Female","NA"),fill=c('#FFC0CB','#CC0000','#808080'),cex=0.5)
graphics.off()



## change clustering paramaters? 
pdf(file="/Users/dilsherdhillon/Box/Dilsher/Andrew/Flow project/Andrew_Flow/Outputs/BCG_CD4_Heatmap_v2_10012018_WardDMethod_DistanceCorrelation.pdf")
dist_samples<-dist(t(bcg_mat),method = "euclidean")
hc_samples<-hclust(dist_samples,method = "complete")

dist_genes<-dist(1-abs(cor(t(bcg_mat),use = "pairwise.complete.obs",method = "spearman"))/2,method = "euclidean")
hc_genes<-hclust(dist_genes,method = "ward.D")

heatmap.plus(bcg_mat, scale="row",col=blues9,cexRow=1,cexCol=0.2, margins = c(20,13), main="BCG Stimulus for CD4",
             ColSideColors=colors,Rowv = as.dendrogram(hc_genes),Colv = as.dendrogram(hc_samples))


## no scaling 
heatmap.plus(bcg_mat,col=blues9,cexRow=1,cexCol=0.2, margins = c(20,13), main="BCG Stimulus for CD4",
             ColSideColors=colors,Rowv = as.dendrogram(hc_genes),Colv = as.dendrogram(hc_samples))
graphics.off()





#################### 11th January 2019 #################
## Is there batch effect? These samples were run on different dates - is there evidence of different batches affecting 
## values ? 
bcg_date<-bcg%>%
  left_join(.,tracking_file)%>%
  modify_at(c(11),as.factor)

bcg_date%>%
  ggplot(.,aes(Experiment_Date,Perforin_fc))+geom_jitter(size=0.5)+geom_boxplot()+theme(axis.text.x = element_text(angle = 45, hjust = 1))

## This is fold change  - what about absolute counts ?
dta_v3<-dta_v2%>%
  left_join(.,tracking_file)%>%
  modify_at(c(23),as.factor)%>%
  mutate(stimulus = case_when(
    grepl("BCG", Sample, ignore.case = TRUE) ~ "BCG",
    grepl("SEB", Sample, ignore.case = TRUE) ~ "SEB",
    grepl("HIV", Sample, ignore.case = TRUE) ~ "HIV",
    grepl("E,2f,C", Sample, ignore.case = TRUE) ~
      "E2FC",
    FALSE ~ "NA"
  ))


#### Write a function where we can filter using stimulus and plot any gene ###
dd_boxplot_function<-function(df,var,cat_var,x){
  ## This function takes in a tidy df, continous variable (var), categorical variable(cat_var)
  ## cat_var should be a factor
  ## x is a character (the stimulus)
  var<-enquo(var)
  cat_var<-enquo(cat_var)
  df%>%
    filter(stimulus==x)%>%
    ggplot(.,aes(!!cat_var,!!var))+geom_jitter(size=0.5)+geom_boxplot()+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    labs(title=paste0(x))
}

pdf("test.pdf")
par(mfcol=c(2,1))
colnames<-names(dta_v3)
dd_boxplot_function(dta_v3,PD1,Experiment_Date,"BCG")
dd_boxplot_function(dta_v3,PD1_fc,Experiment_Date,"BCG")
graphics.off()




#### HOW TO LOOP THROUGH ALL COLUMNS ####
fc_vars <- dta_v3 %>% select(ends_with("fc")) %>% colnames()
bcg_fc_plots<-fc_vars %>%
  syms() %>%
  map(function(var) dd_boxplot_function(dta_v3,!!var,Experiment_Date,"BCG"))
###  WORKED !!!!!!

