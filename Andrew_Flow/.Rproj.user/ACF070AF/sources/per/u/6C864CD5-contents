
#################### 11th January 2019 #################
## Is there batch effect? These samples were run on different dates - is there evidence of different batches affecting 
## values ? 
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



############################ FUNCTION TO PLOT ALL ############################

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
##############################################################################



pdf("test.pdf")
par(mfcol=c(2,1))
colnames<-names(dta_v3)
dd_boxplot_function(dta_v3,PD1,Experiment_Date,"BCG")
dd_boxplot_function(dta_v3,PD1_fc,Experiment_Date,"BCG")
graphics.off()




#### HOW TO LOOP THROUGH ALL COLUMNS ####
fc_vars <- dta_v3 %>% select(ends_with("fc")) %>% colnames()
pdf("BCG_Stimulus_FC_BATCH_EFFECT_Plots.pdf")
fc_vars %>%
  syms() %>%
  map(function(var) dd_boxplot_function(dta_v3,!!var,Experiment_Date,"BCG"))
graphics.off()


###  WORKED !!!!!!

