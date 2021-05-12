

                                            
#------------------------------------------------------------------------------------------------
                                            
# EXAMPLE SCRIPT FOR THE ANALYSIS OF CLUSTERED SURVEY DATA
                                            
#------------------------------------------------------------------------------------------------
                                            
                                            

# Let's assume an anti-TB drug resistance survey (DRS) was to enroll 2181 new bacteriologically confirmed pulmonary TB patients (i.e. sample size = 2181).
# Stratified cluster sampling involving three strata and variable cluster size was chosen for the survey sampling.The required sample was proportionally allocated to
# each stratum based on the total stratum caseload.Stratum 1 was allocated 10 clusters and a total sample size of 650 new TB cases; stratum 2 was allocated 14 clusters 
# and 760 cases;finally stratum 3 was allocated 14 clusters and a sample size of 771 cases. All clusters were to enroll patients over a fix time period in each stratum.
 

#Before loading the libraries, you should install the required R packages if these are not yet installed.
#You may install packages using the RStudio menu (see "Tools" and then see "Install Packages") or by
#typing the following line: install.packages("nameOfRequiredPackage") 
#It is possible to install multiple packages at once with the following script line (this example
#installs 2 packages): install.packages(c("nameOfPackage1", "nameOfPackage2")


#load required libraries

library(data.table) 
library(dplyr)
library(tidyr)
library(ggpubr)
library(kableExtra)
library(tidyverse)
library(broom)
library(survey)
library(mice)
library(miceadds)
library(mitools)
                     
                                            

# open file; replace "." with the file location for "drsSampleData.csv":

#d1<-fread("/......./drsSampleData.csv")           


# view the dataset

names(d1) # view column names
View(d1) # view dataset

#check the size of the dataset

ncol(d1) # number of columns
nrow(d1) # number of observations (rows)

# check that patient DRS IDs are unique

anyDuplicated(d1$drsID) # there are no duplicate DRS IDs

# check for missing values in the dataset

na_count <-sapply(d1, function(y) sum(length(which(is.na(y)))))
data.frame(na_count)

# explore the numbers of observations within each categorical variable level in the dataset:

d1_summary<-select(d1, -c(drsID,clusterID,age))
d1_summary %>%
  map( table )

# there are 3 ntm in this dataset.Check whether pDST results are available for NTM

d1_ntm<- d1 %>% select(cultureResult, rifPdst, inhPdst, lfxPdst, bdqPdst, lzdPdst) %>%
                            filter(d1$cultureResult=="ntm")

(d1_ntm) # there are pDST results for rifPdst, inhPdst, lfxPdst, bdqPdst, lzdPdst. 


# we are only interested in prevalence of resistance among MTB. We will update d1 to dismiss
# pDST results from ntm by recoding the pDST results among the NTM as missing (i.e. "NA")

d1$rifPdst<-ifelse(d1$cultureResult=="ntm",NA,d1$rifPdst)
d1$inhPdst<-ifelse(d1$cultureResult=="ntm",NA,d1$inhPdst)
d1$lfxPdst<-ifelse(d1$cultureResult=="ntm",NA,d1$lfxPdst)
d1$bdqPdst<-ifelse(d1$cultureResult=="ntm",NA,d1$bdqPdst)
d1$lzdPdst<-ifelse(d1$cultureResult=="ntm",NA,d1$lzdPdst)

# Check that there are no pDST results for NTM

d2_ntm<- d1 %>% select(cultureResult, rifPdst, inhPdst, lfxPdst, bdqPdst, lzdPdst) %>%
  filter(d1$cultureResult=="ntm")

(d2_ntm) # there are no pDST results among NTM. 

# review again the number of missing values in the dataset

na_count <-sapply(d1, function(y) sum(length(which(is.na(y)))))
data.frame(na_count)


# create the resistant variables of interest, which will be used in the analyses:

# Begin with creating a new column ("rifXpertPdst") containing the final rifampicin resistance results, taking into
# consideration the Xpert MTB/RIF and phenotypic DST (pDST) results for rifampicin. In this example, a TB case is
#considered resistant to rifampicin, if a resistant result is obtained by at least one of the tests. The new variable
# is coded as a binary outcome, where "1" = rifampicin resistant; "0" rifampicin susceptible. In the script below, "&" = "and"; "|" = "or"


d1<-d1 %>%
  mutate(rifXpertPdst = case_when(
    (endsWith(xpertResult, "RI") & (rifPdst == "r"))|
      (endsWith(xpertResult, "RR") | (rifPdst == "r")) ~ 1,
    (endsWith(xpertResult, "RS") & (rifPdst == "s"))| 
      (endsWith(xpertResult, "RS") & is.na(rifPdst)) |
      (endsWith(xpertResult, "RI") & (rifPdst == "s")) ~ 0,
    endsWith(xpertResult, "RI") & is.na(rifPdst) ~ NA_real_,# NA_real_ is used to specify a numeric missing value
    endsWith(xpertResult, "Negative") ~ NA_real_)) # we would use NA_character_ if the new column was a factor variable

# we can then remove "rifPdst" from the dataset, as this variable will no longer be used:

d1$rifPdst<- NULL

# create a column for MDR-TB (i.e. resistance to rifampicin and isoniazid), that considers "rifXpertPdst" and "inhPdst":

d1<-d1 %>%
  mutate(mdr = case_when(
    rifXpertPdst == 1 & (inhPdst == "r")  ~ 1,
    (rifXpertPdst == 0 & (inhPdst == "s"))|
    (rifXpertPdst == 0 & (inhPdst == "r"))| 
    (rifXpertPdst == 1 & (inhPdst == "s")) ~ 0,
    TRUE~ NA_real_)) # MDR status is coded as missing for criteria other than above (i.e. when either rifXpertPdst or inhPdst are missing)


# create a column for Hr-TB (i.e. suscpetible to rifampicin and resistant to isoniazid), that considers "rifXpertPdst" and "inhPdst":

d1<-d1 %>%
  mutate(hr = case_when(
    rifXpertPdst == 0 & (inhPdst == "r")  ~ 1,
    (rifXpertPdst == 1 & (inhPdst == "r"))| 
    (rifXpertPdst == 0 & (inhPdst == "s"))|
    (rifXpertPdst == 1 & (inhPdst == "s"))~ 0,
    TRUE~ NA_real_)) 

# create a column for pre-XDR TB, which is defined as resistance to FQs among RR-TB

d1<-d1 %>%
  mutate(preXdr = case_when(
    rifXpertPdst == 1 & (lfxPdst == "r")  ~ 1,
    (rifXpertPdst == 1 & (lfxPdst == "s"))| 
    (rifXpertPdst == 0 & (lfxPdst == "s"))|
    (rifXpertPdst == 0 & (lfxPdst == "r"))~ 0,
    TRUE~ NA_real_)) 


# create a column for XDR-TB, which is defined as resistance to FQs plus at least one grup A drug (bdq and/or lzd) among RR-TB

d1<-d1 %>%
  mutate(xdr = case_when(
    preXdr == 1 & (bdqPdst == "r" | lzdPdst == "r")  ~ 1,
    (preXdr == 1 & (bdqPdst == "s" & lzdPdst == "s"))| 
      (preXdr == 0) ~ 0,
    TRUE~ NA_real_)) 


# logistic regression to estimate resistance to individual drugs, requires that drugs are coded as binary variables
# (0/1) as opposite to factors (r/s). We hence re-code the data for the remaining drugs accordingly:

d1$inhPdst<-ifelse(d1$inhPdst=="r",1,0)
d1$lfxPdst<-ifelse(d1$lfxPdst=="r",1,0)
d1$bdqPdst<-ifelse(d1$bdqPdst=="r",1,0)
d1$lzdPdst<-ifelse(d1$lzdPdst=="r",1,0)

# we group age into seven levels and drop the numeric (continuous) "age" variable


d1<-d1%>%
  mutate(ageGroup = case_when(d1$age < 15 ~ '0-14',
                          between(d1$age, 15, 24) ~ '15-24',
                          between(d1$age, 25, 34) ~ '25-34',
                          between(d1$age, 35, 44) ~ '35-44',
                          between(d1$age, 45, 54) ~ '45-54',
                          between(d1$age, 55, 64) ~ '55-64',
                          d1$age >= 65 ~ '65+'))


# view the number of cases by age group:

table(d1$ageGroup) 


## Note that newly created columns appear last in the dataset; move "age" column before "sex" column

d1 <- d1 %>% relocate(ageGroup, .before = sex)


#------------------------------------------------------------------------------------------------

# CALCULATE REVELANT COUNTS TO POPULATE A CONSORT CHART FOR THE SURVEY

#------------------------------------------------------------------------------------------------

# An example patient enrolment flowchart (see powerpoint slide) has been provided with this R code and the Stata code.
# The Stata code summarizes - step by step -  how to calculate all the numbers to populate each of the
# boxes of the consort chart. Here, we provide a brief example illustrating how calculating counts for a consort chart, may be replicated in R.
# We briefly show how to calculate the counts for one of the boxes in the consort chart. Counts for the
# remaining boxes, can easily be calculated by replacing variable names and modifying the filters in the example script
# below.


# In order to be able to count missing data in the consort chart (or patient enrollment chart), we first
# create a new dataset (herein referred to as "d2") where all the missing values are replaced with "unknown":

d2<-d1 %>% replace(is.na(.), "unknown")


# the following script calculates the total number of MTB cases by Xpert MTB/RIF, for which a rifampicin 
# resistance result is available (i.e. excluding rifampicin indeterminates (mtbRI)), and also disaggregates the 
# count by previous treatment history group:

d2 %>%filter(xpertResult != "mtbNegative" & xpertResult != "mtbRI")%>% # this filter excludes all cases where MTB was not confirmed by Xpert MTB/RIF, or where the rifampicin resistance result was indeterminate
  group_by(treatment) %>%
  summarise (n = n()) %>%  # calculate number of cases by history of previous anti-TB treatment
  bind_rows(summarise(.,
                      across(where(is.numeric), sum),
                      across(where(is.character), ~"Total"))) # calculate the total number of MTB cases by Xpert MTB/RIF, with a rifampicin resistance result, and print this as a "Total" row


# the above code shows there are a total of 2363 MTB cases with a rifampicin resistance result by MTB/RIF; of these,
# 2090 are new cases and 265 are previously treated cases. Previous treatment history is unknown for
# 8 MTB cases. 



#-------------------------------------------------------------------------------------------------------------------------------------------------

# CREATE A TABLE DESCRIBING THE CLINICAL AND DEMOGRAPHIC CHARACTERISTICS OF TB PATIENTS, DISSAGGREGATED BY PREVIOUS TREATMENT HISTORY

#-------------------------------------------------------------------------------------------------------------------------------------------------



# begin by creating a data frame (c1), that summarises the frequency and percentage of sex, age, HIV and location
# variable levels among new TB cases

c11<- d1 %>% filter(treatment == "new" & xpertResult != "mtbNegative" & !is.na(sex))%>%
  group_by(sex) %>%
  summarise (n = n()) %>%
  mutate("%" = round((n / sum(n))*100,1)) %>%
  rename(levels = sex) 

c12<- d1 %>% filter(treatment == "new" & xpertResult != "mtbNegative" & !is.na(ageGroup))%>%
  group_by(ageGroup) %>%
  summarise (n = n()) %>%
  mutate("%" = round((n / sum(n))*100,1)) %>%
  rename(levels = ageGroup) 

c13<- d1 %>% filter(treatment == "new" & xpertResult != "mtbNegative" & !is.na(HIV))%>%
  group_by(HIV) %>%
  summarise (n = n()) %>%
  mutate("%" = round((n / sum(n))*100,1)) %>%
  rename(levels = HIV) 


c14<- d1 %>% filter(treatment == "new" & xpertResult != "mtbNegative" & !is.na(location))%>%
  group_by(location) %>%
  summarise (n = n()) %>%
  mutate("%" = round((n / sum(n))*100,1)) %>% 
  rename(levels = location) 

c1<-rbind (c11, c12, c13, c14)

(c1)

# continue by creating a data frame (c2), that summarises the frequency and percentage of sex, age, HIV and location
# variable levels among previously treated TB cases

c21<- d1 %>% filter(treatment == "retreatment" & xpertResult != "mtbNegative" & !is.na(sex))%>%
  group_by(sex) %>%
  summarise (n = n()) %>%
  mutate("%" = round((n / sum(n))*100,1)) %>%
  rename(levels2 = sex) 

c22<- d1 %>% filter(treatment == "retreatment" & xpertResult != "mtbNegative" & !is.na(ageGroup))%>%
  group_by(ageGroup) %>%
  summarise (n = n()) %>%
  mutate("%" = round((n / sum(n))*100,1)) %>%
  rename(levels2 = ageGroup) 

c23<- d1 %>% filter(treatment == "retreatment" & xpertResult != "mtbNegative" & !is.na(HIV))%>%
  group_by(HIV) %>%
  summarise (n = n()) %>%
  mutate("%" = round((n / sum(n))*100,1)) %>%
  rename(levels2 = HIV) 


c24<- d1 %>% filter(treatment == "retreatment" & xpertResult != "mtbNegative" & !is.na(location))%>%
  group_by(location) %>%
  summarise (n = n()) %>%
  mutate("%" = round((n / sum(n))*100,1)) %>% 
  rename(levels2 = location) 

c2<-rbind (c21, c22, c23, c24)

(c2)


# finish by creating a data frame (c3), that summarises the frequency and percentage of sex, age, HIV and location
# variable levels among all TB cases combined

c31<- d1 %>% filter(xpertResult != "mtbNegative" & !is.na(sex))%>%
  group_by(sex) %>%
  summarise (n = n()) %>%
  mutate("%" = round((n / sum(n))*100,1)) %>%
  rename(levels3 = sex) 

c32<- d1 %>% filter(xpertResult != "mtbNegative" & !is.na(ageGroup))%>%
  group_by(ageGroup) %>%
  summarise (n = n()) %>%
  mutate("%" = round((n / sum(n))*100,1)) %>%
  rename(levels3 = ageGroup) 

c33<- d1 %>% filter(xpertResult != "mtbNegative" & !is.na(HIV))%>%
  group_by(HIV) %>%
  summarise (n = n()) %>%
  mutate("%" = round((n / sum(n))*100,1)) %>%
  rename(levels3 = HIV) 


c34<- d1 %>% filter(xpertResult != "mtbNegative" & !is.na(location))%>%
  group_by(location) %>%
  summarise (n = n()) %>%
  mutate("%" = round((n / sum(n))*100,1)) %>% 
  rename(levels3 = location) 

c3<-rbind(c31,c32,c33,c34)

(c3)

# combine the three datasets above and drop unnecessary columns

t1<-cbind(c1,c2,c3)

(t1)

t1$levels2<- NULL
t1$levels3<- NULL

(t1)

# you may want to save an editable copy of the dataset

#write.csv(t1, "/......./Table1.csv", row.names=F) 


# assemble a table using kableExtra features

  kbl(t1[1:14, 1:7], align="lrrrrrr",caption = "Clinical and demographic characteristics of TB cases", col.names = c("variables", "n", "%", "n", "%", "n", "%")) %>%
  kable_classic(full_width = F, html_font = "Cambria", font_size = 18) %>%
  pack_rows("sex", 1, 2) %>%
  pack_rows("age (years)", 3, 9) %>%
  pack_rows("HIV", 10, 11) %>%
  pack_rows("location", 12, 14) %>%
  add_header_above(c(" ", "new" = 2, "previously treated" = 2, "combined" = 2))


  
#------------------------------------------------------------------------------------------------
  
# PLOT AGE-SEX TB POPULATION PYRAMIDS
  
#------------------------------------------------------------------------------------------------
  

#Prepare dataset for age-sex population pyramid: all TB cases

dataPP1<-d1 %>% 
  filter(xpertResult != "mtbNegative" & !is.na(sex) & !is.na(ageGroup)) %>%
  group_by(sex,ageGroup) %>% 
  summarise(population= n())

view(dataPP1)

#look at the total count of males and females; the total count will be added in the legend of the plot

countMf1<-dataPP1 %>%
  group_by(sex) %>%
  summarise(total = sum(population))

(countMf1)

#look at the total number of TB cases; the total count will be added to the title of the plot

allTB<-sum(dataPP1$population)

(allTB)

#costume the plot

plotPP1<-ggplot(data = dataPP1, 
                 mapping = aes(x = ageGroup, fill = sex, 
                               y = ifelse(test = sex == "male", yes = -population, no = population))) +
  theme_bw()+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"),labels = c(paste("females (n=", gsub(" ","",paste(countMf1[1,2],")"))),paste("males (n=", gsub(" ","",paste(countMf1[2,2],")")))))+
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(dataPP1$population) * c(-1,1)) +
  labs(y = "population", x="age (years)") +
  ggtitle("TB population structure")+
  ggtitle(paste("TB cases (n=", gsub(" ","",paste(allTB,")"))))+
  theme(plot.title = element_text(hjust = 0.5))+
  coord_flip()

#plot the age-sex TB population structure

plotPP1


#### next, we will prepare the data to plot the age sex population structure stratified by previous treatment history and
#### by rifampicin resistance status:

##Population structure of new TB cases

dataPP2<-d1 %>% 
  filter(xpertResult != "mtbNegative" & treatment == "new" & !is.na(sex) & !is.na(ageGroup)) %>%
  group_by(sex,ageGroup) %>% 
  summarise(population= n())

view(dataPP2)

countMf2<-dataPP2 %>%
  group_by(sex) %>%
  summarise(total = sum(population))

(countMf2)

newTB<-sum(dataPP2$population)

(newTB)

plotPP2<-ggplot(data = dataPP2, 
                mapping = aes(x = ageGroup, fill = sex, 
                              y = ifelse(test = sex == "male", yes = -population, no = population))) +
  theme_bw()+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"),labels = c(paste("females (n=", gsub(" ","",paste(countMf2[1,2],")"))),paste("males (n=", gsub(" ","",paste(countMf2[2,2],")")))))+
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(dataPP2$population) * c(-1,1)) +
  labs(y = "population", x="age (years)") +
  ggtitle(paste("New (n=", gsub(" ","",paste(newTB,")"))))+
  theme(plot.title = element_text(hjust = 0.5))+
  coord_flip()

##Population structure of previously treated TB cases


dataPP3<-d1 %>% 
  filter(xpertResult != "mtbNegative" & treatment == "retreatment" & !is.na(sex) & !is.na(ageGroup)) %>%
  group_by(sex,ageGroup) %>% 
  summarise(population= n())

view(dataPP3)

countMf3<-dataPP3 %>%
  group_by(sex) %>%
  summarise(total = sum(population))

(countMf3)

retreatmentTB<-sum(dataPP3$population)

(retreatmentTB)

plotPP3<-ggplot(data = dataPP3, 
                mapping = aes(x = ageGroup, fill = sex, 
                              y = ifelse(test = sex == "male", yes = -population, no = population))) +
  theme_bw()+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"),labels = c(paste("females (n=", gsub(" ","",paste(countMf3[1,2],")"))),paste("males (n=", gsub(" ","",paste(countMf3[2,2],")")))))+
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(dataPP3$population) * c(-1,1)) +
  labs(y = "population", x="age (years)") +
  ggtitle(paste("Previously treated (n=", gsub(" ","",paste(retreatmentTB,")"))))+
  theme(plot.title = element_text(hjust = 0.5))+
  coord_flip()


##Population structure of RS-TB cases


dataPP4<-d1 %>% 
  filter(rifXpertPdst == "0" & !is.na(sex) & !is.na(ageGroup)) %>%
  group_by(sex,ageGroup) %>% 
  summarise(population= n())

view(dataPP4)

countMf4<-dataPP4 %>%
  group_by(sex) %>%
  summarise(total = sum(population))

(countMf4)

rsTB<-sum(dataPP4$population)

(rsTB)

plotPP4<-ggplot(data = dataPP4, 
                mapping = aes(x = ageGroup, fill = sex, 
                              y = ifelse(test = sex == "male", yes = -population, no = population))) +
  theme_bw()+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"),labels = c(paste("females (n=", gsub(" ","",paste(countMf4[1,2],")"))),paste("males (n=", gsub(" ","",paste(countMf4[2,2],")")))))+
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(dataPP4$population) * c(-1,1)) +
  labs(y = "population", x="age (years)") +
  ggtitle(paste("Rifampicin susceptible (n=", gsub(" ","",paste(rsTB,")"))))+
  theme(plot.title = element_text(hjust = 0.5))+
  coord_flip()


##Population structure of RR-TB cases


dataPP5<-d1 %>% 
  filter(rifXpertPdst == "1" & !is.na(sex) & !is.na(ageGroup)) %>%
  group_by(sex,ageGroup) %>% 
  summarise(population= n())

view(dataPP5)

countMf5<-dataPP5 %>%
  group_by(sex) %>%
  summarise(total = sum(population))

(countMf5)

rrTB<-sum(dataPP5$population)

(rrTB)

plotPP5<-ggplot(data = dataPP5, 
                mapping = aes(x = ageGroup, fill = sex, 
                              y = ifelse(test = sex == "male", yes = -population, no = population))) +
  theme_bw()+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"),labels = c(paste("females (n=", gsub(" ","",paste(countMf5[1,2],")"))),paste("males (n=", gsub(" ","",paste(countMf5[2,2],")")))))+
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(dataPP5$population) * c(-1,1)) +
  labs(y = "population", x="age (years)") +
  ggtitle(paste("Rifampicin resistant (n=", gsub(" ","",paste(rrTB,")"))))+
  theme(plot.title = element_text(hjust = 0.5))+
  coord_flip()


###########Create a multi-pannel plot for the population structure of new, previously treated, RR-TB and RS-TB cases


# first run the multiplot function:

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# now call the plots; the order of the plots may be changed:


multiplot(plotPP2,plotPP4,plotPP3,plotPP5, cols=2)



#------------------------------------------------------------------------------------------------

# RESISTANCE TO RIFAMPICIN AMONG NEW AND PREVIOUSLY TREATED PULMONARY TB CASES (imputation)

#------------------------------------------------------------------------------------------------


# create a column that contains the weights for each stratum. For this, we first need to count the number of enrolled
# new cases that were confirmed as MTB by Xpert MTB/RIF at the central laboratory within each stratum. We 
# begin by dropping all cases where MTB was not confirmed by Xpert MTB/RIF, before computing the count.
# We create a new column (cs) to store this count. Then, we add a column (ess) that contains the target sample size
# in each stratum. Finally, we create a column of weights (w) that corresponds to the division of
# ess by cs (i.e. ess/cs). The last line of the code drops columns "ess" and "cs". Note that whilst the weights
# are calculated based on new cases only, these same weights are applied to previously treated cases as well. Also
# note that the approach to compute weights may differ in different surveys, and should be considered on a case 
# by case basis.



drif <-d1 %>% filter(xpertResult != "mtbNegative")%>%
  group_by(strata) %>% 
  mutate(cs= sum(treatment=="new", na.rm=TRUE))%>% 
  mutate(ess = case_when(
    strata == 1 ~ 650,
    strata == 2  ~ 760,
    strata == 3  ~ 771))  %>% 
  mutate(w = (ess/cs)) %>% 
  select(-c(ess, cs)) 



#provide function to compute prevalence estimates based on model coefficients and 95% confidence intervals:


invlogit <- function(x) 1/(1+exp(-x)) 



# missing data

(sapply(drif, function(x)sum(is.na(x))))# missing data for treatment history, age and sex. Seven missing values for
#resistance to rifampicin (see variable "rifXpertPdst")


# run a model which does not consider the cluster-based study design, and does not take weights into account

fdrif0 <- glm(rifXpertPdst ~ treatment - 1, family=binomial, data=drif) # no intercept

fdrif0a<-tibble(round((invlogit(coef(fdrif0)))*100,1))             # prevalence by treatment group

fdrif0b<-round((invlogit(confint(fdrif0)))*100,1)        # 95% confidence intervals

fdrif0ab<-cbind(fdrif0a,fdrif0b)%>% 
  rename(prevalence = "round((invlogit(coef(fdrif0))) * 100, 1)") 

(fdrif0ab)


# run a model with sampling weights. Use quasibinomial family to avoid the warning about non-integer successes

fdrif1 <- glm(rifXpertPdst ~  treatment - 1, family=quasibinomial, data=drif, weights=w)


fdrif1a<-tibble(round((invlogit(coef(fdrif1)))*100,1))             # prevalence by treatment group

fdrif1b<-round((invlogit(confint(fdrif1)))*100,1)        # 95% confidence intervals

fdrif1ab<-cbind(fdrif1a,fdrif1b)%>% 
  rename(prevalence = "round((invlogit(coef(fdrif1))) * 100, 1)") 

(fdrif1ab)


# run model adjusted for clustering, with sampling weights (complete data only)

dsdrif <- svydesign(ids=~clusterID, weights=~w, data=drif) # this line provides information on the study design

fdrif2 <- svyglm(rifXpertPdst ~treatment - 1, family=quasibinomial, design=dsdrif)

fdrif2a<-tibble(round((invlogit(coef(fdrif2)))*100,1))             # prevalence by treatment group

fdrif2b<-round((invlogit(confint(fdrif2)))*100,1)        # 95% confidence intervals

fdrif2ab<-cbind(fdrif2a,fdrif2b)%>% 
  rename(prevalence = "round((invlogit(coef(fdrif2))) * 100, 1)") 

(fdrif2ab)

# multiple imputation of missing data, simplify var names, factorize character vars

setDT(drif)

drif2 <- drif[, .(clusterID,
                  
                  treatment=factor(treatment),
                  
                  age=age,
                  
                  sex=factor(sex),
                  
                  rifXpertPdst=factor(rifXpertPdst),
                  
                  w)]

view(drif2)

imp <- mice(drif2, maxit=0)

pred <- imp$pred

meth <- imp$meth

meth['age'] <- 'pmm' #imputation method for numeric data

meth['sex'] <- 'logreg' #imputation method for factor variables with two levels. If the variable has > 2 levels, "polyreg" or "polyr" would be used instead for unordered or ordered variables respectively)

meth['treatment'] <- 'logreg'

meth['rifXpertPdst'] <- 'logreg'

pred[,'clusterID'] <- 0  # do not use cluster as a predictor (fixed effects)



system.time(idrif <- mice(drif2, method=meth, pred=pred, m=20, seed=123)) # takes <2 min 

imp_ <- mids2datlist(idrif)

##the imputation produces 20 imputed datasets. the first rows of the 20 datasets can be viewed with the following lines:

ds20<-datlist_create(datasets = imp_)
print(ds20)
view(ds20)


# run model adjusted for clustering, with sampling weights (imputed data)

idsdrif <- svydesign(ids=~clusterID, weights=~w, data=imputationList(imp_)) # design for the model that uses imputed data

fdrif3 <- with(idsdrif, svyglm(rifXpertPdst ~ treatment - 1, design=.design, family=quasibinomial, weights=w))

out <- summary(MIcombine(fdrif3))

fdrif3ab<-(res <- round((invlogit(out[,c(1,3,4)]))*100,1)) # final results (prevalence by treatment group and 95% confidence intervals)

names(fdrif3ab)[1] <- "prevalence" 
names(fdrif3ab)[2] <- "2.5 %"
names(fdrif3ab)[3] <- "97.5 %" 


(fdrif3ab)


###########################################################
# summarize results on resistance to rifampicin in a table ############################################## 
###########################################################

#create a dataset (r1) with all the rifampicin resistance estimates by treatment group

r1<-rbind(fdrif0ab,fdrif1ab,fdrif2ab,fdrif3ab)
r1<-rownames_to_column(r1, var = "group")

(r1)

# you may want to save an editable copy of the dataset

#write.csv(r1, "/......./Table2.csv", row.names=F) 


r2<-r1 %>% 
  mutate(group = str_replace(group, ".*treatmentnew.*", "new"),
         group = str_replace(group, ".*treatmentretreatment.*", "previously treated"))




kbl(r2[1:8, 1:4], caption = "Resistance to rifampicin", col.names = c("model", "prevalence", "2.5%", "97.5%")) %>%
  kable_classic(full_width = F, html_font = "Cambria", font_size = 16) %>%
  pack_rows("No sampling adjustment; no weights", 1, 2) %>%
  pack_rows("No sampling adjustment; with weights", 3, 4) %>%
  pack_rows("Adjusted for clustering; with weights; complete data only", 5, 6) %>%
  pack_rows("Adjusted for clustering; with weights; imputed data", 7, 8) %>%
  add_header_above(c(" ", " ", "95% confidence interval" = 2)) %>%
  row_spec(7:8, bold = T, color = "white", background = "#D7261E")

#in the final table, we highlight the final rifampicin resistance estimates in red (weighted logistic regression
#adjusted for clustering, using imputed data)


#--------------------------------------------------------------------------------------------------------------------

# RESISTANCE TO OTHER DRUGS AND DRUG COMBINATIONS AMONG NEW AND PREVIOUSLY TREATED PULMONARY TB CASES (no imputation)

#---------------------------------------------------------------------------------------------------------------------

#the function to compute prevalence estimates based on model coefficients and 95% confidence intervals
#has been provided earlier:

#invlogit <- function(x) 1/(1+exp(-x)) 

###################################################################
# DATASET AND DESIGN FOR ANALYSIS OF INH, LFX, MDR, HR, PREXDR, XDR ############################################## 
###################################################################

# for the analysis of these drugs and drug combinations we will use dataset "doth" and design "dsg"

# filter the records with drug susceptibility test results for drugs and drug combinations other than "rifXpertPdst".
# DST is available from all samples that were confirmed as MTB by Xpert MTB/RIF and subsequently by culture:

doth<- drif %>% filter(cultureResult == "mtbc" & !is.na(treatment))

# estimate the percentage resistance using logistic regression adjusted for clustering
# weights are applied to all drugs and drug combinations that were tested among all mtb. 
# this includes inh, lfx, mdr, hr, preXdr, xdr:


dsg <- svydesign(ids=~clusterID, weights=~w, data=doth) # design for analysis of resistance to inh, lfx, mdr, hr, preXdr,xdr

###############################################
# DATASET AND DESIGN FOR ANALYSIS OF BDQ, LZD  ############################################## 
###############################################

# for the analysis of these drugs we will use dataset "dbqlz" and design "dsgbqlz"

# These drugs were not tested among all MTB. Instead, bdq and lzd were only tested among MTB cultures
# that were resistant to rifampicin by either Xpert MTB/RIF or pDST (n=128). In consequence, weights are
# herein not applicable. The design will consider the cluster-based unweighted data. In addition, the  dataset
# is restricted to the relevant subset data:


dbqlz<- drif %>% filter(cultureResult == "mtbc" & rifXpertPdst == 1 & !is.na(treatment))

dsgbqlz<- svydesign(ids=~clusterID, data=dbqlz) # design for analysis of bdq and lzd


###########################################################################################
# fit models to estimate resistance prevalence disaggregated by previous treatment history  ############################################## 
###########################################################################################



#inhPdst

finh1 <- svyglm(inhPdst ~ treatment - 1, family=quasibinomial, design=dsg)

finh1a<-tibble(round((invlogit(coef(finh1)))*100,1))             # prevalence by treatment group

finh1b<-round((invlogit(confint(finh1)))*100,1)        # 95% confidence intervals

finh1ab<-cbind(finh1a,finh1b)%>% 
  rename(prevalence = "round((invlogit(coef(finh1))) * 100, 1)") 

(finh1ab)

#lfxPdst

flfx1 <- svyglm(lfxPdst ~ treatment - 1, family=quasibinomial, design=dsg)

flfx1a<-tibble(round((invlogit(coef(flfx1)))*100,1))             # prevalence by treatment group

flfx1b<-round((invlogit(confint(flfx1)))*100,1)        # 95% confidence intervals

flfx1ab<-cbind(flfx1a,flfx1b)%>% 
  rename(prevalence = "round((invlogit(coef(flfx1))) * 100, 1)") 

(flfx1ab)

#bdqPdst

fbdq1 <- svyglm(bdqPdst ~ treatment - 1, family=binomial, design=dsgbqlz) # note this is calling unweighted design = dsgbqlz

fbdq1a<-tibble(round((invlogit(coef(fbdq1)))*100,1))   # prevalence by treatment group

fbdq1b<-round((invlogit(confint(fbdq1)))*100,1)        # 95% confidence intervals

fbdq1ab<-cbind(fbdq1a,fbdq1b)%>% 
  rename(prevalence = "round((invlogit(coef(fbdq1))) * 100, 1)") 

(fbdq1ab)


#lzdPdst

flzd1 <- svyglm(lzdPdst ~ treatment - 1, family=binomial, design=dsgbqlz) # note this is calling unweighted design = dsgbqlz

flzd1a<-tibble(round((invlogit(coef(flzd1)))*100,1))   # prevalence by treatment group

flzd1b<-round((invlogit(confint(flzd1)))*100,1)        # 95% confidence intervals

flzd1ab<-cbind(flzd1a,flzd1b)%>% 
  rename(prevalence = "round((invlogit(coef(flzd1))) * 100, 1)") 

(flzd1ab)

#mdr

fmdr1 <- svyglm(mdr ~ treatment - 1, family=quasibinomial, design=dsg)

fmdr1a<-tibble(round((invlogit(coef(fmdr1)))*100,1))             # prevalence by treatment group

fmdr1b<-round((invlogit(confint(fmdr1)))*100,1)        # 95% confidence intervals

fmdr1ab<-cbind(fmdr1a,fmdr1b)%>% 
  rename(prevalence = "round((invlogit(coef(fmdr1))) * 100, 1)") 

(fmdr1ab)

#hr

fhr1 <- svyglm(hr ~ treatment - 1, family=quasibinomial, design=dsg)

fhr1a<-tibble(round((invlogit(coef(fhr1)))*100,1))             # prevalence by treatment group

fhr1b<-round((invlogit(confint(fhr1)))*100,1)        # 95% confidence intervals

fhr1ab<-cbind(fhr1a,fhr1b)%>% 
  rename(prevalence = "round((invlogit(coef(fhr1))) * 100, 1)") 

(fhr1ab)


#preXdr

fpreXdr1 <- svyglm(preXdr ~ treatment - 1, family=quasibinomial, design=dsg)

fpreXdr1a<-tibble(round((invlogit(coef(fpreXdr1)))*100,1))             # prevalence by treatment group

fpreXdr1b<-round((invlogit(confint(fpreXdr1)))*100,1)        # 95% confidence intervals

fpreXdr1ab<-cbind(fpreXdr1a,fpreXdr1b)%>% 
  rename(prevalence = "round((invlogit(coef(fpreXdr1))) * 100, 1)") 

(fpreXdr1ab)


#xdr

fxdr1 <- svyglm(xdr ~ treatment - 1, family=quasibinomial, design=dsg)

fxdr1a<-tibble(round((invlogit(coef(fxdr1)))*100,1))             # prevalence by treatment group

fxdr1b<-round((invlogit(confint(fxdr1)))*100,1)        # 95% confidence intervals

fxdr1ab<-cbind(fxdr1a,fxdr1b) %>% 
  rename(prevalence = "round((invlogit(coef(fxdr1))) * 100, 1)") 


(fxdr1ab)


###########################################################
# summarise results on resistance to other drugs in a table ############################################## 
###########################################################

#create a dataset (f1) with all the resistance estimates by treatment group

f1<-rbind(finh1ab, flfx1ab, fbdq1ab, flzd1ab, fmdr1ab, fhr1ab, fpreXdr1ab, fxdr1ab)


f1<-rownames_to_column(f1, var = "group")


(f1)

# you may want to save an editable copy of the dataset

#write.csv(f1, "/......./Table3.csv", row.names=F) 



f1<-f1 %>% 
  mutate(group=recode(group,treatmentnew="treatmentnew0",treatmentretreatment="treatmentretreatment0"))


f1<-f1 %>%
  arrange(parse_number(group))



f1<-f1 %>% 
  mutate(group=recode(group,treatmentnew="treatmentnew0",treatmentretreatment="treatmentretreatment0"))



f2<-f1 %>% 
  mutate(group = str_replace(group, ".*treatmentnew.*", "new"),
         group = str_replace(group, ".*treatmentretreatment.*", "previously treated"))



kbl(f2[1:16, 1:4], caption = "Resistance to other drugs and drug combinations", col.names = c("drug or drug combination", "prevalence", "2.5%", "97.5%")) %>%
  kable_classic(full_width = F, html_font = "Cambria", font_size = 16) %>%
  pack_rows("inh", 1, 2) %>%
  pack_rows("lfx", 3, 4) %>%
  pack_rows("bdq", 5, 6) %>%
  pack_rows("lzd", 7, 8) %>%
  pack_rows("mdr", 9, 10) %>%
  pack_rows("hr", 11, 12) %>%
  pack_rows("pre-xdr", 13, 14) %>%
  pack_rows("xdr", 15, 16) %>%
  add_header_above(c(" ", " ", "95% confidence interval" = 2))


#-------------------------------------------------------------------------------------------------------------

# RESISTANCE TO BDQ, LZD AND PRE-XDR- AND XDR-TB AMONG RIFAMPICIN RESISTANT PULMONARY TB CASES (no imputation)

#-------------------------------------------------------------------------------------------------------------

#the function to compute prevalence estimates based on model coefficients and 95% confidence intervals
#has been provided earlier:

#invlogit <- function(x) 1/(1+exp(-x)) 


#############################################################################
# DATASET AND DESIGN FOR ANALYSIS OF BDQ, LZD, PRE-XDR AND XDR AMONG RR-TB  ############################################## 
#############################################################################

# for the analysis of these drugs and drug combinations we will use dataset "dbrrtb" and design "dsgrrtb"

# Weights are herein not applicable because we are looking at the prevalence of
# selected drugs and drug combinations among RR-TB. The design will consider the cluster-based unweighted data. 
# In addition, the  dataset is restricted to the relevant subset data:


dbrrtb<- drif %>% filter(cultureResult == "mtbc" & rifXpertPdst == 1)

dsgrrtb<- svydesign(ids=~clusterID, data=dbrrtb) 


############################################################
# fit models to estimate resistance prevalence among RR-TB  ############################################## 
############################################################

#bdqPdst

fbdq2 <- svyglm(bdqPdst ~ 1, family=binomial, design=dsgrrtb) 

fbdq2a<-tibble(round((invlogit(coef(fbdq2)))*100,1))             # prevalence among all RR-TB cases

fbdq2b<-round((invlogit(confint(fbdq2)))*100,1)        # 95% confidence intervals

fbdq2ab<-cbind(fbdq2a,fbdq2b)%>% 
  rename(prevalence = "round((invlogit(coef(fbdq2))) * 100, 1)")%>% 
  mutate(varlevel = "bdq")%>% 
  relocate(varlevel) #move column to first position 

(fbdq2ab)

#lzdPdst

flzd2 <- svyglm(lzdPdst ~ 1, family=binomial, design=dsgrrtb) 

flzd2a<-tibble(round((invlogit(coef(flzd2)))*100,1))             # prevalence among all RR-TB cases

flzd2b<-round((invlogit(confint(flzd2)))*100,1)        # 95% confidence intervals

flzd2ab<-cbind(flzd2a,flzd2b)%>% 
  rename(prevalence = "round((invlogit(coef(flzd2))) * 100, 1)")%>% 
  mutate(varlevel = "lzd")%>% 
  relocate(varlevel) #move column to first position 

(flzd2ab)


#preXdr

fpreXdr2 <- svyglm(preXdr ~ 1, family=binomial, design=dsgrrtb)

fpreXdr2a<-tibble(round((invlogit(coef(fpreXdr2)))*100,1))             # prevalence among all RR-TB cases

fpreXdr2b<-round((invlogit(confint(fpreXdr2)))*100,1)        # 95% confidence intervals

fpreXdr2ab<-cbind(fpreXdr2a,fpreXdr2b)%>% 
  rename(prevalence = "round((invlogit(coef(fpreXdr2))) * 100, 1)") %>% 
  mutate(varlevel = "pre-xdr")%>% 
  relocate(varlevel) #move column to first position


(fpreXdr2ab)

#xdr

fxdr2 <- svyglm(xdr ~ 1, family=binomial, design=dsgrrtb)

fxdr2a<-tibble(round((invlogit(coef(fxdr2)))*100,1))             # prevalence among all RR-TB cases

fxdr2b<-round((invlogit(confint(fxdr2)))*100,1)        # 95% confidence intervals

fxdr2ab<-cbind(fxdr2a,fxdr2b) %>% 
  rename(prevalence = "round((invlogit(coef(fxdr2))) * 100, 1)") %>% 
  mutate(varlevel = "xdr")%>% 
  relocate(varlevel) #move column to first position



(fxdr2ab)


###########################################################
# summarise results on resistance among RR-TB in a table ############################################## 
###########################################################

#create a dataset (f2) with all the resistance estimates 

f2<-rbind(fbdq2ab, flzd2ab, fpreXdr2ab, fxdr2ab)
row.names(f2) <- NULL

(f2)

# you may want to save an editable copy of the dataset

#write.csv(f2, "/......./Table4.csv", row.names=F) 



kbl(f2[1:4, 1:4], align = "lccc", caption = "Resistance to selected drugs and combinations among RR-TB", col.names = c("drug or drug combination", "prevalence", "2.5%", "97.5%")) %>%
  kable_classic(full_width = F, html_font = "Cambria", font_size = 16) %>%
  add_header_above(c(" ", " ", "95% confidence interval" = 2))


#--------------------------------------------------------------------------------------

# ANALYSIS OF POTENTIAL PREDICTORS OF RIFAMPICIN RESISTANT TB [RR-TB] (no imputation)
# UNIVARIATE LOGISTIC REGRESSION ANALYSES ADJUSTED BY PREVIOUS TREATMENT HISTORY

#--------------------------------------------------------------------------------------


#The following univariate regression analyses are based on complete data only (i.e. no imputed data)
#and account for the cluster design of the survey plus the weights.The analyses for variables other than
#previous treatment history are adjusted by the history of anti-TB treatment. You should explore the results
#of unadjusted univariate analyses (for an example, see univariate analysis for treatment history below; treatment
#history would be replaced by the explanatory variable of interest). 
#Unadjusted analyses are not presented here.                                   



## treatment history

dt1<-filter(drif, !is.na(drif$treatment) & !is.na(drif$rifXpertPdst)) # remove rows with missing values; models 1 & 2 are then based on the same number of observations
nrow(dt1)
dsdt1 <- svydesign(ids=~clusterID, weights=~w, data=dt1) #specify design
model1 <- svyglm(rifXpertPdst ~treatment, family=quasibinomial, design=dsdt1) # to see output type summary (model1)
model2 <- svyglm(rifXpertPdst ~1, family=quasibinomial, design=dsdt1)
round(exp(coef(model1)),1) # odds ratio for treatment history
round(exp(confint(model1)),1) # 95% confidence intervals for the odds ratio
anova(model1, model2) # we compute the p-value for treatment history, by comparing a model with and without the explanatory variable of interest 




## sex (adjusted by treatment history)

dt2<-filter(drif, !is.na(drif$treatment) & !is.na(drif$rifXpertPdst) & !is.na(drif$sex)) # remove rows with missing values; models 3 & 4 are then based on the same number of observations
nrow(dt2)
dsdt2 <- svydesign(ids=~clusterID, weights=~w, data=dt2) #specify design
model3 <- svyglm(rifXpertPdst ~treatment + sex, family=quasibinomial, design=dsdt2)# to see output type summary (model3)
model4 <- svyglm(rifXpertPdst ~treatment, family=quasibinomial, design=dsdt2)
round(exp(coef(model3)),1) # odds ratio for sex
round(exp(confint(model3)),1) # 95% confidence intervals for the odds ratio
anova(model3,model4) # we compute the p-value for sex, by comparing a model with and without the explanatory variable of interest


## ageGroup (adjusted by treatment history)

dt3<-filter(drif, !is.na(drif$treatment) & !is.na(drif$rifXpertPdst) & !is.na(drif$ageGroup)) # remove rows with missing values; models 5 & 6 are then based on the same number of observations
nrow(dt3)
dsdt3 <- svydesign(ids=~clusterID, weights=~w, data=dt3) #specify design
model5 <- svyglm(rifXpertPdst ~treatment + ageGroup, family=quasibinomial, design=dsdt3)# to see output type summary (model5)
model6 <- svyglm(rifXpertPdst ~treatment, family=quasibinomial, design=dsdt3)
round(exp(coef(model5)),1) # odds ratio for ageGroup
round(exp(confint(model5)),1) # 95% confidence intervals for the odds ratio
anova(model5,model6) # we compute the p-value for ageGroup, by comparing a model with and without the explanatory variable of interest 


## HIV (adjusted by treatment history)

dt4<-filter(drif, !is.na(drif$treatment) & !is.na(drif$rifXpertPdst) & !is.na(drif$HIV)) # remove rows with missing values; models 7 & 8 are then based on the same number of observations
nrow(dt4)
dsdt4 <- svydesign(ids=~clusterID, weights=~w, data=dt4) #specify design
model7 <- svyglm(rifXpertPdst ~treatment + HIV, family=quasibinomial, design=dsdt4)# to see output type summary (model7)
model8 <- svyglm(rifXpertPdst ~treatment, family=quasibinomial, design=dsdt4)
round(exp(coef(model7)),1) # odds ratio for HIV status
round(exp(confint(model7)),1) # 95% confidence intervals for the odds ratio
anova(model7,model8) # we compute the p-value for HIV status, by comparing a model with and without the explanatory variable of interest 




## location (adjusted by treatment history)

dt5<-filter(drif, !is.na(drif$treatment) & !is.na(drif$rifXpertPdst) & !is.na(drif$location)) # remove rows with missing values; models 9 & 10 are then based on the same number of observations
nrow(dt5)
dsdt5 <- svydesign(ids=~clusterID, weights=~w, data=dt5) #specify design
model9 <- svyglm(rifXpertPdst ~treatment + location, family=quasibinomial, design=dsdt5)# to see output type summary (model9)
model10 <- svyglm(rifXpertPdst ~treatment, family=quasibinomial, design=dsdt5)
round(exp(coef(model9)),1) # odds ratio for location
round(exp(confint(model9)),1) # 95% confidence intervals for the odds ratio
anova(model9,model10) # we compute the p-value for location, by comparing a model with and without the explanatory variable of interest 




######################CREATE TABLE FOR ADJUSTED UNIVARIATE ANALYSES###################################


###TREATMENT HISTORY (MODEL1)

fm1a<-as.data.frame(tibble(round(exp(coef(model1)),1))[2,1])             # Odds ratio for previous treatment
fm1b<-as.data.frame(t(round(exp(confint(model1)),1) [2,]))      # 95% confidence intervals of the odds ratio
fm1c<-as.data.frame(round(anova(model1,model2)[[6]][[1]],3))#  we extract p-value for treatment history; to know what column and row of data is needed we can run str(anova(model1,model2))


totm1<-dt1 %>%group_by(treatment) %>%
  summarise (n = n()) %>% rename_at( 2, ~"total observations")%>% rename_at( 1, ~"varlevel")
rrm1<- dt1 %>%
  filter(rifXpertPdst == "1")%>%
  group_by(treatment) %>%
  summarise (n = n()) %>% rename_at( 2, ~"RR-TB")%>% rename_at( 1, ~"varlevel")

countm1<-left_join(totm1,rrm1, by="varlevel")


fm1ab<-cbind(fm1a,fm1b,fm1c)%>% 
  rename_at( 1, ~"odds ratio")%>%
  rename_at( 4, ~"pvalue" ) %>% mutate(varlevel = "retreatment")%>% 
  relocate(varlevel) #move column to first position
(fm1ab)

tm1<-left_join(countm1,fm1ab,by="varlevel")
(tm1)

###SEX (MODEL3)


fm3a<-as.data.frame(tibble(round(exp(coef(model3)),1))[3,1])             # Odds ratio for sex
fm3b<-as.data.frame(t(round(exp(confint(model3)),1) [3,]))      # 95% confidence intervals of the odds ratio
fm3c<-as.data.frame(round(anova(model3,model4)[[6]][[1]],3))#  we extract p-value for sex; to know what column and row of data is needed we can run str(anova(model3,model4))


totm3<-dt2 %>%group_by(sex) %>%
  summarise (n = n()) %>% rename_at( 2, ~"total observations")%>% rename_at( 1, ~"varlevel")
rrm3<- dt2 %>%
  filter(rifXpertPdst == "1")%>%
  group_by(sex) %>%
  summarise (n = n()) %>% rename_at( 2, ~"RR-TB")%>% rename_at( 1, ~"varlevel")

countm3<-left_join(totm3,rrm3, by="varlevel")


fm3ab<-cbind(fm3a,fm3b,fm3c)%>% 
  rename_at( 1, ~"odds ratio")%>%
  rename_at( 4, ~"pvalue" ) %>% mutate(varlevel = "male")%>% 
  relocate(varlevel) #move column to first position
(fm3ab)

tm3<-left_join(countm3,fm3ab,by="varlevel")
(tm3)


###AGE GROUP (MODEL5)


fm5a<-as.data.frame(tibble(round(exp(coef(model5)),1))[3:8,1])             # Odds ratio for age group
fm5b<-as.data.frame(round(exp(confint(model5)),1) [3:8,])     # 95% confidence intervals of the odds ratio
fm5c<-as.data.frame(round(anova(model5,model6)[[6]][[1]],3))#  we extract p-value for sex; to know what column and row of data is needed we can run str(anova(model3,model4))


totm5<-dt3 %>%group_by(ageGroup) %>%
  summarise (n = n()) %>% rename_at( 2, ~"total observations")%>% rename_at( 1, ~"varlevel")
rrm5<- dt3 %>%
  filter(rifXpertPdst == "1")%>%
  group_by(ageGroup) %>%
  summarise (n = n()) %>% rename_at( 2, ~"RR-TB")%>% rename_at( 1, ~"varlevel")

countm5<-left_join(totm5,rrm5, by="varlevel")


fm5ab<-cbind(fm5a,fm5b,fm5c)%>% 
  rename_at( 1, ~"odds ratio")%>%
  rename_at( 4, ~"pvalue" ) %>% mutate(varlevel = c("15-24", "25-34", "35-44","45-54","55-64","65+"))%>% 
  relocate(varlevel) #move column to first position
(fm5ab)

tm5<-left_join(countm5,fm5ab,by="varlevel")
(tm5)



###HIV (MODEL7)

fm7a<-as.data.frame(tibble(round(exp(coef(model7)),1))[3,1])             # Odds ratio for HIV
fm7b<-as.data.frame(t(round(exp(confint(model7)),1) [3,]))      # 95% confidence intervals of the odds ratio
fm7c<-as.data.frame(round(anova(model7,model8)[[6]][[1]],3))#  we extract p-value for HIV; to know what column and row of data is needed we can run str(anova(model7,model8))


totm7<-dt4 %>%group_by(HIV) %>%
  summarise (n = n()) %>% rename_at( 2, ~"total observations")%>% rename_at( 1, ~"varlevel")
rrm7<- dt4 %>%
  filter(rifXpertPdst == "1")%>%
  group_by(HIV) %>%
  summarise (n = n()) %>% rename_at( 2, ~"RR-TB")%>% rename_at( 1, ~"varlevel")

countm7<-left_join(totm7,rrm7, by="varlevel")


fm7ab<-cbind(fm7a,fm7b,fm7c)%>% 
  rename_at( 1, ~"odds ratio")%>%
  rename_at( 4, ~"pvalue" ) %>% mutate(varlevel = "positive")%>% 
  relocate(varlevel) #move column to first position
(fm7ab)

tm7<-left_join(countm7,fm7ab,by="varlevel")
(tm7)


###LOCATION (MODEL9)


fm9a<-as.data.frame(tibble(round(exp(coef(model9)),1))[3:4,1])             # Odds ratio for location
fm9b<-as.data.frame(round(exp(confint(model9)),1) [3:4,])     # 95% confidence intervals of the odds ratio
fm9c<-as.data.frame(round(anova(model9,model10)[[6]][[1]],3))#  we extract p-value for location; to know what column and row of data is needed we can run str(anova(model9,model10))


totm9<-dt5 %>%group_by(location) %>%
  summarise (n = n()) %>% rename_at( 2, ~"total observations")%>% rename_at( 1, ~"varlevel")
rrm9<- dt5 %>%
  filter(rifXpertPdst == "1")%>%
  group_by(location) %>%
  summarise (n = n()) %>% rename_at( 2, ~"RR-TB")%>% rename_at( 1, ~"varlevel")

countm9<-left_join(totm9,rrm9,by="varlevel")


fm9ab<-cbind(fm9a,fm9b,fm9c)%>% 
  rename_at( 1, ~"odds ratio")%>%
  rename_at( 4, ~"pvalue" ) %>% mutate(varlevel = c("rural", "urban"))%>% 
  relocate(varlevel) #move column to first position
(fm9ab)

tm9<-left_join(countm9,fm9ab,by="varlevel")
(tm9)

trfa<-rbind(tm1,tm3,tm5,tm7,tm9)

(trfa)

# ensure that the p-value is also noted in the row that corresponds to the reference level. You can type (trfa) to see the changes in the dataset

trfa <- trfa %>%
  fill("pvalue", .direction = "up")

# we now remove the duplicated p-values, so that the p-value is only noted in the first row of each variable

trfa$pvalue <- replace(trfa$pvalue, duplicated(trfa$pvalue), NA)

# replace p-value "0.000" with "<0.001"

trfa <- trfa %>%                               
  mutate(pvalue = replace(pvalue, pvalue == 0, "<0.001"))


#rename the "pvalue" column to "p-value" for aesthetics purposes

names(trfa)[names(trfa) == "pvalue"] <- "p-value"

#replace NA with blank spaces for the purposes of aesthetics when publishing the table

opts <- options(knitr.kable.NA = "")



# you may want to save an editable copy of the dataset


#write.csv(trfa, "/......./Table5.csv", row.names=F) 


kbl(trfa[1:16, 1:7], align = "lcccccl", caption = "Udjusted univariate analyses: Risk Factors for development of rifampicin-resistant TB", col.names = c("variable", "total observations", "RR-TB", "odds ratio", "2.5%", "97.5%", "p-value")) %>%
  kable_classic(full_width = F, html_font = "Cambria", font_size = 16) %>%
  pack_rows("previous treatment history", 1, 2) %>%
  pack_rows("sex", 3, 4) %>%
  pack_rows("age (years)", 5, 11) %>%
  pack_rows("HIV", 12, 13) %>%
  pack_rows("location", 14, 16) %>%
  add_header_above(c(" ", " ", " ", " ", "95% confidence interval" = 2, ""))
  
  

#--------------------------------------------------------------------------------------

# ANALYSIS OF POTENTIAL PREDICTORS OF RIFAMPICIN RESISTANT TB [RR-TB] (no imputation)
# MULTIPLE REGRESSION ANALYSIS CONSIDERING PREVIOUS TREATMENT HISTORY, SEX AND AGE GROUP

#--------------------------------------------------------------------------------------



dt6<-filter(drif, !is.na(drif$treatment) & !is.na(drif$rifXpertPdst) & !is.na(drif$sex)& !is.na(drif$ageGroup)) # remove rows with missing values
nrow(dt6)
dsdt6 <- svydesign(ids=~clusterID, weights=~w, data=dt6) #specify design
model11 <- svyglm(rifXpertPdst ~treatment + sex + ageGroup, family=quasibinomial, design=dsdt6)# to see output type summary (model11)
round(exp(coef(model11)),1) # odds ratios
round(exp(confint(model11)),1) # 95% confidence intervals for the odds ratios
anova(model11) # p-value for each explanatory variable




###############################CREATE TABLE FOR MULTIPLE REGRESSION ANALYSIS###############################


fm11a<-as.data.frame(tibble(round(exp(coef(model11)),1))[2:9,1])             # Odds ratio for age group
fm11b<-as.data.frame(round(exp(confint(model11)),1) [2:9,])     # 95% confidence intervals of the odds ratio
fm11cT<-as.data.frame(round(anova(model11)[[1]] [[6]],3))  %>% #extract p-value for previous treatment history
  rename_at( 1, ~"pvalue") %>% 
  mutate(varlevel = "new")%>% # the p-value will be annotated next to the reference variable level (i.e. new) in the final table
  relocate(varlevel)   
fm11cS<-as.data.frame(round(anova(model11)[[2]] [[6]],3)) %>% #extract p-value for sex
  rename_at( 1, ~"pvalue") %>% 
  mutate(varlevel = "female")%>% # the p-value will be annotated next to the reference variable level (i.e. female) in the final table
  relocate(varlevel) 
fm11cAG<-as.data.frame(round(anova(model11)[[3]] [[6]],3))%>% #extract p-value for age group
  rename_at( 1, ~"pvalue") %>% 
  mutate(varlevel = "0-14")%>% # the p-value will be annotated next to the reference variable level (i.e. 0-14 years of age) in the final table
  relocate(varlevel) 

fm11c<-rbind(fm11cT,fm11cS,fm11cAG)

totm11a<-dt6 %>%group_by(treatment) %>%
  summarise (n = n()) %>% rename_at( 2, ~"total observations")%>% rename_at( 1, ~"varlevel")
rrm11a<- dt6 %>%
  filter(rifXpertPdst == "1")%>%
  group_by(treatment) %>%
  summarise (n = n()) %>% rename_at( 2, ~"RR-TB")%>% rename_at( 1, ~"varlevel")


totm11b<-dt6 %>%group_by(sex) %>%
  summarise (n = n()) %>% rename_at( 2, ~"total observations")%>% rename_at( 1, ~"varlevel")
rrm11b<- dt6 %>%
  filter(rifXpertPdst == "1")%>%
  group_by(sex) %>%
  summarise (n = n()) %>% rename_at( 2, ~"RR-TB")%>% rename_at( 1, ~"varlevel")


totm11c<-dt6 %>%group_by(ageGroup) %>%
  summarise (n = n()) %>% rename_at( 2, ~"total observations")%>% rename_at( 1, ~"varlevel")
rrm11c<- dt6 %>%
  filter(rifXpertPdst == "1")%>%
  group_by(ageGroup) %>%
  summarise (n = n()) %>% rename_at( 2, ~"RR-TB")%>% rename_at( 1, ~"varlevel")

countm11a<-left_join(totm11a,rrm11a, by="varlevel")
countm11b<-left_join(totm11b,rrm11b, by="varlevel")
countm11c<-left_join(totm11c,rrm11c, by="varlevel")



countm11<-rbind(countm11a, countm11b, countm11c)


fm11ab<-cbind(fm11a,fm11b)%>% 
  rename_at( 1, ~"odds ratio")%>%
  mutate(varlevel = c("retreatment", "male", "15-24", "25-34", "35-44","45-54","55-64","65+"))%>% 
  relocate(varlevel) #move column to first position
(fm11ab)

tm11a<- left_join(countm11, fm11ab, by="varlevel")
tm11b<- left_join(tm11a,fm11c, by="varlevel")


# replace p-value "0.000" with "<0.001"

tm11b <- tm11b %>%                               
  mutate(pvalue = replace(pvalue, pvalue == 0, "<0.001"))


#rename the "pvalue" column to "p-value" for aesthetics purposes

names(tm11b)[names(tm11b) == "pvalue"] <- "p-value"

#replace NA with blank spaces for the purposes of aesthetics when publishing the table


opts <- options(knitr.kable.NA = "")


# you may want to save an editable copy of the dataset


#write.csv(tm11b, "/......./Table6.csv", row.names=F) 


kbl(tm11b[1:11, 1:7], align = "lcccccl", caption = "Multiple regression analysis: Risk Factors for development of rifampicin-resistant TB", col.names = c("variable", "total observations", "RR-TB", "odds ratio", "2.5%", "97.5%", "p-value")) %>%
  kable_classic(full_width = F, html_font = "Cambria", font_size = 16) %>%
  pack_rows("previous treatment history", 1, 2) %>%
  pack_rows("sex", 3, 4) %>%
  pack_rows("age (years)", 5, 11) %>%
  add_header_above(c(" ", " ", " ", " ", "95% confidence interval" = 2, ""))





