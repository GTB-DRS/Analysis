
                                  #################################
################################### VISUALISATION OF CLUSTER SIZE ##############################################################
                                  #################################
                                  
#The following script produces an example plot to help visualise the size of clusters for cluster-based surveys. The script and
#plot may be adapted to accommodate for other forms of cluster size visualisation compared to expected targets.


#Let's assume that a survey used cluster sampling. Thirty-eight clusters were selected using probability
#proportional to size sampling (PPS). All clusters were asked to enroll 60 bacteriologically confirmed
# new pulmonary TB cases (fixed cluster size). In order to assess the need for weights in the regression analysis
#fitted to estimate the prevalence of resistance to rifampicin, it is useful to visualise and compare
#the expected cluster size (n=60 cases per cluster) with the observed cluster size (i.e. total number of new
#cases ultimately enrolled in each cluster)

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


# open file; replace "." with the file location for "dsrSampleData.csv":
                                  

d1<-fread("/......./dsrSampleData.csv")                          
                                  
# view the dataset

names(d1) # view column names
View(d1) # view dataset

#check the size of the dataset

ncol(d1) # number of columns
nrow(d1) # number of observations

# check that patient DRS IDs are unique

anyDuplicated(d1$drsID) # there are no duplicate DRS IDs

# check for missing values in the dataset

na_count <-sapply(d1, function(y) sum(length(which(is.na(y)))))
data.frame(na_count) # there are missing valued for treatment history

#replace missing values in treatment column with "unknown"

d1$treatment[is.na(d1$treatment)] = "unknown"
table(d1$treatment)



# Plot the size of the clusters (i.e. number of new bacteriologically confirmed pulmonary TB cases enrolled in each cluster, as informed by the initial Xpert MTB/RIF test) compared to the expected [target] size. 
# Begin by creating a new dataset (dataCS), that summarises the number of enrolled (and confirmed) TB cases by cluster (new and overall), and arranges the data in the
# format required for the plot:

dataCS<-d1 %>% filter(xpertResult != "mtbNegative")%>%
  group_by(clusterID) %>% 
  summarise(observedAll= n(), observedNew= sum(treatment=="new"))%>% 
  mutate (targetNew=60) %>% 
  gather(group, size, observedAll:targetNew)


View(dataCS)

# Produce the plot. For this, ensure clusterID is recognised as a factor rather than a number:

dataCS$clusterID<-as.factor(dataCS$clusterID)

csFigure<-ggdotchart(dataCS, x = "clusterID", y ="size", 
           color = "group",                 
           palette = c("#E7B80090","#00AFBB20", "#FC4E0720"), 
           sorting = "descending",                       
           add = "segments",                             
           rotate = TRUE,                                
           group = "group",                                
           dot.size = 6, 
           font.xtickslab= 7,   
           font.ytickslab= 7, 
           label = round(dataCS$size),                        
           font.label = list(color = "black", size = 8, 
                             vjust = 0.5),               
           ggtheme = theme_pubr())

plot(csFigure)

#If you wish to store a high resolution tiff picture of the
#plot, set up a working directory where you want to store the figure, and run the following script:

setwd("/......./.....")
tiff("csFigure.tiff", width= 3858.333333, height= 3229.166667, units="px", res=300, compression = 'lzw')
plot(csFigure)
dev.off()

#Note you can modify the properties of the picture as needed.






