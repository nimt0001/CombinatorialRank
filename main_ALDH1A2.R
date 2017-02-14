########################## KM Analysis new for Frontier Oncology paper - ALDH1A2 in prostate cancer
library(survival)
library(MASS)
library(cgdsr)
library(HapEstXXR)

rm(list=ls(all=TRUE)) 
mycgds = CGDS("http://www.cbioportal.org/public-portal/")

#### Get available case lists - PRAD ALL - 
### Note that the TCGA database is updated regularly and results may change from time to time 
mycancerstudy = getCancerStudies(mycgds)[118,1] #118 18
mycaselist = getCaseLists(mycgds,mycancerstudy)[2,1]  #[2,1] [1,1]

# Get available genetic profiles
mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy) #All
mygeneticprofile_cnv = mygeneticprofile[2,1] #Copy number variation
mygeneticprofile_exp = getGeneticProfiles(mycgds,mycancerstudy)[4,1]  
mygeneticprofile_exp_z = getGeneticProfiles(mycgds,mycancerstudy)[3,1] #3 
mygeneticprofile_mut = getGeneticProfiles(mycgds,mycancerstudy)[6,1]  
myRNAdata_z = getProfileData(mycgds,c('ALDH1A2','CYP26C1','CYP26A1', 'CYP26B1', 'RDH10', 'ADH5', 'DHRS3', 'ADH4', 'ADH7', 'ADH1B', 'ADH1A'),mygeneticprofile_exp_z,mycaselist)
myRNAdata = getProfileData(mycgds,c('ALDH1A2','CYP26C1','CYP26A1', 'CYP26B1', 'RDH10', 'ADH5', 'DHRS3', 'ADH4', 'ADH7', 'ADH1B', 'ADH1A'),mygeneticprofile_exp,mycaselist)
all_event_rna <- apply(as.data.frame(myRNAdata_z), 1, function(x) ifelse(abs(x) > 1.96,1,0))   
sum(all_event_rna,na.rm=1)

myclinicaldata = getClinicalData(mycgds,mycaselist)
clinical_IDs <- rownames(myclinicaldata)
all_clin <- data.frame(myclinicaldata$DFS_MONTHS) #DFS_MONTHS OS_MONTHS
colnames(all_clin) <- c("new_dfsMonths")
rownames(all_clin) <- clinical_IDs
all_clin$dfs_events <- ifelse(myclinicaldata$DFS_STATUS == "DiseaseFree", 0,1)#DFS_STATUS == "DiseaseFree" OS_STATUS == "LIVING"

ps <- powerset(1:11)
my_event_rna_list <- list()
s_list <- list()
s1_list <- list()
pv_list <- list()

for (i in 1:length(ps)){
  my_event_rna <- rep(0, 499) 
  for (j in 1:length(ps[[i]])){
    my_event_rna <- bitwOr(my_event_rna, all_event_rna[ps[[i]][j],]) 
  }
  my_event_rna_list[[length(my_event_rna_list)+1]] <- my_event_rna
  s_temp <- Surv(as.numeric(as.character(all_clin$new_dfsMonths))[],all_clin$dfs_events[])
  s_new <- survfit(s_temp~my_event_rna)  
  s1_new <- tryCatch(survdiff(s_temp~my_event_rna), error = function(e) return(NA))
  pv_new <- ifelse(is.na(s1_new),next,(round(1 - pchisq(s1_new$chisq, length(s1_new$n) - 1),3)))[[1]]
  s_list[[length(s_list)+1]] <- s_new
  s1_list[[length(s1_list)+1]] <- s1_new
  pv_list[[length(pv_list)+1]] <- pv_new
}

plot(1:length(ps),unlist(pv_list)[1:length(ps)], log="y", xlab="Gene set ID", ylab="p-value",
     lwd=2, main="Gene set ID versus p-values")

# Plot the KM curve
my_event_rna = my_event_rna_list[[id=257]] 
s_temp <- Surv(as.numeric(as.character(all_clin$new_dfsMonths)),all_clin$dfs_events)
s_new <- survfit(s_temp~my_event_rna)  
s1_new <- tryCatch(survdiff(s_temp~my_event_rna), error = function(e) return(NA))
pv_new <- ifelse(is.na(s1_new),next,(round(1 - pchisq(s1_new$chisq, length(s1_new$n) - 1),3)))[[1]]
my_event_rna_name <- paste("") 
plot(s_new,col=c(1:3), frame=F, mark.time=TRUE, lwd=2,main=paste("Kaplan-Meier survival analysis",my_event_rna_name,sep="\n"))

# Cox regression 
prostate <- c()
prostate$time <- as.numeric(as.character(all_clin$new_dfsMonths))
prostate$status <- all_clin$dfs_events
prostate$SurvObj <- with(prostate, Surv(time, status == 1))
prostate$rna <- my_event_rna
res.cox1 <- coxph(SurvObj ~ rna, data =  prostate) #age + rna + gleason + tstage
summary(res.cox1)

########################## Combinatorial Ranking based on KM analysis - ALDH1A2 in Breast Cancer
library(survival)
library(MASS)
library(cgdsr)
library(HapEstXXR)

rm(list=ls(all=TRUE)) 
mycgds = CGDS("http://www.cbioportal.org/public-portal/")

#### Get available case lists - Breast Cancer
### Note that the TCGA database is updated regularly and results may change from time to time 
mycancerstudy = getCancerStudies(mycgds)[18,1] #118 18
mycaselist = getCaseLists(mycgds,mycancerstudy)[1,1]  #[2,1] [1,1]

# Get available genetic profiles
mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy) #All
mygeneticprofile_cnv = mygeneticprofile[6,1] #Copy number variation
mygeneticprofile_exp = getGeneticProfiles(mycgds,mycancerstudy)[3,1]  
mygeneticprofile_exp_z = getGeneticProfiles(mycgds,mycancerstudy)[2,1] 
mygeneticprofile_mut = getGeneticProfiles(mycgds,mycancerstudy)[4,1]  
myRNAdata_z = getProfileData(mycgds,c('ALDH1A2','CYP26C1','CYP26A1', 'CYP26B1', 'RDH10', 'ADH5', 'DHRS3', 'ADH4', 'ADH7', 'ADH1A'),mygeneticprofile_exp_z,mycaselist)
myRNAdata = getProfileData(mycgds,c('ALDH1A2','CYP26C1','CYP26A1', 'CYP26B1', 'RDH10', 'ADH5', 'DHRS3', 'ADH4', 'ADH7', 'ADH1A'),mygeneticprofile_exp,mycaselist)
all_event_rna <- apply(as.data.frame(myRNAdata_z), 1, function(x) ifelse(abs(x) > 1.96,1,0))   
sum(all_event_rna,na.rm=1)

myclinicaldata = getClinicalData(mycgds,mycaselist)
clinical_IDs <- rownames(myclinicaldata)
all_clin <- data.frame(myclinicaldata$OS_MONTHS) #DFS_MONTHS OS_MONTHS
colnames(all_clin) <- c("new_dfsMonths")
rownames(all_clin) <- clinical_IDs
all_clin$dfs_events <- ifelse(myclinicaldata$OS_STATUS == "LIVING", 0,1)#DFS_STATUS == "DiseaseFree" OS_STATUS == "LIVING"

ps <- powerset(1:dim(myRNAdata)[2])
my_event_rna_list <- list()
s_list <- list()
s1_list <- list()
pv_list <- list()

for (i in 1:length(ps)){
  my_event_rna <- rep(0, dim(all_event_rna)[2]) 
  for (j in 1:length(ps[[i]])){
    my_event_rna <- bitwOr(my_event_rna, all_event_rna[ps[[i]][j],]) 
  }
  my_event_rna_list[[length(my_event_rna_list)+1]] <- my_event_rna
  s_temp <- Surv(as.numeric(as.character(all_clin$new_dfsMonths))[],all_clin$dfs_events[])
  s_new <- survfit(s_temp~my_event_rna)  
  s1_new <- tryCatch(survdiff(s_temp~my_event_rna), error = function(e) return(NA))
  pv_new <- ifelse(is.na(s1_new),next,(round(1 - pchisq(s1_new$chisq, length(s1_new$n) - 1),3)))[[1]]
  s_list[[length(s_list)+1]] <- s_new
  s1_list[[length(s1_list)+1]] <- s1_new
  pv_list[[length(pv_list)+1]] <- pv_new
}

plot(1:length(ps),unlist(pv_list)[1:length(ps)], log="y", xlab="Gene set ID", ylab="p-value",
     lwd=2, main="Gene set ID versus p-values")

# Plot the KM curve for the optimal gene set (OGS)
my_event_rna = my_event_rna_list[[which.min(unlist(pv_list))]]
s_temp <- Surv(as.numeric(as.character(all_clin$new_dfsMonths)),all_clin$dfs_events)
s_new <- survfit(s_temp~my_event_rna)  
s1_new <- tryCatch(survdiff(s_temp~my_event_rna), error = function(e) return(NA))
pv_new <- ifelse(is.na(s1_new),next,(round(1 - pchisq(s1_new$chisq, length(s1_new$n) - 1),3)))[[1]]
my_event_rna_name <- paste("") 
plot(s_new,col=c(1:3), frame=F, mark.time=TRUE, lwd=2,main=paste("Kaplan-Meier survival analysis",my_event_rna_name,sep="\n"))

# Cox regression 
brca <- c()
brca$time <- as.numeric(as.character(all_clin$new_dfsMonths))
brca$status <- all_clin$dfs_events
brca$SurvObj <- with(brca, Surv(time, status == 1))
brca$rna <- my_event_rna
res.cox1 <- coxph(SurvObj ~ rna, data =  brca) 
summary(res.cox1)






