RequestID <- "DR3616"
FolderName <- paste(RequestID, "Ryan White Conference Abstract")
ProgramName <- "Ryan White Viral Load Analysis"  # If the DR has just one program, then that ProgramNam is usually the same as FolderName. IF YOU CHANGE IT, make sure to start it with RequestID.
DataSource <- "eHARS"
ContactName <- "Sam Parmar"
DRdir <- paste(".../Data Requests/", FolderName, "/", sep = "")

# Created by Sam Parmar on 2018-11-16
# Validated by _name_ on 2015-0x-xx
# Purpose: Arithmetic Mean and Geometric Mean Viral Load Program for 2019 RW Nat'l Conference Workshop
# Production schedule: (if program is to be run on a regular basis)
# Limitations and Warnings
# Program derived from "DR3495 Viral Load Disparity Abstract.SAS"
# Program Flow Description (high level review of the steps of the program)---- 
#  1) Define key constants, formats, etc.
#  2) Standardize eHARS_calc_VLstart and eHARS_calc_VLend 
#  3) Calculate Log10 of most recent individual viral load during Year.
#  4) Functions to calculate geometric and arithmetic community viral loads and 95% CIs
#  5) Modify relevant comparison categories after VL variable (example: year, race_rec)
#  X) Clean up. 

# When this program is ready for the final run (i.e., it is what you want and error free), put its full path and filename as the argument of final_run(), and run that within R. That will save a log and output file in the same folder as the program.
# e.g., final_run("Program_name.R")
#source(".../final_run.R")

# Message string for clients to know who to contact about program output. 
ContactMe <- paste(RequestID, ContactName, "EPI", Sys.Date(), "Data Source:", DataSource)

# Add the following to all plots created by ggplot2: +labs(caption = ContactMe)

# Leave mark at beginning of program log see when the run started. 
print("###############################################")
paste("Run started at", Sys.time())
print("###############################################")

# Note that the above info is for internal purposes

#  1) Define key constants, formats, etc. ; ----

#this function will install the necessary packages (if not installed already) and then load libraries
list.of.packages <- c("tidyverse","gmodels")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) {
  install.packages(new.packages)
}
lapply(list.of.packages, require, character.only = TRUE)

#declare and import HIV data:

#declare file directory here (make sure to use forward slash / )
file_dir=".../DR3254 RWSP Viral Load Protocol for Providers/"
#include data file name here (program currently requires csv datafile):
file_name="CombinedDat.csv"

full_fPath<-paste(file_dir,file_name,sep="")

#read csv data into R 
Combined0<-read.csv(full_fPath, header=TRUE)

#declare variables of interest
list_c<-c("ehars_uid_rwg","Year","eHARS_calc_VLend","RWGclient","SrvInYr","MedVstGap",
"Retained","cur_county_name2","cur_zip_cd","race_rec","birth_sex","sextgcnd","cur_age_num",
"Risk","rwg_FacilityID","linktocare","prevdx_status", "pct_poverty")

#the rest of the program is based on the variables above and a few others
#that you will need to rename accordingly 
#you can remove or add variables as needed as long as they are included in the csv data and declared in the program correctly

#eHARS_calc_VLend is the last reported viral load for PLWHA
#SrvInYr is Service in Year (0 or 1)
#MedVstGap is gap in medical services during the year (0 or 1)
#retained is retained in care (0 or 1)
#cur_county_name2 is current county name
#cur_zip_cd is current zip code
#race_rec is recent race
#birth_sex is sex at birth
#sextgcnd is gender
#cur_age_num is current age
#Risk is HIV risk category
#rwg_FacilityID is the ryan white facility id
#link to care is linked to care status
#prevdx_status is HIV status
#pct_poverty is percentage of poverty

#result_range_lower is the lower limit of detection for VL test
#result_range_higher is the higher limit of detection for VL test
#eHARS_calc_VLstart is the first viral load reported in the time period
#eHARS_calc_VLend is the last viral load reported in the time period
#eHARS_calc_VLcount is the number of viral load tests conducted in the time period

#2) Standardize eHARS_calc_VLstart and eHARS_calc_VLend ----
# In recent years, the lower limit of detection (LLD) of HIV viral load has decreased to <20 RNA copies/mL. When
# reporting undetectable results (results below the LLD), laboratories use a variety of conventions. Typically,
# they report the results as either 0 or <50 copies/mL. To standardize these results and facilitate a more accurate
# statistical analysis, VL test results reported as 0 or <50 copies/mL will be analyzed using a number equal to half
# the reported LLD. For example, if a VL was reported as <50 and the respective LLD is reported to be 20, then the VL
# result is replaced with 10. This is the CDC-recommended process for analyzing HIV viral load (see website below).
# 
# The lower/upper detectable limits of many viral load results reported in eHARS is not reported, especially for
# early years. Where these values are missing, they were replaced by the most commonly reported values for the year
# in which the viral load test sample was drawn. These values were set for 2011-2015 as follows:


Combined <-Combined0 %>% mutate(
  result_range_lower=ifelse(is.na(result_range_lower), 50, result_range_lower),
  result_range_upper=ifelse(is.na(result_range_upper), 750000, result_range_upper)) %>% 
  mutate(result_range_upper = as.numeric(result_range_upper),
         result_range_lower = as.numeric(result_range_lower),
         eHARS_calc_VLcount = as.numeric(eHARS_calc_VLcount),
         eHARS_calc_VLend = as.numeric(eHARS_calc_VLend),
         eHARS_calc_VLstart = as.numeric(eHARS_calc_VLstart)) 

Combined <- Combined %>% 
  mutate(eHARS_calc_VLstart=case_when(
    eHARS_calc_VLcount>0 & !is.na(result_range_lower) & eHARS_calc_VLstart >= result_range_lower & eHARS_calc_VLstart <= result_range_upper ~ eHARS_calc_VLstart,
    eHARS_calc_VLcount>0 & !is.na(result_range_lower) & eHARS_calc_VLstart < result_range_lower ~ 0.5*result_range_lower,
    eHARS_calc_VLcount>0 & !is.na(result_range_lower) & eHARS_calc_VLstart > result_range_lower ~ result_range_upper+1
    )
,eHARS_calc_VLend=case_when(
    eHARS_calc_VLcount>0 & !is.na(result_range_lower) &  eHARS_calc_VLend >= result_range_lower & eHARS_calc_VLend <= result_range_upper ~ eHARS_calc_VLend,
    eHARS_calc_VLcount>0 & !is.na(result_range_lower) & eHARS_calc_VLend < result_range_lower ~ 0.5*result_range_lower,
    eHARS_calc_VLcount>0 & !is.na(result_range_lower) & eHARS_calc_VLend > result_range_lower ~ result_range_upper+1
))

Combined <- Combined %>%
  mutate(
    eHARS_calc_VLend= case_when(
      eHARS_calc_VLcount == 1 ~ eHARS_calc_VLstart,
      eHARS_calc_VLcount > 1 & is.na(eHARS_calc_VLend) ~ eHARS_calc_VLstart,
      eHARS_calc_VLcount > 1 & !is.na(eHARS_calc_VLend) ~ eHARS_calc_VLend
))

#select only variables included in list above called list_c
Combined <- Combined %>% select(list_c) 

# 3) Calculate Log10 of most recent individual viral load during Year. ** -----
  
# The CDC recommends the following for statistical analysis of community viral load:
#   
# "The rationale for this logarithmic transformation is that it helps to
# normalize the distribution of viral load values and reduces the influence of
# outlying measurements for persons having extreme viremia (due to acute infection,
# advanced HIV disease, concurrent sexually transmitted infections, or random
# variability). Since the interpretation of viral load measurements is often more
# intuitive on a linear scale, we recommend calculation of geometric mean (GM) for
# viral load. The GM is thus calculated through log transformation by averaging
# the log transformed values and transforming the average back to the original
# (linear) scale. The base used for the log transformation has no effect on the
# final GM estimate. However, using log base 10 has an advantage by its
# relationship to the value on the original scale - for example, a value of 2 on
# the log10 scale is 100 on the original scale, 3 corresponding to 1000, 4
# corresponding 10000, and so forth."
# 
# Centers for Disease Control and Prevention. (2011). Guidance on community viral
# load: A family of measures, definitions, and method for calculation. Retrieved
# from http://www.ct.gov/dph/lib/dph/aids_and_chronic/surveillance/statewide/
#   community_viralload_guidance.pdf
# 
# ** CD4/VL values chosen because, "Data suggest that during a period of up to 2 years,
#    patients with viral rebound to levels below 10,000 copies/mL may be able to
#    maintain their CD4 counts despite a failing regimen."
# - http://www.hivguidelines.org/wp-content/uploads/viral-load-report.pdf ;
# 
# ** Geometric mean is often used to evaluate data covering several orders of magnitude ;

#create dataset called FinDat where there are no missing Viral load values and there is a 
FinDat <- Combined %>% filter(!is.na(eHARS_calc_VLend)) %>% ## Where last VL during period not missing 
  mutate(
    LogVLend= case_when(
      eHARS_calc_VLend > 0  ~ log10(eHARS_calc_VLend)
))

# 4) Functions to calculate geometric and arithmetic community viral loads and 95% CIs ----
#inputs: raw viral load values for arithmetic mean calc or log base 10 VL values of viral load 
#for geometric viral load calc
#outputs: mean, low 95% CI, and high 95% CI

#geometric mean calculation function
calc_geoCVL<-function(dataset, VL_variable,...) {
  confidence=0.95
  VL_variable <-enquo(VL_variable)
  group_col<-quos(...)
  suppressWarnings(
  dataset %>% group_by(!!!group_col) %>% summarise(mean = ci(!!VL_variable, confidence)[1],
                                                     lowCI = ci(!!VL_variable, confidence)[2],
                                                     hiCI = ci(!!VL_variable, confidence)[3]
                                                     #,sd = ci(!!VL_variable, confidence)[4]
                                                   ) %>% mutate(
                                                     mean=case_when(!is.na(mean) ~ 10**mean),
                                                     lowCI=case_when( 
                                                       !is.na(lowCI) & lowCI >= 0 ~ 10**lowCI,
                                                       !is.na(lowCI) & lowCI < 0 ~ 0),
                                                     
                                                     hiCI=case_when( 
                                                       !is.na(hiCI) & hiCI >= 0 ~ 10**hiCI,
                                                       !is.na(hiCI) & hiCI < 0 ~ 0))  
  )
}

#arithmetic mean calculation function
calc_ariCVL<-function(dataset, VL_variable,...) {
  confidence=0.95
  VL_variable <-enquo(VL_variable)
  group_col<-quos(...)
  suppressWarnings(
    dataset %>% group_by(!!!group_col) %>% summarise(mean = ci(!!VL_variable, confidence)[1],
                                                     lowCI = ci(!!VL_variable, confidence)[2],
                                                     hiCI = ci(!!VL_variable, confidence)[3]
                                                     #,sd = ci(!!VL_variable, confidence)[4]
    ) %>% mutate(
      lowCI=case_when( 
        !is.na(lowCI) & lowCI >= 0 ~ lowCI,
        !is.na(lowCI) & lowCI < 0 ~ 0),
      hiCI=case_when( 
        !is.na(hiCI) & hiCI >= 0 ~ hiCI,
        !is.na(hiCI) & hiCI < 0 ~ 0))  
  )
}


# 5) modify relevant comparison categories after VL variable (example: year, race_rec)  ----
#do this in the in calc_geoCVL() and grepl() functions 
#for both geometric mean grouping, arithmetic mean grouping, and merge

#geometric mean groupings
GeoMeanYear <- FinDat %>% calc_geoCVL(., LogVLend, Year) %>% rename_if(!grepl("Year", names(.)), funs(sprintf("geo_%s", .)))
GeoMeanRace <- FinDat %>% calc_geoCVL(., LogVLend, Year, race_rec) %>%  rename_if(!grepl("Year|race_rec", names(.)), funs(sprintf("geo_%s", .)))
GeoMeanAge <- FinDat %>% calc_geoCVL(., LogVLend, Year, cur_age_num) %>%  rename_if(!grepl("Year|cur_age_num", names(.)), funs(sprintf("geo_%s", .)))
GeoMeanGender <- FinDat %>% calc_geoCVL(., LogVLend, Year,sextgcnd) %>%  rename_if(!grepl("Year|sextgcnd", names(.)), funs(sprintf("geo_%s", .)))
GeoMeanCount <- FinDat %>% calc_geoCVL(., LogVLend, Year, cur_county_name2) %>%  rename_if(!grepl("Year|cur_county_name2", names(.)), funs(sprintf("geo_%s", .)))
GeoMeanRisk <- FinDat %>% calc_geoCVL(., LogVLend, Year, Risk) %>%  rename_if(!grepl("Year|Risk", names(.)), funs(sprintf("geo_%s", .)))
GeoMeanRetained <- FinDat %>% calc_geoCVL(., LogVLend, Year, Retained) %>%  rename_if(!grepl("Year|Retained", names(.)), funs(sprintf("geo_%s", .)))
GeoMeanRWGStat <- FinDat %>% calc_geoCVL(., LogVLend, Year, RWGclient) %>%  rename_if(!grepl("Year|RWGclient", names(.)), funs(sprintf("geo_%s", .)))
GeoMeanZip <- FinDat %>% calc_geoCVL(., LogVLend, Year, cur_zip_cd) %>%  rename_if(!grepl("Year|cur_zip_cd", names(.)), funs(sprintf("geo_%s", .)))

#arithmetic mean groupings
AriMeanYear <- FinDat %>% calc_ariCVL(., eHARS_calc_VLend, Year) %>% rename_if(!grepl("Year", names(.)), funs(sprintf("ari_%s", .)))
AriMeanRace <- FinDat %>% calc_ariCVL(., eHARS_calc_VLend, Year, race_rec) %>% rename_if(!grepl("Year|race_rec", names(.)), funs(sprintf("ari_%s", .)))
AriMeanAge <- FinDat %>% calc_ariCVL(., eHARS_calc_VLend, Year, cur_age_num) %>% rename_if(!grepl("Year|cur_age_num", names(.)), funs(sprintf("ari_%s", .)))
AriMeanGender <- FinDat %>% calc_ariCVL(., eHARS_calc_VLend, Year, sextgcnd) %>% rename_if(!grepl("Year|sextgcnd", names(.)), funs(sprintf("ari_%s", .)))
AriMeanCount <- FinDat %>% calc_ariCVL(., eHARS_calc_VLend, Year, cur_county_name2) %>% rename_if(!grepl("Year|cur_county_name2", names(.)), funs(sprintf("ari_%s", .)))
AriMeanRisk <- FinDat %>% calc_ariCVL(., eHARS_calc_VLend, Year, Risk) %>% rename_if(!grepl("Year|Risk", names(.)), funs(sprintf("ari_%s", .)))
AriMeanRetained <- FinDat %>% calc_ariCVL(., eHARS_calc_VLend, Year, Retained) %>% rename_if(!grepl("Year|Retained", names(.)), funs(sprintf("ari_%s", .)))
AriMeanRWGStat <- FinDat %>% calc_ariCVL(., eHARS_calc_VLend, Year, RWGclient) %>% rename_if(!grepl("Year|RWGclient", names(.)), funs(sprintf("ari_%s", .)))
AriMeanZip <- FinDat %>% calc_ariCVL(., eHARS_calc_VLend, Year, cur_zip_cd) %>% rename_if(!grepl("Year|cur_zip_cd", names(.)), funs(sprintf("ari_%s", .)))

#merge geometric mean results side by side with arithmetic mean results
MeanYear <-merge(GeoMeanYear, AriMeanYear, by="Year")
MeanRace <-merge(GeoMeanRace, AriMeanRace, by=c("Year","race_rec"))
MeanAge <-merge(GeoMeanAge, AriMeanAge, by=c("Year","cur_age_num"))
MeanGender <-merge(GeoMeanGender, AriMeanGender, by=c("Year","sextgcnd"))
MeanCount <-merge(GeoMeanCount, AriMeanCount, by=c("Year","cur_county_name2"))
MeanRisk <-merge(GeoMeanRisk, AriMeanRisk, by=c("Year","Risk"))
MeanRetained <-merge(GeoMeanRetained, AriMeanRetained, by=c("Year","Retained"))
MeanRWGStat <-merge(GeoMeanRWGStat, AriMeanRWGStat, by=c("Year","RWGclient"))
MeanZip <-merge(GeoMeanZip, AriMeanZip, by=c("Year","cur_zip_cd"))

#  X) Clean up. ----
# Delete all objects that were created.
rm(RequestID, FolderName, ProgramName, DataSource, ContactName, DRdir, final_run, ContactMe)

# Leave mark at end of program to see when the run ended.
print("###############################################")
paste("Run ended at", Sys.time())
print("###############################################")
