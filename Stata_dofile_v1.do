
////////////////////////* SAMPLE STATA DOFILE FOR ANALYSIS OF DATA FROM ANTI-TB DRUG RESISTANCE SURVEYS, MAY 2021 *//////////////////////

/* This analysis is based on a cluster-based survey with variable cluster size. The clusters are distributed across three different strata. Each strata has its own target sample size. */
/* Clusters in the same strata enrol consecutive patients over the same period of time until the strata sample size is met. */

/* Set working drive */
cd "C:\Users\deanan\OneDrive - World Health Organization\Documents\dr surveys\data analysis tools\sample data and dofile\may 2021" /* change working drive to match yours */


/******************************************************Data cleaning*************************************************************************/

/* import data */
clear
insheet using "Stata_sampleDataFinal.csv", comma names 


/* check that patient ID is unique */
duplicates list drsid


/******** Patient characteristics *******/

/* generate binary treatment history variable */
tab treatment, m
gen treat_hist = .
replace treat_hist = 1 if treatment=="retreatment"
replace treat_hist = 0 if treatment=="new"
label variable treat_hist "previous treatment history"
label define treat_hist 0 "new" 1 "retreatment" 
label values treat_hist treat_hist
tab treat_hist, m
drop treatment


/* generate age groups */
hist age, freq
label variable age "age"
gen agegrp=age
label variable agegrp "age group"
recode agegrp min/14=1 15/24=2 25/34=3 35/44=4 45/54=5 55/64=6 65/max=7
label define agegrp 1 "0-14 yrs" 2 "15-24 yrs" 3 "25-34 yrs" 4 "35-44 yrs" 5 "45-54 yrs" 6 "55-64 yrs" 7 "65+ yrs"
label values agegrp agegrp
codebook agegrp
order agegrp, after(age)


/* review and label sex variable */
tab sex, m
gen sex2 = .
replace sex2 = 1 if sex=="female"
replace sex2 = 0 if sex=="male"
label variable sex2 "sex"
label define sex2 0 "Male" 1 "Female" 
label values sex2 sex2
tab sex2, m
drop sex
rename sex2 sex


/* review and label sex variable */
tab hiv, m
gen hiv2 = .
replace hiv2 = 1 if hiv=="positive"
replace hiv2 = 0 if hiv=="negative"
label variable hiv2 "hiv"
label define hiv2 0 "negative" 1 "positive" 
label values hiv2 hiv2
tab hiv2, m
drop hiv
rename hiv2 hiv



/* histogram by age and sex - the script below must be run all together, not line by line */
/* see output in pdng file "Stata_age_sex_distribution" */
preserve
drop if agegrp==. | sex==.
bysort sex agegrp: gen totsize=_N 
collapse (mean) totsize, by(agegrp sex)
reshape wide totsize, i(agegrp) j(sex) 
rename totsize1 female
rename totsize0 male
replace male = -male
gen minus_100 = -100 
twoway bar male agegrp, horizontal xvarlab(Male)|| bar  female agegrp, horizontal xvarlab(Female) || sc  agegrp minus_100, mlabel(agegrp) mlabcolor(black) msymbol(i) ||, ///
xtitle("Number of patients", margin(0 0 0 5)) ytitle("") plotregion(style(none)) ysca(noline) ylabel(none) xsca(noline titlegap(-3.5)) xlabel(-100 "100" -200 "200" -300 "300" 100(100)300 , ///
tlength(0) grid gmin gmax) legend(label(1 Male) label(2 Female)) legend(order(1 2)) title("Age distribution")
restore


/* generate residence variable */
tab location, m
gen residence = .
replace residence = 3 if location=="rural"
replace residence = 2 if location=="periUrban"
replace residence = 1 if location=="urban"
label define residence 1 "urban" 2 "periUrban" 3 "rural"
label values residence residence
tab residence, m
drop location



/****** Laboratory results****/

/* generate Xpert variable */
tab xpertresult, m
gen xpert=.
replace xpert=0 if xpertresult=="mtbRS" 
replace xpert=1 if xpertresult=="mtbRR" 
replace xpert=2 if xpertresult=="mtbNegative" 
label variable xpert "Xpert"
label define xpert 0 "MTB pos, RIF SUS"  1 "MTB pos, RIF RES" 2 "MTB neg" 
label values xpert xpert
tab xpert, m
drop xpertresult
codebook xpert


/* generate binary phenotypic DST variables */       

tab cultureresult rifpdst, m /* note that some with phenotypic DST for rifampicin are actually NTMs, not MTBC and so pheno DST should not be used */

gen rifPheno=.
replace rifPheno=1 if rifpdst=="r" & cultureresult=="mtbc"
replace rifPheno=0 if rifpdst=="s" & cultureresult=="mtbc"
label variable rifPheno "rif pheno"
label define rifPheno 0 "susceptible" 1 "resistant"
label values rifPheno rifPheno
tab rifPheno, m
drop rifpdst

gen inhPheno=.
replace inhPheno=1 if inhpdst=="r" & cultureresult=="mtbc"
replace inhPheno=0 if inhpdst=="s" & cultureresult=="mtbc"
label variable inhPheno "inh pheno"
label define inhPheno 0 "susceptible" 1 "resistant"
label values inhPheno inhPheno
tab inhPheno, m
drop inhpdst

gen lfxPheno=.
replace lfxPheno=1 if lfxpdst=="r" & cultureresult=="mtbc"
replace lfxPheno=0 if lfxpdst=="s" & cultureresult=="mtbc"
label variable lfxPheno "lfx pheno"
label define lfxPheno 0 "susceptible" 1 "resistant"
label values lfxPheno lfxPheno
tab lfxPheno, m
drop lfxpdst

gen bdqPheno=.
replace bdqPheno=1 if bdqpdst=="r" & cultureresult=="mtbc"
replace bdqPheno=0 if bdqpdst=="s" & cultureresult=="mtbc"
label variable bdqPheno "bdq pheno"
label define bdqPheno 0 "susceptible" 1 "resistant"
label values bdqPheno bdqPheno
tab bdqPheno, m
drop bdqpdst

gen lzdPheno=.
replace lzdPheno=1 if lzdpdst=="r" & cultureresult=="mtbc"
replace lzdPheno=0 if lzdpdst=="s" & cultureresult=="mtbc"
label variable lzdPheno "bdq pheno"
label define lzdPheno 0 "susceptible" 1 "resistant"
label values lzdPheno lzdPheno
tab lzdPheno, m
drop lzdpdst


/* generate binary combined DST variables */

gen rifXpertPdst=.
replace rifXpertPdst = 0 if (rifPheno==0 | xpert==0)
replace rifXpertPdst = 1 if (rifPheno==1 | xpert==1) /* resistant by either Xpert or pheno DST is classified as resistant */
label variable rifXpertPdst "combined RR-TB"
label define rifXpertPdst 0 "susceptible" 1 "resistant"
label values rifXpertPdst rifXpertPdst
tab rifXpertPdst, m

gen mdr=.
replace mdr=0 if (inhPheno==0 | rifXpertPdst==0)  
replace mdr=1 if inhPheno==1 & rifXpertPdst==1  
label variable mdr "MDR-TB"
label define mdr 0 "not MDR" 1 "MDR"
label values mdr mdr
tab mdr, m

gen preXdr=. /* pre-XDR-TB is defined as RR-TB plus resistance to fluoroquinolones */
replace preXdr = 0 if rifXpertPdst==0
replace preXdr = 0 if rifXpertPdst==1 & lfxPheno==0
replace preXdr = 1 if rifXpertPdst==1 & lfxPheno==1
label variable preXdr "pre-XDR-TB"
label define preXdr 0 "not pre-XDR" 1 "pre-XDR"
label values preXdr preXdr

gen xdr=. /* XDR-TB is defined as RR-TB plus resistance to fluoroquinolones plus resistance to any other group A drug */
replace xdr = 0 if preXdr==0
replace xdr = 0 if preXdr==1 & bdqPheno==0 & lzdPheno==0
replace xdr = 1 if preXdr==1 & (bdqPheno==1 | lzdPheno==1)
label variable xdr "XDR-TB"
label define xdr 0 "not XDR" 1 "XDR"
label values xdr xdr

gen hr=. /* Hr-TB is defined as resistance to isoniazid but susceptibility to rifampicin */
replace hr = 0 if inhPheno==0
replace hr = 0 if mdr==1
replace hr = 1 if inhPheno==1 & rifXpertPdst==0
label variable hr "Hr-TB"
label define hr 0 "not Hr-TB" 1 "Hr-TB"
label values hr hr

save data_clean, replace


/************************************************** Flowchart of eligible patients and outcomes *******************************************************/
/* see powerpoint slide "Powerpoint_patient_enrolment_flowchart" */

use data_clean, clear

tab treat_hist, m /* eligible patients by treatment history */
tab treat_hist xpert, m /* Xpert results by treatment history */
bysort treat_hist: tab xpert rifPheno if (xpert!=2 & xpert!=.), m /* missing phenotypic DST for RIF by treatment history among MTB pos */
tab treat_hist rifXpertPdst, m /* combined RIF result (mixed Xpert and pheno DST) by treatment history - resistance detected by either test is classified as resistant */
tab treat_hist preXdr if rifXpertPdst==1, m /* second-line DST results among RR-TB */
tab treat_hist xdr if preXdr==1, m


/***************************************************** Assess patterns of missing data *******************************************************************/
/* The key consideration in performing multiple imputation in the context of a DRS is to investigate if the prevalence of RR-TB 
is biased due to missing data for rif DST. It is important to understand exactly how missing data for RR-TB are split across the
categories of other key variables in the dataset. Results from this analysis are given in the Excel file, "Stata_analysis_results" and will inform the multiple imputation model. */

use data_clean, clear 
count /*2,413*/

/* First drop the MTB negative by Xpert that are assumed to not be TB */
drop if xpert==2

/* Table showing patterns of missing data: 1 means complete data for this variable, 0 means some patients are missing data for this variable*/
misstable patterns age sex treat_hist rifXpertPdst, freq 

/* Assess whether key variables are associated with RR-TB. */
tab treat_hist rifXpertPdst, m
mhodds rifXpertPdst treat_hist

tab sex rifXpertPdst, m
mhodds rifXpertPdst sex
mhodds rifXpertPdst sex, by(treat_hist)

tab hiv rifXpertPdst, m
mhodds rifXpertPdst hiv
mhodds rifXpertPdst hiv, by(treat_hist)

tab residence rifXpertPdst, m
mhodds rifXpertPdst residence, c (3, 1) 
mhodds rifXpertPdst residence, c (3, 1) by(treat_hist)
mhodds rifXpertPdst residence, c (2, 1) 
mhodds rifXpertPdst residence, c (2, 1) by(treat_hist)

tab agegrp rifXpertPdst, m
mhodds rifXpertPdst agegrp, c(1,3) /* compares odds of RR-TB in agegrp 1 compared to agegrp 3 */
mhodds rifXpertPdst agegrp, c(1,3)  by(treat_hist)
mhodds rifXpertPdst agegrp, c(2,3) /* compares odds of RR-TB in agegrp 2 compared to agegrp 3 */
mhodds rifXpertPdst agegrp, c(2,3)  by(treat_hist)
mhodds rifXpertPdst agegrp, c(4,3) 
mhodds rifXpertPdst agegrp, c(4,3)  by(treat_hist)
mhodds rifXpertPdst agegrp, c(5,3) 
mhodds rifXpertPdst agegrp, c(5,3)  by(treat_hist)
mhodds rifXpertPdst agegrp, c(6,3) 
mhodds rifXpertPdst agegrp, c(6,3)  by(treat_hist)
mhodds rifXpertPdst agegrp, c(7,3) 
mhodds rifXpertPdst agegrp, c(7,3)  by(treat_hist)

/* Assess whether key variables are associated with having a missing result for RR-TB. First create binary variable to show if RR-TB result is missing. */
gen rifXpertPdst_miss=1 if rifXpertPdst==.
replace rifXpertPdst_miss=0 if rifXpertPdst!=.
tab rifXpertPdst_miss, m

mhodds rifXpertPdst_miss treat_hist
mhodds rifXpertPdst_miss sex
mhodds rifXpertPdst_miss hiv
mhodds rifXpertPdst_miss residence, c(2,1)
mhodds rifXpertPdst_miss residence, c(3,1)
mhodds rifXpertPdst_miss agegrp, c(1,3) 
mhodds rifXpertPdst_miss agegrp, c(2,3) 
mhodds rifXpertPdst_miss agegrp, c(4,3) 
mhodds rifXpertPdst_miss agegrp, c(5,3) 
mhodds rifXpertPdst_miss agegrp, c(6,3) 
mhodds rifXpertPdst_miss agegrp, c(7,3) 


/*******************************************************Sampling weights for new patients - two different types*****************************************/

/* Sampling weights are calculated based on the number of new cases. 
This is because there is only a target sample size for new cases, while previously treated cases are enrolled opportunistically.
The weights calculated for new cases are then also applied to previously treated cases, assuming that the ratio of new:previously treated cases in each strata is similar. */

*The sample size of each strata is:
*Strata 1 = 650
*Stata 2 = 760
*Strata 3 = 771


/****** WEIGHTS 1: Restricted dataset containing just patients with rif results. This is the dataset for the non-imputed analysis. ******/
use data_clean, clear
keep if rifXpertPdst!=. /* just keep patients with a rifampicin result */

bysort strata treat_hist: gen totsize=_N 
by strata (treat_hist), sort: replace totsize=totsize[1] 
gen pweight_restricted = 650/totsize if strata==1 /* you must replace 650 by the target sample size for new patients in the given strata */
replace pweight_restricted = 760/totsize if strata==2 
replace pweight_restricted = 771/totsize if strata==3 
svyset cluster [pw=pweight_restricted]
drop totsize

save data_weights_restricted.dta, replace


/****** WEIGHTS 2: Dataset including those patients without rif results. This is the dataset for imputation ******/
use data_clean, clear
drop if xpert==2 /* drop MTB negative by Xpert as these cases are assumed to not be TB */
tab xpert, m
bysort strata treat_hist: gen totsize=_N 
by strata (treat_hist), sort: replace totsize=totsize[1] 
gen pweight_imp = 650/totsize if strata==1 /* you must replace 650 by the target sample size for new patients in the given strata */
replace pweight_imp = 760/totsize if strata==2 
replace pweight_imp = 771/totsize if strata==3 
svyset cluster [pw=pweight_imp]
drop totsize

save data_weights_impute.dta, replace /* this dataset will be used for multiple imputation */



/********************************************************* RR-TB analyses without imputation********************************************************/

/* This section compares different approaches to calculating proportions and 95% CIs:
 (i)  Crude analysis (individual random sampling without clustering), with and without sampling weights to adjust for under- or over-enrolment of new cases
 (ii) Logistic regression, with and without sampling weights
 (iii) Logistic regression with sampling weights and robust standard errors to account for clustering
 (iv) Random effects model to account for clustering in both point estimate and 95% CIs

 Results are summarised in provided Excel spreadsheet, "Analysis results", and should be compared to assess the impact of different approaches */


/* to estimate design effect among new cases */
use data_weights_restricted, clear
svyset cluster
svy: logit rifXpertPdst if treat_hist==0
estat effects /* design effect is 1.7 */



/***** (i) Crude analysis ****/
use data_weights_restricted, clear 
svyset cluster [pw=pweight_restricted]
tab rifXpertPdst if treat_hist==0
*cii #total #rifXpertPdst (to calculate 95% CIs, need to copy the total number of patients and total number with RR-TB from table, and replace in line below)
cii prop 2090 88
tab rifXpertPdst if treat_hist==1
cii prop 265 91


/* Can also apply sampling weights to account for under- or over-enrolment of cases*/
svy: tab rifXpertPdst if treat_hist==0, ci /* for new cases */
svy: tab rifXpertPdst if treat_hist==1, ci /* for previously treated cases */


/* (ii) Logistic regression without taking into account clustering */
logit rifXpertPdst if treat_hist==0 /* for new cases */
di invlogit(_b[_cons])
di invlogit(_b[_cons] - invnormal(0.975)*_se[_cons])
di invlogit(_b[_cons] + invnormal(0.975)*_se[_cons])

logit rifXpertPdst if treat_hist==1 /* for previously treated cases */
di invlogit(_b[_cons])
di invlogit(_b[_cons] - invnormal(0.975)*_se[_cons])
di invlogit(_b[_cons] + invnormal(0.975)*_se[_cons])


/* Now add sampling weights to account for under- or over-enrolment */
logit rifXpertPdst if treat_hist==0 [pw=pweight_restricted] /* for new cases */
di invlogit(_b[_cons])
di invlogit(_b[_cons] - invnormal(0.975)*_se[_cons])
di invlogit(_b[_cons] + invnormal(0.975)*_se[_cons])

logit rifXpertPdst if treat_hist==1 [pw=pweight_restricted] /* for previously treated cases */
di invlogit(_b[_cons])
di invlogit(_b[_cons] - invnormal(0.975)*_se[_cons])
di invlogit(_b[_cons] + invnormal(0.975)*_se[_cons])


/* (iii) Now add robust standard errors to account for clustering */
logit rifXpertPdst if treat_hist ==0 [pw=pweight_restricted], robust cluster(cluster) /* for new cases */
di invlogit(_b[_cons])
di invlogit(_b[_cons] - invnormal(0.975)*_se[_cons])
di invlogit(_b[_cons] + invnormal(0.975)*_se[_cons])

logit rifXpertPdst if treat_hist ==1 [pw=pweight_restricted], robust cluster(cluster) /* for previously treated cases */
di invlogit(_b[_cons])
di invlogit(_b[_cons] - invnormal(0.975)*_se[_cons])
di invlogit(_b[_cons] + invnormal(0.975)*_se[_cons])


/* (iv) Random effects model */
xtlogit rifXpertPdst if treat_hist==0, re i(cluster) /* for new cases */
di invlogit(_b[_cons])
di invlogit(_b[_cons] - invnormal(0.975)*_se[_cons])
di invlogit(_b[_cons] + invnormal(0.975)*_se[_cons])

xtlogit rifXpertPdst if treat_hist==1, re i(cluster) /* for previously treated cases */
di invlogit(_b[_cons])
di invlogit(_b[_cons] - invnormal(0.975)*_se[_cons])
di invlogit(_b[_cons] + invnormal(0.975)*_se[_cons])



/********************************************************* RR-TB analyses with imputation**************************************************/

/* Try different models using variables that are either associated with RR-TB or are associated with having a missing RR-TB result (if any are identified). */

use data_weights_impute, clear
svyset cluster [pw=pweight_imp]

mi set flong
mi register imputed rifXpertPdst treat_hist sex agegrp

mi impute chained ///
(ologit, augment) agegrp (logit, augment) rifXpertPdst (logit, augment) sex (logit, augment) treat_hist, add(10) rseed(101) burnin(20) chaindots dryrun
	
mi impute chained ///
(ologit, augment) agegrp (logit, augment) rifXpertPdst (logit, augment) sex (logit, augment) treat_hist, add(10) rseed(101) burnin(20) chaindots report

gen _mj=_mi_m
gen _mi=_mi_id

save model1, replace 

/* You have now created 10 imputed datasets (add(10)), each of which was saved after the
imputation of chained equations process completed 20 cycles (burnin(20)). These ten datasets
including the original one (with the missing data) and are stacked one on top of the other and saved
in a file called model1.dta (save (model1, replace)). You first set a random number seed (seed(101)) in order to be able to replicate
exactly the same dataset at a later stage. */

/* to perform the analysis below, you need to first install an add-on module called mim by typing: ssc install mim */
use model1, clear
mi unset
svyset cluster [pw=pweight_imp]
bysort _mj: tab rifXpertPdst, m /* compare the number of imputed cases across datasets */

mim, category(fit) storebv: logit rifXpertPdst  if treat_hist==0 [pw=pweight_imp] /* for new cases */
di invlogit(_b[_cons])
di invlogit(_b[_cons] - invnormal(0.975)*_se[_cons])
di invlogit(_b[_cons] + invnormal(0.975)*_se[_cons])

mim, category(fit) storebv: logit rifXpertPdst  if treat_hist==1 [pw=pweight_imp] /* for previously treated cases */
di invlogit(_b[_cons])
di invlogit(_b[_cons] - invnormal(0.975)*_se[_cons])
di invlogit(_b[_cons] + invnormal(0.975)*_se[_cons])



/*****************************************Resistance to other drugs ***************************************************************************/

/* Any INH resistance by treatment history */
use data_weights_restricted, clear  

svyset cluster [pw=pweight_restricted] /* we are applying the same weights applied for the non-imputed analysis of RR-TB plus robust standard errors to account for clustering */
logit inhPheno if treat_hist==0 [pw=pweight_restricted], robust cluster(cluster)
di invlogit(_b[_cons])
di invlogit(_b[_cons] - invnormal(0.975)*_se[_cons])
di invlogit(_b[_cons] + invnormal(0.975)*_se[_cons])

logit inhPheno if treat_hist==1 [pw=pweight_restricted], robust cluster(cluster)
di invlogit(_b[_cons])
di invlogit(_b[_cons] - invnormal(0.975)*_se[_cons])
di invlogit(_b[_cons] + invnormal(0.975)*_se[_cons])


/* Hr-TB by treatment history */
use data_weights_restricted, clear  

svyset cluster [pw=pweight_restricted] /* we are applying the same weights applied for the non-imputed analysis of RR-TB plus robust standard errors to account for clustering */
logit hr if treat_hist==0 [pw=pweight_restricted], robust cluster(cluster)
di invlogit(_b[_cons])
di invlogit(_b[_cons] - invnormal(0.975)*_se[_cons])
di invlogit(_b[_cons] + invnormal(0.975)*_se[_cons])

logit hr if treat_hist==1 [pw=pweight_restricted], robust cluster(cluster)
di invlogit(_b[_cons])
di invlogit(_b[_cons] - invnormal(0.975)*_se[_cons])
di invlogit(_b[_cons] + invnormal(0.975)*_se[_cons])


/* MDR-TB by treatment history */
use data_weights_restricted, clear 

svyset cluster [pw=pweight_restricted] /* we are applying the same weights applied for the non-imputed analysis of RR-TB plus robust standard errors to account for clustering */
logit mdr if treat_hist==0, robust cluster(cluster)
di invlogit(_b[_cons])
di invlogit(_b[_cons] - invnormal(0.975)*_se[_cons])
di invlogit(_b[_cons] + invnormal(0.975)*_se[_cons])

logit mdr if treat_hist==1 , robust cluster(cluster)
di invlogit(_b[_cons])
di invlogit(_b[_cons] - invnormal(0.975)*_se[_cons])
di invlogit(_b[_cons] + invnormal(0.975)*_se[_cons])


/* LFX by treatment history */
use data_weights_restricted, clear 

svyset cluster [pw=pweight_restricted] /* we are applying the same weights applied for the non-imputed analysis of RR-TB plus robust standard errors to account for clustering */
logit lfxPheno if treat_hist==0, robust cluster(cluster)
di invlogit(_b[_cons])
di invlogit(_b[_cons] - invnormal(0.975)*_se[_cons])
di invlogit(_b[_cons] + invnormal(0.975)*_se[_cons])

logit lfxPheno if treat_hist==1, robust cluster(cluster)
di invlogit(_b[_cons])
di invlogit(_b[_cons] - invnormal(0.975)*_se[_cons])
di invlogit(_b[_cons] + invnormal(0.975)*_se[_cons])


/* pre-XDR-TB and XDR-TB by treatment history */
use data_weights_restricted 

svyset cluster [pw=pweight_restricted] /* we are applying the same weights applied for the non-imputed analysis of RR-TB plus robust standard errors to account for clustering */
logit preXdr if treat_hist==0, robust cluster(cluster)
di invlogit(_b[_cons])
di invlogit(_b[_cons] - invnormal(0.975)*_se[_cons])
di invlogit(_b[_cons] + invnormal(0.975)*_se[_cons])

logit preXdr if treat_hist==1, robust cluster(cluster)
di invlogit(_b[_cons])
di invlogit(_b[_cons] - invnormal(0.975)*_se[_cons])
di invlogit(_b[_cons] + invnormal(0.975)*_se[_cons])

logit xdr if treat_hist==0, robust cluster(cluster)
di invlogit(_b[_cons])
di invlogit(_b[_cons] - invnormal(0.975)*_se[_cons])
di invlogit(_b[_cons] + invnormal(0.975)*_se[_cons])

logit xdr if treat_hist==1, robust cluster(cluster)
di invlogit(_b[_cons])
di invlogit(_b[_cons] - invnormal(0.975)*_se[_cons])
di invlogit(_b[_cons] + invnormal(0.975)*_se[_cons])


/* pre-XDR-TB and XDR-TB among RR-TB */
use data_weights_restricted /* we are using the same dataset as above but not applying weights nor robust standard errors */
keep if rifXpertPdst==1 & preXdr!=. /* and restricting dataset to only RR-TB cases with a pre-XDR-TB result */

logit preXdr 
di invlogit(_b[_cons])
di invlogit(_b[_cons] - invnormal(0.975)*_se[_cons])
di invlogit(_b[_cons] + invnormal(0.975)*_se[_cons])

logit xdr
di invlogit(_b[_cons])
di invlogit(_b[_cons] - invnormal(0.975)*_se[_cons])
di invlogit(_b[_cons] + invnormal(0.975)*_se[_cons])



/**************************************************Risk factor analysis for RR-TB *************************************************************/
/* Surveys are not designed to investigate possible risk factors. With the low numbers of RR-TB found in most surveys, a significant association
is rarely found, except for a history of previous treatment. */


use data_weights_restricted, clear /* this uses the dataset containing only patients with a rifampicin result, not the imputed dataset */


/**** Univariate analysis *****/
/* you can also use the mhodds command used earlier under "Assessing patters of missing data", or as per below */

/* treatment history */
tab rifXpertPdst treat_hist, col
logistic rifXpertPdst treat_hist

/* sex */
tab rifXpertPdst sex, m col
logistic rifXpertPdst sex  

/* hiv */
tab rifXpertPdst hiv, m col
logistic hiv rifXpertPdst hiv 



/******Multivariate analysis ******/
/* Choose variables with p <0.2 in univariate analysis for inclusion in multivariate model. 
However, treatment history is likely to be the only variable significantly associated with RR-TB in most surveys.
Re-assess this relationship when adjusted for age and sex. 
To compare nested models, perform likelihood ratio test. */

xi: logistic rifXpertPdst treat_hist sex ib3.agegrp /* this specifies category 3 as the ref category for agegrp */


