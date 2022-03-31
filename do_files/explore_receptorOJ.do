**// Exploring results of Mirai Validation output
**//16/12/21
**// Olasubomi J. Omoleye
** D Huo, modified on Nov/16/2021 
** Mount to directory in huo-lab/Image
** Review accepted by Olasubomi O


cd /Volumes/huo-lab/Image/ojomoleye/projects/mirai_validation/


import delimited "/Volumes/huo-lab/Image/ojomoleye/projects/mirai_validation/data/Mirai_results_hrneg1116.csv", clear 

// count if patient_id !=study_id
drop v1 study_id

rename _year_risk one_year_risk
rename v4 two_year_risk
rename v5 three_year_risk
rename v6 four_year_risk
rename v7 five_year_risk

** [0] check the data 
codebook patient_id exam_id
tab years_to_cancer case
tab years_to_last_followup case
count if years_to_cancer== years_to_last_followup
count if years_to_cancer== years_to_last_followup & case=="True"
**graph matrix *_year_risk
** the estimated risk at years 1-5 are highly correlated. The ranking by any of them are similar: 
pwcorr *_year_risk
spearman *_year_risk 
sum *_year_risk 
loneway one_year_risk patient_id

// search SJ-6-3  snp15_6
** install somersd package; 
gen case_num=0 if case=="False" //numeric representation of case status
replace case_num=1 if case=="True"


*** [Section 1]
** DH: suggest to revise as follows because a control with 2 years follow-up could develop cancer later. She contribed the 2-year risk evaluation but not beyond 2 years. 
** Not clear how "years_to_last_followup" was rounded and defined. It seems it is the year at which the mammograph is taken? Here I assume "years_to_last_followup" means follow-up at least to the recorded year. 


duplicates tag patient_id , gen(dup)
sort patient_id years_to_last_followup exam_id
// tab dup
// list patient_id exam_id years_to_last_followup years_to_cancer if dup==9 & case=="True"
// list patient_id exam_id years_to_last_followup years_to_cancer if dup==9 & case=="False" , sepby(patient_id)
** take the last (max) year of follow-up: this is the time lapsed from enrollment (the first exam) to the last date of no cancer or cancer diagnosis. 
by patient_id: egen year_FU = max(years_to_last_followup)
// tab year_FU case 


** assume that year_FU==1 for cases means the cases have a cancer after 1 year but before 2 years. 
** year_FU==1 for controls mean that the controls have been followed up at least 1 year but < 2 years. 
** so to evaluate the risk by 1 year, all cases without cancer in the first year and all controls follwed at 1 years are "controls"
gen case_1=0 if (case=="True" & year_FU>=1 ) | (case=="False" & year_FU>=1) 
gen case_2=0 if (case=="True" & year_FU>=2 ) | (case=="False" & year_FU>=2) 
gen case_3=0 if (case=="True" & year_FU>=3 ) | (case=="False" & year_FU>=3) 
gen case_4=0 if (case=="True" & year_FU>=4 ) | (case=="False" & year_FU>=4) 
gen case_5=0 if (case=="True" & year_FU>=5 ) | (case=="False" & year_FU>=5) 
replace case_1=1 if case=="True" & year_FU < 1
replace case_2=1 if case=="True" & year_FU < 2
replace case_3=1 if case=="True" & year_FU < 3
replace case_4=1 if case=="True" & year_FU < 4
replace case_5=1 if case=="True" & year_FU < 5
// tab1 case_1 - case_5


// AUC at year 1, ignore the clustering within patient using Delong et al's method
roctab case_1 one_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen censor_1 = 1-case_1 
gen one_year_surv = 1-one_year_risk 
tab year_FU case_1
gen year_FU1 = 2 if case_1 ==0 
replace year_FU1 = 1 if case_1 ==1 
** C-index ignores the clustering: this is similar to roctab 
somersd year_FU one_year_surv, tr(c) cenind(censor_1)
** C-index allows the clustering 
somersd year_FU one_year_surv, tr(c) cluster(patient_id) cenind(censor_1)
** C-index with the 2 dummy times: case has shorter time than control. This is identical to the model just above 
tab year_FU1 case_1 
somersd year_FU1 one_year_surv, tr(c) cenind(censor_1)
 
// AUC at year 2, ignore the clustering within patient using Delong et al's method
roctab case_2 two_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen censor_2 = 1-case_2
gen two_year_surv = 1- two_year_risk 
tab year_FU case_2
gen year_FU2 = 2 if case_2 ==0 
replace year_FU2 = 1 if case_2 ==1 
** C-index ignores the clustering: this is similar to roctab ?  
somersd year_FU2 two_year_surv, tr(c) cenind(censor_2)
** C-index allows the clustering 
somersd year_FU2 two_year_surv, tr(c) cluster(patient_id) cenind(censor_2)
 
 // AUC at year 3, ignore the clustering within patient using Delong et al's method
roctab case_3 three_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen censor_3 = 1-case_3
gen three_year_surv = 1- three_year_risk 
tab year_FU case_3
gen year_FU3 = 2 if case_3 ==0 
replace year_FU3 = 1 if case_3 ==1 
** C-index ignores the clustering: this is similar to roctab ?  
somersd year_FU3 three_year_surv, tr(c) cenind(censor_3)
** C-index allows the clustering 
somersd year_FU3 three_year_surv, tr(c) cluster(patient_id) cenind(censor_3)

// AUC at year 4, ignore the clustering within patient using Delong et al's method
roctab case_4 four_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen censor_4 = 1-case_4
gen four_year_surv = 1- four_year_risk 
tab year_FU case_4
gen year_FU4 = 2 if case_4 ==0 
replace year_FU4 = 1 if case_4 ==1 
** C-index ignores the clustering: this is similar to roctab ?  
somersd year_FU4 four_year_surv, tr(c) cenind(censor_4)
** C-index allows the clustering 
somersd year_FU4 four_year_surv, tr(c) cluster(patient_id) cenind(censor_4)
 
// AUC at year 5, ignore the clustering within patient using Delong et al's method
roctab case_5 five_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen censor_5 = 1-case_5
gen five_year_surv = 1- five_year_risk 
tab year_FU case_5
gen year_FU5 = 2 if case_5 ==0 
replace year_FU5 = 1 if case_5 ==1 
** C-index ignores the clustering: this is similar to roctab ?  
somersd year_FU5 five_year_surv, tr(c) cenind(censor_5)
** C-index allows the clustering 
somersd year_FU5 five_year_surv, tr(c) cluster(patient_id) cenind(censor_5)


// C-index Using 5-year risk as the predictor of breast cancer risk over time. Using all mortality data even those beyond year 5: 
gen surv_param = 1 - five_year_risk  
gen censind = 1 - case_num 
** somersd survival-time ~ predictor: here we use all available exams, but exams long time ago from case diagnosis might have weak prediction. 
somersd year_FU surv_param , tr(c) cluster(patient_id) ce(censind) 

// C-index Using 5-year risk as the predictor of breast cancer risk over time. Using data up to year 5
gen censind5 = 1 - case_5 
** somersd survival-time ~ predictor 
somersd year_FU surv_param , tr(c) cluster(patient_id) cenind(censind5) 



**=============================================================================
** HORMONE RECEPTOR NEGATIVE COHORT
**=============================================================================
*** [Section 2]
*** DH guessed the "years_to_last_followup" means years from the current exam to the last follow-up. 
*** MIRAI's result probabily take the exams exact 1 year ago to predict 1-year risk; exams of exact 2 years to predict 2-year risk.... 
** if years_to_last_followup>=1, this exam provides info on 1-year risk, regardless case/control 
gen cancer_1=0 if (case=="True" & years_to_last_followup>=1 ) | (case=="False" & years_to_last_followup>=1) 
gen cancer_2=0 if (case=="True" & years_to_last_followup>=2 ) | (case=="False" & years_to_last_followup>=2) 
gen cancer_3=0 if (case=="True" & years_to_last_followup>=3 ) | (case=="False" & years_to_last_followup>=3) 
gen cancer_4=0 if (case=="True" & years_to_last_followup>=4 ) | (case=="False" & years_to_last_followup>=4) 
** if years_to_last_followup <5 years, this exam has no info for 5-year risk, i.e. censored
gen cancer_5=0 if (case=="True" & years_to_last_followup>=5 ) | (case=="False" & years_to_last_followup>=5)  
** if years_to_last_followup<1, this exam provides info on 1-year risk for case; 
replace cancer_1=1 if case=="True" & years_to_last_followup < 1
replace cancer_2=1 if case=="True" & years_to_last_followup < 2
** Any cancers occured within the last 3 years before the last follow-up should be included in the 3-year evaluation; 
replace cancer_3=1 if case=="True" & years_to_last_followup < 3
replace cancer_4=1 if case=="True" & years_to_last_followup < 4
replace cancer_5=1 if case=="True" & years_to_last_followup < 5
tab1 cancer_1 - cancer_5, m 


// AUC at year 1, ignore the clustering within patient using Delong et al's method
roctab cancer_1 one_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen cenind_1 = 1-cancer_1 
*gen one_year_surv = 1-one_year_risk 
tab years_to_last_followup cancer_1
gen surv_casecontrol1 = 2- cancer_1 
** C-index ignores the clustering: this is similar to roctab 
somersd surv_casecontrol1 one_year_surv, tr(c) cenind(cenind_1)
** C-index with the 2 dummy times: case has shorter time than control. The final result 
tab surv_casecontrol1 cancer_1 
somersd surv_casecontrol1 one_year_surv, tr(c) cluster(patient_id) cenind(cenind_1)
 
// AUC at year 2, ignore the clustering within patient using Delong et al's method
roctab cancer_2 two_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen cenind_2 = 1-cancer_2 
tab years_to_last_followup cancer_2
gen surv_casecontrol2 = 2- cancer_2 
** C-index ignores the clustering: this is similar to roctab 
somersd surv_casecontrol2 two_year_surv, tr(c) cenind(cenind_2)
** C-index with the 2 dummy times: case has shorter time than control. The final result 
tab surv_casecontrol2 cancer_2 
somersd surv_casecontrol2 two_year_surv, tr(c) cluster(patient_id) cenind(cenind_2)

// AUC at year 3, ignore the clustering within patient using Delong et al's method
roctab cancer_3 three_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen cenind_3 = 1-cancer_3 
tab years_to_last_followup cancer_3
gen surv_casecontrol3 = 2- cancer_3 
** C-index ignores the clustering: this is similar to roctab 
somersd surv_casecontrol3 three_year_surv, tr(c) cenind(cenind_3)
** C-index with the 2 dummy times: case has shorter time than control. The final result
tab surv_casecontrol3 cancer_3 
somersd surv_casecontrol3 three_year_surv, tr(c) cluster(patient_id) cenind(cenind_3)

// AUC at year 4, ignore the clustering within patient using Delong et al's method
roctab cancer_4 four_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen cenind_4 = 1-cancer_4 
tab years_to_last_followup cancer_4
gen surv_casecontrol4 = 2- cancer_4 
** C-index ignores the clustering: this is similar to roctab 
somersd surv_casecontrol4 four_year_surv, tr(c) cenind(cenind_4)
** C-index with the 2 dummy times: case has shorter time than control. The final result 
tab surv_casecontrol4 cancer_4 
somersd surv_casecontrol4 four_year_surv, tr(c) cluster(patient_id) cenind(cenind_4)

// AUC at year 5, ignore the clustering within patient using Delong et al's method
roctab cancer_5 five_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen cenind_5 = 1-cancer_5 
tab years_to_last_followup cancer_5
gen surv_casecontrol5 = 2- cancer_5 
** C-index ignores the clustering: this is similar to roctab 
somersd surv_casecontrol5 five_year_surv, tr(c) cenind(cenind_5)
** C-index with the 2 dummy times: case has shorter time than control. The final result 
tab surv_casecontrol5 cancer_5 
somersd surv_casecontrol5 five_year_surv, tr(c) cluster(patient_id) cenind(cenind_5)


// unclear how the Harrell's C-index from MIRAI's package is calculated. The risk estimates from the most recent exams? 
// let me try: 
somersd years_to_last_followup five_year_surv , tr(c) cluster(patient_id) ce(censind) 

// revised c-index values for full cohort
gen censind5A = 1 - cancer_5 
replace censind5A = 1 if case=="False"
somersd years_to_last_followup five_year_surv , tr(c) cluster(patient_id) ce(censind5A) 

end


**=================================================================================



import delimited "/Volumes/huo-lab/Image/ojomoleye/projects/mirai_validation/data/Mirai_results_hrpos1116.csv", clear 

// count if patient_id !=study_id
drop v1 study_id

rename _year_risk one_year_risk
rename v4 two_year_risk
rename v5 three_year_risk
rename v6 four_year_risk
rename v7 five_year_risk

** [0] check the data 
codebook patient_id exam_id
tab years_to_cancer case
tab years_to_last_followup case
count if years_to_cancer== years_to_last_followup
count if years_to_cancer== years_to_last_followup & case=="True"
**graph matrix *_year_risk
** the estimated risk at years 1-5 are highly correlated. The ranking by any of them are similar: 
pwcorr *_year_risk
spearman *_year_risk 
sum *_year_risk 
loneway one_year_risk patient_id

// search SJ-6-3  snp15_6
** install somersd package; 
gen case_num=0 if case=="False" //numeric representation of case status
replace case_num=1 if case=="True"


*** [Section 1]
** DH: suggest to revise as follows because a control with 2 years follow-up could develop cancer later. She contribed the 2-year risk evaluation but not beyond 2 years. 
** Not clear how "years_to_last_followup" was rounded and defined. It seems it is the year at which the mammograph is taken? Here I assume "years_to_last_followup" means follow-up at least to the recorded year. 


duplicates tag patient_id , gen(dup)
sort patient_id years_to_last_followup exam_id
// tab dup
// list patient_id exam_id years_to_last_followup years_to_cancer if dup==9 & case=="True"
// list patient_id exam_id years_to_last_followup years_to_cancer if dup==9 & case=="False" , sepby(patient_id)
** take the last (max) year of follow-up: this is the time lapsed from enrollment (the first exam) to the last date of no cancer or cancer diagnosis. 
by patient_id: egen year_FU = max(years_to_last_followup)
// tab year_FU case 


** assume that year_FU==1 for cases means the cases have a cancer after 1 year but before 2 years. 
** year_FU==1 for controls mean that the controls have been followed up at least 1 year but < 2 years. 
** so to evaluate the risk by 1 year, all cases without cancer in the first year and all controls follwed at 1 years are "controls"
gen case_1=0 if (case=="True" & year_FU>=1 ) | (case=="False" & year_FU>=1) 
gen case_2=0 if (case=="True" & year_FU>=2 ) | (case=="False" & year_FU>=2) 
gen case_3=0 if (case=="True" & year_FU>=3 ) | (case=="False" & year_FU>=3) 
gen case_4=0 if (case=="True" & year_FU>=4 ) | (case=="False" & year_FU>=4) 
gen case_5=0 if (case=="True" & year_FU>=5 ) | (case=="False" & year_FU>=5) 
replace case_1=1 if case=="True" & year_FU < 1
replace case_2=1 if case=="True" & year_FU < 2
replace case_3=1 if case=="True" & year_FU < 3
replace case_4=1 if case=="True" & year_FU < 4
replace case_5=1 if case=="True" & year_FU < 5
// tab1 case_1 - case_5


// AUC at year 1, ignore the clustering within patient using Delong et al's method
roctab case_1 one_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen censor_1 = 1-case_1 
gen one_year_surv = 1-one_year_risk 
tab year_FU case_1
gen year_FU1 = 2 if case_1 ==0 
replace year_FU1 = 1 if case_1 ==1 
** C-index ignores the clustering: this is similar to roctab 
somersd year_FU one_year_surv, tr(c) cenind(censor_1)
** C-index allows the clustering 
somersd year_FU one_year_surv, tr(c) cluster(patient_id) cenind(censor_1)
** C-index with the 2 dummy times: case has shorter time than control. This is identical to the model just above 
tab year_FU1 case_1 
somersd year_FU1 one_year_surv, tr(c) cenind(censor_1)
 
// AUC at year 2, ignore the clustering within patient using Delong et al's method
roctab case_2 two_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen censor_2 = 1-case_2
gen two_year_surv = 1- two_year_risk 
tab year_FU case_2
gen year_FU2 = 2 if case_2 ==0 
replace year_FU2 = 1 if case_2 ==1 
** C-index ignores the clustering: this is similar to roctab ?  
somersd year_FU2 two_year_surv, tr(c) cenind(censor_2)
** C-index allows the clustering 
somersd year_FU2 two_year_surv, tr(c) cluster(patient_id) cenind(censor_2)
 
 // AUC at year 3, ignore the clustering within patient using Delong et al's method
roctab case_3 three_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen censor_3 = 1-case_3
gen three_year_surv = 1- three_year_risk 
tab year_FU case_3
gen year_FU3 = 2 if case_3 ==0 
replace year_FU3 = 1 if case_3 ==1 
** C-index ignores the clustering: this is similar to roctab ?  
somersd year_FU3 three_year_surv, tr(c) cenind(censor_3)
** C-index allows the clustering 
somersd year_FU3 three_year_surv, tr(c) cluster(patient_id) cenind(censor_3)

// AUC at year 4, ignore the clustering within patient using Delong et al's method
roctab case_4 four_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen censor_4 = 1-case_4
gen four_year_surv = 1- four_year_risk 
tab year_FU case_4
gen year_FU4 = 2 if case_4 ==0 
replace year_FU4 = 1 if case_4 ==1 
** C-index ignores the clustering: this is similar to roctab ?  
somersd year_FU4 four_year_surv, tr(c) cenind(censor_4)
** C-index allows the clustering 
somersd year_FU4 four_year_surv, tr(c) cluster(patient_id) cenind(censor_4)
 
// AUC at year 5, ignore the clustering within patient using Delong et al's method
roctab case_5 five_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen censor_5 = 1-case_5
gen five_year_surv = 1- five_year_risk 
tab year_FU case_5
gen year_FU5 = 2 if case_5 ==0 
replace year_FU5 = 1 if case_5 ==1 
** C-index ignores the clustering: this is similar to roctab ?  
somersd year_FU5 five_year_surv, tr(c) cenind(censor_5)
** C-index allows the clustering 
somersd year_FU5 five_year_surv, tr(c) cluster(patient_id) cenind(censor_5)


// C-index Using 5-year risk as the predictor of breast cancer risk over time. Using all mortality data even those beyond year 5: 
gen surv_param = 1 - five_year_risk  
gen censind = 1 - case_num 
** somersd survival-time ~ predictor: here we use all available exams, but exams long time ago from case diagnosis might have weak prediction. 
somersd year_FU surv_param , tr(c) cluster(patient_id) ce(censind) 

// C-index Using 5-year risk as the predictor of breast cancer risk over time. Using data up to year 5
gen censind5 = 1 - case_5 
** somersd survival-time ~ predictor 
somersd year_FU surv_param , tr(c) cluster(patient_id) cenind(censind5) 



**=============================================================================
** HORMONE RECEPTOR POSITIVE COHORT
**=============================================================================
*** [Section 2]
*** DH guessed the "years_to_last_followup" means years from the current exam to the last follow-up. 
*** MIRAI's result probabily take the exams exact 1 year ago to predict 1-year risk; exams of exact 2 years to predict 2-year risk.... 
** if years_to_last_followup>=1, this exam provides info on 1-year risk, regardless case/control 
gen cancer_1=0 if (case=="True" & years_to_last_followup>=1 ) | (case=="False" & years_to_last_followup>=1) 
gen cancer_2=0 if (case=="True" & years_to_last_followup>=2 ) | (case=="False" & years_to_last_followup>=2) 
gen cancer_3=0 if (case=="True" & years_to_last_followup>=3 ) | (case=="False" & years_to_last_followup>=3) 
gen cancer_4=0 if (case=="True" & years_to_last_followup>=4 ) | (case=="False" & years_to_last_followup>=4) 
** if years_to_last_followup <5 years, this exam has no info for 5-year risk, i.e. censored
gen cancer_5=0 if (case=="True" & years_to_last_followup>=5 ) | (case=="False" & years_to_last_followup>=5)  
** if years_to_last_followup<1, this exam provides info on 1-year risk for case; 
replace cancer_1=1 if case=="True" & years_to_last_followup < 1
replace cancer_2=1 if case=="True" & years_to_last_followup < 2
** Any cancers occured within the last 3 years before the last follow-up should be included in the 3-year evaluation; 
replace cancer_3=1 if case=="True" & years_to_last_followup < 3
replace cancer_4=1 if case=="True" & years_to_last_followup < 4
replace cancer_5=1 if case=="True" & years_to_last_followup < 5
tab1 cancer_1 - cancer_5, m 


// AUC at year 1, ignore the clustering within patient using Delong et al's method
roctab cancer_1 one_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen cenind_1 = 1-cancer_1 
*gen one_year_surv = 1-one_year_risk 
tab years_to_last_followup cancer_1
gen surv_casecontrol1 = 2- cancer_1 
** C-index ignores the clustering: this is similar to roctab 
somersd surv_casecontrol1 one_year_surv, tr(c) cenind(cenind_1)
** C-index with the 2 dummy times: case has shorter time than control. The final result 
tab surv_casecontrol1 cancer_1 
somersd surv_casecontrol1 one_year_surv, tr(c) cluster(patient_id) cenind(cenind_1)
 
// AUC at year 2, ignore the clustering within patient using Delong et al's method
roctab cancer_2 two_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen cenind_2 = 1-cancer_2 
tab years_to_last_followup cancer_2
gen surv_casecontrol2 = 2- cancer_2 
** C-index ignores the clustering: this is similar to roctab 
somersd surv_casecontrol2 two_year_surv, tr(c) cenind(cenind_2)
** C-index with the 2 dummy times: case has shorter time than control. The final result 
tab surv_casecontrol2 cancer_2 
somersd surv_casecontrol2 two_year_surv, tr(c) cluster(patient_id) cenind(cenind_2)

// AUC at year 3, ignore the clustering within patient using Delong et al's method
roctab cancer_3 three_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen cenind_3 = 1-cancer_3 
tab years_to_last_followup cancer_3
gen surv_casecontrol3 = 2- cancer_3 
** C-index ignores the clustering: this is similar to roctab 
somersd surv_casecontrol3 three_year_surv, tr(c) cenind(cenind_3)
** C-index with the 2 dummy times: case has shorter time than control. The final result
tab surv_casecontrol3 cancer_3 
somersd surv_casecontrol3 three_year_surv, tr(c) cluster(patient_id) cenind(cenind_3)

// AUC at year 4, ignore the clustering within patient using Delong et al's method
roctab cancer_4 four_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen cenind_4 = 1-cancer_4 
tab years_to_last_followup cancer_4
gen surv_casecontrol4 = 2- cancer_4 
** C-index ignores the clustering: this is similar to roctab 
somersd surv_casecontrol4 four_year_surv, tr(c) cenind(cenind_4)
** C-index with the 2 dummy times: case has shorter time than control. The final result 
tab surv_casecontrol4 cancer_4 
somersd surv_casecontrol4 four_year_surv, tr(c) cluster(patient_id) cenind(cenind_4)

// AUC at year 5, ignore the clustering within patient using Delong et al's method
roctab cancer_5 five_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen cenind_5 = 1-cancer_5 
tab years_to_last_followup cancer_5
gen surv_casecontrol5 = 2- cancer_5 
** C-index ignores the clustering: this is similar to roctab 
somersd surv_casecontrol5 five_year_surv, tr(c) cenind(cenind_5)
** C-index with the 2 dummy times: case has shorter time than control. The final result 
tab surv_casecontrol5 cancer_5 
somersd surv_casecontrol5 five_year_surv, tr(c) cluster(patient_id) cenind(cenind_5)


// unclear how the Harrell's C-index from MIRAI's package is calculated. The risk estimates from the most recent exams? 
// let me try: 
somersd years_to_last_followup five_year_surv , tr(c) cluster(patient_id) ce(censind) 

// revised c-index values for full cohort
gen censind5A = 1 - cancer_5 
replace censind5A = 1 if case=="False"
somersd years_to_last_followup five_year_surv , tr(c) cluster(patient_id) ce(censind5A) 

end




**=================================================================================



import delimited "/Volumes/huo-lab/Image/ojomoleye/projects/mirai_validation/data/Mirai_results_her2pos1116.csv", clear 

// count if patient_id !=study_id
drop v1 study_id

rename _year_risk one_year_risk
rename v4 two_year_risk
rename v5 three_year_risk
rename v6 four_year_risk
rename v7 five_year_risk

** [0] check the data 
codebook patient_id exam_id
tab years_to_cancer case
tab years_to_last_followup case
count if years_to_cancer== years_to_last_followup
count if years_to_cancer== years_to_last_followup & case=="True"
**graph matrix *_year_risk
** the estimated risk at years 1-5 are highly correlated. The ranking by any of them are similar: 
pwcorr *_year_risk
spearman *_year_risk 
sum *_year_risk 
loneway one_year_risk patient_id

// search SJ-6-3  snp15_6
** install somersd package; 
gen case_num=0 if case=="False" //numeric representation of case status
replace case_num=1 if case=="True"


*** [Section 1]
** DH: suggest to revise as follows because a control with 2 years follow-up could develop cancer later. She contribed the 2-year risk evaluation but not beyond 2 years. 
** Not clear how "years_to_last_followup" was rounded and defined. It seems it is the year at which the mammograph is taken? Here I assume "years_to_last_followup" means follow-up at least to the recorded year. 


duplicates tag patient_id , gen(dup)
sort patient_id years_to_last_followup exam_id
// tab dup
// list patient_id exam_id years_to_last_followup years_to_cancer if dup==9 & case=="True"
// list patient_id exam_id years_to_last_followup years_to_cancer if dup==9 & case=="False" , sepby(patient_id)
** take the last (max) year of follow-up: this is the time lapsed from enrollment (the first exam) to the last date of no cancer or cancer diagnosis. 
by patient_id: egen year_FU = max(years_to_last_followup)
// tab year_FU case 


** assume that year_FU==1 for cases means the cases have a cancer after 1 year but before 2 years. 
** year_FU==1 for controls mean that the controls have been followed up at least 1 year but < 2 years. 
** so to evaluate the risk by 1 year, all cases without cancer in the first year and all controls follwed at 1 years are "controls"
gen case_1=0 if (case=="True" & year_FU>=1 ) | (case=="False" & year_FU>=1) 
gen case_2=0 if (case=="True" & year_FU>=2 ) | (case=="False" & year_FU>=2) 
gen case_3=0 if (case=="True" & year_FU>=3 ) | (case=="False" & year_FU>=3) 
gen case_4=0 if (case=="True" & year_FU>=4 ) | (case=="False" & year_FU>=4) 
gen case_5=0 if (case=="True" & year_FU>=5 ) | (case=="False" & year_FU>=5) 
replace case_1=1 if case=="True" & year_FU < 1
replace case_2=1 if case=="True" & year_FU < 2
replace case_3=1 if case=="True" & year_FU < 3
replace case_4=1 if case=="True" & year_FU < 4
replace case_5=1 if case=="True" & year_FU < 5
// tab1 case_1 - case_5


// AUC at year 1, ignore the clustering within patient using Delong et al's method
roctab case_1 one_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen censor_1 = 1-case_1 
gen one_year_surv = 1-one_year_risk 
tab year_FU case_1
gen year_FU1 = 2 if case_1 ==0 
replace year_FU1 = 1 if case_1 ==1 
** C-index ignores the clustering: this is similar to roctab 
somersd year_FU one_year_surv, tr(c) cenind(censor_1)
** C-index allows the clustering 
somersd year_FU one_year_surv, tr(c) cluster(patient_id) cenind(censor_1)
** C-index with the 2 dummy times: case has shorter time than control. This is identical to the model just above 
tab year_FU1 case_1 
somersd year_FU1 one_year_surv, tr(c) cenind(censor_1)
 
// AUC at year 2, ignore the clustering within patient using Delong et al's method
roctab case_2 two_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen censor_2 = 1-case_2
gen two_year_surv = 1- two_year_risk 
tab year_FU case_2
gen year_FU2 = 2 if case_2 ==0 
replace year_FU2 = 1 if case_2 ==1 
** C-index ignores the clustering: this is similar to roctab ?  
somersd year_FU2 two_year_surv, tr(c) cenind(censor_2)
** C-index allows the clustering 
somersd year_FU2 two_year_surv, tr(c) cluster(patient_id) cenind(censor_2)
 
 // AUC at year 3, ignore the clustering within patient using Delong et al's method
roctab case_3 three_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen censor_3 = 1-case_3
gen three_year_surv = 1- three_year_risk 
tab year_FU case_3
gen year_FU3 = 2 if case_3 ==0 
replace year_FU3 = 1 if case_3 ==1 
** C-index ignores the clustering: this is similar to roctab ?  
somersd year_FU3 three_year_surv, tr(c) cenind(censor_3)
** C-index allows the clustering 
somersd year_FU3 three_year_surv, tr(c) cluster(patient_id) cenind(censor_3)

// AUC at year 4, ignore the clustering within patient using Delong et al's method
roctab case_4 four_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen censor_4 = 1-case_4
gen four_year_surv = 1- four_year_risk 
tab year_FU case_4
gen year_FU4 = 2 if case_4 ==0 
replace year_FU4 = 1 if case_4 ==1 
** C-index ignores the clustering: this is similar to roctab ?  
somersd year_FU4 four_year_surv, tr(c) cenind(censor_4)
** C-index allows the clustering 
somersd year_FU4 four_year_surv, tr(c) cluster(patient_id) cenind(censor_4)
 
// AUC at year 5, ignore the clustering within patient using Delong et al's method
roctab case_5 five_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen censor_5 = 1-case_5
gen five_year_surv = 1- five_year_risk 
tab year_FU case_5
gen year_FU5 = 2 if case_5 ==0 
replace year_FU5 = 1 if case_5 ==1 
** C-index ignores the clustering: this is similar to roctab ?  
somersd year_FU5 five_year_surv, tr(c) cenind(censor_5)
** C-index allows the clustering 
somersd year_FU5 five_year_surv, tr(c) cluster(patient_id) cenind(censor_5)


// C-index Using 5-year risk as the predictor of breast cancer risk over time. Using all mortality data even those beyond year 5: 
gen surv_param = 1 - five_year_risk  
gen censind = 1 - case_num 
** somersd survival-time ~ predictor: here we use all available exams, but exams long time ago from case diagnosis might have weak prediction. 
somersd year_FU surv_param , tr(c) cluster(patient_id) ce(censind) 

// C-index Using 5-year risk as the predictor of breast cancer risk over time. Using data up to year 5
gen censind5 = 1 - case_5 
** somersd survival-time ~ predictor 
somersd year_FU surv_param , tr(c) cluster(patient_id) cenind(censind5) 



**=============================================================================
** HER 2 RECEPTOR POSITIVE COHORT
**=============================================================================
*** [Section 2]
*** DH guessed the "years_to_last_followup" means years from the current exam to the last follow-up. 
*** MIRAI's result probabily take the exams exact 1 year ago to predict 1-year risk; exams of exact 2 years to predict 2-year risk.... 
** if years_to_last_followup>=1, this exam provides info on 1-year risk, regardless case/control 
gen cancer_1=0 if (case=="True" & years_to_last_followup>=1 ) | (case=="False" & years_to_last_followup>=1) 
gen cancer_2=0 if (case=="True" & years_to_last_followup>=2 ) | (case=="False" & years_to_last_followup>=2) 
gen cancer_3=0 if (case=="True" & years_to_last_followup>=3 ) | (case=="False" & years_to_last_followup>=3) 
gen cancer_4=0 if (case=="True" & years_to_last_followup>=4 ) | (case=="False" & years_to_last_followup>=4) 
** if years_to_last_followup <5 years, this exam has no info for 5-year risk, i.e. censored
gen cancer_5=0 if (case=="True" & years_to_last_followup>=5 ) | (case=="False" & years_to_last_followup>=5)  
** if years_to_last_followup<1, this exam provides info on 1-year risk for case; 
replace cancer_1=1 if case=="True" & years_to_last_followup < 1
replace cancer_2=1 if case=="True" & years_to_last_followup < 2
** Any cancers occured within the last 3 years before the last follow-up should be included in the 3-year evaluation; 
replace cancer_3=1 if case=="True" & years_to_last_followup < 3
replace cancer_4=1 if case=="True" & years_to_last_followup < 4
replace cancer_5=1 if case=="True" & years_to_last_followup < 5
tab1 cancer_1 - cancer_5, m 


// AUC at year 1, ignore the clustering within patient using Delong et al's method
roctab cancer_1 one_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen cenind_1 = 1-cancer_1 
*gen one_year_surv = 1-one_year_risk 
tab years_to_last_followup cancer_1
gen surv_casecontrol1 = 2- cancer_1 
** C-index ignores the clustering: this is similar to roctab 
somersd surv_casecontrol1 one_year_surv, tr(c) cenind(cenind_1)
** C-index with the 2 dummy times: case has shorter time than control. The final result 
tab surv_casecontrol1 cancer_1 
somersd surv_casecontrol1 one_year_surv, tr(c) cluster(patient_id) cenind(cenind_1)
 
// AUC at year 2, ignore the clustering within patient using Delong et al's method
roctab cancer_2 two_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen cenind_2 = 1-cancer_2 
tab years_to_last_followup cancer_2
gen surv_casecontrol2 = 2- cancer_2 
** C-index ignores the clustering: this is similar to roctab 
somersd surv_casecontrol2 two_year_surv, tr(c) cenind(cenind_2)
** C-index with the 2 dummy times: case has shorter time than control. The final result 
tab surv_casecontrol2 cancer_2 
somersd surv_casecontrol2 two_year_surv, tr(c) cluster(patient_id) cenind(cenind_2)

// AUC at year 3, ignore the clustering within patient using Delong et al's method
roctab cancer_3 three_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen cenind_3 = 1-cancer_3 
tab years_to_last_followup cancer_3
gen surv_casecontrol3 = 2- cancer_3 
** C-index ignores the clustering: this is similar to roctab 
somersd surv_casecontrol3 three_year_surv, tr(c) cenind(cenind_3)
** C-index with the 2 dummy times: case has shorter time than control. The final result
tab surv_casecontrol3 cancer_3 
somersd surv_casecontrol3 three_year_surv, tr(c) cluster(patient_id) cenind(cenind_3)

// AUC at year 4, ignore the clustering within patient using Delong et al's method
roctab cancer_4 four_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen cenind_4 = 1-cancer_4 
tab years_to_last_followup cancer_4
gen surv_casecontrol4 = 2- cancer_4 
** C-index ignores the clustering: this is similar to roctab 
somersd surv_casecontrol4 four_year_surv, tr(c) cenind(cenind_4)
** C-index with the 2 dummy times: case has shorter time than control. The final result 
tab surv_casecontrol4 cancer_4 
somersd surv_casecontrol4 four_year_surv, tr(c) cluster(patient_id) cenind(cenind_4)

// AUC at year 5, ignore the clustering within patient using Delong et al's method
roctab cancer_5 five_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen cenind_5 = 1-cancer_5 
tab years_to_last_followup cancer_5
gen surv_casecontrol5 = 2- cancer_5 
** C-index ignores the clustering: this is similar to roctab 
somersd surv_casecontrol5 five_year_surv, tr(c) cenind(cenind_5)
** C-index with the 2 dummy times: case has shorter time than control. The final result 
tab surv_casecontrol5 cancer_5 
somersd surv_casecontrol5 five_year_surv, tr(c) cluster(patient_id) cenind(cenind_5)


// unclear how the Harrell's C-index from MIRAI's package is calculated. The risk estimates from the most recent exams? 
// let me try: 
somersd years_to_last_followup five_year_surv , tr(c) cluster(patient_id) ce(censind) 

// revised c-index values for full cohort
gen censind5A = 1 - cancer_5 
replace censind5A = 1 if case=="False"
somersd years_to_last_followup five_year_surv , tr(c) cluster(patient_id) ce(censind5A) 

end


**=================================================================================



import delimited "/Volumes/huo-lab/Image/ojomoleye/projects/mirai_validation/data/Mirai_results_her2neg1116.csv", clear 

// count if patient_id !=study_id
drop v1 study_id

rename _year_risk one_year_risk
rename v4 two_year_risk
rename v5 three_year_risk
rename v6 four_year_risk
rename v7 five_year_risk

** [0] check the data 
codebook patient_id exam_id
tab years_to_cancer case
tab years_to_last_followup case
count if years_to_cancer== years_to_last_followup
count if years_to_cancer== years_to_last_followup & case=="True"
**graph matrix *_year_risk
** the estimated risk at years 1-5 are highly correlated. The ranking by any of them are similar: 
pwcorr *_year_risk
spearman *_year_risk 
sum *_year_risk 
loneway one_year_risk patient_id

// search SJ-6-3  snp15_6
** install somersd package; 
gen case_num=0 if case=="False" //numeric representation of case status
replace case_num=1 if case=="True"


*** [Section 1]
** DH: suggest to revise as follows because a control with 2 years follow-up could develop cancer later. She contribed the 2-year risk evaluation but not beyond 2 years. 
** Not clear how "years_to_last_followup" was rounded and defined. It seems it is the year at which the mammograph is taken? Here I assume "years_to_last_followup" means follow-up at least to the recorded year. 


duplicates tag patient_id , gen(dup)
sort patient_id years_to_last_followup exam_id
// tab dup
// list patient_id exam_id years_to_last_followup years_to_cancer if dup==9 & case=="True"
// list patient_id exam_id years_to_last_followup years_to_cancer if dup==9 & case=="False" , sepby(patient_id)
** take the last (max) year of follow-up: this is the time lapsed from enrollment (the first exam) to the last date of no cancer or cancer diagnosis. 
by patient_id: egen year_FU = max(years_to_last_followup)
// tab year_FU case 


** assume that year_FU==1 for cases means the cases have a cancer after 1 year but before 2 years. 
** year_FU==1 for controls mean that the controls have been followed up at least 1 year but < 2 years. 
** so to evaluate the risk by 1 year, all cases without cancer in the first year and all controls follwed at 1 years are "controls"
gen case_1=0 if (case=="True" & year_FU>=1 ) | (case=="False" & year_FU>=1) 
gen case_2=0 if (case=="True" & year_FU>=2 ) | (case=="False" & year_FU>=2) 
gen case_3=0 if (case=="True" & year_FU>=3 ) | (case=="False" & year_FU>=3) 
gen case_4=0 if (case=="True" & year_FU>=4 ) | (case=="False" & year_FU>=4) 
gen case_5=0 if (case=="True" & year_FU>=5 ) | (case=="False" & year_FU>=5) 
replace case_1=1 if case=="True" & year_FU < 1
replace case_2=1 if case=="True" & year_FU < 2
replace case_3=1 if case=="True" & year_FU < 3
replace case_4=1 if case=="True" & year_FU < 4
replace case_5=1 if case=="True" & year_FU < 5
// tab1 case_1 - case_5


// AUC at year 1, ignore the clustering within patient using Delong et al's method
roctab case_1 one_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen censor_1 = 1-case_1 
gen one_year_surv = 1-one_year_risk 
tab year_FU case_1
gen year_FU1 = 2 if case_1 ==0 
replace year_FU1 = 1 if case_1 ==1 
** C-index ignores the clustering: this is similar to roctab 
somersd year_FU one_year_surv, tr(c) cenind(censor_1)
** C-index allows the clustering 
somersd year_FU one_year_surv, tr(c) cluster(patient_id) cenind(censor_1)
** C-index with the 2 dummy times: case has shorter time than control. This is identical to the model just above 
tab year_FU1 case_1 
somersd year_FU1 one_year_surv, tr(c) cenind(censor_1)
 
// AUC at year 2, ignore the clustering within patient using Delong et al's method
roctab case_2 two_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen censor_2 = 1-case_2
gen two_year_surv = 1- two_year_risk 
tab year_FU case_2
gen year_FU2 = 2 if case_2 ==0 
replace year_FU2 = 1 if case_2 ==1 
** C-index ignores the clustering: this is similar to roctab ?  
somersd year_FU2 two_year_surv, tr(c) cenind(censor_2)
** C-index allows the clustering 
somersd year_FU2 two_year_surv, tr(c) cluster(patient_id) cenind(censor_2)
 
 // AUC at year 3, ignore the clustering within patient using Delong et al's method
roctab case_3 three_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen censor_3 = 1-case_3
gen three_year_surv = 1- three_year_risk 
tab year_FU case_3
gen year_FU3 = 2 if case_3 ==0 
replace year_FU3 = 1 if case_3 ==1 
** C-index ignores the clustering: this is similar to roctab ?  
somersd year_FU3 three_year_surv, tr(c) cenind(censor_3)
** C-index allows the clustering 
somersd year_FU3 three_year_surv, tr(c) cluster(patient_id) cenind(censor_3)

// AUC at year 4, ignore the clustering within patient using Delong et al's method
roctab case_4 four_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen censor_4 = 1-case_4
gen four_year_surv = 1- four_year_risk 
tab year_FU case_4
gen year_FU4 = 2 if case_4 ==0 
replace year_FU4 = 1 if case_4 ==1 
** C-index ignores the clustering: this is similar to roctab ?  
somersd year_FU4 four_year_surv, tr(c) cenind(censor_4)
** C-index allows the clustering 
somersd year_FU4 four_year_surv, tr(c) cluster(patient_id) cenind(censor_4)
 
// AUC at year 5, ignore the clustering within patient using Delong et al's method
roctab case_5 five_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen censor_5 = 1-case_5
gen five_year_surv = 1- five_year_risk 
tab year_FU case_5
gen year_FU5 = 2 if case_5 ==0 
replace year_FU5 = 1 if case_5 ==1 
** C-index ignores the clustering: this is similar to roctab ?  
somersd year_FU5 five_year_surv, tr(c) cenind(censor_5)
** C-index allows the clustering 
somersd year_FU5 five_year_surv, tr(c) cluster(patient_id) cenind(censor_5)


// C-index Using 5-year risk as the predictor of breast cancer risk over time. Using all mortality data even those beyond year 5: 
gen surv_param = 1 - five_year_risk  
gen censind = 1 - case_num 
** somersd survival-time ~ predictor: here we use all available exams, but exams long time ago from case diagnosis might have weak prediction. 
somersd year_FU surv_param , tr(c) cluster(patient_id) ce(censind) 

// C-index Using 5-year risk as the predictor of breast cancer risk over time. Using data up to year 5
gen censind5 = 1 - case_5 
** somersd survival-time ~ predictor 
somersd year_FU surv_param , tr(c) cluster(patient_id) cenind(censind5) 



**=============================================================================
** HER 2 RECEPTOR NEGATIVE COHORT
**=============================================================================
*** [Section 2]
*** DH guessed the "years_to_last_followup" means years from the current exam to the last follow-up. 
*** MIRAI's result probabily take the exams exact 1 year ago to predict 1-year risk; exams of exact 2 years to predict 2-year risk.... 
** if years_to_last_followup>=1, this exam provides info on 1-year risk, regardless case/control 
gen cancer_1=0 if (case=="True" & years_to_last_followup>=1 ) | (case=="False" & years_to_last_followup>=1) 
gen cancer_2=0 if (case=="True" & years_to_last_followup>=2 ) | (case=="False" & years_to_last_followup>=2) 
gen cancer_3=0 if (case=="True" & years_to_last_followup>=3 ) | (case=="False" & years_to_last_followup>=3) 
gen cancer_4=0 if (case=="True" & years_to_last_followup>=4 ) | (case=="False" & years_to_last_followup>=4) 
** if years_to_last_followup <5 years, this exam has no info for 5-year risk, i.e. censored
gen cancer_5=0 if (case=="True" & years_to_last_followup>=5 ) | (case=="False" & years_to_last_followup>=5)  
** if years_to_last_followup<1, this exam provides info on 1-year risk for case; 
replace cancer_1=1 if case=="True" & years_to_last_followup < 1
replace cancer_2=1 if case=="True" & years_to_last_followup < 2
** Any cancers occured within the last 3 years before the last follow-up should be included in the 3-year evaluation; 
replace cancer_3=1 if case=="True" & years_to_last_followup < 3
replace cancer_4=1 if case=="True" & years_to_last_followup < 4
replace cancer_5=1 if case=="True" & years_to_last_followup < 5
tab1 cancer_1 - cancer_5, m 


// AUC at year 1, ignore the clustering within patient using Delong et al's method
roctab cancer_1 one_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen cenind_1 = 1-cancer_1 
*gen one_year_surv = 1-one_year_risk 
tab years_to_last_followup cancer_1
gen surv_casecontrol1 = 2- cancer_1 
** C-index ignores the clustering: this is similar to roctab 
somersd surv_casecontrol1 one_year_surv, tr(c) cenind(cenind_1)
** C-index with the 2 dummy times: case has shorter time than control. The final result 
tab surv_casecontrol1 cancer_1 
somersd surv_casecontrol1 one_year_surv, tr(c) cluster(patient_id) cenind(cenind_1)
 
// AUC at year 2, ignore the clustering within patient using Delong et al's method
roctab cancer_2 two_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen cenind_2 = 1-cancer_2 
tab years_to_last_followup cancer_2
gen surv_casecontrol2 = 2- cancer_2 
** C-index ignores the clustering: this is similar to roctab 
somersd surv_casecontrol2 two_year_surv, tr(c) cenind(cenind_2)
** C-index with the 2 dummy times: case has shorter time than control. The final result 
tab surv_casecontrol2 cancer_2 
somersd surv_casecontrol2 two_year_surv, tr(c) cluster(patient_id) cenind(cenind_2)

// AUC at year 3, ignore the clustering within patient using Delong et al's method
roctab cancer_3 three_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen cenind_3 = 1-cancer_3 
tab years_to_last_followup cancer_3
gen surv_casecontrol3 = 2- cancer_3 
** C-index ignores the clustering: this is similar to roctab 
somersd surv_casecontrol3 three_year_surv, tr(c) cenind(cenind_3)
** C-index with the 2 dummy times: case has shorter time than control. The final result
tab surv_casecontrol3 cancer_3 
somersd surv_casecontrol3 three_year_surv, tr(c) cluster(patient_id) cenind(cenind_3)

// AUC at year 4, ignore the clustering within patient using Delong et al's method
roctab cancer_4 four_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen cenind_4 = 1-cancer_4 
tab years_to_last_followup cancer_4
gen surv_casecontrol4 = 2- cancer_4 
** C-index ignores the clustering: this is similar to roctab 
somersd surv_casecontrol4 four_year_surv, tr(c) cenind(cenind_4)
** C-index with the 2 dummy times: case has shorter time than control. The final result 
tab surv_casecontrol4 cancer_4 
somersd surv_casecontrol4 four_year_surv, tr(c) cluster(patient_id) cenind(cenind_4)

// AUC at year 5, ignore the clustering within patient using Delong et al's method
roctab cancer_5 five_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen cenind_5 = 1-cancer_5 
tab years_to_last_followup cancer_5
gen surv_casecontrol5 = 2- cancer_5 
** C-index ignores the clustering: this is similar to roctab 
somersd surv_casecontrol5 five_year_surv, tr(c) cenind(cenind_5)
** C-index with the 2 dummy times: case has shorter time than control. The final result 
tab surv_casecontrol5 cancer_5 
somersd surv_casecontrol5 five_year_surv, tr(c) cluster(patient_id) cenind(cenind_5)


// unclear how the Harrell's C-index from MIRAI's package is calculated. The risk estimates from the most recent exams? 
// let me try: 
somersd years_to_last_followup five_year_surv , tr(c) cluster(patient_id) ce(censind) 

// revised c-index values for full cohort
gen censind5A = 1 - cancer_5 
replace censind5A = 1 if case=="False"
somersd years_to_last_followup five_year_surv , tr(c) cluster(patient_id) ce(censind5A) 

end




import delimited "/Volumes/huo-lab/Image/ojomoleye/projects/mirai_validation/data/Mirai_results_tripleneg1116.csv", clear 

// count if patient_id !=study_id
drop v1 study_id

rename _year_risk one_year_risk
rename v4 two_year_risk
rename v5 three_year_risk
rename v6 four_year_risk
rename v7 five_year_risk

** [0] check the data 
codebook patient_id exam_id
tab years_to_cancer case
tab years_to_last_followup case
count if years_to_cancer== years_to_last_followup
count if years_to_cancer== years_to_last_followup & case=="True"
**graph matrix *_year_risk
** the estimated risk at years 1-5 are highly correlated. The ranking by any of them are similar: 
pwcorr *_year_risk
spearman *_year_risk 
sum *_year_risk 
loneway one_year_risk patient_id

// search SJ-6-3  snp15_6
** install somersd package; 
gen case_num=0 if case=="False" //numeric representation of case status
replace case_num=1 if case=="True"


*** [Section 1]
** DH: suggest to revise as follows because a control with 2 years follow-up could develop cancer later. She contribed the 2-year risk evaluation but not beyond 2 years. 
** Not clear how "years_to_last_followup" was rounded and defined. It seems it is the year at which the mammograph is taken? Here I assume "years_to_last_followup" means follow-up at least to the recorded year. 


duplicates tag patient_id , gen(dup)
sort patient_id years_to_last_followup exam_id
// tab dup
// list patient_id exam_id years_to_last_followup years_to_cancer if dup==9 & case=="True"
// list patient_id exam_id years_to_last_followup years_to_cancer if dup==9 & case=="False" , sepby(patient_id)
** take the last (max) year of follow-up: this is the time lapsed from enrollment (the first exam) to the last date of no cancer or cancer diagnosis. 
by patient_id: egen year_FU = max(years_to_last_followup)
// tab year_FU case 


** assume that year_FU==1 for cases means the cases have a cancer after 1 year but before 2 years. 
** year_FU==1 for controls mean that the controls have been followed up at least 1 year but < 2 years. 
** so to evaluate the risk by 1 year, all cases without cancer in the first year and all controls follwed at 1 years are "controls"
gen case_1=0 if (case=="True" & year_FU>=1 ) | (case=="False" & year_FU>=1) 
gen case_2=0 if (case=="True" & year_FU>=2 ) | (case=="False" & year_FU>=2) 
gen case_3=0 if (case=="True" & year_FU>=3 ) | (case=="False" & year_FU>=3) 
gen case_4=0 if (case=="True" & year_FU>=4 ) | (case=="False" & year_FU>=4) 
gen case_5=0 if (case=="True" & year_FU>=5 ) | (case=="False" & year_FU>=5) 
replace case_1=1 if case=="True" & year_FU < 1
replace case_2=1 if case=="True" & year_FU < 2
replace case_3=1 if case=="True" & year_FU < 3
replace case_4=1 if case=="True" & year_FU < 4
replace case_5=1 if case=="True" & year_FU < 5
// tab1 case_1 - case_5


// AUC at year 1, ignore the clustering within patient using Delong et al's method
roctab case_1 one_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen censor_1 = 1-case_1 
gen one_year_surv = 1-one_year_risk 
tab year_FU case_1
gen year_FU1 = 2 if case_1 ==0 
replace year_FU1 = 1 if case_1 ==1 
** C-index ignores the clustering: this is similar to roctab 
somersd year_FU one_year_surv, tr(c) cenind(censor_1)
** C-index allows the clustering 
somersd year_FU one_year_surv, tr(c) cluster(patient_id) cenind(censor_1)
** C-index with the 2 dummy times: case has shorter time than control. This is identical to the model just above 
tab year_FU1 case_1 
somersd year_FU1 one_year_surv, tr(c) cenind(censor_1)
 
// AUC at year 2, ignore the clustering within patient using Delong et al's method
roctab case_2 two_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen censor_2 = 1-case_2
gen two_year_surv = 1- two_year_risk 
tab year_FU case_2
gen year_FU2 = 2 if case_2 ==0 
replace year_FU2 = 1 if case_2 ==1 
** C-index ignores the clustering: this is similar to roctab ?  
somersd year_FU2 two_year_surv, tr(c) cenind(censor_2)
** C-index allows the clustering 
somersd year_FU2 two_year_surv, tr(c) cluster(patient_id) cenind(censor_2)
 
 // AUC at year 3, ignore the clustering within patient using Delong et al's method
roctab case_3 three_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen censor_3 = 1-case_3
gen three_year_surv = 1- three_year_risk 
tab year_FU case_3
gen year_FU3 = 2 if case_3 ==0 
replace year_FU3 = 1 if case_3 ==1 
** C-index ignores the clustering: this is similar to roctab ?  
somersd year_FU3 three_year_surv, tr(c) cenind(censor_3)
** C-index allows the clustering 
somersd year_FU3 three_year_surv, tr(c) cluster(patient_id) cenind(censor_3)

// AUC at year 4, ignore the clustering within patient using Delong et al's method
roctab case_4 four_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen censor_4 = 1-case_4
gen four_year_surv = 1- four_year_risk 
tab year_FU case_4
gen year_FU4 = 2 if case_4 ==0 
replace year_FU4 = 1 if case_4 ==1 
** C-index ignores the clustering: this is similar to roctab ?  
somersd year_FU4 four_year_surv, tr(c) cenind(censor_4)
** C-index allows the clustering 
somersd year_FU4 four_year_surv, tr(c) cluster(patient_id) cenind(censor_4)
 
// AUC at year 5, ignore the clustering within patient using Delong et al's method
roctab case_5 five_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen censor_5 = 1-case_5
gen five_year_surv = 1- five_year_risk 
tab year_FU case_5
gen year_FU5 = 2 if case_5 ==0 
replace year_FU5 = 1 if case_5 ==1 
** C-index ignores the clustering: this is similar to roctab ?  
somersd year_FU5 five_year_surv, tr(c) cenind(censor_5)
** C-index allows the clustering 
somersd year_FU5 five_year_surv, tr(c) cluster(patient_id) cenind(censor_5)


// C-index Using 5-year risk as the predictor of breast cancer risk over time. Using all mortality data even those beyond year 5: 
gen surv_param = 1 - five_year_risk  
gen censind = 1 - case_num 
** somersd survival-time ~ predictor: here we use all available exams, but exams long time ago from case diagnosis might have weak prediction. 
somersd year_FU surv_param , tr(c) cluster(patient_id) ce(censind) 

// C-index Using 5-year risk as the predictor of breast cancer risk over time. Using data up to year 5
gen censind5 = 1 - case_5 
** somersd survival-time ~ predictor 
somersd year_FU surv_param , tr(c) cluster(patient_id) cenind(censind5) 



**=============================================================================
** TRIPLE NEGATIVE COHORT
**=============================================================================
*** [Section 2]
*** DH guessed the "years_to_last_followup" means years from the current exam to the last follow-up. 
*** MIRAI's result probabily take the exams exact 1 year ago to predict 1-year risk; exams of exact 2 years to predict 2-year risk.... 
** if years_to_last_followup>=1, this exam provides info on 1-year risk, regardless case/control 
gen cancer_1=0 if (case=="True" & years_to_last_followup>=1 ) | (case=="False" & years_to_last_followup>=1) 
gen cancer_2=0 if (case=="True" & years_to_last_followup>=2 ) | (case=="False" & years_to_last_followup>=2) 
gen cancer_3=0 if (case=="True" & years_to_last_followup>=3 ) | (case=="False" & years_to_last_followup>=3) 
gen cancer_4=0 if (case=="True" & years_to_last_followup>=4 ) | (case=="False" & years_to_last_followup>=4) 
** if years_to_last_followup <5 years, this exam has no info for 5-year risk, i.e. censored
gen cancer_5=0 if (case=="True" & years_to_last_followup>=5 ) | (case=="False" & years_to_last_followup>=5)  
** if years_to_last_followup<1, this exam provides info on 1-year risk for case; 
replace cancer_1=1 if case=="True" & years_to_last_followup < 1
replace cancer_2=1 if case=="True" & years_to_last_followup < 2
** Any cancers occured within the last 3 years before the last follow-up should be included in the 3-year evaluation; 
replace cancer_3=1 if case=="True" & years_to_last_followup < 3
replace cancer_4=1 if case=="True" & years_to_last_followup < 4
replace cancer_5=1 if case=="True" & years_to_last_followup < 5
tab1 cancer_1 - cancer_5, m 


// AUC at year 1, ignore the clustering within patient using Delong et al's method
roctab cancer_1 one_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen cenind_1 = 1-cancer_1 
*gen one_year_surv = 1-one_year_risk 
tab years_to_last_followup cancer_1
gen surv_casecontrol1 = 2- cancer_1 
** C-index ignores the clustering: this is similar to roctab 
somersd surv_casecontrol1 one_year_surv, tr(c) cenind(cenind_1)
** C-index with the 2 dummy times: case has shorter time than control. The final result 
tab surv_casecontrol1 cancer_1 
somersd surv_casecontrol1 one_year_surv, tr(c) cluster(patient_id) cenind(cenind_1)
 
// AUC at year 2, ignore the clustering within patient using Delong et al's method
roctab cancer_2 two_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen cenind_2 = 1-cancer_2 
tab years_to_last_followup cancer_2
gen surv_casecontrol2 = 2- cancer_2 
** C-index ignores the clustering: this is similar to roctab 
somersd surv_casecontrol2 two_year_surv, tr(c) cenind(cenind_2)
** C-index with the 2 dummy times: case has shorter time than control. The final result 
tab surv_casecontrol2 cancer_2 
somersd surv_casecontrol2 two_year_surv, tr(c) cluster(patient_id) cenind(cenind_2)

// AUC at year 3, ignore the clustering within patient using Delong et al's method
roctab cancer_3 three_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen cenind_3 = 1-cancer_3 
tab years_to_last_followup cancer_3
gen surv_casecontrol3 = 2- cancer_3 
** C-index ignores the clustering: this is similar to roctab 
somersd surv_casecontrol3 three_year_surv, tr(c) cenind(cenind_3)
** C-index with the 2 dummy times: case has shorter time than control. The final result
tab surv_casecontrol3 cancer_3 
somersd surv_casecontrol3 three_year_surv, tr(c) cluster(patient_id) cenind(cenind_3)

// AUC at year 4, ignore the clustering within patient using Delong et al's method
roctab cancer_4 four_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen cenind_4 = 1-cancer_4 
tab years_to_last_followup cancer_4
gen surv_casecontrol4 = 2- cancer_4 
** C-index ignores the clustering: this is similar to roctab 
somersd surv_casecontrol4 four_year_surv, tr(c) cenind(cenind_4)
** C-index with the 2 dummy times: case has shorter time than control. The final result 
tab surv_casecontrol4 cancer_4 
somersd surv_casecontrol4 four_year_surv, tr(c) cluster(patient_id) cenind(cenind_4)

// AUC at year 5, ignore the clustering within patient using Delong et al's method
roctab cancer_5 five_year_risk
// trick the somersd to get C-index as a proxy of AUC. 
gen cenind_5 = 1-cancer_5 
tab years_to_last_followup cancer_5
gen surv_casecontrol5 = 2- cancer_5 
** C-index ignores the clustering: this is similar to roctab 
somersd surv_casecontrol5 five_year_surv, tr(c) cenind(cenind_5)
** C-index with the 2 dummy times: case has shorter time than control. The final result 
tab surv_casecontrol5 cancer_5 
somersd surv_casecontrol5 five_year_surv, tr(c) cluster(patient_id) cenind(cenind_5)


// unclear how the Harrell's C-index from MIRAI's package is calculated. The risk estimates from the most recent exams? 
// let me try: 
somersd years_to_last_followup five_year_surv , tr(c) cluster(patient_id) ce(censind) 

// revised c-index values for full cohort
gen censind5A = 1 - cancer_5 
replace censind5A = 1 if case=="False"
somersd years_to_last_followup five_year_surv , tr(c) cluster(patient_id) ce(censind5A) 

end












