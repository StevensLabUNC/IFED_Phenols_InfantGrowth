/*************************************************************************************************************
  Programmer: Danielle Stevens
  Purpose of syntax: Clean and examine IFED data
  Last updated on: 11/17/2026
**************************************************************************************************************/

libname ifed "J:\StevensLab\IFED Study\Original Data";
options nofmterr;


/*********************************************/
/*Read long back into SAS after imputing in R*/
/*This is so we can run the LME models       */
/*********************************************/
proc import out=ifed_long
datafile="J:\StevensLab\IFED Study\IFED Phenols & Growth_DRS\Analysis\Derived data\Imputed dataset with outcomes long.csv"
replace;
run;
/*Get rid of the original data*/
data ifed_long; set ifed_long;
where _imp ne 0;
run;
proc freq data=ifed_long; tables ifedid*_imp / missing out=out; run;
proc freq data=ifed_long; tables ifedid / missing out=out; run;

/*Log-transform visit-specific exposures*/
%let exposure_list1 = bpa_obrien_v68 bp_3_obrien_v68 dcp_24_obrien_v68 dcp_25_obrien_v68 m_pb_obrien_v68 p_pb_obrien_v68 tcs_obrien_v68 bpa_obrien_v12 bp_3_obrien_v12 dcp_24_obrien_v12 dcp_25_obrien_v12 m_pb_obrien_v12 p_pb_obrien_v12 tcs_obrien_v12;
%macro loopexp;
%local i next_exposure;
%let i=1;
%do %while (%scan(&exposure_list1,&i) ne );
%let next_exposure = %scan(&exposure_list1,&i);

data ifed_long; set ifed_long;
LN_&next_exposure = log(&next_exposure);
run;

%let i = %eval(&i + 1);
%end;
%mend;
%loopexp;


/*Macro setup to combine results from type 3 tests*/
%macro combchi(df=,ddf=,chi=);
/******************************************************************************
Version 1.0, 4-27-00

This is a BETA version. USE AT YOUR OWN RISK!. Report any errors, problems
or questions to the e-mail address below. The latest version of this macro
can be found at the web site below:

Author:  Paul D. Allison, University of Pennsylvania
         allison@ssc.upenn.edu
         http://www.ssc.upenn.edu/~allison/

MACRO COMBCHI combines chi-square statistics from an analysis of several data
sets created by multiple imputation, using the method described on p. 115 of

Schafer, J.L. (1997) Analysis of Incomplete Multivariate Data.  London: 
   Chapman and Hall. 

The chi-square statistics can be either Wald statistics or likelihood ratio
statistics.  All that's needed are the several chi-square values and the
segrees of freedom. 

COMBCHI requires the installation of IML.

Example of usage:

Suppose a 3 d.f. test on four data sets produces chi-squares of 5.8, 7.2,
6.1 and 8.5.  Submit the statement:
%combchi(df=3, chi=5.8 7.2 6.1 8.5)

The following output is printed

                 F        DF       DDF

         2.1201613         3 342.22381


                      P

                   0.0974097

The macro calculates an F-statistic of 2.12 with 3 and 342 degrees of freedom.  
The associated p-value is .097.

*****************************************************************************/
proc iml;
  df=&df;
  g2={&chi};
  m=ncol(g2);
  g=sqrt(g2);
  mg2=sum(g2)/m;
  r=(1+1/m)*(ssq(g)-(sum(g)**2)/m)/(m-1);
  f=(mg2/df - r*(m-1)/(m+1))/(1+r);
  ddf=&ddf;
  p=1-probf(f,df,ddf);
  print f df ddf;
  print p;
run;
%mend combchi;

/*******************************/
/*Spline models for the outcome*/
/*******************************/
proc sort data=ifed_long; by _imp ifedid ageinweeks; run;


/*Loop to run in overall dataframe*/ 
%macro loopmodels(exposure, iqr, ddf, int);

/*IQR-standardize*/
data ifed_long; set ifed_long;
iqr_&exposure = &exposure / &iqr;
run;

/*glimmix with cubic b-spline & random intercepts: No interaction*/
proc GLIMMIX data=ifed_long;
 class ifedid FeedingGroup mage meduc childsgender crace gestationalage smkpreg alcpregany;
 effect sp_age = spline(ageinweeks / basis=bspline knotmethod=list(3.42857, 7.85714, 15.8571, 24.2857) degree=3);
 model bmiz =  iqr_&exposure sp_age 
               FeedingGroup mage meduc childsgender crace gestationalage smkpreg alcpregany/ SOLUTION;
 random intercept / SUBJECT=ifedid S VCORR G V TYPE=UN; 
 by _imp;
 estimate 'bmiz' iqr_&exposure 1 / cl;
 ods output estimates=est;
run; 
quit;
proc sort data=est; by label _imp; run;
ods output parameterestimates=est2;
proc mianalyze data=est;
	by label;
	modeleffects Estimate;
    stderr StdErr;
    run;
data noint_&exposure; retain label exposure estimate lclmean uclmean betaCI1 betaCI2 Probt; set est2; 
exposure = "&exposure";
betaCI1 = TRIM(LEFT(PUT(Estimate,9.1)))||" ("||TRIM(LEFT(PUT(LCLMean,9.1)))||", "||TRIM(LEFT(PUT(UCLMean,9.1)))||")";
betaCI2 = TRIM(LEFT(PUT(Estimate,9.2)))||" ("||TRIM(LEFT(PUT(LCLMean,9.2)))||", "||TRIM(LEFT(PUT(UCLMean,9.2)))||")";
keep label exposure estimate lclmean uclmean betaCI1 betaCI2 Probt;
run;


/*glimmix with cubic b-spline & random intercepts: Interaction*/
proc GLIMMIX data=ifed_long;
 class ifedid FeedingGroup mage meduc childsgender crace gestationalage smkpreg alcpregany;
 effect sp_age = spline(ageinweeks / basis=bspline knotmethod=list(3.42857, 7.85714, 15.8571, 24.2857) degree=3);
 model bmiz =  iqr_&exposure sp_age iqr_&exposure*sp_age 
               FeedingGroup mage meduc childsgender crace gestationalage smkpreg alcpregany/ SOLUTION;
 random intercept / SUBJECT=ifedid S VCORR G V TYPE=UN; 
 by _imp;
 estimate 'bmiz at 3 months' iqr_&exposure 1 iqr_&exposure*sp_age [0,1 12][1,1 12] / cl;
 estimate 'bmiz at 6 months' iqr_&exposure 1 iqr_&exposure*sp_age [0,1 24][1,1 24] / cl;
 estimate 'bmiz at 9 months' iqr_&exposure 1 iqr_&exposure*sp_age [0,1 36][1,1 36] / cl;
 ods output estimates=est Tests3=Tests3;
run; 
quit;

proc sort data=est; by label _imp; run;
ods output parameterestimates=est2;
proc mianalyze data=est;
	by label;
	modeleffects Estimate;
    stderr StdErr;
    run;
data int_&exposure; retain label exposure estimate lclmean uclmean betaCI1 betaCI2 Probt; set est2; 
exposure = "&exposure";
betaCI1 = TRIM(LEFT(PUT(Estimate,9.1)))||" ("||TRIM(LEFT(PUT(LCLMean,9.1)))||", "||TRIM(LEFT(PUT(UCLMean,9.1)))||")";
betaCI2 = TRIM(LEFT(PUT(Estimate,9.2)))||" ("||TRIM(LEFT(PUT(LCLMean,9.2)))||", "||TRIM(LEFT(PUT(UCLMean,9.2)))||")";
keep label exposure estimate lclmean uclmean betaCI1 betaCI2 Probt;
run;

%let varlist =;
data _null_; set tests3; where effect = &int; call symput('varlist', trim(resolve('&varlist'))||" "||trim(FValue)); run;
%put &varlist;
%combchi(df=7, ddf=&ddf, chi=&varlist)
%mend;
%loopmodels(bpa_obrien_gm_ln, 0.92, 1751, "iqr_bpa_obrien_gm_ln*sp_age"); /*p-int=0.9848364*/
%loopmodels(bp_3_obrien_gm_ln, 1.67, 1759, "iqr_bp_3_obrien_gm_ln*sp_age"); /*p-int=0.9888482*/
%loopmodels(dcp_24_obrien_gm_ln, 1.09,1751, "iqr_dcp_24_obrien_gm_ln*sp_age");  /*p-int=0.9799976*/
%loopmodels(dcp_25_obrien_gm_ln, 2.08, 1739, "iqr_dcp_25_obrien_gm_ln*sp_age");  /*p-int=0.4524687*/
%loopmodels(m_pb_obrien_gm_ln, 2.37, 1784, "iqr_m_pb_obrien_gm_ln*sp_age"); /*p-int=0.979603*/
%loopmodels(p_pb_obrien_gm_ln, 2.66, 1773, "iqr_p_pb_obrien_gm_ln*sp_age"); /*p-int=0.9848791*/
%loopmodels(tcs_obrien_gm_ln, 1.35, 1645, "iqr_tcs_obrien_gm_ln*sp_age"); /*p-int=0.965*/

/*No interaction terms*/
data beta_estimates; merge noint_dcp_24_obrien_gm_ln noint_bpa_obrien_gm_ln noint_bp_3_obrien_gm_ln noint_dcp_25_obrien_gm_ln noint_m_pb_obrien_gm_ln noint_p_pb_obrien_gm_ln noint_tcs_obrien_gm_ln; by exposure; run;
proc transpose data=beta_estimates out=flip_beta_estimates;
id label;
var betaCI2;
by exposure;
run;
ods csv file="J:\StevensLab\IFED Study\IFED Phenols & Growth_DRS\Analysis\Results\BMI z-score estimates.csv";
proc print data=beta_estimates; run;
ods csv close;
ods csv file="J:\StevensLab\IFED Study\IFED Phenols & Growth_DRS\Analysis\Results\Table of BMI z-score Estimates.csv";
proc print data=flip_beta_estimates; run;
ods csv close;

/*Interaction terms*/
data beta_estimates; merge int_dcp_24_obrien_gm_ln int_bpa_obrien_gm_ln int_bp_3_obrien_gm_ln int_dcp_25_obrien_gm_ln int_m_pb_obrien_gm_ln int_p_pb_obrien_gm_ln int_tcs_obrien_gm_ln; by exposure; run;
proc transpose data=beta_estimates out=flip_beta_estimates;
id label;
var betaCI2;
by exposure;
run;
ods csv file="J:\StevensLab\IFED Study\IFED Phenols & Growth_DRS\Analysis\Results\BMI z-score estimates by age.csv";
proc print data=beta_estimates; run;
ods csv close;
ods csv file="J:\StevensLab\IFED Study\IFED Phenols & Growth_DRS\Analysis\Results\Table of BMI z-score Estimates by Age.csv";
proc print data=flip_beta_estimates; run;
ods csv close;



/*Run unadjusted models*/ 
%macro loopmodels(exposure);
proc GLIMMIX data=ifed_long;
 class ifedid FeedingGroup mage meduc childsgender crace gestationalage smkpreg alcpregany;
 effect sp_age = spline(ageinweeks / basis=bspline knotmethod=list(3.42857, 7.85714, 15.8571, 24.2857) degree=3);
 model bmiz =  iqr_&exposure sp_age /solution;
 random intercept / SUBJECT=ifedid S VCORR G V TYPE=UN; 
 by _imp;
 estimate 'bmiz' iqr_&exposure 1 / cl;
 ods output estimates=est;
run; 
quit;
proc sort data=est; by label _imp; run;
ods output parameterestimates=est2;
proc mianalyze data=est;
	by label;
	modeleffects Estimate;
    stderr StdErr;
    run;
data &exposure; retain label exposure estimate lclmean uclmean betaCI1 betaCI2 Probt; set est2; 
exposure = "&exposure";
betaCI1 = TRIM(LEFT(PUT(Estimate,9.1)))||" ("||TRIM(LEFT(PUT(LCLMean,9.1)))||", "||TRIM(LEFT(PUT(UCLMean,9.1)))||")";
betaCI2 = TRIM(LEFT(PUT(Estimate,9.2)))||" ("||TRIM(LEFT(PUT(LCLMean,9.2)))||", "||TRIM(LEFT(PUT(UCLMean,9.2)))||")";
keep label exposure estimate lclmean uclmean betaCI1 betaCI2 Probt;
run;
%mend;
%loopmodels(bpa_obrien_gm_ln);
%loopmodels(bp_3_obrien_gm_ln); 
%loopmodels(dcp_24_obrien_gm_ln); 
%loopmodels(dcp_25_obrien_gm_ln); 
%loopmodels(m_pb_obrien_gm_ln);
%loopmodels(p_pb_obrien_gm_ln);
%loopmodels(tcs_obrien_gm_ln);

data beta_estimates; merge dcp_24_obrien_gm_ln bpa_obrien_gm_ln bp_3_obrien_gm_ln dcp_25_obrien_gm_ln m_pb_obrien_gm_ln p_pb_obrien_gm_ln tcs_obrien_gm_ln; by exposure; run;
proc transpose data=beta_estimates out=flip_beta_estimates;
id label;
var betaCI2;
by exposure;
run;
ods csv file="J:\StevensLab\IFED Study\IFED Phenols & Growth_DRS\Analysis\Results\Unadjusted BMI z-score estimates.csv";
proc print data=beta_estimates; run;
ods csv close;
ods csv file="J:\StevensLab\IFED Study\IFED Phenols & Growth_DRS\Analysis\Results\Table of Unadjusted BMI z-score Estimates.csv";
proc print data=flip_beta_estimates; run;
ods csv close;



/**************************/
/*By sex assigned at birth*/
/**************************/
%macro loopmodels(exposure, gender);
proc GLIMMIX data=ifed_long;
 class ifedid FeedingGroup mage meduc crace gestationalage smkpreg alcpregany;
 effect sp_age = spline(ageinweeks / basis=bspline knotmethod=list(3.42857, 7.85714, 15.8571, 24.2857) degree=3);
 model bmiz =  iqr_&exposure sp_age childsgender childsgender*iqr_&exposure
               FeedingGroup mage meduc crace gestationalage smkpreg alcpregany/ SOLUTION;
 random intercept / SUBJECT=ifedid S VCORR G V TYPE=UN; 
 by _imp;
 ods output parameterestimates=parms;
run; 
quit;
proc sort data=parms; by effect _imp; run;
ods output parameterestimates=parms2;
proc mianalyze data=parms;
	by effect;
	modeleffects Estimate;
    stderr StdErr;
	where effect contains "iqr";
    run;
data parms_&exposure; retain effect exposure Probt; set parms2; 
exposure = "&exposure";
keep effect exposure Probt;
where effect contains "child";
run;


proc GLIMMIX data=ifed_long;
 class ifedid FeedingGroup mage meduc childsgender crace gestationalage smkpreg alcpregany;
 effect sp_age = spline(ageinweeks / basis=bspline knotmethod=list(3.42857, 7.85714, 15.8571, 24.2857) degree=3);
 model bmiz =  iqr_&exposure sp_age
               FeedingGroup mage meduc crace gestationalage smkpreg alcpregany/ SOLUTION;
 random intercept / SUBJECT=ifedid S VCORR G V TYPE=UN; 
 by _imp;
 estimate 'bmiz' iqr_&exposure 1 / cl;
 ods output estimates=est;
 where childsgender = &gender;
run; 
quit;

proc sort data=est; by label _imp; run;
ods output parameterestimates=est2;
proc mianalyze data=est;
	by label;
	modeleffects Estimate;
    stderr StdErr;
    run;

data &exposure; retain label gender exposure estimate lclmean uclmean betaCI1 betaCI2 Probt; set est2; 
exposure = "&exposure";
gender = &gender;
betaCI1 = TRIM(LEFT(PUT(Estimate,9.1)))||" ("||TRIM(LEFT(PUT(LCLMean,9.1)))||", "||TRIM(LEFT(PUT(UCLMean,9.1)))||")";
betaCI2 = TRIM(LEFT(PUT(Estimate,9.2)))||" ("||TRIM(LEFT(PUT(LCLMean,9.2)))||", "||TRIM(LEFT(PUT(UCLMean,9.2)))||")";
keep label gender exposure estimate lclmean uclmean betaCI1 betaCI2 Probt;
run;
%mend;
%loopmodels(bpa_obrien_gm_ln, 2);
%loopmodels(bp_3_obrien_gm_ln, 2); 
%loopmodels(dcp_24_obrien_gm_ln, 2); 
%loopmodels(dcp_25_obrien_gm_ln, 2); 
%loopmodels(m_pb_obrien_gm_ln, 2);
%loopmodels(p_pb_obrien_gm_ln, 2);
%loopmodels(tcs_obrien_gm_ln, 2);
data male_beta_estimates; merge dcp_24_obrien_gm_ln bpa_obrien_gm_ln bp_3_obrien_gm_ln dcp_25_obrien_gm_ln m_pb_obrien_gm_ln p_pb_obrien_gm_ln tcs_obrien_gm_ln; by exposure; run;
data int_estimates; merge parms_dcp_24_obrien_gm_ln parms_bpa_obrien_gm_ln parms_bp_3_obrien_gm_ln parms_dcp_25_obrien_gm_ln parms_m_pb_obrien_gm_ln parms_p_pb_obrien_gm_ln parms_tcs_obrien_gm_ln; by exposure; run;
ods csv file="J:\StevensLab\IFED Study\IFED Phenols & Growth_DRS\Analysis\Results\pvals from sex stratified analyses.csv";
proc print data=int_estimates; run;
ods csv close;

%loopmodels(bpa_obrien_gm_ln, 1);
%loopmodels(bp_3_obrien_gm_ln, 1); 
%loopmodels(dcp_24_obrien_gm_ln, 1); 
%loopmodels(dcp_25_obrien_gm_ln, 1); 
%loopmodels(m_pb_obrien_gm_ln, 1);
%loopmodels(p_pb_obrien_gm_ln, 1);
%loopmodels(tcs_obrien_gm_ln, 1);
data female_beta_estimates; merge dcp_24_obrien_gm_ln bpa_obrien_gm_ln bp_3_obrien_gm_ln dcp_25_obrien_gm_ln m_pb_obrien_gm_ln p_pb_obrien_gm_ln tcs_obrien_gm_ln; by exposure; run;

data beta_estimates; merge male_beta_estimates female_beta_estimates; by gender; run;
proc transpose data=beta_estimates out=flip_beta_estimates;
id label;
var betaCI2;
by gender exposure;
run;

ods csv file="J:\StevensLab\IFED Study\IFED Phenols & Growth_DRS\Analysis\Results\BMI z-score estimates by sex.csv";
proc print data=beta_estimates; run;
ods csv close;
ods csv file="J:\StevensLab\IFED Study\IFED Phenols & Growth_DRS\Analysis\Results\Table of BMI z-score Estimates by Sex.csv";
proc print data=flip_beta_estimates; run;
ods csv close;
/*NOTE: ESTIMATES ARE THE SAME AFTER MODIFYING KNOT POINTS FOR MALES TO BETTER RESEMBLE THE AGE DISTRIBUTION IN THAT SUBSET*/

/*Sensitivity analysis in females restricted to before 28 weeks*/
%macro loopmodels(exposure, gender);
data newdat; set ifed_long;
if childsgender = 1 AND ageinweeks > 28 then DELETE;
run; 
proc GLIMMIX data=newdat;
 class ifedid FeedingGroup mage meduc crace gestationalage smkpreg alcpregany;
 effect sp_age = spline(ageinweeks / basis=bspline knotmethod=list(3.42857, 7.85714, 15.8571, 24.2857) degree=3);
 model bmiz =  iqr_&exposure sp_age childsgender childsgender*iqr_&exposure
               FeedingGroup mage meduc crace gestationalage smkpreg alcpregany/ SOLUTION;
 random intercept / SUBJECT=ifedid S VCORR G V TYPE=UN; 
 by _imp;
 ods output parameterestimates=parms;
run; 
quit;
proc sort data=parms; by effect _imp; run;
ods output parameterestimates=parms2;
proc mianalyze data=parms;
	by effect;
	modeleffects Estimate;
    stderr StdErr;
	where effect contains "iqr";
    run;
data parms_&exposure; retain effect exposure Probt; set parms2; 
exposure = "&exposure";
keep effect exposure Probt;
where effect contains "child";
run;


proc GLIMMIX data=ifed_long;
 class ifedid FeedingGroup mage meduc childsgender crace gestationalage smkpreg alcpregany;
 effect sp_age = spline(ageinweeks / basis=bspline knotmethod=list(3.42857, 7.85714, 15.8571, 24.2857) degree=3);
 model bmiz =  iqr_&exposure sp_age
               FeedingGroup mage meduc crace gestationalage smkpreg alcpregany/ SOLUTION;
 random intercept / SUBJECT=ifedid S VCORR G V TYPE=UN; 
 by _imp;
 estimate 'bmiz' iqr_&exposure 1 / cl;
 ods output estimates=est;
 where ageinweeks LE 28 AND childsgender = &gender;
run; 
quit;

proc sort data=est; by label _imp; run;
ods output parameterestimates=est2;
proc mianalyze data=est;
	by label;
	modeleffects Estimate;
    stderr StdErr;
    run;

data &exposure; retain label gender exposure estimate lclmean uclmean betaCI1 betaCI2 Probt; set est2; 
exposure = "&exposure";
gender = &gender;
betaCI1 = TRIM(LEFT(PUT(Estimate,9.1)))||" ("||TRIM(LEFT(PUT(LCLMean,9.1)))||", "||TRIM(LEFT(PUT(UCLMean,9.1)))||")";
betaCI2 = TRIM(LEFT(PUT(Estimate,9.2)))||" ("||TRIM(LEFT(PUT(LCLMean,9.2)))||", "||TRIM(LEFT(PUT(UCLMean,9.2)))||")";
keep label gender exposure estimate lclmean uclmean betaCI1 betaCI2 Probt;
run;
%mend;
%loopmodels(bpa_obrien_gm_ln, 1);
%loopmodels(bp_3_obrien_gm_ln, 1); 
%loopmodels(dcp_24_obrien_gm_ln, 1); 
%loopmodels(dcp_25_obrien_gm_ln, 1); 
%loopmodels(m_pb_obrien_gm_ln, 1);
%loopmodels(p_pb_obrien_gm_ln, 1);
%loopmodels(tcs_obrien_gm_ln, 1);
data int_estimates; merge parms_dcp_24_obrien_gm_ln parms_bpa_obrien_gm_ln parms_bp_3_obrien_gm_ln parms_dcp_25_obrien_gm_ln parms_m_pb_obrien_gm_ln parms_p_pb_obrien_gm_ln parms_tcs_obrien_gm_ln; by exposure; run;
ods csv file="J:\StevensLab\IFED Study\IFED Phenols & Growth_DRS\Analysis\Results\pvals from sex stratified analyses before 28 weeks.csv";
proc print data=int_estimates; run;
ods csv close;

data female_beta_estimates; merge dcp_24_obrien_gm_ln bpa_obrien_gm_ln bp_3_obrien_gm_ln dcp_25_obrien_gm_ln m_pb_obrien_gm_ln p_pb_obrien_gm_ln tcs_obrien_gm_ln; by exposure; run;
proc transpose data=female_beta_estimates out=flip_beta_estimates;
id label;
var betaCI2;
by exposure;
run;

ods csv file="J:\StevensLab\IFED Study\IFED Phenols & Growth_DRS\Analysis\Results\BMI z-score estimates for females before 28 weeks.csv";
proc print data=beta_estimates; run;
ods csv close;
ods csv file="J:\StevensLab\IFED Study\IFED Phenols & Growth_DRS\Analysis\Results\Table of BMI z-score Estimates for Females Before 28 Weeks.csv";
proc print data=flip_beta_estimates; run;
ods csv close;



/****************/
/*By infant diet*/
/****************/
%macro loopmodels(exposure, feedinggroup);
proc GLIMMIX data=ifed_long;
 class ifedid childsgender mage meduc crace gestationalage smkpreg alcpregany;
 effect sp_age = spline(ageinweeks / basis=bspline knotmethod=list(3.42857, 7.85714, 15.8571, 24.2857) degree=3);
 model bmiz =  iqr_&exposure sp_age FeedingGroup feedinggroup*iqr_&exposure
               childsgender mage meduc crace gestationalage smkpreg alcpregany/ SOLUTION;
 random intercept / SUBJECT=ifedid S VCORR G V TYPE=UN; 
 by _imp;
 ods output parameterestimates=parms;
run; 
quit;
proc sort data=parms; by effect _imp; run;
ods output parameterestimates=parms2;
proc mianalyze data=parms;
	by effect;
	modeleffects Estimate;
    stderr StdErr;
	where effect contains "iqr";
    run;
data parms_&exposure; retain effect exposure Probt; set parms2; 
exposure = "&exposure";
keep effect exposure Probt;
where effect contains "eedin";
run;


proc GLIMMIX data=ifed_long;
 class ifedid FeedingGroup mage meduc childsgender crace gestationalage smkpreg alcpregany;
 effect sp_age = spline(ageinweeks / basis=bspline knotmethod=list(3.42857, 7.85714, 15.8571, 24.2857) degree=3);
 model bmiz =  iqr_&exposure sp_age
               mage meduc childsgender crace gestationalage smkpreg alcpregany/ SOLUTION;
 random intercept / SUBJECT=ifedid S VCORR G V TYPE=UN; 
 by _imp;
 estimate 'bmiz' iqr_&exposure 1/ cl;
 ods output estimates=est;
 where feedinggroup = &feedinggroup;
run; 
quit;

proc sort data=est; by label _imp; run;
ods output parameterestimates=est2;
proc mianalyze data=est;
	by label;
	modeleffects Estimate;
    stderr StdErr;
    run;

data &exposure; retain label feedinggroup exposure estimate lclmean uclmean betaCI1 betaCI2 Probt; set est2; 
exposure = "&exposure";
feedinggroup = &feedinggroup;
betaCI1 = TRIM(LEFT(PUT(Estimate,9.1)))||" ("||TRIM(LEFT(PUT(LCLMean,9.1)))||", "||TRIM(LEFT(PUT(UCLMean,9.1)))||")";
betaCI2 = TRIM(LEFT(PUT(Estimate,9.2)))||" ("||TRIM(LEFT(PUT(LCLMean,9.2)))||", "||TRIM(LEFT(PUT(UCLMean,9.2)))||")";
keep label feedinggroup exposure estimate lclmean uclmean betaCI1 betaCI2 Probt;
run;
%mend;
%loopmodels(bpa_obrien_gm_ln, 2);
%loopmodels(bp_3_obrien_gm_ln, 2); 
%loopmodels(dcp_24_obrien_gm_ln, 2); 
%loopmodels(dcp_25_obrien_gm_ln, 2); 
%loopmodels(m_pb_obrien_gm_ln, 2);
%loopmodels(p_pb_obrien_gm_ln, 2);
%loopmodels(tcs_obrien_gm_ln, 2);
data cowsmilk_beta_estimates; merge dcp_24_obrien_gm_ln bpa_obrien_gm_ln bp_3_obrien_gm_ln dcp_25_obrien_gm_ln m_pb_obrien_gm_ln p_pb_obrien_gm_ln tcs_obrien_gm_ln; by exposure; run;
data int_estimates; merge parms_dcp_24_obrien_gm_ln parms_bpa_obrien_gm_ln parms_bp_3_obrien_gm_ln parms_dcp_25_obrien_gm_ln parms_m_pb_obrien_gm_ln parms_p_pb_obrien_gm_ln parms_tcs_obrien_gm_ln; by exposure; run;
ods csv file="J:\StevensLab\IFED Study\IFED Phenols & Growth_DRS\Analysis\Results\pvals from diet stratified analyses.csv";
proc print data=int_estimates; run;
ods csv close;


%loopmodels(bpa_obrien_gm_ln, 1);
%loopmodels(bp_3_obrien_gm_ln, 1); 
%loopmodels(dcp_24_obrien_gm_ln, 1); 
%loopmodels(dcp_25_obrien_gm_ln, 1); 
%loopmodels(m_pb_obrien_gm_ln, 1);
%loopmodels(p_pb_obrien_gm_ln, 1);
%loopmodels(tcs_obrien_gm_ln, 1);
data breastmilk_beta_estimates; merge dcp_24_obrien_gm_ln bpa_obrien_gm_ln bp_3_obrien_gm_ln dcp_25_obrien_gm_ln m_pb_obrien_gm_ln p_pb_obrien_gm_ln tcs_obrien_gm_ln; by exposure; run;

data beta_estimates; merge cowsmilk_beta_estimates breastmilk_beta_estimates; by feedinggroup; run;
proc transpose data=beta_estimates out=flip_beta_estimates;
id label;
var betaCI2;
by feedinggroup exposure;
run;

ods csv file="J:\StevensLab\IFED Study\IFED Phenols & Growth_DRS\Analysis\Results\BMI z-score estimates by diet.csv";
proc print data=beta_estimates; run;
ods csv close;
ods csv file="J:\StevensLab\IFED Study\IFED Phenols & Growth_DRS\Analysis\Results\Table of BMI z-score Estimates by Diet.csv";
proc print data=flip_beta_estimates; run;
ods csv close;




/*********************************/
/*Run for weight-for-age z-score */
/*********************************/

%macro loopmodels(exposure, iqr);
proc GLIMMIX data=ifed_long;
 class ifedid FeedingGroup mage meduc childsgender crace gestationalage smkpreg alcpregany;
 effect sp_age = spline(ageinweeks / basis=bspline knotmethod=list(3.42857, 7.85714, 15.8571, 24.2857) degree=3);
 model waz =  iqr_&exposure sp_age 
               FeedingGroup mage meduc childsgender crace gestationalage smkpreg alcpregany/ SOLUTION;
 random intercept / SUBJECT=ifedid S VCORR G V TYPE=UN; 
 by _imp;
 estimate 'waz' iqr_&exposure 1 / cl;
 ods output estimates=est;
run; 
quit;

proc sort data=est; by label _imp; run;
ods output parameterestimates=est2;
proc mianalyze data=est;
	by label;
	modeleffects Estimate;
    stderr StdErr;
    run;

data &exposure; retain label exposure estimate lclmean uclmean betaCI1 betaCI2 Probt; set est2; 
exposure = "&exposure";
betaCI1 = TRIM(LEFT(PUT(Estimate,9.1)))||" ("||TRIM(LEFT(PUT(LCLMean,9.1)))||", "||TRIM(LEFT(PUT(UCLMean,9.1)))||")";
betaCI2 = TRIM(LEFT(PUT(Estimate,9.2)))||" ("||TRIM(LEFT(PUT(LCLMean,9.2)))||", "||TRIM(LEFT(PUT(UCLMean,9.2)))||")";
keep label exposure estimate lclmean uclmean betaCI1 betaCI2 Probt;
run;
%mend;
%loopmodels(bpa_obrien_gm_ln); 
%loopmodels(bp_3_obrien_gm_ln);
%loopmodels(dcp_24_obrien_gm_ln);
%loopmodels(dcp_25_obrien_gm_ln);
%loopmodels(m_pb_obrien_gm_ln);
%loopmodels(p_pb_obrien_gm_ln);
%loopmodels(tcs_obrien_gm_ln);

data beta_estimates; merge dcp_24_obrien_gm_ln bpa_obrien_gm_ln bp_3_obrien_gm_ln dcp_25_obrien_gm_ln m_pb_obrien_gm_ln p_pb_obrien_gm_ln tcs_obrien_gm_ln; by exposure; run;
proc transpose data=beta_estimates out=flip_beta_estimates;
id label;
var betaCI2;
by exposure;
run;

ods csv file="J:\StevensLab\IFED Study\IFED Phenols & Growth_DRS\Analysis\Results\Weight-for-Age Z-score estimates for waz.csv";
proc print data=beta_estimates; run;
ods csv close;
ods csv file="J:\StevensLab\IFED Study\IFED Phenols & Growth_DRS\Analysis\Results\Table of Weight-for-Age Z-score Estimates.csv";
proc print data=flip_beta_estimates; run;
ods csv close;





/*****************************************/
/*Run for BMI z-score after 12 weeks only*/
/*****************************************/
data ifed_long_12plus; set ifed_long;
where ageinweeks GE 12;
run;

%macro loopmodels(exposure, iqr);
proc GLIMMIX data=ifed_long_12plus;
 class ifedid FeedingGroup mage meduc childsgender crace gestationalage smkpreg alcpregany;
 effect sp_age = spline(ageinweeks / basis=bspline knotmethod=list(15.11429, 19.51429, 24, 28.42857) degree=3);
 model bmiz =  iqr_&exposure sp_age 
               FeedingGroup mage meduc childsgender crace gestationalage smkpreg alcpregany/ SOLUTION;
 random intercept / SUBJECT=ifedid S VCORR G V TYPE=UN; 
 by _imp;
 estimate 'bmiz' iqr_&exposure 1 / cl;
 ods output estimates=est;
run; 
quit;

proc sort data=est; by label _imp; run;
ods output parameterestimates=est2;
proc mianalyze data=est;
	by label;
	modeleffects Estimate;
    stderr StdErr;
    run;

data &exposure; retain label exposure estimate lclmean uclmean betaCI1 betaCI2 Probt; set est2; 
exposure = "&exposure";
betaCI1 = TRIM(LEFT(PUT(Estimate,9.1)))||" ("||TRIM(LEFT(PUT(LCLMean,9.1)))||", "||TRIM(LEFT(PUT(UCLMean,9.1)))||")";
betaCI2 = TRIM(LEFT(PUT(Estimate,9.2)))||" ("||TRIM(LEFT(PUT(LCLMean,9.2)))||", "||TRIM(LEFT(PUT(UCLMean,9.2)))||")";
keep label exposure estimate lclmean uclmean betaCI1 betaCI2 Probt;
run;
%mend;
%loopmodels(bpa_obrien_gm_ln); 
%loopmodels(bp_3_obrien_gm_ln);
%loopmodels(dcp_24_obrien_gm_ln);
%loopmodels(dcp_25_obrien_gm_ln);
%loopmodels(m_pb_obrien_gm_ln);
%loopmodels(p_pb_obrien_gm_ln);
%loopmodels(tcs_obrien_gm_ln);


data beta_estimates; merge dcp_24_obrien_gm_ln bpa_obrien_gm_ln bp_3_obrien_gm_ln dcp_25_obrien_gm_ln m_pb_obrien_gm_ln p_pb_obrien_gm_ln tcs_obrien_gm_ln; by exposure; run;
proc transpose data=beta_estimates out=flip_beta_estimates;
id label;
var betaCI2;
by exposure;
run;

ods csv file="J:\StevensLab\IFED Study\IFED Phenols & Growth_DRS\Analysis\Results\BMI z-score estimates after 12 weeks.csv";
proc print data=beta_estimates; run;
ods csv close;
ods csv file="J:\StevensLab\IFED Study\IFED Phenols & Growth_DRS\Analysis\Results\Table of BMI z-score Estimates for BMI after 12 weeks.csv";
proc print data=flip_beta_estimates; run;
ods csv close;




/***********************************/
/*Exposure at 6-8 and then 12 weeks*/
/***********************************/
proc sort data=ifed_long; by _imp ifedid ageinweeks; run;

/*Loop to run in overall dataframe*/ 
%macro loopmodels(exposure, iqr);

/*IQR-standardize*/
data ifed_long; set ifed_long;
iqr_&exposure = &exposure / &iqr;
run;

/*glimmix with cubic b-spline & random intercepts: No interaction*/
proc GLIMMIX data=ifed_long;
 class ifedid FeedingGroup mage meduc childsgender crace gestationalage smkpreg alcpregany;
 effect sp_age = spline(ageinweeks / basis=bspline knotmethod=list(3.42857, 7.85714, 15.8571, 24.2857) degree=3);
 model bmiz =  iqr_&exposure sp_age 
               FeedingGroup mage meduc childsgender crace gestationalage smkpreg alcpregany/ SOLUTION;
 random intercept / SUBJECT=ifedid S VCORR G V TYPE=UN; 
 by _imp;
 estimate 'bmiz' iqr_&exposure 1 / cl;
 ods output estimates=est;
run; 
quit;
proc sort data=est; by label _imp; run;
ods output parameterestimates=est2;
proc mianalyze data=est;
	by label;
	modeleffects Estimate;
    stderr StdErr;
    run;
data &exposure; retain label exposure estimate lclmean uclmean betaCI1 betaCI2 Probt; set est2; 
exposure = "&exposure";
betaCI1 = TRIM(LEFT(PUT(Estimate,9.1)))||" ("||TRIM(LEFT(PUT(LCLMean,9.1)))||", "||TRIM(LEFT(PUT(UCLMean,9.1)))||")";
betaCI2 = TRIM(LEFT(PUT(Estimate,9.2)))||" ("||TRIM(LEFT(PUT(LCLMean,9.2)))||", "||TRIM(LEFT(PUT(UCLMean,9.2)))||")";
keep label exposure estimate lclmean uclmean betaCI1 betaCI2 Probt;
run;
%mend;

/*6-8 weeks*/
%loopmodels(ln_bpa_obrien_v68, 0.92); 
%loopmodels(ln_bp_3_obrien_v68, 1.67);
%loopmodels(ln_dcp_24_obrien_v68, 1.09);
%loopmodels(ln_dcp_25_obrien_v68, 2.08);
%loopmodels(ln_m_pb_obrien_v68, 2.37);
%loopmodels(ln_p_pb_obrien_v68, 2.66);
%loopmodels(ln_tcs_obrien_v68, 1.35);

data beta_estimates; merge ln_dcp_24_obrien_v68 ln_bpa_obrien_v68 ln_bp_3_obrien_v68 ln_dcp_25_obrien_v68 ln_m_pb_obrien_v68 ln_p_pb_obrien_v68 ln_tcs_obrien_v68; by exposure; run;
proc transpose data=beta_estimates out=flip_beta_estimates;
id label;
var betaCI2;
by exposure;
run;

ods csv file="J:\StevensLab\IFED Study\IFED Phenols & Growth_DRS\Analysis\Results\BMI z-score estimates at 6-8 weeks.csv";
proc print data=beta_estimates; run;
ods csv close;
ods csv file="J:\StevensLab\IFED Study\IFED Phenols & Growth_DRS\Analysis\Results\Table of BMI z-score Estimates at 6-8 Weeks.csv";
proc print data=flip_beta_estimates; run;
ods csv close;


/*12 weeks*/
%loopmodels(ln_bpa_obrien_v12, 0.92); 
%loopmodels(ln_bp_3_obrien_v12, 1.67);
%loopmodels(ln_dcp_24_obrien_v12, 1.09);
%loopmodels(ln_dcp_25_obrien_v12, 2.08);
%loopmodels(ln_m_pb_obrien_v12, 2.37);
%loopmodels(ln_p_pb_obrien_v12, 2.66);
%loopmodels(ln_tcs_obrien_v12, 1.35);

data beta_estimates; merge ln_dcp_24_obrien_v12 ln_bpa_obrien_v12 ln_bp_3_obrien_v12 ln_dcp_25_obrien_v12 ln_m_pb_obrien_v12 ln_p_pb_obrien_v12 ln_tcs_obrien_v12; by exposure; run;
proc transpose data=beta_estimates out=flip_beta_estimates;
id label;
var betaCI2;
by exposure;
run;

ods csv file="J:\StevensLab\IFED Study\IFED Phenols & Growth_DRS\Analysis\Results\BMI z-score estimates at 12 weeks.csv";
proc print data=beta_estimates; run;
ods csv close;
ods csv file="J:\StevensLab\IFED Study\IFED Phenols & Growth_DRS\Analysis\Results\Table of BMI z-score Estimates at 12 Weeks.csv";
proc print data=flip_beta_estimates; run;
ods csv close;






/**********************************/
/*Non-Creatinine Adjusted Exposure*/
/**********************************/
proc sort data=ifed_long; by _imp ifedid ageinweeks; run;

/*Loop to run in overall dataframe*/ 
%macro loopmodels(exposure, iqr);

/*IQR-standardize*/
data ifed_long; set ifed_long;
iqr_&exposure = &exposure / &iqr;
run;

/*glimmix with cubic b-spline & random intercepts: No interaction*/
proc GLIMMIX data=ifed_long;
 class ifedid FeedingGroup mage meduc childsgender crace gestationalage smkpreg alcpregany;
 effect sp_age = spline(ageinweeks / basis=bspline knotmethod=list(3.42857, 7.85714, 15.8571, 24.2857) degree=3);
 model bmiz =  iqr_&exposure sp_age 
               FeedingGroup mage meduc childsgender crace gestationalage smkpreg alcpregany/ SOLUTION;
 random intercept / SUBJECT=ifedid S VCORR G V TYPE=UN; 
 by _imp;
 estimate 'bmiz' iqr_&exposure 1 / cl;
 ods output estimates=est;
run; 
quit;
proc sort data=est; by label _imp; run;
ods output parameterestimates=est2;
proc mianalyze data=est;
	by label;
	modeleffects Estimate;
    stderr StdErr;
    run;
data &exposure; retain label exposure estimate lclmean uclmean betaCI1 betaCI2 Probt; set est2; 
exposure = "&exposure";
betaCI1 = TRIM(LEFT(PUT(Estimate,9.1)))||" ("||TRIM(LEFT(PUT(LCLMean,9.1)))||", "||TRIM(LEFT(PUT(UCLMean,9.1)))||")";
betaCI2 = TRIM(LEFT(PUT(Estimate,9.2)))||" ("||TRIM(LEFT(PUT(LCLMean,9.2)))||", "||TRIM(LEFT(PUT(UCLMean,9.2)))||")";
keep label exposure estimate lclmean uclmean betaCI1 betaCI2 Probt;
run;
%mend;

/*6-8 weeks*/
%loopmodels(bpa_gm_ln, 0.92); 
%loopmodels(bp_3_gm_ln, 1.67);
%loopmodels(dcp_24_gm_ln, 1.09);
%loopmodels(dcp_25_gm_ln, 2.08);
%loopmodels(m_pb_gm_ln, 2.37);
%loopmodels(p_pb_gm_ln, 2.66);
%loopmodels(tcs_gm_ln, 1.35);

data beta_estimates; merge dcp_24_gm_ln bpa_gm_ln bp_3_gm_ln dcp_25_gm_ln m_pb_gm_ln p_pb_gm_ln tcs_gm_ln; by exposure; run;
proc transpose data=beta_estimates out=flip_beta_estimates;
id label;
var betaCI2;
by exposure;
run;

ods csv file="J:\StevensLab\IFED Study\IFED Phenols & Growth_DRS\Analysis\Results\BMI z-score estimates for non-creatinine adjusted exposures.csv";
proc print data=beta_estimates; run;
ods csv close;
ods csv file="J:\StevensLab\IFED Study\IFED Phenols & Growth_DRS\Analysis\Results\Table of BMI z-score Estimates for Non-Creatinine Adjusted Exposures.csv";
proc print data=flip_beta_estimates; run;
ods csv close;


