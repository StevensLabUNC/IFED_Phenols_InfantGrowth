/*************************************************************************************************************
  Programmer: Danielle Stevens
  Purpose of syntax: Clean and examine IFED data
  Last updated on: 11/17/2026
**************************************************************************************************************/

libname ifed "J:\StevensLab\IFED Study\Original Data";
libname ifed2 "J:\StevensLab\IFED Study\IFED Phenols & Growth_DRS\Analysis\Derived data";
options nofmterr;

/* Formats for derived variables;
          value AMPMf    1 = "1. AM"
                         2 = "2. PM";

          value Educf    1 = "1. <HS"
                         2 = "2. HS or GED"
                         3 = "3. Some college"
                         4 = "4. College";

          value ETHNF    1 = "Non-Hispanic/Latino-a"
                         2 = "Hispanic/Latino-a"
                         3 = "Unknown or Not Reported";

          value FeedGF   1 = "Breast Milk"
                         2 = "Cow-Milk Formula"
                         3 = "Soy-Milk Formula";

          value racef    1 = "1. Black"
                         2 = "2. White"
                         3 = "3. Other/multi/unknown";

value Race2f
  135 = "135. American Indian or Alaskan Native, Black or African American, White"
  235 = "235. Asian, Black or African American, White"
  13 = "13. American Indian or Alaskan Native, Black or African American"
  23 = "23. Asian, Black or African American"
  25 = "25. Asian, White"
  34 = "34. Black or African American, Native Hawaiian or Other Pacific Islander"
  35 = "35. Black or African American, White"
  1 = "1. American Indian or Alaskan Native"
  2 = "2. Asian"
  3 = "3. Black or African American"
  4 = "4. Native Hawaiian or Other Pacific Islander"
  5 = "5. White"
  6 = "6. Unknown or Not Reported";

        value seasonf   1 = '1. Winter'
			    2 = '2. Spring'
			    3 = '3. Summer'
			    4 = '4. Fall';

           value sexf   1 = "Female" --> NOTE: THIS WAS INCORRECT ORIGINALLY AND FIXED
                        2 = "Male";

        value specmisf  .S = "s. skipped"
                        .M = "m. missing";

            value YNMf  1 = "1. YES"
                        2 = "2. NO"
                       .R = "r. refused"
                       .D = "d. don't know"
                       .M = "m. missing"
                       .V = "v. missing due to version"
                       .S = "s. skipped"
                       .N = "n. no visit";
*/



/****************************/
/*Get together the exposures*/
/****************************/
data ifed1; set ifed.Ifd124001_phenol_20220421; 
/*Eliminate n=1 with soy milk consumption*/
if FeedingGroup = 3 then DELETE;

/*Weight as a variable often reads in incorrectly - let's fix for analysis*/
infant_weight = Weight; 

/*Age variables*/
AgeinWeeks=AgeInDays/7; 
AgeinMonths=AgeInDays/30; 

/*Set creatinine to missing so you can impute all concentrations together*/
creat = Creat_I;
if Creat_I = 1/sqrt(2) then creat = .;

/*Create flag for missing exposures: 1 if missing all exposures and creatinine, not due to < LOD*/
flag_exp = 0;
if (BPA=. AND BPA_Com = " ") AND 
   (BP_3=. AND BP_3_Com = " ") AND
   (DCP_24=. AND DCP_24_Com = " ") AND
   (DCP_25=. AND DCP_24_Com = " ") AND 
   (M_PB=. AND M_PB_Com = " ") AND 
   (P_PB=. AND P_PB_Com = " ") AND
   (TCS=. AND TCS_Com = " ") AND 
    Creat_I = . THEN flag_exp = 1;
run;

proc means data=ifed1 max;
var AgeinWeeks;
run; *40 weeks or 10 months;


/*Creat_I --> imputed if below LOD using lod/sqrt(2). This occurred for only 1 sample with <LOD.

Season batch AMPM --> expected to be predictors of exposure

BPA BP_3 DCP_24 DCP_25 M_PB P_PB TCS

TCS detected in 64% of samples, but the rest detected in >90%*/

/*Count number of exposure assessments/participant*/
data exposures; set ifed1;
where flag_exp = 0;
by ifedid;
   if first.ifedid then num_exp = 1;
   else num_exp + 1;
   if last.ifedid then output;
retain ifedid num_exp;
keep ifedid num_exp;
run;
proc freq data=exposures;
tables num_exp;
run; /*n=47 with just 1 assessmnet and n=153 with 2 assessments*/

data ifed2; merge ifed1 exposures;
by ifedid;
run;
proc freq data=ifed2;
tables visit*flag_exp;
run; /*113 samples at the 6 week, 73 at the 8 week, 166 at the 12 week visit*/
proc means data=ifed2;
class flag_exp;
var ageinweeks;
run; /*5-13 weeks so exposure was assessed within the allotted timeframe*/


data exposures; set ifed2;
where flag_exp = 0;
by ifedid;
   if first.ifedid then count_exp = 1;
   else count_exp + 1;
retain ifedid count_exp;
keep ifedid visit count_exp;
run;
data ifed3; merge ifed2 exposures;
by ifedid visit;
run;
%macro exposureoutput(exposure);
proc means data=ifed3 n q1 median q3;
var &exposure;
output out=out n=n q1=q1 median=q2 q3=q3;
where visit in (6 8) and flag_exp = 0;
run;
%mend;
%exposureoutput(bpa);
%exposureoutput(bp_3);
%exposureoutput(dcp_24);
%exposureoutput(dcp_25);
%exposureoutput(m_pb);
%exposureoutput(p_pb);
%exposureoutput(tcs);

/*LOD levels*/
proc freq data=ifed3;
tables bpa_com bp_3_com dcp_24_com dcp_25_com m_pb_com p_pb_com TCS_Com / missing;
where flag_exp = 0 AND visit in (6 8);
run;
proc freq data=ifed3;
tables bpa_com bp_3_com dcp_24_com dcp_25_com m_pb_com p_pb_com TCS_Com / missing;
where flag_exp = 0 AND visit in (12);
run;
proc freq data=ifed3;
tables bpa_com bp_3_com dcp_24_com dcp_25_com m_pb_com p_pb_com TCS_Com / missing;
where flag_exp = 0;
run;


/*Potential covariates:

  FeedingGroup 
  ChildsGender
  Crace or Childsraces and childsethnicity
  GestationalAge
  parity
  smkpreg
  SmkOthers12Wk and SmkInside12Wk
  AlcPregAny
  Edu_cat4
  MothersAgeAtChildsDOB

*/

data first; set ifed3;
by ifedid;
if first.ifedid then output;
run;
proc freq data=first;
tables FeedingGroup ChildsGender / missing;
tables Crace Childsraces crace*childsethnicity / missing;
tables parity / missing;
tables smkpreg smkothers12wk*smkinside12wk / missing;
tables alcpregany / missing;
tables edu_cat4 / missing;
run;
proc univariate data=first plots;
var mothersageatchildsdob gestationalage;
run;

data first; set first;
smkpost = .;
  if SmkOthers12Wk = 2 OR SmkInside12Wk = 2 then smkpost=0; *no smoking if any no;
  if SmkOthers12Wk = 1 OR SmkInside12Wk = 1 then smkpost=1; *overwrite if any yes;

mage = .;
   if . < mothersageatchildsdob LE 20 then mage = 1;
   if 20 < mothersageatchildsdob LE 25 then mage = 2;
   if 25 < mothersageatchildsdob LE 30 then mage = 3;
   if 30 < mothersageatchildsdob then mage = 4;

meduc = .;
  if . < edu_cat4 LE 2 then meduc = 1;
  if 2 < edu_cat4 LE 4 then meduc = 2;
run;
proc freq data=first;
tables smkpreg*smkpost  / missing;
tables FeedingGroup / missing;
tables mage gestationalage meduc/ missing;
run; /*Plenty of overlap but not enough smoking postpartum for me to want to lump together. Let's just do smkpreg for analyses*/

/*Too many skipped parity - let's not include as covariate*/
data ifed3; set ifed3;
mage = .;
   if . < mothersageatchildsdob LE 20 then mage = 1;
   if 20 < mothersageatchildsdob LE 25 then mage = 2;
   if 25 < mothersageatchildsdob LE 30 then mage = 3;
   if 30 < mothersageatchildsdob then mage = 4;

meduc = .;
  if . < edu_cat4 LE 2 then meduc = 1;
  if 2 < edu_cat4 LE 4 then meduc = 2;
run;


/******************/
/*Examine outcomes*/
/******************/
%macro plotoutcomes(outcome);
proc sgplot data=ifed3;
 series x=AgeInDays y=&outcome / group=IFEDID;
 label AgeInDays="Age (days)";
 run;

proc means data=ifed3;
class AgeInMonths;
var &outcome;
output out=out mean=mean;
run;
proc sgplot data=out;
 series x=AgeInMonths y=mean;
 where _TYPE_ ne 0;
 label AgeInMonths="Age (months)";
 run;

proc means data=ifed3;
class visit;
var &outcome;
output out=&outcome n=n mean=mean;
run;
data &outcome; set &outcome;
where visit ne .;
n_&outcome = TRIM(LEFT(PUT(visit, 9.0)))||" (n="||TRIM(LEFT(PUT(n, 9.0)))||")";
keep visit n_&outcome;
run;
%mend;

%plotoutcomes(infant_weight);
%plotoutcomes(height);
%plotoutcomes(waz);
%plotoutcomes(haz);
%plotoutcomes(bmiz);

/*Merge into main dataset*/
proc sort data=ifed3; by visit ifedid; run;
data ifed4; merge ifed3(in=aaa) infant_weight height waz haz bmiz; by visit; if aaa; run;


/*Do for select groups*/
%macro plotoutcomes(wherevar, wherecond, output, outcome);
data data; set ifed4;
where &wherevar = &wherecond;
run;

proc means data=data;
class visit;
var &outcome;
output out=&output n=n mean=mean;
run;
data &output; set &output;
where visit ne .;
n_&output = TRIM(LEFT(PUT(visit, 9.0)))||" (n="||TRIM(LEFT(PUT(n, 9.0)))||")";
keep visit n_&output;
run;
%mend;

%plotoutcomes(FeedingGroup, 1, breastmilk_weight, infant_weight);
%plotoutcomes(FeedingGroup, 2, cowsmilk_weight, infant_weight);
%plotoutcomes(FeedingGroup, 1, breastmilk_height, height);
%plotoutcomes(FeedingGroup, 2, cowsmilk_height, height);
%plotoutcomes(FeedingGroup, 1, breastmilk_waz, waz);
%plotoutcomes(FeedingGroup, 2, cowsmilk_waz, waz);
%plotoutcomes(FeedingGroup, 1, breastmilk_haz, haz);
%plotoutcomes(FeedingGroup, 2, cowsmilk_haz, haz);
%plotoutcomes(FeedingGroup, 1, breastmilk_bmiz, bmiz);
%plotoutcomes(FeedingGroup, 2, cowsmilk_bmiz, bmiz);

%plotoutcomes(childsgender, 1, female_weight, infant_weight);
%plotoutcomes(childsgender, 2, male_weight, infant_weight);
%plotoutcomes(childsgender, 1, female_height, height);
%plotoutcomes(childsgender, 2, male_height, height);
%plotoutcomes(childsgender, 1, female_waz, waz);
%plotoutcomes(childsgender, 2, male_waz, waz);
%plotoutcomes(childsgender, 1, female_haz, haz);
%plotoutcomes(childsgender, 2, male_haz, haz);
%plotoutcomes(childsgender, 1, female_bmiz, bmiz);
%plotoutcomes(childsgender, 2, male_bmiz, bmiz);

/*Merge into main dataset*/
proc sort data=ifed4; by visit; run;
data ifed4; merge ifed4(in=aaa)
 breastmilk_weight breastmilk_height breastmilk_waz breastmilk_haz breastmilk_bmiz
 cowsmilk_weight cowsmilk_height cowsmilk_waz cowsmilk_haz cowsmilk_bmiz
 male_weight male_height male_waz male_haz male_bmiz
 female_weight female_height female_waz female_haz female_bmiz; 
 by visit;
 if aaa;
 run;


/*How many visits?*/
proc sort data=ifed4;
by ifedid visit;
run;
data ifed4; set ifed4;
drop visnum;
run;
data ifed4; set ifed4;
by ifedid visit;
if first.ifedid then visnum=1;
else visnum+1;
run;
proc freq data=ifed4;
tables visnum visit / missing;
run; /*12 visits total, all 199 infants have at least 4 visits*/
data last_vis; set ifed4;
by ifedid visit;
if first.ifedid then max_vis=1;
else max_vis+1;
retain ifedid max_vis;
if last.ifedid then output;
run;
proc freq data=last_vis;
tables max_vis visit / missing;
run; /*51 infants have 12 visits, 32 have 11 visits, 59 have 10 visits, 31 have 9 visits, and the other 26 infants have <9 visits*/
proc means data=last_vis median q1 q3;
var max_vis;
*where BMIz ne .; *same regardless of having this statement or not;
run; 

/*Do any infants have > 1 visit / week?*/
proc freq data=ifed4;
tables ifedid*ageinweeks / out=out;
run;
proc sort data=out;
by count;
run; /*None with >1 visit/week*/



/******************************/
/*Flip long to wide: outcomes */
/******************************/

data ifed_wide; set ifed4;
by ifedid;
retain ifedid weight_1-weight_12 height_1-height_12 waz_1-waz_12 haz_1-haz_12 whz_1-whz_12 bmiz_1-bmiz_12 age_dys_1-age_dys_12;
array a{1:12} weight_1-weight_12;
array b{1:12} height_1-height_12;
array c{1:12} waz_1-waz_12;
array d{1:12} haz_1-haz_12;
array e{1:12} whz_1-whz_12;
array f{1:12} bmiz_1-bmiz_12;
array g(1:12) age_dys_1-age_dys_12;
if first.ifedid then do;
do z=1 to 12 by 1;
   a(z)=.;
   b(z)=.;
   c(z)=.;
   d(z)=.;
   e(z)=.;
   f(z)=.;
   g(z)=.;
end;
end;
a(visnum)=weight;
b(visnum)=height;
c(visnum)=waz;
d(visnum)=haz;
e(visnum)=whz;
f(visnum)=bmiz;
g(visnum)=ageindays;
if last.ifedid then output;
keep ifedid feedinggroup mage meduc childsgender crace gestationalage smkpreg alcpregany 
     weight_1-weight_12 height_1-height_12 waz_1-waz_12 haz_1-haz_12 whz_1-whz_12 bmiz_1-bmiz_12 age_dys_1-age_dys_12;
run;

/*Heggeseth & Jewell 2018 discusses subtracting the mean level for each individual in order to normalize trajectories so that the pattern/trajectory of growth 
  matters more than the level of growth, which tends to generally dominate group assignment for LCGM

  We will do so and re-run the LCGA to see what kinds of trajectories we are getting*/

data ifed_wide2; set ifed_wide;
mean_bmiz = mean(of bmiz_1-bmiz_12);
mean_waz = mean(of waz_1-waz_12);
mean_haz = mean(of haz_1-haz_12);

array a(*) norm_bmiz_1-norm_bmiz_12;
array b(*) base_bmiz_1-base_bmiz_12;
array c(*) bmiz_1-bmiz_12;

array d(*) norm_waz_1-norm_waz_12;
array e(*) base_waz_1-base_waz_12;
array f(*) waz_1-waz_12;

array g(*) norm_haz_1-norm_haz_12;
array h(*) base_haz_1-base_haz_12;
array j(*) haz_1-haz_12;

do z=1 to 12;
   a(z) = c(z)-mean_bmiz; *normalize to the mean value;
   b(z) = c(z)-bmiz_1; *normalize to the baseline value - we didn't end up using this in the analysis;

   d(z) = f(z)-mean_waz; *normalize to the mean value;
   e(z) = f(z)-waz_1; *normalize to the baseline value - we didn't end up using this in the analysis;

   g(z) = j(z)-mean_haz; *normalize to the mean value;
   h(z) = j(z)-haz_1; *normalize to the baseline value - we didn't end up using this in the analysis;
end;
run;

/*Create splines for age*/
proc sort data=ifed4 nodup;
by ifedid ageinweeks;
run;
proc univariate data=ifed4 ;
var ageinweeks;
output out=percentiles pctlpts = 25 50 75 pctlpre=p_; /*calculate the 25th and 75th percentiles*/
run; /* 4, 11.4286, 23 */
proc print data=percentiles; run;
proc univariate data=ifed4 ;
var ageinweeks;
output out=percentiles pctlpts = 33 66 pctlpre=p_; /*calculate the 33rd and 66th percentiles*/
run; 
proc print data=percentiles; run;
proc univariate data=ifed4 ;
var ageinweeks;
output out=percentiles pctlpts = 20 40 60 80 pctlpre=p_; /*calculate the 33rd and 66th percentiles*/
run; 
proc print data=percentiles; run;
/*Manually code splines based on percentiles*/
data ifed4; set ifed4;
spline1 = ageinweeks;
spline2 = (ageinweeks > 4)*(ageinweeks - 4);
spline3 = (ageinweeks > 11.4286)*(ageinweeks - 11.4286);
spline4 = (ageinweeks > 23)*(ageinweeks-23);
run;


/*******************************/
/*Flip long to wide: exposures */
/*******************************/
data exp_v68; set ifed4;
where visit in (6 8) and flag_exp = 0;
keep ifedid visit season ampm ageinweeks batch creat bpa bp_3 dcp_24 dcp_25 m_pb p_pb tcs;
run; /*n=186*/
proc sql noprint;
   select cats(name,'=',name,'_v68')
          into :list
          separated by ' '
          from dictionary.columns
          where libname = 'WORK' and memname = 'EXP_V68';
quit;
proc datasets library = work nolist;
   modify exp_v68;
   rename &list;
quit;
data exp_v68; set exp_v68;
ifedid = ifedid_v68; drop ifedid_v68; run;

data exp_v12; set ifed4;
where visit in (12) and flag_exp = 0;
keep ifedid season ampm ageinweeks batch creat bpa bp_3 dcp_24 dcp_25 m_pb p_pb tcs;
run; /*n=166*/
proc sql noprint;
   select cats(name,'=',name,'_v12')
          into :list
          separated by ' '
          from dictionary.columns
          where libname = 'WORK' and memname = 'EXP_V12';
quit;
proc datasets library = work nolist;
   modify EXP_v12;
   rename &list;
quit;
data exp_v12; set exp_v12;
ifedid = ifedid_V12; drop ifedid_v12; run;
DATA exposures; merge exp_v68 exp_v12; by ifedid; run;


data ifed.ifed_long; set ifed4; run;
data ifed.ifed_wide; set ifed_wide2; run;
data ifed.exposures; set exposures; run;

