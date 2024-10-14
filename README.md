# SGLT-2-inhibitor-and-dementia-risk-propensity-score-stratification
/*** Section 1: exclusion criteria ***/

* Extract records from the raw data;
%macro records(i=);
%do n = 1 %to &i;
proc sql;
create table medication as
select a.id, a.product_code, a.date, b.class
from source.therapy&n as a, patient as b
where a.product_code = b.prodcodeid;
quit;
%end;
%mend;
%records(i=number);

* Subset based on exclusion criterion; 
proc sql;
create table patient1 as 
select a.* 
from patient as a inner join medication as b
on a.id = b.id where condition
order by id;
quit;


/*** Section 2: propensity score fine stratification ***/

* Logistic regression to obtain the propensity scores;
proc logistic data= data; 
class categorical covariates;
model treatment = continuous and categorical covariates;
output out = data_ps prob = ps;
run;

* Check overlap of propensity score between treatment groups;
proc univariate data = data_ps;
var ps;
class treatment;
histogram /outhistogram = h1 midpoints = 0 to 1 by 0.0025;
run;

* Remove subjects with non-overlapping propensity scores;
proc sql;
create table psrange as
select *, max(ps) as maxps, min(ps) as minps
from data_ps
group by treatment;
create table pstrim as
select *, max(minps) as trimminps, min(maxps) as trimmaxps
from psrange having ps >= calculated trimminps and ps <= calculated trimmaxps
order by id;
quit;

* Create 50 strata by using exposure group alone;
data exposed; set pstrim; where treatment = "exposure"; run;
proc rank data=exposed groups = 50 out = strata_exp; ranks strata_ps; var ps; run;

data fs_exp (drop = strata_ps); set strata_exp; strata = 1+strata_ps; run;

proc means data = fs_exp noprint; 
class strata;
var ps;
output out = bounds_for_fs n = n
min(ps) = PS_min max(ps) = PS_max;
run;

data bounds_for_fs; set bounds_for_fs; where strata ne .; run;
proc transpose data=bounds_for_fs out = bounds prefix = psrange; var ps_min; run;

* Find comparison subjects with propensity scores falling within each stratum;
data unexposed; set pstrim; where treatment = "reference"; run;
data unexp1; merge unexposed bounds; run;

proc stdize out = unexp2 reponly missing = mean; run;

data unexp_strata; 
set unexp2;
array abc(*) psrange: ;
do i = 3 to dim(abc) while (strata = .);
if 0 <= ps < psrange2 then strata = 1; 
if abc(i-1) <= ps < abc(i) then strata = (i-1); 
end;
run;

data fs_unexp (drop = i psrange:); set unexp_strata; if strata = . then strata = 50; run;

proc sort data = fs_exp; by id; run;
proc sort data = fs_unexp; by id; run;
data stratified; set fs_exp fs_unexp; by id; run;

* Calculate the stratum-specific weights;
proc sort data = stratified out = matched_ps1; by strata; run;
data total_exposed; set matched_ps1; where treatment = "exposure"; run;

data total_exposed1; 
set total_exposed; 
by strata; 
if first.strata then total_Exp = 0; total_Exp + 1;
if last.strata then output;
keep strata total_exp;
run;

data total_unexposed; set matched_ps1; where treatment = "reference"; run;

data total_unexposed1; 
set total_unexposed; 
by strata;
if first.strata then total_unExp = 0; total_unExp + 1;
if last.strata then output;
keep strata total_unexp;
run;

data strata_n;
merge total_exposed1 total_unexposed1;
by strata;
if total_exp = . then total_exp = 0; 
if total_unexp = . then total_unexp = 0;
run;

data strata_n; set strata_n; if total_exp = 0 or total_unexp = 0 then delete; run;

proc sql noprint; select sum(total_exp) into: sum_exp from strata_n; quit;
proc sql noprint; select sum(total_unexp) into: sum_unexp from strata_n; quit;

data strata_n; 
set strata_n;
weight=(total_exp/&sum_exp)/(total_unexp/&sum_unexp);
weighted_n_unexp=total_unexp*((total_exp/&sum_exp)/(total_unexp/&sum_unexp));
run;

proc sql noprint; select sum (weighted_n_unexp) into: sum_unexp_weighted  from strata_n; quit;

data numbers;
length n $50;
n='After dropping uninformative strata';
total_exposed= &sum_exp; 
total_unexposed=&sum_unexp; 
total_weighted_unexp=&sum_unexp_weighted;
run;

proc sql;
create table with_weights as
select stratified.*, strata_n.weight
from stratified left join strata_n
on stratified.strata = strata_n.strata
where weight ne .; 
quit;

data PS_weights; 
set with_weights;
psweight=1; 
if treatment = "reference" then psweight = weight; 
drop weight;
run;

proc sort data = PS_weights; by id; run;

/*** Section 3: survival analysis ***/

* KM curve;
proc lifetest data = data_model plot = survival(atrisk);
time time* status(0);
strata treatment;
run;
* Crude hazard ratio;
proc freq data = data_model; table treatment*status / nopercent norow; run;
proc means data = data_model sum; class treatment; var time; run;
proc phreg data = data_model; 
class treatment; 
model time*status(0) = treatment / rl; 
run;

* Adjusted hazard ratio;
proc phreg data = model_weighted;
class treatment; 
model time*status(0) = treatment / rl; 
weight psweight; 
run;


