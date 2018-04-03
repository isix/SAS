# Goodness-of-Fit Statistics for Quantile Regression
# (https://www.federalreserve.gov/econresdata/feds/2016/files/2016002r1pap.pdf) 

%macro quant_gof(data = , y = , x = , tau = 0.5);
***********************************************************;
* THE MACRO CALCULATES GOODNESS-OF-FIT STATISTICS FOR     *;
* QUANTILE REGRESSION                                     *;
* ------------------------------------------------------- *;
* REFERENCE:                                              *;
*  GOODNESS OF FIT AND RELATED INFERENCE PROCESSES FOR    *;
*  QUANTILE REGRESSION, KOENKER AND MACHADO, 1999         *;
***********************************************************;

options nodate nocenter;
title;

* UNRESTRICTED QUANTILE REGRESSION *;
ods select ParameterEstimates ObjFunction;
ods output ParameterEstimates = _est;
proc quantreg data = &data ci = resampling(nrep = 500);
  model &y = &x / quantile = &tau nosummary nodiag seed = 1;
  output out = _full p = _p;
run;

* RESTRICTED QUANTILE REGRESSION *;
ods select none;
proc quantreg data = &data ci = none;
  model &y = / quantile = &tau nosummary nodiag;
  output out = _null p = _p;
run;
ods select all; 

proc sql noprint;
  select sum(df) into :p from _est;
quit;

proc iml;
  use _full;
  read all var {&y _p} into A;
  close _full;

  use _null;
  read all var {&y _p} into B;
  close _null;

  * DEFINE A FUNCTION CALCULATING THE SUM OF ABSOLUTE DEVIATIONS *;
  start loss(x);
    r = x[, 1] - x[, 2];
    z = j(nrow(r), 1, 0);
    l = sum(&tau * (r <> z) + (1 - &tau) * (-r <> z));
    return(l);
  finish;
  
  r1 = 1 - loss(A) / loss(B);
  adj_r1 = 1 - ((nrow(A) - 1) * loss(A)) / ((nrow(A) - &p) * loss(B));
  aic = 2 * nrow(A) * log(loss(A) / nrow(A)) + 2 * &p;
  aicc = 2 * nrow(A) * log(loss(A) / nrow(A)) + 2 * &p * nrow(A) / (nrow(A) - &p - 1);
  bic = 2 * nrow(A) * log(loss(A) / nrow(A)) + &p * log(nrow(A));
  
  l = {"R1" "ADJUSTED R1" "AIC" "AICC" "BIC"};
  v = r1 // adj_r1 // aic // aicc // bic;
  print v[rowname = l format = 20.8 label = "Fit Statistics"];
quit;

%mend quant_gof;
