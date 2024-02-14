{smcl}
{* *! version 1.0.2  May 1, 2023 @ 00:00:00}{...}
{cmd:help fect}
{hline}

{title:Title}

{p2colset 5 20 22 2}{...}  {p2col :{hi:fect} {hline 2}}
 A Practical Guide to Counterfactual Estimators for Causal Inference with Time-Series Cross-Sectional Data{p_end} {p2colreset}{...}


{title:Syntax}

{p 8 17 2}
{cmdab:fect}
{it:{help varname:outcome}}
[if]
{it:,{help varname:treat(varname)}}
{it:{help varname:unit(varname)}}
{it:{help varname:time(varname)}}
[{it:options}]



{synoptset 23 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt treat(varlist)}} specifying the treatment indicator. {p_end}
{synopt:{opt unit(varlist)}} specifying the unit (group) indicator. {p_end}
{synopt:{opt time(varlist)}} specifying the time indicator. {p_end}
{synopt:{opt cov(varlist)}} specifying time-varying covariates. {p_end}
{synopt:{opt weight(varlist)}} specifying weight indicator, will be used to adjust various treatment effects. {p_end}
{synopt:{opt force(string)}} a string indicating whether unit or time fixed effects will be imposed. Must be one of the following, {opt two-way} (default), {opt unit}, {opt time} and {opt none}.{p_end}
{synopt:{opt degree(integer)}} an integer specifying the order of the polynomial trend term, should be larger than 1, default to 3. {p_end}
{synopt:{opt nknots(integer)}} an integer specifying the number of knots used for a restricted cubic spline, must be between 3 and 7, default to 3.{p_end}
{synopt:{opt r(integer)}} an integer specifying the number of factors. If {opt cv} is on, the cross validation procedure will select the optimal number of factors from 1 to r.{p_end}
{synopt:{opt lambda(numlist)}} a single or sequence of positive numbers specifying the hyper-parameter sequence for matrix completion method. If {opt lambda} is a sequence, cross-validation will be performed. {p_end}
{synopt:{opt nlambda(integer)}} an integer specifying the length of hyper-parameter sequence for matrix completion method, default to 10.{p_end}
{synopt:{opt cv}} a flag indicating whether cross-validation will be performed to select the optimal number of factors or hyper-parameter in matrix completion algorithm.{p_end}
{synopt:{opt kfold(integer)}} an integer specifying number of cross-validation rounds, default to 10. {p_end}
{synopt:{opt cvtreat}} a flag indicating  whether to only use observations of treated units as testing set in cross validation.{p_end}
{synopt:{opt cvobs(integer)}} an integer specifying the length of continuous observations within a unit in the testing set, default to 3.{p_end}
{synopt:{opt method(string)}} a string specifying which matrix completion algorithm will be used. {opt fe}(default) for two-way fixed effects model, {opt ife} for interactive fixed effects model, {opt mc} for matrix copletion method, {opt polynomial} for polynomial trend terms, {opt bspline} for regression splines and {opt both} for both of the IFE model and the MC model. {p_end} 
{synopt:{opt alpha(real)}} significant level for hypothesis test and CIs, default to 0.05. {p_end}
{synopt:{opt se}} a flag indicating whether uncertainty estimates will be produced, needed for the equivalence test and the placebo test. {p_end}
{synopt:{opt vartype(string)}} a string specifying the type of variance estimator. Must be one of the following, {opt bootstrap}(default) or {opt jackknife}. {p_end}
{synopt:{opt nboots(integer)}} an integer specifying the number of bootstrap runs. {p_end}
{synopt:{opt tol(real)}} a positive number indicating the tolerance level, default to 1e-5, a small {opt tol} can yield a more accurate result but may make the programe time-consuming. {p_end}
{synopt:{opt maxiterations(integer)}} an integer specifying the upper limit of number of interactions, default to 2000. {p_end}
{synopt:{opt preperiod(integer)}} an negative integer specifying the range of pre-treatment period used for goodness-of-fit test(Wald test and Equivalence test) and ATT plot. If left blank, all pre-treatment periods will be used. {p_end}
{synopt:{opt offperiod(integer)}} an positive integer specifying the range of post-treatment period used for ATT plot. If left blank, all post-treatment periods will be used. {p_end}
{synopt:{opt proportion(real)}} a positive real value specifying pre-treatment periods that have observations larger than the proportion of observations at period 0. Default is 0.3 {p_end}
{synopt:{opt wald}} a flag indicating whether to perform wald test for pre-treatment fitting check. {p_end}
{synopt:{opt placeboTest}} a flag indicating whether to perform placebo test. {p_end}
{synopt:{opt placeboperiod(integer)}} an positive integer specifying the range of pre-treatment period that will be assigned as "placebo" treatment period. {p_end}
{synopt:{opt equiTest}} a flag indicating whether to perform equivalence test. {p_end}
{synopt:{opt permutation}} a flag indicating whether to perform permutation test. {p_end}


{syntab:Advanced}
{synopt:{opt seed(integer)}} an integer that sets the seed in random number generation, default to 100. {p_end}
{synopt:{opt minT0(integer)}} an integer specifying the minimum value of observed periods that a unit is under control, default to 5. {p_end}
{synopt:{opt maxmissing(integer)}} an integer. Units with number of missing values greater than it will be removed, default to NULL(ignored). {p_end}
{synopt:{opt title(string)}} set the title of the graph {p_end}
{synopt:{opt xlabel(string)}} set the label for the x-axis in the graph {p_end}
{synopt:{opt ylabel(string)}} set the label for the y-axis in the graph {p_end}
{synopt:{opt saving(string)}} specify a filename with which the graph will be stored {p_end}
{synopt:{opt counterfactual}} a flag indicating whether to generate the counterfactual after the estimation. {p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
{p_end}


{title:Description}

{pstd} {opt fect} implements counterfactual estimators in TSCS data analysis. These estimators first impute counterfactuals for 
  each treated observation in a TSCS dataset by fitting an outcome model (fixed effects model, interactive fixed effects model, or
  matrix completion) using the untreated observations. They then estimate the individualistic treatment effect for each treated 
  observation by subtracting the predicted counterfactual outcome from its observed outcome. Finally, the average treatment effect
  on the treated (ATT) or period-specific ATTs are calculated. A placebo test and an equivalence test are included to evaluate the
  validity of identification assumptions behind these estimators. Data must be with a dichotomous treatment. {p_end}



{title:Examples}

We use simdata1.dta to illustrate how {cmd:fect} works. 

{pstd}Load simdata1 {p_end}
{phang2}{cmd:. use "https://raw.githubusercontent.com/xuyiqing/fect_stata/master/simdata1.dta", clear}{p_end}

{pstd}Estimate ATT using the fixed effects model(FE) {p_end}
{phang2}{cmd:. fect Y, treat(D) unit(id) time(time) cov(X1 X2) method("fe")}{p_end}

{pstd}Estimate ATT using the interactive fixed effects model(IFE) {p_end}
{phang2}{cmd:. fect Y, treat(D) unit(id) time(time) cov(X1 X2) method("ife") r(2)}{p_end}

{pstd}Estimate ATT using the matrix completion model(MC){p_end}
{phang2}{cmd:. fect Y, treat(D) unit(id) time(time) cov(X1 X2) method("mc") lambda(0.003)}{p_end}

{pstd}Estimate ATT using the IFE model after the cross-validation{p_end}
{phang2}{cmd:. fect Y, treat(D) unit(id) time(time) cov(X1 X2) method("ife") r(4) cv}{p_end}

{pstd}Estimate ATT using the MC model after the cross-validation{p_end}
{phang2}{cmd:. fect Y, treat(D) unit(id) time(time) cov(X1 X2) method("mc") lambda(0.002 0.003 0.004 0.005)}{p_end}

{pstd}Estimate ATT and its confidence intervals using the IFE model{p_end}
{phang2}{cmd:. fect Y, treat(D) unit(id) time(time) cov(X1 X2) method("ife") r(2) se}{p_end}

{pstd}Estimate ATT and its confidence intervals using the MC model and jackknife variance{p_end}
{phang2}{cmd:. fect Y, treat(D) unit(id) time(time) cov(X1 X2) method("mc") lambda(0.003) se vartype("jackknife")}{p_end}

{pstd}Estimate ATT and its confidence intervals using the IFE model and keep pre-treatment periods that have observations larger than the 50 percent of observations at period 0 {p_end}
{phang2}{cmd:. fect Y, treat(D) unit(id) time(time) cov(X1 X2) method("ife") proportion(0.5) r(2) se}{p_end}

{pstd}Conduct the Wald Test {p_end}
{phang2}{cmd:. fect Y,  treat(D) unit(id) time(time) cov(X1 X2) se method("ife") r(2) preperiod(-14) offperiod(5) wald nboots(100)}{p_end}

{pstd}Conduct the Equivalence Test {p_end}
{phang2}{cmd:. fect Y,  treat(D) unit(id) time(time) cov(X1 X2) se method("ife") r(2) preperiod(-14) offperiod(0) equiTest}{p_end}

{pstd}Conduct the Placebo Test {p_end}
{phang2}{cmd:. fect Y,  treat(D) unit(id) time(time) cov(X1 X2) se method("ife") r(2) placeboTest}{p_end}

{pstd}Conduct the Permutation Test {p_end}
{phang2}{cmd:. fect Y,  treat(D) unit(id) time(time) cov(X1 X2) method("ife") r(2) permutation}{p_end}


{title:Saved results}

{pstd}{cmd:fect} returns the estimated average treatment effects(ATT), by-periods ATTs and their confidence intervals {p_end}

{synoptset 24 tabbed}{...}
{syntab: Scalars}
{synopt:{cmd:e(optimal_parameter)}} saves the optimal hyper-parameter in the {opt ife} model or the {opt mc} {p_end}
{synopt:{cmd:e(min_mspe)}} saves the least mean squared prediction error in cross-validation {p_end}
{synopt:{cmd:e(permutation_pvalue)}} saves the p-value in the permutation test {p_end}
{synopt:{cmd:e(placebo_pvalue)}} saves the p-value in the placeboTest test {p_end}

{syntab: Matrixs}
{synopt:{cmd:e(CV)}} saves mean squared prediction errors in cross-validation when {opt cv} is on {p_end}
{synopt:{cmd:e(coefs)}} saves the coefficients of time-varying observables {p_end}
{synopt:{cmd:e(ATT)}} saves estimated average treatment effect(ATT) {p_end}
{synopt:{cmd:e(ATTs)}} saves the estimated by-period average treatment effect(ATTs) {p_end}
{synopt:{cmd:e(placebo_ATT)}} saves estimated pseudo average treatment effect in the placebo region when {opt placeboTest} is on {p_end}


{title:Reference}
{p 4 8 2}

{pstd}Jushan Bai. 2009. "Panel Data Models with Interactive Fixed
  Effects." Econometrica 77:1229--1279.{p_end}

{pstd} Yiqing Xu. 2017. "Generalized Synthetic Control Method: Causal Inference
  with Interactive Fixed Effects Models." Political Analysis, Vol. 25, 
  Iss. 1, January 2017, pp. 57-76. Available at: https://doi.org/10.1017/pan.2016.2. {p_end}

{pstd} Athey, Susan, et al. 2018 "Matrix completion methods for causal panel data models." arXiv preprint arXiv:1710.10251. Available 
  at: https://https://arxiv.org/abs/1710.10251. {p_end}

{pstd} Licheng Liu, et al. 2020. "A Practical Guide to Counterfactual Estimators for Causal Inference with Time-Series Cross-Sectional 
  Data." Working paper. Available at: https://polmeth.mit.edu/sites/default/files/documents/Yiqing_Xu.pdf. {p_end}

{pstd} For more details about the matrix completion method, see https://github.com/susanathey/MCPanel. {p_end}

  
{title:Authors}

  Licheng Liu (MIT); Ye Wang (NYU); Yiqing Xu(Stanford); Ziyi Liu(UChicago); Shijian Liu(Duke)



 

  
  
 


 
