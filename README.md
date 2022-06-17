# Power analysis of Poisson Regression with Bernoulli Distribution
  
# Aim of study

Power analysis is essential in the design of clinical or experimental study. G-power, PASS, and WebPower from R package are three main software for computing power, sample size, and effect size. In this study, we focus on the power analysis of Poisson Regression with Bernoulli Distribution. We find a slight difference for Poisson regression among the results from those software. We plan to discuss and find the reason for the difference and give some suggestions about power analysis for Poisson regression. We hope this study could help to understand these software's pros and cons and be careful when selecting software to calculate power and sample size. 

# Study definitions

In a simple Poisson regression model with only one covariate of Bernoulli distribution, the outcome variable Y is a discrete count variable, which assumes to follow Poisson distribution. The probability of Y with parameter of Poisson distribution "$\lambda$" and during the exposure time "t" is,    

$$P(Y=y|\lambda,t)=\frac{e^{-\lambda t}(\lambda t)^{y}}{y!}$$
The log rate of Poisson regression model is, 

$$log(\lambda)=b_0+b_1X$$

The definitions of this study as below:

Alternative: two-sided or one-side
$e^{b_0}$: the mean event rate assumed under $H_0$  
$e^{b_1}$: the relative increase of the event rate over the mean rate assumed under $H_0$  
Mean exposure: denote by t, the time unit during which the events are counted  
Parameter of Bernoulli distribution: denote by p, the probability when X=1  
Type I error: denote by $\alpha$, the probability that the test falsely rejects the null hypothesis  
Sample size: denote by N, the smallest sample size for designed power
Power: the probability that the test correctly rejects the null hypothesis
Effect size: the ratio of $\lambda$ under null hypothesis to alternative hypothesis.


# Methods for different software 

Based on the large sample approximation procedure for logistic regression by Whittemore(1982), Signorini(1991) gives a general approach to calculate power for Poisson regression with one covariate. This method considers $b_1$ as the parameter of interest, and uses Wald test $Z=\frac{\hat b_1}{Var(\hat b_1)}$. For a two-sided hypothesis, $H_0: b_1=0 \  vs. H_1: b_1\neq0$. Given the significant level of $\alpha$, power of $1-\beta$, assuming N is large enough to apply the asymptotic results, the formula for power analysis is,

$$N=\frac{\bigg(Z_{1-\alpha/2}\sqrt{V(b_1|H_0)}+Z_{1-\beta}\sqrt{V(b_1|H_1)} \bigg)^2}{te^{b_0}b_1^2}(1)$$
For Bernoulli distribution, the variance under the null hypothesis is,  
$$V(b_1|H_0)=\frac{1}{p(1-p)}$$
The variance under the alternative hypothesis is,  
$$V(\hat b_1|H_1)=\frac{1}{1-p}+\frac{1}{pe^{b_1}}$$

PASS uses this method. From the formula, we could derive each type of power analysis with R if we inputs other parameters. 

Unlike formula (1) uses the variance at the null hypothesis $b_1=0$, Demidenko(2007) gives a general formula for logistic regression derived under the assumption that the variance is evaluated at the MLE. Demidenko states that Whittemore(1982) does not entirely reflect the way Wald tests are used in practice. Specifically, the formula from Signorini(1991) is not accurate since the null hypothesis is the unknown parameter takes value of 0, assuming $b_1=0$. Demidenko shows the test statistics is $Z=\frac{\sqrt{n}(\hat b_1-0)}{\sqrt{\hat V}}$, $V$ is the asymptotic variance of $\sqrt{n}\hat b_1$ under the alternative hypothesis $H_1: b_1 \neq0$. So the formula for Poisson regression is,

$$N=\frac{\big(Z_{1-\alpha/2}+Z_{1-\beta}\big)^2V}{tb_1^2}(2)$$
For Bernoulli distribution, 
$$V=\frac{1}{(1-p)e^{b_0}}+\frac{1}{pe^{b_0}e^{b_1}}$$

Webpower from R package use this method with function wp.Poisson for calculate power and sample size for Poisson regression. In this study, we expand the range of this function and could calculate effect size with the same statistic approximation.   

Moreover, based on method of Demidenko(2007), G-power adds a variance correction to compensate the variance distortions from skewed X distribution.   
With variance correction, instead of using null hypothesis $b_1=0$, we use the null hypothesis with $b_1=b^*$. The formula for Poisson regression is, 

$$N=\frac{\bigg(Z_{1-\alpha/2}\sqrt{V(\hat 
b_1|H_1)}+Z_{1-\beta}\sqrt{V(\hat b_1|H_0^*)} \bigg)^2}{tb_1^2}(3)$$
For Bernoulli distribution,
$$b^*=log(\int f(x)e^{b_0+b_1x}dx)=log(pe^{b_0}e^{b_1}+(1-p)e^{b_0})$$

$$V(\hat b_1|H_0^*)=\frac{1}{p(1-p)e^{b^*}}$$
$$var(\hat b_1|H_1)=\frac{1}{(1-p)e^{b_0}}+\frac{1}{pe^{b_1}e^{b_0}}$$

Besides the methods above, we also develop Monte Carlo simulation to estimate power of Poisson regression. Based on generating random data with specified sample size and parameter of distribution, we use R command "glm" on those randomly generated data and obtain the p-value for the test. We reject the null hypothesis when the p-value is less than our significance level. We repeat these steps for n times. Our estimate of power is the proportion of times that the null hypothesis is rejected. 


# R and Shiny App

For this study, we use R for programming. We didn't use specific R package for power analysis since it is hard to find an accuracy R package for Poisson regression. Because the code of G-power and PASS is not public, and WebPower can only use for power and sample size calculation. To get the results from simulation and different software, we code by ourselves based on approximation methods they use. 

For calculating effect size, we use uniroot function from R to solve nonlinear equation. For this function, we need to specify a finite interval for effect size. We set the effect size between 0-1000. If effect size exceeds this range, the uniroot may can not get a solution.  

In order to more intuitively show the running results of different software，we use Shinny APP to create an interactive interface for power analysis. The App could display power analysis by different software. With different parameter for each type of power analysis, we could get a table and a plot based on results from different software. Our code is public on Github (https://github.com/jessieliuyan/PSS_Poisson), run it with R code "runGitHub('PSS_Poisson','jessieliuyan', ref="main")".  

# Discussion  
In this study, we compare different software results of Poisson Regression with Bernoulli Distribution. We derive the formulas behind these software and implement the simulation method. We found the result from G-power and WebPower is very close, and the accuracy of power analysis is similar. The result from PASS is the most inaccurate and should be used carefully. We create Shiny App to display the results, which could help researchers analyze power when designing the study. 

This study is part of Power and Sample Size project. We only do the power analysis with a simple Poisson regression model with Bernoulli distribution for time limitation. We could explore more after this study. First, we could compare the difference between different distributions. We could consider Normal, Uniform, Exponential, Lognormal, etc. Second, we could consider more covariates and do a power analysis of multiple Poisson regression. Third, considering overdispersion of discrete count outcome, we could do the power analysis with quasi-Poisson regression or a negative binomial regression model. 

# Reference
Whittemore, A. S. (1981). Sample size for logistic regression with small response probabilities. Journal of the 84 American Statistical Association, 76, 27-32.    
Signorini, D. F. (1991). Sample size for Poisson regression.
Biometrika, 78, 446-450.    
Demidenko, E. (2007). Sample size determination for logistic regression revisited. Statistics in Medicine, 26, 3385-3397.  
G* Power 3.1 manual, January 21, 2021.  
R Package 'WebPower' manual, May 18, 2021.  
PASS sample size software, Chapter 870, Poisson Regression.  



