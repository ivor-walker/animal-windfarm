library(MuMIn)
library(lmtest)
library(nlme)
library(dplyr)
library(callr)

#Analysis of EIA data
eia_data <- read.csv("EIA.csv")

#Create a new column for the density and turn tidestate into a factor.
eia_data$density <- eia_data$count / eia_data$area
eia_data$tidestate <- as.factor(eia_data$tidestate)



#Fitting Multiple Covariate Linear Models
#Begin by fitting the two multiple covariate linear models below.
#month as continuous
fit.full <- lm(density ~ tidestate + observationhour + DayOfMonth + MonthOfYear + impact + Year + x.pos + y.pos, data = eia_data)
#month as a factor
fit.full.fac <- lm(density ~ tidestate + observationhour + DayOfMonth + as.factor(MonthOfYear) + impact + Year + x.pos + y.pos, data = eia_data)


#2.How many coefficients does changing month to a factor add to the model?
#10
no_extra_coefficients <- length(coef(fit.full.fac)) - length(coef(fit.full))


#Use AIC to determine whether month is preferred as a continuous or discrete term.
#Define a VIF function courtesy of chatgpt because I can't import car on CS lab
vif <- function(model) {
  # Get the design matrix of the model (excluding the intercept)
  X <- model.matrix(model)[, -1]  # Remove intercept
  
  # Calculate VIF for each predictor
  vif_values <- numeric(ncol(X))
  names(vif_values) <- colnames(X)
  
  for (i in seq_along(vif_values)) {
    # Regress each predictor on all other predictors
    predictor <- X[, i]
    other_predictors <- X[, -i]
    
    # Fit a linear model of the predictor against the other predictors
    r_squared <- summary(lm(predictor ~ other_predictors))$r.squared
    
    # Calculate VIF: VIF = 1 / (1 - R^2)
    vif_values[i] <- 1 / (1 - r_squared)
  }
  
  return(vif_values)
}

#Use the VIF function on the fit.full.fac model, to assess collinearity
#Collinear: impact, year
vif(fit.full.fac)

#Use AIC to determine whether month is preferred as a continuous or discrete term.
#AIC: factor > continuous
AIC(fit.full, fit.full.fac)


#3. Which of the following statements about collinearity is false?
#It is appropriate to assess models with and without the collinear variables and use AIC score to choose the best model. In this instance the preferred model is one with impact removed.
#TRUE: fit.fullfac.noimp has the LOWEST aic which makes it the best model
fit.fullfac.noyr <- update(fit.full.fac, .~. - Year)
fit.fullfac.noimp <- update(fit.full.fac, .~. - impact)
fit.fullfac.nocol <- update(fit.fullfac.noyr, .~. - impact)
AIC(fit.fullfac.noimp, fit.full, fit.full.fac, fit.fullfac.noyr, fit.fullfac.nocol)

#A confidence interval for the Year coefficient is over 48 times wider than it would be for a model with no collinear variables.
#FALSE: confidence interval ratio is 6.546917, <48
confint_full <- confint(fit.full.fac, "Year")
confint_reduced <- confint(fit.fullfac.noimp, "Year")

width_full <- confint_full[2] - confint_full[1]
width_reduced <- confint_reduced[2] - confint_reduced[1]

ratio <- width_full / width_reduced



#Adding Interaction Terms
#We would like to assess how density changes with each year of the study. For this reason, we ignore the suggestion of the VIF analysis and instead remove the impact term from the model. This is actually the preferred model if AIC is used to decide.

#Update the new reduced model fit.fullfac.noimp to include interaction terms between year and the x and y coordinates and name it fit.interac.
fit.interac <- update(fit.fullfac.noimp, . ~ . + Year:x.pos + Year:y.pos)

#Model Selection


#4. True or False? The null hypothesis for the F-test is that a model with a particular covariate included is no better than a model with that covariate removed. In this case, all variables except DayOfMonth and Year:y.pos appear to have significant relationships (at the 5% level) with density.
#FALSE: x:pos is insignificant too
drop1(fit.interac, test = "F")


#Perform a stepwise automated selection using AIC on the model fit.interac. Use the step function with direction = 'both'
fit.interac.step.both <- step(fit.interac, direction = "both")

#Now do all possible subsets selection using the dredge function and the default, AICc.
#Note that the dredge function is in the MuMIn library and that you will need to run the code options(na.action = 'na.fail') before using the dredge function.
options(na.action = "na.fail")
fit.interac.dredge <- dredge(fit.interac)
fit.interac.dredge.best <- get.models(fit.interac.dredge, 1)[[1]]


#5. Which of the following statements about model selection is FALSE?
#The best stepwise model is the same as that from all possible subsets
#TRUE
formula(fit.interac.step.both)
formula(fit.interac.dredge.best)

#Owing to the large sample size, it makes no difference to covariate selection whether we do all possible subsets selection using AIC or AICc.
#TRUE: AICc difference is minimal (<2 points)
AICc(fit.interac, fit.interac.dredge.best)

#The best all possible subsets model contains only the observation hour, x position, y position, and year covariates.
#FALSE: as.factor(MonthOfYear), tideState and x.pos:Year also included
formula(fit.interac.dredge.best)

#The model with 'y.pos:Year' included is within 2 AIC points of the model with it not included and the weights of the two models are almost equal (~25% weight compared with ~25% weight).
#TRUE: no idea how to compare AIC weights from this, but AIC is within two points
fit.interac.noypos <- update(fit.interac, . ~ . - Year:y.pos)
AIC(fit.interac, fit.interac.noypos)


#6. True or False? Using BIC for stepwise selection does not change the covariates selected in the final model.
#TRUE: both includes dayOfMonth and x.pos
n <- nrow(eia_data)
fit.interac.step.bic <- step(fit.interac, direction="both", k=log(n))
formula(fit.interac.step.bic)
formula(fit.interac.step.both)


#Using the model with interaction terms, conduct F-tests to assess covariate significance and compare this with the output from all possible subsets (default settings) and forwards and backwards stepwise selection (AIC, direction = both) on the same model.
drop1(fit.interac.step.both)


#7. Which of the following statements about model selection is FALSE?
#Hypothesis testing (F-test) suggests that the model without the year:y.pos interaction is preferred to the full model.
#FALSE
anova(fit.interac.noypos, fit.interac)

#Forwards and backwards stepwise selection using AIC chooses the full model with the year:x.pos interaction term retained.
#TRUE
formula(fit.interac.step.both)

#Dredge, using AICc, suggests that there is little to choose between a model a) without either interaction term, b) with both interaction terms and c) with only the x.pos interaction term retained.
#TRUE
head(fit.interac.dredge, n=10)

#Forwards and backwards stepwise selection using BIC returns the same model as backwards selection using hypothesis testing.
#FALSE: hypothesis testing includes tidestate, day of month, month of year and pos:year interactions whereas BIC includes pos
drop1(fit.interac, test = "F")
formula(fit.interac.step.bic)


#9. What is the estimate of the error variance for this model? Give your answer to two decimal places.
fit.interac.step.summary <- summary(fit.interac.step.both)
RSS <- sum(fit.interac.step.summary$residuals ^ 2)
df <- fit.interac.step.summary$df[2]
RSS / df



#Linear Model Assumptions
#11. Which of the following, about the validity of assumptions for the AIC-based stepwise selection model with interaction terms is TRUE?
#Can't import ncvTest, but ncvtest is just a wrapper around the breusch-pagan test so use that instead 
bptest(fit.interac.step.both)
#Small p-value, which is evidence of non-constant residual variance



#Generalised Least Squares Modelling
#Fit two GLS models with an exponential mean-variance relationship and one with a power based mean-variance relationship.
eia_data$sqrtdensity <- sqrt(eia_data$density)
fit.gls <- gls(sqrtdensity ~ tidestate + observationhour + impact + x.pos + y.pos + MonthOfYear + impact:x.pos + impact:y.pos, 
               data = eia_data, 
               method = 'ML')
fit.gls.exp <- gls(sqrtdensity ~ tidestate + observationhour + impact + x.pos + y.pos + MonthOfYear + impact:x.pos + impact:y.pos,
                   data = eia_data, 
                   method = 'ML',
                   weights=varExp())


#12. True or False? The exponential model is a better representation of the error variance than a constant.
#TRUE
AIC(fit.gls, fit.gls.exp)
anova(fit.gls, fit.gls.exp)


#13. Save your mean-variance plot and upload to Moodle. Use sensible axis labels and give your plot a title. Your file should be one of jpeg, png or pdf. 
fit.gls.exp.fitted <- fitted(fit.gls.exp)
fit.gls.exp.residuals <- residuals(fit.gls.exp)
fit.gls.exp.residuals_response <- residuals(fit.gls.exp, type = "response")

#Bin the fitted values
cut.fitted <- cut(fit.gls.exp.fitted, breaks = quantile(fit.gls.exp.fitted, probs = seq(0, 1, length = 20)))
#Calculate the variance of the residuals in each bin
cut.fitted.variance <- cut.variance <- tapply(fit.gls.exp.residuals_response, cut.fitted, var)
#Plot the mean fitted value in each bin on the x-axis and the variance of residuals in each bin on the y-axis
cut.fitted.mean <- tapply(fit.gls.exp.fitted, cut.fitted, mean)
plot(
  cut.fitted.mean, cut.fitted.variance,
  xlab = "Mean Fitted Value",
  ylab = "Variance of Residuals",
  main = "Mean-Variance Relationship",
  pch = 16, col = "blue"
)

#Estimate the mean-variance relationship itself using an exponential GLS model
exp.mean.variance <- gls(
  log(cut.fitted.variance) ~ cut.fitted.mean
)
coef_a <- exp(coef(exp.mean.variance)[1])
coef_b <- coef(exp.mean.variance)[2]

#Bin the fitted values of this relationship
mean_seq <- seq(min(cut.fitted.mean), max(cut.fitted.mean), length.out = 100)
predicted_variance <- coef_a * exp(coef_b * mean_seq)
lines(mean_seq, predicted_variance, col = "red", lwd = 2)


#14. Which of the following about the mean-variance relationship is FALSE?
#The exponential model slightly underestimates the variance for the smaller predicted root-density values.
#FALSE:



#Dealing With Correlated Errors
#We have updated the mean-variance relationship but still find some correlation in model residuals. 
#Use an acf plot to visualise this:
par(mfrow = c(1,2))
acf(residuals(fit.gls.exp, type = 'response'))
acf(residuals(fit.gls.exp, type = 'normalized'))
par(mfrow = c(1,1))
#These acf plots should look identical as we have not dealt with any correlation yet.

#Update your GLS model to include an AR(1) correlation matrix with gridcode/day as a blocking structure. 
#Note that you will need to use the new dataset (created below) 
library(dplyr)
eia_data$block <- paste(eia_data$Year, eia_data$MonthOfYear, eia_data$DayOfMonth, eia_data$GridCode, sep = '')
eia_data_2 <- arrange(eia_data, block, Year, MonthOfYear, DayOfMonth, GridCode)

fit_ar_result <- function(p, data) {
  library(nlme)
  
  return (gls(
    sqrtdensity ~ tidestate + observationhour + impact + x.pos + y.pos + MonthOfYear + impact:x.pos + impact:y.pos,
    data = data,
    correlation = corARMA(p = p, q = 0, form = ~ 1 | block),
    method = "ML"
  ))
}

fit_ar <- function(p, data) {
    job <- r_bg(fit_ar_result, args = list(p = p, data = data))
    
    timer = 0;
    while (job$is_alive()) {
      Sys.sleep(1)
      timer = timer + 1;
      if(timer %% 60 == 0){
        print(paste(timer/60, "minutes passed."))
      }
    }
    
    print(paste("p = ", p, " finished."))
    return(job$get_result())
}

fit.glsexp.ar1 <- fit_ar(p = 1, data = eia_data_2)

#Having fitted an AR(1) model, also try an AR(2).
fit.glsexp.ar2 <- fit_ar(p = 2, data = eia_data_2)



#15.Which of the following about GLS models is FALSE?
#The AR(2) model is the best model for the errors since the AIC score is the lowest.
#TRUE: AR(2) has lowest AIC
AIC(fit.glsexp.ar1, fit.glsexp.ar2, fit.gls.exp, fit.interac.dredge.best)

#The normalised residual acf plot for the AR(2) confirms the AIC result, that the AR(2) model fits better the correlation in the residuals than the AR(1) model.
#TRUE: ACF looks smaller across all lags in AR(2) vs AR(1)
par(mfrow = c(1,3))
acf(residuals(fit.gls.exp, type = 'normalized'), main = "No AR")
acf(residuals(fit.glsexp.ar1, type = 'normalized'), main = "AR(1)")
acf(residuals(fit.glsexp.ar2, type = 'normalized'), main = "AR(2)")
par(mfrow = c(1,1))

#The AR(2) model reduces the correlation at most lags to near zero and so has dealt with the issues regarding correlation in model residuals. 
#We can expect that the standard errors for the estimated coefficients have been reduced accordingly and we can trust any model selection results that use hypothesis testing.
#FALSE: Standard errors do increase, but lots of coefficients non zero
summary(fit.glsexp.ar1)
summary(fit.glsexp.ar2)


#Use hypothesis testing (F-tests) for backwards model selection and answer the following question.
#16. After backwards selection using hypothesis testing, select all the variables that remain in your model.
#intercept, tidestate, observationhour, impact, xpos
anova(fit.glsexp.ar2)
fit.ar2.reduced <- gls(
  sqrtdensity ~ tidestate + observationhour + x.pos + impact,
  correlation = corARMA(p = 2, q = 0, form = ~ 1 | block),
  method = "ML",
  data = eia_data_2
)

AIC(fit.glsexp.ar2, fit.ar2.reduced)
#17. Make a prediction from your best model for both before and after impact using the relevant covariate values given below. 
#Give your answers in density and to 2 decimal places.
finalmodel <- fit.glsexp.ar2

tidestate <- "SLACK"
tidestate_levels <- c("EBB", "FLOOD", "SLACK")
observationhour <- 10
MonthOfYear <- 6
x.pos <- 1500
y.pos <- 1000

preImpact <- data.frame(
  tidestate = factor(tidestate, levels = tidestate_levels),
  observationhour = observationhour,
  MonthOfYear = MonthOfYear,
  x.pos = x.pos,
  y.pos = y.pos,
  impact = 0
)

postImpact <- data.frame(
  tidestate = factor(tidestate, levels = tidestate_levels),
  observationhour = observationhour,
  MonthOfYear = MonthOfYear,
  x.pos = x.pos,
  y.pos = y.pos,
  impact = 1
)

preImpact_prediction <- MuMIn:::predict.gls(finalmodel, newdata = preImpact, se.fit = TRUE)$fit
postImpact_prediction <- MuMIn:::predict.gls(finalmodel, newdata = postImpact, se.fit = TRUE)$fit

#Squared to adjust for sqrt transformation
preImpact_prediction <- preImpact_prediction ^ 2
postImpact_prediction <- postImpact_prediction ^ 2

print(preImpact_prediction)
print(postImpact_prediction)


#Model diagnostics for report
finalmodel.residuals <- residuals(finalmodel)
finalmodel.residuals.normalised <- residuals(finalmodel, type="normalized")

acf(finalmodel.residuals.normalised)

plot(fitted(finalmodel), finalmodel.residuals,
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted Values")
abline(h = 0, col = "red")


shapiro.test(sample(finalmodel.residuals, size=100, replace=FALSE))

qqnorm(finalmodel.residuals)
qqline(finalmodel.residuals)


total_count <- nrow(eia_data_2)
zero_count <- sum(eia_data_2$density == 0)
zero_percentage <- (zero_count / total_count) * 100