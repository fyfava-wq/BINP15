# BIOS15 – Chapter 2: The Linear Model I — Linear Regression

**Course:** Processing and Analysis of Biological Data  
**Lecturer:** Øystein H. Opedal, Lund University  
**Lecture:** 2 · 6 November 2026

## 1. What is a linear model?

A linear model describes how one or more predictor variables are related to a response variable:

\[
y_i = \beta_0 + \sum_j \beta_j x_{ij} + \varepsilon_i
\]

| Term | Meaning |
|---|---|
| \(y_i\) | Observed response for observation \(i\) |
| \(x_{ij}\) | Value of predictor \(j\) for observation \(i\) |
| \(\beta_0\) | Intercept: expected response when predictors are zero |
| \(\beta_j\) | Coefficient/slope for predictor \(j\) |
| \(\varepsilon_i\) | Residual: difference between observed and model-predicted response |

For simple linear regression with one predictor:

\[
y_i = \beta_0 + \beta_1 x_i + \varepsilon_i
\]

- **Response variable (y):** dependent variable being explained or predicted.
- **Predictor variable (x):** independent variable used to explain or predict the response.
- **Linear regression** is used when the predictor is continuous.
- Regression can be used to investigate a potentially causal relationship or to predict unknown response values. A regression by itself does not establish causation.

## 2. Ordinary least squares (OLS)

Ordinary least-squares regression estimates the intercept and slope by minimising the **sum of squared residuals** (the squared vertical differences between observations and the fitted regression line).

\[
\hat{y}_i = \hat{\beta}_0 + \hat{\beta}_1 x_i
\]

\[
e_i = y_i - \hat{y}_i
\]

- A **predicted value** \(\hat y_i\) is the response expected from the fitted model.
- A **residual** \(e_i\) is the observed value minus the predicted value.
- Positive residual: observation is above the fitted line.
- Negative residual: observation is below the fitted line.

The slope in simple regression can also be expressed as:

\[
\hat{\beta}_1 = \frac{\operatorname{Cov}(x,y)}{\operatorname{Var}(x)}
\]

### Interpreting the coefficients

- **Intercept:** predicted value of \(y\) when \(x=0\). It may have little biological meaning if zero is outside the relevant range of the predictor.
- **Slope:** expected change in \(y\) for a one-unit increase in \(x\).
- The slope has units of **response units / predictor units**. For example, if leaf width and leaf length are both measured in mm, the slope is in mm/mm.

A useful biological interpretation can describe the predicted response to a realistic change in the predictor, such as a one-standard-deviation increase.

## 3. Assumptions and model diagnostics

The model assumes that the **residuals**, not the raw response or predictor data, are normally distributed:

\[
\varepsilon_i \sim N(0,\sigma^2)
\]

Minor deviations from normality are generally not a problem.

Useful diagnostic approaches shown in the lecture:

- Plot the observations and fitted regression line.
- Inspect residuals against fitted values or predictors.
- Inspect a histogram of residuals.
- Use `plot(model)` to obtain standard diagnostic plots.
- Look for residuals centred around zero and without strong patterns.

A residual histogram alone is not a complete assessment of model fit. Diagnostics should be considered alongside the study design and biological context.

## 4. Fitting a linear regression in R

The lecture uses simulated data where the underlying relationship is known:

```r
set.seed(85)

x <- rnorm(n = 200, mean = 10, sd = 2)
y <- 0.4 * x + rnorm(200, 0, 1)

plot(x, y,
     las = 1,
     xlab = "Leaf length (mm)",
     ylab = "Leaf width (mm)")

model <- lm(y ~ x)

# Inspect the model object and coefficients
str(model)
coef(model)

# Model summary
summary(model)

# Predicted values and residuals
predicted <- predict(model)
residuals(model)

# Plot fitted line and residuals
plot(x, y, las = 1,
     xlab = "Leaf length (mm)",
     ylab = "Leaf width (mm)")
abline(model)
segments(x, y, x, predicted)

# Histogram of residuals
hist(residuals(model), xlab = "Residuals", las = 1)

# Standard diagnostic plots
par(mfrow = c(2, 2))
plot(model)
par(mfrow = c(1, 1))
```

### Important R functions

| Function | Purpose |
|---|---|
| `lm(y ~ x)` | Fit a linear model |
| `summary(model)` | Display coefficient estimates, standard errors, test statistics, p-values and model fit |
| `coef(model)` or `model$coef` | Extract estimated coefficients |
| `str(model)` | Inspect the structure of an R object |
| `predict(model)` | Calculate fitted/predicted values |
| `residuals(model)` | Extract residuals |
| `abline(model)` | Add the fitted regression line to a plot |
| `segments()` | Draw line segments, e.g. from observations to fitted values |
| `hist()` | Plot a histogram |
| `seq()` | Generate a sequence of values |
| `lines()` | Add a line to an existing plot |

To draw a fitted line only over the observed range of x:

```r
newx <- seq(min(x), max(x), length.out = 200)
predy <- predict(model, newdata = data.frame(x = newx))

plot(x, y, las = 1,
     xlab = "Leaf length (mm)",
     ylab = "Leaf width (mm)")
lines(newx, predy)
```

## 5. Understanding `summary(lm())`

A model summary includes several types of information:

| Output | Interpretation |
|---|---|
| Residual quantiles | Distribution of residuals (minimum, quartiles, median, maximum) |
| Estimate | Estimated intercept or slope |
| Std. Error | Uncertainty/precision of the coefficient estimate |
| t value | Test statistic, calculated from estimate divided by its standard error |
| Pr(\>|t\|) | P-value for the coefficient test |
| Residual standard error | Typical residual variation, in response-variable units |
| Multiple R-squared | Proportion of variation in the response explained by the model |
| Adjusted R-squared | R-squared adjusted for the number of predictors |
| F-statistic | Test of the model against a model with no predictors (in the simple regression example, the slope-zero model) |

### Standard error and sample size

The standard error describes the precision of an estimated coefficient. In the lecture's explanation, increasing sample size reduces the standard error, all else being equal, and can therefore increase the test statistic.

Do not confuse:
- **Residual variation:** variation among observations around the fitted relationship.
- **Standard error:** uncertainty in an estimated parameter.

## 6. P-values: statistical support is not biological importance

A p-value is calculated under a null hypothesis (for example, that the true slope is zero). It describes how extreme the observed test statistic would be if that null hypothesis were true.

Important points from the lecture:

- P-values depend on effect strength, data variation and sample size.
- Larger samples can produce smaller p-values, even for small effects.
- With moderate sample sizes, biologically important effects may fail to reach statistical significance.
- Effects that are detectable in small samples can be upwardly biased.
- A p-value does not describe the size or biological importance of an effect.

**Main principle:** interpret parameter estimates, their standard errors and units, and their biological consequences. Do not focus on the p-value alone.

## 7. Interpreting and reporting a slope

In the simulated leaf example, the fitted slope was approximately:

\[
\hat{\beta}_1 = 0.433 \pm 0.035\ \text{mm/mm}
\]

Interpretation: predicted leaf width increased by approximately 0.43 mm for each 1 mm increase in leaf length.

The lecture also illustrates interpreting the effect over a one-standard-deviation increase in the predictor. In that example, predicted leaf width increased by approximately 0.86 mm per standard-deviation increase in leaf length.

When reporting a slope:
1. State the estimated effect.
2. Include its standard error.
3. Include units.
4. Explain what the estimate means biologically.
5. Where useful, translate it into a realistic change in the predictor.

## 8. Coefficient of determination: \(R^2\)

For a simple linear regression:

\[
R^2 = r^2
\]

where \(r\) is the Pearson correlation between x and y.

\(R^2\) describes the proportion of variation in the response variable explained by the fitted model:

\[
R^2 = \frac{\operatorname{Var}(\hat y)}{\operatorname{Var}(y)}
\]

For example, \(R^2 = 0.431\) means that **43.1% of the variation in y is explained by x** in this simple regression.

A high \(R^2\) does not, by itself, tell us that the slope is biologically large or important.

## 9. Log-transformed regression and allometry

Biological scaling relationships are often described by a power law:

\[
y = ax^b
\]

Taking logarithms gives a linear relationship:

\[
\log(y) = \log(a) + b\log(x)
\]

When both response and predictor are natural-log transformed:
- The intercept is \(\log(a)\).
- The slope \(b\) is an **allometric slope**.
- The slope can be interpreted as an elasticity: the percentage change in y associated with a percentage change in x.

For example, a slope of 0.55 means that a 10% increase in body mass is associated with approximately a 5.5% increase in brain mass.

### Example from the lecture: bird brain and body mass

The lecture analyses brain mass against body mass separately for male and female birds:

```r
birds <- read.csv("datasets/allometry/bird_allometry.csv")

males <- birds[birds$Sex == "m", ]
females <- birds[birds$Sex == "f", ]

model_males <- lm(log(brain_mass) ~ log(body_mass), data = males)
model_females <- lm(log(brain_mass) ~ log(body_mass), data = females)

summary(model_males)
summary(model_females)

hist(residuals(model_males))
```

Reported results in the lecture:

| Group | Allometric slope | Standard error | R-squared |
|---|---:|---:|---:|
| Males | 0.550 | 0.0099 | 0.8345 |
| Females | 0.575 | 0.0118 | 0.8307 |

Interpretation: brain mass increased less than proportionally with body mass. The female slope was slightly steeper in the separate models. The lecture notes that a model including sex and an interaction would be a way to estimate the slope difference directly.

## 10. Measurement error in the predictor

The standard linear model assumes that the predictor variable is measured without error. Measurement error in x can bias the estimated slope **downwards** (attenuation bias).

The lecture describes:
- Simulating data with a known relationship.
- Fitting the model without measurement error.
- Adding increasing amounts of error to x.
- Observing how the estimated slope changes.
- Comparing the estimated slope with a correction based on the reliability ratio.

The lecture gives the reliability ratio as:

\[
K = 1 - \frac{\sigma^2_{me}}{\sigma^2_x}
\]

and the corrected slope as:

\[
\beta' = \frac{\hat{\beta}}{K}
\]

where \(\sigma^2_{me}\) is the measurement-error variance and \(\sigma^2_x\) is the variance in the predictor.

The notes also ask what happens when error is added to y instead. This is an exercise to investigate rather than a result to assume.

## 11. Programming strategy

The lecture recommends building code step by step:

1. **Verbal code:** describe the task in ordinary language.
2. **Pseudocode:** list the functions and operations needed.
3. **R code:** implement the steps.

Start by performing the relevant operation once. Only then add indexing, loops or functions. This makes it easier to identify and debug errors.

## 12. Writing methods, results and figures

### Methods

- Start with the aim of the analysis before technical details.
- Describe the model terms in words.
- Include the model equation where useful.
- State transformations and how groups or subsets were analysed.

### Results

- Start with summary statistics or patterns of variation.
- Focus on biological meaning rather than statistical terminology.
- Quantify effects and give examples in meaningful units.
- Avoid unnecessary introductory sentences.
- Refer to tables and figures after making the relevant statement.
- Normally, do not include discussion or explanations of why the result occurred in the Results section.

### Figures

- Label axes and include units.
- Keep regression lines within the observed data range.
- Include a legend when several symbols or colours are used.
- Use clear, reasonably square plotting panels.

### Tables

- Avoid vertical lines.
- Explain all parameters in the table legend.

---

## Key takeaways

- A linear regression estimates the relationship between a response and one or more predictors.
- OLS estimates coefficients by minimising squared residuals.
- Residuals are observed minus predicted values; the normality assumption concerns residuals, not the raw data.
- Interpret slopes in their units and translate them into biological changes.
- Standard errors describe coefficient precision; p-values alone do not measure biological importance.
- \(R^2\) is the proportion of response variation explained by the model.
- Log-log regression turns a power-law relationship into a linear model; its slope is an elasticity.
- Measurement error in the predictor can bias the estimated slope downwards.
- Reproducible, readable code and biologically meaningful reporting are part of the analysis.

## Exercises from the chapter and lecture

- Draw a fitted regression line using the estimated intercept and slope.
- Identify the intercept, slope, predicted values and residuals on a scatterplot.
- Use non-parametric bootstrapping to estimate the standard error of a regression slope.
- Fit a linear regression to a real dataset with continuous response and predictor variables.
- Analyse bird brain–body allometry using log-transformed variables.
- Investigate how measurement error in x and y affects the estimated slope.
- Practise describing the analysis in Methods and communicating the biological result in Results.

---

*Study notes based on the Chapter 2 reading and Lecture 2 slides for BIOS15. The R code examples are adapted into consistent `<-` assignment style for readability.*
