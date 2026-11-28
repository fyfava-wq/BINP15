# BIOS15 — Multiple Regression and ANCOVA

**Course:** Processing and Analysis of Biological Data  
**Lecture material:** Chapter 4, *The Linear Model III: Multiple Regression and Analysis of Covariance* (Opedal, 2026); Lecture 3 slides, *ANCOVA and Multiple Regression* (Tsuboi & Opedal, 2026)  
**Reading:** Graham (2003), *Confronting Multicollinearity in Ecological Multiple Regression*

---

## 1. Where these models fit in the linear-model framework

The linear-model framework can accommodate different types and numbers of predictors:

| Analysis | Model | Predictors |
|---|---|---|
| Simple regression | `y ~ x` | One continuous |
| One-way ANOVA | `y ~ group` | One categorical |
| Two-way ANOVA | `y ~ A * B` | Two categorical |
| ANCOVA | `y ~ x * group` | One continuous and one categorical |
| Multiple regression | `y ~ x1 + x2` | Several continuous |

The response variable (`y`) is continuous in these examples.

## 2. Multiple regression

### What is it?

Multiple regression is a linear model with several continuous explanatory variables:

```math
y = \beta_0 + \beta_1x_1 + \beta_2x_2 + \dots + \beta_px_p + \varepsilon
```

- `y`: response variable.
- `x1`, `x2`, ...: explanatory variables (predictors).
- `β0`: intercept.
- `β1`, `β2`, ...: regression coefficients (slopes).
- `ε`: residual error.

Each coefficient estimates the **marginal effect** of its predictor: the expected change in the response for a one-unit change in that predictor, while the other predictors are held constant.

This makes it possible to estimate the independent association of several predictors with a response, or to control for other variables (for example, sampling effort).

```r
m <- lm(y ~ x1 + x2, data = df)
summary(m)
```

### Interpreting the model output

`summary(m)` reports:

- **Estimate:** estimated coefficient.
- **Std. Error:** uncertainty of the coefficient estimate.
- **t value and p-value:** test of whether an individual coefficient differs from zero, conditional on the other predictors in the model.
- **Residual standard error:** typical size of residuals, in response-variable units.
- **Multiple R-squared:** proportion of variation in the response explained by the model.
- **Adjusted R-squared:** R-squared adjusted for the number of predictors.
- **F-statistic:** test of the model against an intercept-only model.

A significant coefficient indicates evidence for an association after accounting for the other predictors in the model. Statistical significance does not by itself establish biological importance.

### Explained variation and shared contributions

The model's fitted values (`y_hat`) represent the variation explained by all predictors together. When predictors are correlated, their explained contributions overlap.

For two predictors, the variance of the combined fitted values includes a covariance term:

```math
\mathrm{Var}(\hat y_1+\hat y_2)
=
\mathrm{Var}(\hat y_1)
+
\mathrm{Var}(\hat y_2)
+
2\mathrm{Cov}(\hat y_1,\hat y_2)
```

Therefore, the variation explained by the predictors separately does not necessarily sum to the total variation explained by the full model.

The lecture demonstrates this by holding the other predictor at its mean when calculating the fitted values associated with one predictor:

```r
coefs <- summary(m)$coef

y_hat <- coefs[1, 1] +
  coefs[2, 1] * x1 +
  coefs[3, 1] * x2

y_hat1 <- coefs[1, 1] +
  coefs[2, 1] * x1 +
  coefs[3, 1] * mean(x2)

y_hat2 <- coefs[1, 1] +
  coefs[2, 1] * mean(x1) +
  coefs[3, 1] * x2

var(y_hat) / var(y)
var(y_hat1) / var(y)
var(y_hat2) / var(y)
```

The individual contributions can overlap because the predictors covary.

### Standardising predictors

Predictors measured in different units or with different variances can make raw slopes difficult to compare. A z-transformation gives a variable a mean of zero and a standard deviation of one:

```math
z_x = \frac{x-\bar{x}}{SD(x)}
```

```r
x1_z <- (x1 - mean(x1)) / sd(x1)
x2_z <- (x2 - mean(x2)) / sd(x2)

m_z <- lm(y ~ x1_z + x2_z)
summary(m_z)
```

- The model fit (including R-squared) does not change just because predictors are standardised.
- The slopes now describe the expected change in `y` per one-standard-deviation increase in the predictor.
- The intercept represents the expected response when the standardised predictors are zero, i.e. when they are at their means.

The lecture also shows mean-scaling:

```r
x1_m <- (x1 - mean(x1)) / mean(x1)
x2_m <- (x2 - mean(x2)) / mean(x2)

lm(y ~ x1_m + x2_m)
```

This expresses changes relative to the predictor means. In the lecture, these proportional slopes are described as **elasticities**.

## 3. Multicollinearity

### What is it?

Multicollinearity occurs when two or more explanatory variables are strongly correlated. It becomes difficult to separate their independent effects because they contain overlapping information.

Consequences include:

- Less precise coefficient estimates (larger standard errors).
- Reduced power to detect individual predictor effects.
- Coefficients that can change substantially when predictors are added or removed.
- Difficulty deciding which correlated predictor is biologically responsible for a shared association with the response.

A predictor can be associated with the response in a simple regression but appear non-significant in a multiple regression because the model estimates its effect after accounting for the other correlated predictors.

### Checking multicollinearity

The lecture gives a rule of thumb that correlations above approximately **0.6–0.7** may indicate a potential problem.

A common measure is the **variance inflation factor (VIF)**:

```math
VIF_i = \frac{1}{1-R_i^2}
```

Here, $R_i^2$ comes from a regression in which predictor $i$ is the response and all the other predictors are explanatory variables.

```r
# Example with two predictors
m1 <- lm(x1 ~ x2)
r2 <- summary(m1)$r.squared
1 / (1 - r2)
```

The course material notes that rules of thumb for concerning VIF values vary, from greater than 3 to greater than 10. Graham (2003) demonstrates that even VIF values around 2 can affect inference in ecological data. Treat thresholds as guidelines, not absolute rules.

### What can be done?

The course material discusses several possibilities:

- **Avoid redundant predictors:** simplify the model if several variables measure essentially the same property.
- **Use a priori models:** define the model based on biological hypotheses rather than selecting variables only by p-values.
- **Compare alternative models:** if correlated predictors are both biologically interesting, fit and compare models containing each candidate predictor.
- **Use model-selection approaches:** for example, all-subsets comparisons using adjusted R-squared, AIC, Mallows' Cp, PRESS, or MSE.
- **Use cross-validation:** assess whether a model predicts independent data well, rather than relying only on fit to the data used to build it.
- **Consider methods discussed by Graham (2003):** residual/sequential regression, principal components regression, and structural equation modelling.

Do not automatically remove a predictor just because it is correlated with another. Removing it can discard biologically meaningful information. The choice should be guided by the research question and biological reasoning.

### Overfitting

Overfitting occurs when a model includes too many predictors relative to the information in the data. It may explain a large amount of variation in the data used to fit it but perform poorly when predicting independent data.

A useful aim is to find a parsimonious model: a model that is sufficiently simple while still addressing the biological question.

The lecture contrasts:
- A statistical selection approach, such as backward selection from a full model.
- An a priori approach, where a biological hypothesis is used to specify the model before testing.

Backward selection starts with a full model and sequentially removes the least significant terms. The lecture cautions that this focuses on hypothesis testing and does not necessarily prioritise biological interpretation.

## 4. Analysis of Covariance (ANCOVA)

### What is ANCOVA?

ANCOVA is a linear model containing both:
- a **continuous predictor** (covariate), and
- a **categorical predictor** (group).

It combines ideas from regression and ANOVA.

A simple ANCOVA model with an interaction is:

```r
m <- lm(y ~ x * group, data = df)
```

The `*` includes both main effects and the interaction:

```r
y ~ x * group
# equivalent to:
y ~ x + group + x:group
```

### What do the terms test?

| Model term | Biological/statistical question | In a plot |
|---|---|---|
| `x` | Does the response change with the continuous predictor? | Slope |
| `group` | Do groups differ in intercept (at x = 0)? | Vertical position of lines |
| `x:group` | Does the relationship (slope) differ between groups? | Whether lines are parallel |

- **No interaction:** groups have parallel regression lines (common slope), but may have different intercepts.
- **Interaction:** slopes differ between groups; the effect of the continuous predictor depends on group.

The interpretation of a group difference at `x = 0` depends on whether zero is meaningful and within the observed data range. If zero is far outside the data, the intercept contrast may not be biologically useful. Centre the continuous predictor when a more meaningful reference value is needed.

### ANCOVA example from the lecture

The deer example asks whether male antler length is associated with female skull length and whether this relationship differs among breeding group sizes:

- Response: `logAL_x_M` (log mean male antler length).
- Continuous predictor: `logPSL_x_F` (log mean female skull length).
- Categorical predictor: `BGS` (breeding group size: A = 1–2, B = 3–5, C = 6 or more individuals).

The model sequence is:

```r
# Step 1: common relationship, without group
m0 <- lm(logAL_x_M ~ logPSL_x_F, data = deer)

# Step 2: add group differences in intercept
m1 <- lm(logAL_x_M ~ logPSL_x_F + BGS, data = deer)

# Step 3: allow slopes to differ between groups
m2 <- lm(logAL_x_M ~ logPSL_x_F * BGS, data = deer)

summary(m0)
summary(m1)
summary(m2)
anova(m2)
```

The interaction asks whether the allometric slope differs among breeding group sizes.

In the lecture example, the joint interaction test was:

```math
F_{2,30}=3.01,\quad p=0.0645
```

This is weak evidence for differences in slopes, but it is not statistically significant at the conventional 0.05 level.

### Interpreting ANCOVA coefficients

With the usual treatment contrasts, one group is the reference group:

- The intercept is the expected response for the reference group when `x = 0`.
- The coefficient for `x` is the slope for the reference group.
- A group coefficient is the difference in intercept from the reference group.
- An interaction coefficient is the difference in slope from the reference group.

For another group, its slope is the reference-group slope plus its interaction coefficient.

For example, in the lecture's model:

```r
m2 <- lm(logAL_x_M ~ logPSL_x_F * BGS, data = deer)
```

the slope for group B is the `logPSL_x_F` coefficient plus the `logPSL_x_F:BGSB` coefficient.

### Reparameterising the model to show group-specific lines

The same model can be written without a global intercept, giving a separate intercept and slope for each group:

```r
m3 <- lm(logAL_x_M ~ -1 + BGS + BGS:logPSL_x_F,
         data = deer)

summary(m3)
```

This is an alternative parameterisation of the same fitted model. The lecture checks that the log-likelihood is unchanged:

```r
logLik(m2)
logLik(m3)
```

This form can make group-specific intercepts and slopes easier to read directly.

## 5. Model summaries in R: `anova()` vs `summary()`

For a fitted linear model:

```r
m <- lm(y ~ x * group, data = df)

anova(m)
summary(m)
```

- `anova(m)` gives an analysis-of-variance table with degrees of freedom, sums of squares, mean squares, F-statistics and p-values for model terms.
- `summary(m)` gives coefficient estimates, standard errors, t-statistics and p-values for individual coefficients, plus residual standard error and R-squared.

For a model with interactions, the coefficient tests in `summary()` compare individual coefficients to the reference group. The ANOVA test for the interaction assesses the interaction term as a whole. These answer related but not identical questions.

## 6. Assumptions and diagnostics

The lecture identifies the following assumptions for ANCOVA:
- Residuals are approximately normally distributed within each group.
- Variances are equal between groups.

As with other linear models, inspect residuals and diagnostic plots rather than relying on a single test alone.

```r
par(mfrow = c(2, 2))
plot(m)
```

Useful plots include residuals versus fitted values and a normal Q–Q plot. Look for systematic patterns, strong changes in spread, and substantial departures from the expected residual distribution.

## 7. Exploring biological trait data

The flower-trait exercise uses measurements from nine natural populations in Costa Rica. The measured traits are:

| Abbreviation | Trait | Units |
|---|---|---|
| ASD | Anther–stigma distance | mm |
| GAD | Gland–anther distance | mm |
| GSD | Gland–stigma distance | mm |
| LBL | Lower bract length | mm |
| LBW | Lower bract width | mm |
| UBL | Upper bract length | mm |
| UBW | Upper bract width | mm |
| GW | Gland width | mm |
| GA | Gland area | mm² |

The course notes describe the assumed functions:
- ASD is relevant to the ability to self-pollinate.
- GAD and GSD affect the fit of flowers to pollinators.
- Upper and lower bracts act as advertisements (similar to petals in other flowers).
- The gland produces the reward for pollinators.

The exercise recommends exploring distributions, outliers, correlations among traits, proportional variability, and differences between populations before fitting models.

```r
blossoms <- read.csv("datasets/blossoms/blossoms.csv")
names(blossoms)

# Mean upper bract width by population
tapply(blossoms$UBW, blossoms$pop, mean, na.rm = TRUE)
```

The exercise also suggests summarising traits by population and fitting linear models for relationships between traits, on both arithmetic and log scales.

## 8. Writing up results

A clear Methods section should state:
1. The biological question.
2. The response variable.
3. The predictors and their types.
4. The model(s) fitted and why the terms were included.

A clear Results section should report:
- The direction and size of the estimated effects.
- The relevant model test and uncertainty (for example, F-statistic, degrees of freedom and p-value).
- Model fit where relevant (for example, R-squared).
- Whether the interaction is supported, and what that means biologically.

Avoid interpreting a p-value as a measure of biological importance. Interpret coefficient estimates in the context of the variables and their units.

---

## 9. Quick R syntax reference

```r
# Multiple regression
lm(y ~ x1 + x2, data = df)

# ANCOVA without interaction (common slope, group intercepts)
lm(y ~ x + group, data = df)

# ANCOVA with interaction (group-specific slopes)
lm(y ~ x * group, data = df)

# Same interaction expanded
lm(y ~ x + group + x:group, data = df)

# Group-specific intercepts and slopes
lm(y ~ -1 + group + group:x, data = df)

# Model coefficient summary
summary(m)

# ANOVA table for model terms
anova(m)

# Residual diagnostic plots
par(mfrow = c(2, 2))
plot(m)

# R-squared
summary(m)$r.squared

# Variance inflation factor for x1 with x2 as the only other predictor
m1 <- lm(x1 ~ x2, data = df)
1 / (1 - summary(m1)$r.squared)

# Standardise a predictor
x_z <- (x - mean(x)) / sd(x)

# Mean-scale a predictor
x_m <- (x - mean(x)) / mean(x)
```

## 10. Key takeaways

- Multiple regression estimates the marginal association of each predictor while holding the other predictors constant.
- Correlated predictors share explanatory information; their individual effects can be difficult to separate.
- Multicollinearity can inflate standard errors and make coefficient estimates less reliable.
- Overfitting can produce a model that fits the current data but predicts poorly in new data.
- ANCOVA combines continuous and categorical predictors.
- In ANCOVA, the interaction tests whether slopes differ between groups.
- `anova()` tests model terms; `summary()` reports individual coefficient estimates and tests.
- Model choice should be guided by the biological question, not only by p-values.

## References and source notes

- Opedal, Ø. H. (2026). *Processing and Analysis of Biological Data: Chapter 4. The Linear Model III: Multiple regression and Analysis of Covariance*. Department of Biology, Lund University.
- Tsuboi, M. & Opedal, Ø. H. (2026). *Processing and Analysis of Biological Data: Lecture 3. ANCOVA and Multiple Regression*. BIOS15.
- Graham, M. H. (2003). Confronting multicollinearity in ecological multiple regression. *Ecology*, 84(11), 2809–2815.

**Source note:** This summary follows the uploaded course chapter and lecture slides. The multicollinearity section also summarises Graham (2003). The flower-trait section follows the uploaded course exercise and its trait diagram.
