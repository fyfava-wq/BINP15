# BIOS15 — Generalized Linear Models I: Logistic Regression

**Course:** Processing and Analysis of Biological Data  
**Lecture:** Lecture 5 — GLM I: Logistic regression  
**Main course source:** Opedal, Ø. H. (2026), *Chapter 5. Generalized Linear Models (GLM)*  
**Additional sources:** Lecture 5 slides; Warton & Hui (2011), *The arcsine is asinine: the analysis of proportions in ecology*

---

## 1. Generalized linear models

Generalized linear models (GLMs) extend ordinary linear models by allowing response variables to follow distributions other than the normal distribution.

A GLM has two central components:

1. A **linear predictor**, combining explanatory variables and their coefficients.
2. A **link function**, which connects the linear predictor to the expected value of the response.

The linear predictor is commonly written as:

\[
\eta = \beta_0 + \sum_j \beta_j x_j
\]

The expected response is connected to the linear predictor through:

\[
g(\mu) = \eta
\qquad\text{or equivalently}\qquad
\mu = g^{-1}(\eta)
\]

The distribution and link function are chosen to suit the response data. The course introduces binomial errors for binary/proportional data and Poisson or negative-binomial errors for count data.

## 2. Binomial data

The binomial distribution describes the number of successes in a set of \(n\) Bernoulli trials. Each trial has two possible outcomes:

- Success, coded as 1.
- Failure, coded as 0.

The binomial distribution is described by:

- \(n\): number of trials.
- \(p\): probability of success.

For the number of successes, the expected value and variance are:

\[
E(X)=np
\qquad
\mathrm{Var}(X)=np(1-p)
\]

For the observed proportion \(\hat p=X/n\):

\[
E(\hat p)=p
\qquad
\mathrm{Var}(\hat p)=\frac{p(1-p)}{n}
\]

The variance depends on the probability and the number of trials. It is largest around \(p=0.5\), and smaller near 0 and 1. Therefore, proportional/binomial data do not generally meet the constant-variance assumption of an ordinary Gaussian linear model.

The number of trials matters: an observed proportion of 0.8 based on 50 trials contains more information than 0.8 based on only a few trials.

## 3. Logit transformation and inverse-logit

The **logit** (log odds) transformation is:

\[
\mathrm{logit}(p)=\log\left(\frac{p}{1-p}\right)
\]

It transforms probabilities between 0 and 1 to values spanning the whole real line.

The inverse-logit transforms a linear predictor back to a probability:

\[
p=\mathrm{logit}^{-1}(\eta)
=\frac{1}{1+\exp(-\eta)}
\]

Useful reference points:

| Linear predictor \(\eta\) | Probability \(p\) |
|---:|---:|
| Very negative | Approaches 0 |
| 0 | 0.5 |
| Very positive | Approaches 1 |

The logistic curve is S-shaped. It is steepest around \(p=0.5\) and flattens near 0 and 1.

```r
logit <- function(p) log(p / (1 - p))
invlogit <- function(eta) 1 / (1 + exp(-eta))

eta <- seq(-5, 5, by = 0.01)
plot(eta, invlogit(eta), type = "l",
     xlab = "Log odds", ylab = "Probability")
```

## 4. Logistic regression

Logistic regression is a GLM with **binomial errors** and usually the **logit link**. It is suitable for binary outcomes and binomial counts/proportions.

For one continuous predictor:

\[
\mathrm{logit}(p)=\beta_0+\beta_1x
\]

The model is linear on the log-odds scale, not on the probability scale.

In R:

```r
m <- glm(y ~ x, family = binomial(link = "logit"), data = df)
summary(m)
```

For individual binary observations, `y` contains 0s and 1s.

### Interpreting coefficients

The coefficients returned by `summary()` are on the **link (logit) scale**.

- \(\beta_0\) is the log odds when \(x=0\).
- \(\beta_1\) is the change in log odds for a one-unit increase in \(x\).
- A positive \(\beta_1\) means that increasing \(x\) increases the odds, and therefore the probability, of success.
- A negative \(\beta_1\) means that increasing \(x\) decreases the odds and probability of success.

The coefficient is not a constant change in probability. Because the inverse-logit is nonlinear, the same change in log odds can correspond to different absolute probability changes depending on the starting probability.

Exponentiating a coefficient gives an odds ratio:

\[
\exp(\beta_1)
\]

For a one-unit increase in \(x\), the odds are multiplied by \(\exp(\beta_1)\).

### Back-transforming predictions

Predicted probabilities are often easier to interpret biologically than log odds.

```r
coefs <- summary(m)$coef

x_pred <- seq(min(df$x), max(df$x), length.out = 100)
eta_hat <- coefs[1, 1] + coefs[2, 1] * x_pred
p_hat <- invlogit(eta_hat)

plot(df$x, df$y, xlab = "x", ylab = "Observed outcome")
lines(x_pred, p_hat)
```

For a model with multiple predictors, use `predict()` to obtain predictions:

```r
pred <- predict(m, newdata = new_data, type = "response")
```

`type = "response"` returns predictions on the response scale (probabilities for a binomial-logit model).

### Predictor value corresponding to a probability of 0.5

A log odds of zero corresponds to \(p=0.5\). For a one-predictor model:

\[
0=\beta_0+\beta_1x
\qquad\Longrightarrow\qquad
x=-\frac{\beta_0}{\beta_1}
\]

This value can be a useful biological benchmark. In toxicology, for example, the dose corresponding to a 50% probability of death is called the LD50.

## 5. Binary observations and probabilities

A predicted probability of 0.70 does not mean an individual observation should equal 0.70. Individual binary observations remain 0 or 1.

For 100 similar individuals with predicted probability 0.70, approximately 70 successes would be expected on average. The observed number will vary because the model describes a probability, not a deterministic outcome.

## 6. Model fit and discrimination

The ordinary \(R^2\) from Gaussian linear regression is not directly provided or automatically appropriate for logistic regression. Logistic regression is fitted using likelihood, and its response distribution and variance structure differ from Gaussian least squares.

The course chapter introduces two approaches:

### Pseudo-\(R^2\)

Several pseudo-\(R^2\) measures compare the fitted model with a null model (typically an intercept-only model). Different definitions can give different values, so state which measure is being reported.

The course demonstrates `r.squaredGLMM()` from the `MuMIn` package:

```r
library(MuMIn)
r.squaredGLMM(m)
```

### Tjur's coefficient of discrimination

Tjur's \(D\) compares the mean predicted probability for observations that were actually successes (1s) with the mean predicted probability for observations that were failures (0s):

\[
D=\overline{\hat p}_{\,y=1}-\overline{\hat p}_{\,y=0}
\]

A value near zero indicates little separation between the two groups. Larger values indicate better discrimination.

```r
p_hat <- predict(m, type = "response")

D <- mean(p_hat[y == 1]) - mean(p_hat[y == 0])
D
```

Discrimination is about how well predicted probabilities separate observed 1s from observed 0s. It is not the same as calibration (whether predicted probabilities agree with observed frequencies).

## 7. Binomial data with multiple trials per observation

When each row represents several trials, the number of trials must be included in the model.

The course presents two equivalent ways to fit a binomial GLM when the data are represented as a proportion and a trial count.

### Proportion response with weights

```r
glm(proportion ~ x,
    family = binomial(link = "logit"),
    weights = n,
    data = df)
```

Here `proportion` is the proportion of successes and `n` is the number of trials.

### Successes and failures as a two-column response

```r
glm(cbind(successes, failures) ~ x,
    family = binomial(link = "logit"),
    data = df)
```

The first column contains the number of successes and the second the number of failures.

Both formulations describe the same binomial information when specified consistently.

## 8. Proportions and the arcsine transformation

Warton and Hui (2011) argue that the arcsine square-root transformation should not be used for either binomial or non-binomial proportional data.

For binomial data, they recommend logistic regression because it:

- Ensures fitted proportions remain between 0 and 1.
- Gives coefficients with a natural interpretation on the odds scale.
- Models the binomial mean–variance relationship directly.
- Generally provides greater interpretability and statistical power than an arcsine-transformed linear model.

The arcsine transformation only approximately stabilises variance under particular conditions, including equal sample sizes, and can behave poorly near proportions of 0 and 1.

For binomial data, use a binomial model rather than transforming proportions and fitting an ordinary Gaussian model.

### Overdispersion in binomial models

Overdispersion means that the observed variation is greater than expected under the assumed binomial distribution.

Warton and Hui (2011) discuss checking for additional unexplained variation and accounting for it, for example with random effects in a mixed-effects logistic regression. The course chapter also notes that a random-intercept term can be used to model extra variation.

A mixed-effects logistic model can be fitted using `glmer()` from the `lme4` package, for example:

```r
library(lme4)

m_mixed <- glmer(
  cbind(successes, failures) ~ x + (1 | group),
  family = binomial(link = "logit"),
  data = df
)
```

The random effect shown here is an illustrative grouping structure; choose it based on the experimental design.

## 9. Course data exercise: seed dormancy and germination

The course dataset describes seed germination after different durations of dry after-ripening. Once seeds are ripe, they normally germinate under favourable moist conditions.

The data include four populations and the following variables:

| Variable | Meaning |
|---|---|
| `pop` | Population ID |
| `mother` | Maternal plant ID |
| `crossID` | Unique cross identifier |
| `blocktray` | Sowing tray / experimental block |
| `timetosowing` | Days from seed dispersal to watering |
| `MCseed` | Population-mean-centred seed mass (mg) |
| `nseed` | Number of seeds sown |
| `germ2` | Proportion of seeds germinated |

The exercise asks students to investigate germination success in relation to after-ripening duration, compare patterns among populations, and consider other factors that may affect germination.

```r
dat <- read.csv("dormancy.csv")
names(dat)

# Select one population
subdat <- dat[dat$pop == "CC", ]

# Calculate successes and failures
germ <- subdat$germ2 * subdat$nseed
notgerm <- subdat$nseed - germ

# Fit using successes and failures
mod1 <- glm(
  cbind(germ, notgerm) ~ timetosowing,
  family = binomial,
  data = subdat
)

# Fit using proportions and trial-count weights
mod2 <- glm(
  germ2 ~ timetosowing,
  family = binomial,
  weights = nseed,
  data = subdat
)

logLik(mod1)
logLik(mod2)
```

The two model formulations should give the same fit when the data are represented consistently.

The exercise asks for the after-ripening duration corresponding to an expected germination rate of 0.5. For a one-predictor logistic model, this is found by solving the linear predictor for a logit of zero.

## 10. Count data: Poisson regression

The course chapter also introduces count responses. Count data are often skewed, and their variance often increases with their mean.

The Poisson distribution has one parameter, \(\lambda\), with:

\[
E(X)=\lambda
\qquad
\mathrm{Var}(X)=\lambda
\]

A Poisson GLM commonly uses a log link:

\[
\log(\lambda)=\beta_0+\beta_1x
\]

In R:

```r
m_pois <- glm(y ~ x, family = "poisson", data = df)
summary(m_pois)
```

The coefficients are on the log scale. A one-unit increase in \(x\) changes the expected count multiplicatively by \(\exp(\beta_1)\).

To obtain fitted values on the response scale:

```r
predict(m_pois, newdata = data.frame(x = x_new),
        type = "response")
```

The course notes that log-transforming count data and fitting a Gaussian model can be problematic when there are zeros. Adding a constant before taking logs is generally not recommended as a default solution.

### Poisson overdispersion

The Poisson distribution assumes that the mean and variance are equal. Real count data may have a variance greater than the mean; this is called **overdispersion**.

The course suggests checking the ratio:

\[
\frac{\text{Residual deviance}}{\text{Residual degrees of freedom}}
\]

A value substantially greater than one can indicate overdispersion. When overdispersion is serious, standard Poisson model estimates and standard errors may not be reliable.

### Negative-binomial regression

The negative-binomial distribution adds a parameter that allows variance to increase more than the mean. It can be a suitable alternative when count data are overdispersed relative to a Poisson model.

```r
library(MASS)

m_nb <- glm.nb(y ~ x, data = df)
summary(m_nb)
```

## 11. Hurdle models

When count data contain many zeros, one approach described in the course is to separate the process into two components:

1. A binary model for whether the outcome is zero or non-zero.
2. A count model for the positive counts, conditional on a non-zero observation.

The two predictions can then be combined to obtain an expected count.

The course illustrates this with a binomial model for the occurrence component and a Poisson model for the positive counts. The choice of model should reflect the biological process that generated the zeros and counts.

---

## 12. Quick R reference

```r
# Logit and inverse-logit
logit <- function(p) log(p / (1 - p))
invlogit <- function(eta) 1 / (1 + exp(-eta))

# Logistic regression for individual 0/1 observations
glm(y ~ x, family = binomial(link = "logit"), data = df)

# Binomial model for proportions with trial-count weights
glm(prop ~ x, family = binomial, weights = n, data = df)

# Binomial model using successes and failures
glm(cbind(successes, failures) ~ x,
    family = binomial, data = df)

# Predicted probabilities
predict(m, newdata = new_df, type = "response")

# Tjur's coefficient of discrimination
p_hat <- predict(m, type = "response")
mean(p_hat[y == 1]) - mean(p_hat[y == 0])

# Pseudo-R2 measures demonstrated in the course
MuMIn::r.squaredGLMM(m)

# Poisson regression for counts
glm(count ~ x, family = "poisson", data = df)

# Negative-binomial regression for overdispersed counts
MASS::glm.nb(count ~ x, data = df)
```

## 13. Key takeaways

- GLMs connect a linear predictor to the expected response through a link function.
- Binomial models are appropriate for binary outcomes and counts of successes out of a known number of trials.
- Logistic regression uses binomial errors and the logit link.
- Logistic-regression coefficients describe changes in log odds, not constant changes in probability.
- Back-transformed predicted probabilities are often easier to communicate biologically.
- The number of trials matters when modelling observed proportions.
- Ordinary Gaussian \(R^2\) is not directly applicable; pseudo-\(R^2\) and Tjur's \(D\) are possible model-performance summaries.
- Check for overdispersion in binomial and count models.
- Poisson regression is used for counts; negative-binomial regression can accommodate overdispersion.
- Warton and Hui (2011) recommend logistic regression rather than arcsine-transformed linear models for binomial data.

## References

- Opedal, Ø. H. (2026). *Processing and Analysis of Biological Data: Chapter 5. Generalized Linear Models (GLM).* Department of Biology, Lund University.
- Opedal, Ø. H. (2026). *Processing and Analysis of Biological Data: Lecture 5. GLM I: Logistic regression.* BIOS15, Lund University.
- Warton, D. I. & Hui, F. K. C. (2011). The arcsine is asinine: the analysis of proportions in ecology. *Ecology*, 92(1), 3–10. https://doi.org/10.1890/10-0340.1

**Source note:** This summary follows the uploaded BIOS15 chapter and lecture slides. The section on arcsine transformations and proportional data summarises Warton & Hui (2011). Code examples are adapted from the course materials; illustrative variable names are identified in context.
