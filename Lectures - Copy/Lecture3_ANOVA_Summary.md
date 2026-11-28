# Chapter 3: Analysis of Variance (ANOVA)

> Study notes for BIOS15, *Processing and Analysis of Biological Data*. Based on the Chapter 3 reading and Lecture 3 slides.

## 1. What ANOVA does

Analysis of variance (ANOVA) is a linear-model approach for questions where one or more predictors are categorical variables (factors). Instead of estimating regression slopes for continuous predictors, the model estimates group means and contrasts between factor levels.

ANOVA partitions variation in a response into additive sums of squares (SS). In a one-way ANOVA, total variation is divided into variation among groups and residual variation within groups:

\[
SS_T = SS_{among} + SS_{within}
\]

The analysis asks whether the among-group variation is large relative to the within-group variation expected under a common-mean model. A statistically supported factor effect means that at least one group mean differs; it does not by itself identify which groups differ.

The course frames ANOVA within the linear model. Residuals are assumed to be normally distributed with mean zero and a common variance (\(\sigma^2\)); for ANOVA this means residuals are approximately normal within groups and group variances are equal.

## 2. One-way ANOVA: sums of squares and the table

For a factor with \(k\) levels and \(n\) observations:

- **Total SS (SST)**: sum of squared deviations of observations from the grand mean; \(SST/(n-1)\) is the sample variance.
- **Among-groups SS**: sum of squared deviations of group means from the grand mean, weighted by group size. This is the variation attributed to the factor.
- **Residual/within-groups SS (SSE)**: sum of squared deviations of observations from their own group means.
- **Degrees of freedom**: factor \(k-1\); residual \(n-k\); total \(n-1\).
- **Mean square (MS)**: SS divided by its degrees of freedom.
- **F statistic**: factor MS divided by residual MS. It compares among-group variation with residual variation.
- **P-value**: evidence against the null hypothesis that all group means are equal, under the model assumptions.

A typical one-way table has one row for the factor and one for residuals. The factor row's F-test is an omnibus test: a small p-value indicates evidence that not all means are equal, but follow-up contrasts or other comparisons are needed to describe specific differences.

For a one-factor model, the proportion of total variation explained by the factor is the model \(R^2\):

\[
R^2 = SS_{among}/SS_T = 1 - SSE/SST.
\]

## 3. Fitting and interpreting a one-way model in R

The chapter simulates body-size measurements for three food-level groups. The code below follows its example:

```r
set.seed(100)
groups <- as.factor(rep(c("Low", "Medium", "High"), each = 50))
x <- c(rnorm(50, 10, 3),
       rnorm(50, 13, 3),
       rnorm(50, 14, 3))

# Inspect distributions across factor levels
boxplot(x ~ groups, xlab = "", ylab = "Body size (g)")

# Fit the linear model and obtain the ANOVA table
m <- lm(x ~ groups)
anova(m)
summary(m)
```

In the chapter's output, the factor has 2 df (3 groups minus 1), and residuals have 147 df (150 observations minus 3 group means). The example reports factor SS = 319.97 and residual SS = 1200.43, giving total SS = 1520.40 and \(R^2\) about 0.210. Its F statistic is 19.591 (p = 2.866e-08), evidence that at least one mean differs.

R's default treatment contrasts use a reference level. In this example, the alphabetical first level, `High`, is the reference. The intercept estimates its mean; `groupsLow` and `groupsMedium` estimate differences from `High`. Their t-tests test each contrast against the reference, not every possible pairwise comparison. The factor-level order can be set to choose a different reference:

```r
groups <- factor(groups, levels = c("Low", "Medium", "High"))
m <- lm(x ~ groups)
summary(m)
```

To estimate a mean for each level directly (without an intercept), the chapter uses `-1`. `confint()` then gives confidence intervals for those estimated means:

```r
m_means <- lm(x ~ groups - 1)
summary(m_means)$coef
confint(m_means)
```

Contrasts also help quantify effect size. For example, the chapter compares the Low and High means: \(-3.456/13.701 \approx -0.252\), so the estimated Low-group mean is about 25.2% lower than the High-group mean. State the comparison and denominator when reporting a percentage.

## 4. Two-way ANOVA and interactions

With two categorical predictors, the model can estimate both main effects and their interaction. In R formula syntax:

```r
m2 <- lm(response ~ factor1 * factor2, data = dat)
anova(m2)
summary(m2)
```

`factor1 * factor2` expands to `factor1 + factor2 + factor1:factor2`: both main effects and the interaction. A colon (`factor1:factor2`) denotes only the interaction term.

An interaction means the effect of one factor depends on the level of the other factor (and equivalently, the effect of the other depends on the first). For example, if temperature and nitrogen supply are factors, an interaction means the temperature effect differs across nitrogen levels. This dependence is central to interpretation; a single main-effect contrast may not describe the pattern adequately when an interaction is present.

With treatment contrasts, the intercept in a two-way model is the estimated mean at the reference level of each factor. Main-effect coefficients are contrasts from that reference combination. The interaction coefficient is the additional effect when moving across both factors together. The fitted mean for the combination where both factors are non-reference is obtained by adding the intercept, the two relevant main-effect estimates, and their interaction estimate.

## 5. Course factorial example: butterfly hosts

The chapter's exercise uses an experiment in which female butterflies reared on two host plants were allowed to oviposit on the same two host plants. The dataset includes larval and maternal host, sex, mother identity, developmental time, adult weight, and growth rate.

For example, developmental time can be analysed by maternal and larval host as a two-factor model:

```r
dat <- read.csv("datasets/butterflies.csv")
dat$MaternalHost <- paste0(dat$MaternalHost, "M")
dat$LarvalHost <- paste0(dat$LarvalHost, "L")

m_dev <- lm(DevelopmentTime ~ MaternalHost * LarvalHost, data = dat)
anova(m_dev)
summary(m_dev)
```

The chapter first inspects cell means, for example with `tapply()`, and plots developmental time for each larval-host and maternal-host combination. The printed means are 21.69608, 27.00000, 23.51282, and 31.01923 days across the four combinations. The exercise asks students to assess maternal- and larval-host effects on one or more responses, interpret the results, and make a plot that illustrates the pattern. Choose the response and model terms to match the research question and experimental design.

## 6. Assumptions, design, and checks

- **Normal residuals within groups:** inspect residuals or group distributions. Boxplots can help assess distributions and potential outliers; normality of the errors is the model assumption.
- **Equal variances:** the course model assumes a common residual variance across groups.
- **Independent observations and appropriate experimental units:** randomized experiments can provide strong evidence for causality, but design matters. The lecture notes that analyses should help reduce pseudoreplication; repeated or clustered observations may require methods covered later in the course, such as mixed modelling.
- **Randomization and other sources of variation:** experiments aim to minimize variation in factors other than those under study. Statistical output cannot repair a design that confounds effects or treats non-independent observations as independent.
- **Interaction interpretation:** in a factorial model, inspect and explain the interaction before summarizing main effects as if they were constant across the other factor.

Plots should show the data where possible. Boxplots summarize quartiles, median, whiskers (by default extending to 1.5 times the interquartile range), and outlying points; a jittered point display can reveal the individual observations that boxplots hide.

## 7. Writing Methods and Results

### Methods

- Begin with the analytical aim, then give technical details.
- Describe model terms in words; include the model equation when useful.
- State the response, factors, and interaction considered. R syntax is increasingly acceptable as part of the description.

### Results

- Start with summary statistics and patterns of variation.
- Explain the result in biological terms and quantify the differences, rather than focusing only on test statistics.
- Refer to tables and figures after stating the finding.
- Report enough detail to make the claim clear: group or cell means, the relevant contrast or interaction, and the ANOVA test statistics as appropriate.
- Keep interpretation distinct from discussion; the lecture advises that Results normally do not include discussion.

### Figures and tables

- Give axes labels and units; choose a plot that shows the data and the comparisons clearly.
- Use a legend when multiple symbols or colours appear, and keep it inside the graph when practical.
- Avoid very wide or narrow panels when a squarer panel works.
- For tables, avoid vertical rules and explain parameters in the legend.

## 8. Quick interpretation checklist

1. What is the response, and which predictors are factors?
2. What are the experimental units, and are observations independent?
3. Do residual distributions and group variances look consistent with the assumptions?
4. For a one-way ANOVA, what are factor and residual SS, df, MS, F, and p?
5. Is the omnibus test supported? If so, which planned contrasts or comparisons answer the biological question?
6. What reference level is used for coefficient interpretation?
7. For two-way ANOVA, is there an interaction, and how does the effect vary across the other factor?
8. What are the size and direction of the differences in biological units?
9. Does the figure show the pattern and data clearly?

## Sources

- *Chapter 3. The Linear Model II: Analysis of Variance*, Øystein H. Opedal, Department of Biology, Lund University.
- *Lecture 3. ANOVA*, BIOS15, Øystein H. Opedal.
