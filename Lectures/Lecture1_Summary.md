# BINP15 – Lecture 1: Introduction, measurement and summary statistics

**Course:** Processing and Analysis of Biological Data (BIOS15), Lund University  
**Lecture:** 1 – Introduction  
**Purpose of these notes:** A compact revision guide based on the Lecture 1 slides, Chapter 1 notes, the Appendix and the conceptual worksheet.

## 1. What is this course about?

The course introduces quantitative analysis for biology, including:

- Scientific programming and data handling in R.
- Statistical modelling and quantitative interpretation.
- Summary statistics and graphical presentation.
- Reproducible workflows and version control.

The course emphasizes **quantification, parameter estimation and biological interpretation**, rather than treating statistical hypothesis testing as the main goal. A result should describe the size and biological meaning of an effect, not only whether a p-value crosses a threshold.

## 2. Measurement and meaning in biology

Measurement connects the biological question to the numbers used in an analysis.

- A **biological concept/attribute** is what we want to understand (for example, attractiveness, fitness or stress).
- An **operational measurement** is what we actually record (for example, body length, offspring number or hormone concentration).
- A measurement is useful only to the extent that it represents the biological attribute relevant to the question.

A measurement can be precise and repeatable but still be a poor representation of the concept. Always ask:

1. What biological attribute am I interested in?
2. What exactly did I measure?
3. Why should this measurement represent that attribute?
4. What assumptions connect the measurement to the biological interpretation?

This is the central message of Houle et al. (2011), *Measurement and Meaning in Biology*.

## 3. Scale types: what comparisons and calculations are meaningful?

The scale type determines which transformations preserve meaning. The table below summarizes the four scale types introduced in the lecture notes.

| Scale type | What the values mean | Biological example | Meaningful comparisons | Usually meaningful summaries/operations | Not meaningful / avoid |
|---|---|---|---|---|---|
| **Nominal** | Different, unordered categories | Species identity; blood group; colour | Same/different; counts and proportions by category | Frequencies, proportions, mode; bar plots | Ranking categories; differences, means, SD, ratios |
| **Ordinal** | Ordered categories, but gaps between ranks are not known to be equal | Very rare / rare / common / dominant; severity score | Greater/lesser; ranking | Frequencies, proportions, mode; median or rank-based summaries where appropriate | Treating category codes as measured quantities; interpreting differences between ranks; arithmetic mean and SD as if intervals were equal |
| **Interval** | Numeric values with equal-sized differences, but an arbitrary zero | Day of year / days since January 1; temperature in °C | Order and differences | Differences; mean, variance and SD on the chosen interval scale | Ratios such as “twice as late”; CV or other ratio-to-mean interpretations; conclusions that depend on the arbitrary zero |
| **Ratio** | Numeric values with equal-sized differences and a natural zero | Length, mass, duration, count | Order, differences and ratios | Mean, variance, SD, ratios, CV; proportional comparisons | Negative values are not possible for the usual positive ratio-scale examples; transformations must preserve the measurement's meaning |

### Permissible transformations

- **Nominal:** any one-to-one relabelling of categories.
- **Ordinal:** any monotonically increasing transformation that preserves order.
- **Interval:** a positive linear transformation, `x_new = a * x + b` (`a > 0`); differences are preserved up to a constant factor, but the zero can move.
- **Ratio:** multiplication by a positive constant, `x_new = a * x` (`a > 0`); the natural zero remains fixed.

### Examples to remember

- Plant height: 40 cm is twice 20 cm (ratio scale).
- Day of year: day 40 is **not** meaningfully twice as late as day 20 (interval scale; zero is arbitrary).
- An abundance score coded rare = 1, common = 2, dominant = 3 is ordinal. A mean of 2.2 does not describe a meaningful abundance category or a measured distance between categories.
- Changing the origin of a date variable can change its mean and CV, even though the actual dates and differences between dates have not changed.

**Rule:** Do not choose a statistical operation just because software allows it. First establish what the numbers mean.

## 4. Describing the centre of data

| Statistic | Meaning | Useful when | Limitation |
|---|---|---|---|
| **Mean** | Arithmetic average | Values are quantitative and the distribution is reasonably symmetric | Sensitive to skew and extreme observations |
| **Median** | Middle value after sorting | Data are skewed or contain extreme observations | Does not use the magnitude of every observation in the same way as the mean |
| **Mode** | Most frequent value/category | Especially useful for nominal or ordinal categories | May be non-unique or uninformative for some continuous data |
| **Geometric mean** | Exponentiated mean of log values | Positive multiplicative quantities; one example discussed is long-term fitness | Requires positive values and a suitable multiplicative interpretation |

For a roughly symmetric distribution, mean and median are often similar. In a right-skewed distribution, the long upper tail tends to pull the mean above the median.

A genuine extreme observation is not automatically an error and should not be removed solely because it is extreme. Investigate its origin and retain it if it is a valid observation unless there is a justified analysis reason to handle it differently.

## 5. Describing variation and uncertainty

| Quantity | Formula (sample notation) | What it describes | Units |
|---|---|---|---|
| Variance | `s^2` | Dispersion around the mean; squared deviations | Squared units |
| Standard deviation (SD) | `s = sqrt(s^2)` | Spread of individual observations around the mean | Same units as the data |
| Standard error of the mean (SE) | `SE = SD / sqrt(n)` | Precision/uncertainty of the estimated mean | Same units as the data |
| Coefficient of variation (CV) | `CV = SD / mean` | SD relative to the mean; proportional variation | Dimensionless (often reported as %) |

Important distinctions:

- SD describes **variation among observations**, not the precision of the estimated mean.
- SE describes the uncertainty/precision of an estimate; it generally decreases as sample size increases.
- A population can have substantial biological variation but a precisely estimated mean if the sample is sufficiently large.
- Do not write `mean ± SD` if the intention is to communicate uncertainty in the mean. State clearly whether the quantity is SD or SE.
- Variances are additive in variance partitioning; standard deviations are not.

### Comparing variation: absolute vs proportional

A larger mean often comes with a larger SD. To compare variation relative to the mean, the lecture introduces:

`CV = SD / mean`

Example from the conceptual worksheet:

- Species A: mean = 10 mm, SD = 2 mm → CV = 0.20 (20%).
- Species B: mean = 100 mm, SD = 5 mm → CV = 0.05 (5%).

Species B has the larger absolute SD (5 mm), but Species A has the larger proportional variation (CV).

**Caution:** CV depends on a meaningful natural zero and is therefore not appropriate for interval-scale variables such as day of year. Shifting the arbitrary origin changes the mean and hence the CV.

## 6. Log transformations and effect sizes

For positive values, log transformations can express multiplicative differences as additive differences:

- `log(a / b) = log(a) - log(b)`.
- When two values are relatively close, a difference on the natural-log scale multiplied by 100 is approximately a percentage difference.
- The lecture notes that SD on a log scale is similar to SD divided by the mean (approximately CV) when variance is small relative to the mean.

Use a log transformation only when it is appropriate for the variable and biological question; it is not a universal fix for skewness.

## 7. P-values and interpretation

The lecture encourages attention to effect sizes, estimation and biological interpretation rather than a binary “significant / non-significant” conclusion.

Key points from the ASA statement (Wasserstein & Lazar, 2016) and Berner & Amrhein (2022):

- A p-value is not the probability that the null hypothesis is true.
- A p-value alone does not measure the size or biological importance of an effect.
- `p < 0.05` is not a magical boundary separating true from false findings.
- A non-significant result does not establish that there is no effect.
- Interpret estimates in context and consider their uncertainty (for example, confidence/compatibility intervals), study design and biological relevance.

Berner & Amrhein (2022) argue that p-values can vary substantially between studies and that dichotomous conclusions based on a default threshold can be misleading. Their recommendation is to describe and discuss a range of effect sizes compatible with the data, rather than reducing inference to a threshold decision.

## 8. Reproducible analysis and Git

A reproducible analysis should make it possible for another person (or your future self) to follow the steps from raw data to results.

Good habits introduced in the lecture and appendix:

- Use an RStudio **Project** for each analysis/course project. The project folder becomes the working directory.
- Keep data, scripts, results and figures in clearly named folders.
- Use descriptive file names; avoid spaces in variable names and use underscores (for example, `body_size_mm`).
- Keep raw data separate from processed data; document transformations.
- Prefer code-based imports (`read.csv()`, `read.table()`, or `readxl::read_excel()`) over repeated manual importing.
- Make missing values explicit as `NA`; understand `sep`, `dec` and `header` when importing text files.
- Separate data formatting, analysis/model fitting, plotting and reusable functions when a workflow becomes substantial.
- Use Git/GitHub for version history, backup and sharing code. A sequence of commits records what changed and when; it is more informative than files named `final`, `final2`, `really_final`.

## 9. R programming concepts from the Appendix

| Concept | Purpose | Example |
|---|---|---|
| Sequence | Create a sequence of values | `seq(1, 20, 2)` |
| Random normal data | Simulate data | `rnorm(n = 100, mean = 5, sd = 1)` |
| Reproducible randomness | Set random-number seed | `set.seed(1)` |
| `for` loop | Repeat an operation for each index | `for (i in seq_along(x)) { ... }` |
| `if` / `else` | Choose an operation based on a condition | `if (x[i] > 10) ... else ...` |
| `ifelse()` | Vectorized conditional operation | `ifelse(x > 10, 1, 0)` |
| `while` loop | Repeat until a condition is no longer true | `while (out < 21) { ... }` |
| Function | Package a repeatable operation | `computeSE <- function(x) { ... }` |

For a standard error that handles missing values:

```r
computeSE <- function(x) {
  n <- sum(!is.na(x))
  sqrt(var(x, na.rm = TRUE) / n)
}
```

## 10. Lecture 1 checklist

- [ ] I can distinguish a biological concept from the measurement used to represent it.
- [ ] I can classify variables as nominal, ordinal, interval or ratio.
- [ ] I can explain why some transformations and comparisons are not meaningful for a given scale.
- [ ] I can choose between mean, median and mode based on the variable and distribution.
- [ ] I can distinguish SD (variation) from SE (precision of an estimate).
- [ ] I can calculate and interpret CV, and explain why it is unsuitable for interval-scale variables.
- [ ] I can explain why p-values should not be interpreted as a binary verdict.
- [ ] I can describe how an R Project and Git support reproducibility.

## Sources

Course materials supplied for Lecture 1:

- Opedal, Ø. H. **BIOS15 Lecture 1: Introduction** (2026).
- Opedal, Ø. H. **Chapter 1: Preface and summary statistics** (2026).
- Opedal, Ø. H. **Appendix** (2026).
- BIOS15 **Lecture 1 conceptual worksheet: Measurement, variation, and reproducible analysis**.
- Houle, D., Pélabon, C., Wagner, G. P. & Hansen, T. F. (2011). Measurement and meaning in biology. *The Quarterly Review of Biology*, 86, 3–34.
- Wasserstein, R. L. & Lazar, N. A. (2016). The ASA Statement on p-Values: Context, Process, and Purpose. *The American Statistician*, 70, 129–133.
- Berner, D. & Amrhein, V. (2022). Why and how we should join the shift from significance testing to estimation. *Journal of Evolutionary Biology*, 35, 777–787.
