# sdmTMB model description

## Other vignettes available

**If this vignette is being viewed on CRAN, note that many other
vignettes describing how to use sdmTMB are available on the
[documentation site](https://pbs-assess.github.io/sdmTMB/) under
[Articles](https://pbs-assess.github.io/sdmTMB/articles/).**

## Notation conventions

This appendix uses the following notation conventions, which generally
follow the guidance in Edwards & Auger-Méthé (2019):

- Greek symbols for parameters,

- the Latin/Roman alphabet for data (except $\mathbf{Q}$ and
  $\mathbf{H}$, which are used by convention),

- bold symbols for vectors or matrices (e.g., $\mathbf{ω}$ is a vector
  and $\omega_{\mathbf{s}}$ is the value of $\mathbf{ω}$ at point in
  space $\mathbf{s}$),

- $\phi$ for all distribution dispersion parameters for consistency with
  the code,

- ${\mathbb{E}}\lbrack y\rbrack$ to define the expected value (mean) of
  variable $y$,

- ${Var}\lbrack y\rbrack$ to define the expected variance of the
  variable $y$,

- a $^{*}$ superscript represents interpolated or projected values as
  opposed to values at knot locations (e.g., $\mathbf{ω}$
  vs. ${\mathbf{ω}}^{*}$), and

- where possible, notation has been chosen to match VAST (Thorson 2019)
  to maintain consistency (e.g., $\mathbf{ω}$ for spatial fields and
  ${\mathbf{ϵ}}_{t}$ for spatiotemporal fields).

We include tables of all major indices (Table 1) and symbols (Table 2).

| Symbol       | Description                                      |
|:-------------|:-------------------------------------------------|
| $\mathbf{s}$ | Index for space; a vector of x and y coordinates |
| $t$          | Index for time                                   |
| $g$          | Group                                            |

Table 1: Subscript notation

| Symbol                        | Code                 | Description                                                                     |
|:------------------------------|:---------------------|:--------------------------------------------------------------------------------|
| $y$                           | `y_i`                | Observed response data                                                          |
| $\mu$                         | `mu_i`               | Mean                                                                            |
| $\phi$                        | `phi`                | A dispersion parameter for a distribution                                       |
| $f$                           | `fit$family$link`    | Link function                                                                   |
| $f^{- 1}$                     | `fit$family$linkinv` | Inverse link function                                                           |
| $\mathbf{β}$                  | `b_j`                | Parameter vector                                                                |
| $\mathbf{X}$                  | `X_ij`               | A predictor model matrix                                                        |
| $O_{\mathbf{s},t}$            | `offset`             | An offset variable at point $\mathbf{s}$ and time $t$                           |
| $\omega_{\mathbf{s}}$         | `omega_s`            | Spatial random field at point $\mathbf{s}$ (knot)                               |
| $\omega_{\mathbf{s}}^{*}$     | `omega_s_A`          | Spatial random field at point $\mathbf{s}$ (interpolated)                       |
| $\zeta_{\mathbf{s}}$          | `zeta_s`             | Spatially varying coefficient random field at point $\mathbf{s}$ (knot)         |
| $\zeta_{\mathbf{s}}^{*}$      | `zeta_s_A`           | Spatially varying coefficient random field at point $\mathbf{s}$ (interpolated) |
| $\epsilon_{\mathbf{s},t}$     | `epsilon_st`         | Spatiotemporal random field at point $\mathbf{s}$ and time $t$ (knot)           |
| $\epsilon_{\mathbf{s},t}^{*}$ | `epsilon_st_A`       | Spatiotemporal random field at point $\mathbf{s}$ and time $t$ (interpolated)   |
| $\delta_{\mathbf{s},t}$       | `b_t`                | AR(1) or random walk spatiotemporal deviations (knot)                           |
| $\alpha_{g}$                  | `RE`                 | IID random intercept deviation for group $g$                                    |
| $\mathbf{\Sigma}_{\omega}$    | `-`                  | Spatial random field covariance matrix                                          |
| $\mathbf{\Sigma}_{\zeta}$     | `-`                  | Spatially varying coefficient random field covariance matrix                    |
| $\mathbf{\Sigma}_{\epsilon}$  | `-`                  | Spatiotemporal random field covariance matrix                                   |
| $\mathbf{Q}_{\omega}$         | `Q_s`                | Spatial random field precision matrix                                           |
| $\mathbf{Q}_{\zeta}$          | `Q_s`                | Spatially varying coefficient random field precision matrix                     |
| $\mathbf{Q}_{\epsilon}$       | `Q_st`               | Spatiotemporal random field precision matrix                                    |
| $\sigma_{\alpha}^{2}$         | `sigma_G`            | IID random intercept variance                                                   |
| $\sigma_{\epsilon}^{2}$       | `sigma_E`            | Spatiotemporal random field marginal variance                                   |
| $\sigma_{\omega}^{2}$         | `sigma_O`            | Spatial random field marginal variance                                          |
| $\sigma_{\zeta}^{2}$          | `sigma_Z`            | Spatially varying coefficient random field marginal variance                    |
| $\kappa_{\omega}$             | `kappa(0)`           | Spatial decorrelation rate                                                      |
| $\kappa_{\epsilon}$           | `kappa(1)`           | Spatiotemporal decorrelation rate                                               |
| $\rho$                        | `ar1_rho`            | Correlation between random fields in subsequent time steps                      |
| $\rho_{\gamma}$               | `rho_time`           | Correlation between time-varying coefficients in subsequent time steps          |
| $\mathbf{A}$                  | `A`                  | Sparse projection matrix to interpolate between knot and data locations         |
| $\mathbf{H}$                  | `H`                  | 2-parameter rotation matrix used to define anisotropy                           |

Table 2: Symbol notation, code representation (in model output or in
model template code), and descriptions.

## sdmTMB model structure

The complete sdmTMB model can be written as

$$\begin{aligned}
{{\mathbb{E}}\left\lbrack y_{\mathbf{s},t} \right\rbrack} & {= \mu_{\mathbf{s},t},} \\
\mu_{\mathbf{s},t} & {= f^{- 1}\left( \mathbf{X}_{\mathbf{s},t}^{main}{\mathbf{β}} + O_{\mathbf{s},t} + \alpha_{g} + \mathbf{X}_{\mathbf{s},t}^{tvc}{\mathbf{γ}}_{\mathbf{t}} + \mathbf{X}_{\mathbf{s},t}^{svc}\zeta_{\mathbf{s}} + \omega_{\mathbf{s}} + \epsilon_{\mathbf{s},t} \right),}
\end{aligned}$$

where

- $y_{\mathbf{s},t}$ represents the response data at point $\mathbf{s}$
  and time $t$;
- $\mu$ represents the mean;
- $f$ represents a link function (e.g., log or logit) and $f^{- 1}$
  represents its inverse;
- $\mathbf{X}^{main}$, $\mathbf{X}^{tvc}$, and $\mathbf{X}^{svc}$
  represent design matrices (the superscript identifiers ‘main’ = main
  effects, ‘tvc’ = time varying coefficients, and ‘svc’ = spatially
  varying coefficients);
- $\mathbf{β}$ represents a vector of fixed-effect coefficients;
- $O_{\mathbf{s},t}$ represents an offset: a covariate (usually log
  transformed) with a coefficient fixed at one;
- $\alpha_{g}$ represents random intercepts by group $g$,
  $\alpha_{g} \sim N\left( 0,\sigma_{\alpha}^{2} \right)$;
- $\gamma_{t}$ represents time-varying coefficients (a random walk),
  $\gamma_{t} \sim N\left( \gamma_{t - 1},\sigma_{\gamma}^{2} \right)$;
- $\zeta_{\mathbf{s}}$ represents spatially varying coefficients (a
  random field),
  $\zeta_{\mathbf{s}} \sim {MVN}\left( \mathbf{0},\mathbf{\Sigma}_{\zeta} \right)$;
- $\omega_{\mathbf{s}}$ represents a spatial component (a random field),
  $\omega_{\mathbf{s}} \sim {MVN}\left( \mathbf{0},\mathbf{\Sigma}_{\omega} \right)$;
  and
- $\epsilon_{\mathbf{s},t}$ represents a spatiotemporal component (a
  random field),
  $\epsilon_{\mathbf{s},t} \sim {MVN}\left( \mathbf{0},\mathbf{\Sigma}_{\epsilon} \right)$.

A single sdmTMB model will rarely, if ever, contain all of the above
components. Next, we will split the model to describe the various parts
in more detail using ‘$\ldots$’ to represent the other optional
components.

### Main effects

$$\begin{aligned}
\mu_{\mathbf{s},t} & {= f^{- 1}\left( \mathbf{X}_{\mathbf{s},t}^{main}{\mathbf{β}}\ldots \right)}
\end{aligned}$$

Within
[`sdmTMB()`](https://pbs-assess.github.io/sdmTMB/reference/sdmTMB.md),
$\mathbf{X}_{\mathbf{s},t}^{main}{\mathbf{β}}$ is defined by the
`formula` argument and represents the main-effect model matrix and a
corresponding vector of coefficients. This main effect formula can
contain optional penalized smoothers or non-linear functions as defined
below.

#### Smoothers

Smoothers in sdmTMB are implemented with the same formula syntax
familiar to mgcv (Wood 2017) users fitting GAMs (generalized additive
models). Smooths are implemented in the formula using `+ s(x)`, which
implements a smooth from
[`mgcv::s()`](https://rdrr.io/pkg/mgcv/man/s.html). Within these
smooths, the same syntax commonly used in
[`mgcv::s()`](https://rdrr.io/pkg/mgcv/man/s.html) can be applied,
e.g. 2-dimensional smooths may be constructed with `+ s(x, y)`; smooths
can be specific to various factor levels, `+ s(x, by = group)`; smooths
can vary according to a continuous variable, `+ s(x, by = x2)`; the
basis function dimensions may be specified, e.g. `+ s(x, k = 4)` (see
[`?mgcv::choose.k`](https://rdrr.io/pkg/mgcv/man/choose.k.html)); and
various types of splines may be constructed such as cyclic splines to
model seasonality, e.g. `+ s(month, bs = "cc", k = 12)`.

While mgcv can fit unpenalized (e.g., B-splines) or penalized splines
(P-splines), sdmTMB only implements penalized splines. The penalized
splines are constructed in sdmTMB using the function
[`mgcv::smooth2random()`](https://rdrr.io/pkg/mgcv/man/smooth2random.html),
which transforms splines into random effects (and associated design
matrices) that are estimable in a mixed-effects modelling framework.
This is the same approach as is implemented in the R packages gamm4
(Wood & Scheipl 2020) and brms (Bürkner 2017).

#### Linear break-point threshold models

The linear break-point or “hockey stick” model can be used to describe
threshold or asymptotic responses. This function consists of two pieces,
so that for $x < b_{1}$, $s(x) = x \cdot b_{0}$, and for $x > b_{1}$,
$s(x) = b_{1} \cdot b_{0}$. In both cases, $b_{0}$ represents the slope
of the function up to a threshold, and the product $b_{1} \cdot b_{0}$
represents the value at the asymptote. No constraints are placed on
parameters $b_{0}$ or $b_{1}$.

These models can be fit by including `+ breakpt(x)` in the model
formula, where `x` is a covariate. The formula can contain a single
break-point covariate.

#### Logistic threshold models

Models with logistic threshold relationships between a predictor and the
response can be fit with the form

$$s(x) = \tau + \psi\ \left\lbrack 1 + e^{- \ln{(19)} \cdot {(x - s50)}/{(s95 - s50)}} \right\rbrack^{- 1},$$

where $s$ represents the logistic function, $\psi$ is a scaling
parameter (controlling the height of the y-axis for the response;
unconstrained), $\tau$ is an intercept, $s50$ is a parameter controlling
the point at which the function reaches 50% of the maximum ($\psi$), and
$s95$ is a parameter controlling the point at which the function reaches
95% of the maximum. The parameter $s50$ is unconstrained but $s95$ is
constrained to be larger than $s50$.

These models can be fit by including `+ logistic(x)` in the model
formula, where `x` is a covariate. The formula can contain a single
logistic covariate.

### Spatial random fields

Spatial random fields, $\omega_{\mathbf{s}}$, are included if
`spatial = 'on'` (or `TRUE`) and omitted if `spatial = 'off'` (or
`FALSE`).

$$\begin{aligned}
\mu_{\mathbf{s},t} & {= f^{- 1}\left( \ldots + \omega_{\mathbf{s}} + \ldots \right),} \\
{\mathbf{ω}} & {\sim \operatorname{MVNormal}\left( \mathbf{0},\mathbf{\Sigma}_{\omega} \right),} \\
 & 
\end{aligned}$$ The marginal standard deviation of $\mathbf{ω}$ is
indicated by `Spatial SD` in the printed model output or as `sigma_O` in
the output of `sdmTMB::tidy(fit, "ran_pars")`. The ‘O’ is for ‘omega’
($\omega$).

Internally, the random fields follow a Gaussian Markov random field
(GMRF)

$${\mathbf{ω}} \sim {MVNormal}\left( \mathbf{0},\sigma_{\omega}^{2}\mathbf{Q}_{\omega}^{- 1} \right),$$
where $\mathbf{Q}_{\omega}$ is a sparse precision matrix and
$\sigma_{\omega}^{2}$ is the marginal variance.

### Spatiotemporal random fields

Spatiotemporal random fields are included by default if there are
multiple time elements (`time` argument is not `NULL`) and can be set to
IID (independent and identically distributed, `'iid'`; default), AR(1)
(`'ar1'`), random walk (`'rw'`), or off (`'off'`) via the
`spatiotemporal` argument. These text values are case insensitive.

Spatiotemporal random fields are represented by ${\mathbf{ϵ}}_{t}$
within sdmTMB. This has been chosen to match the representation in VAST
(Thorson 2019). The marginal standard deviation of ${\mathbf{ϵ}}_{t}$ is
indicated by `Spatiotemporal SD` in the printed model output or as
`sigma_E` in the output of `sdmTMB::tidy(fit, "ran_pars")`. The ‘E’ is
for ‘epsilon’ ($\epsilon$).

#### IID spatiotemporal random fields

IID spatiotemporal random fields (`spatiotemporal = 'iid'`) can be
represented as

$$\begin{aligned}
\mu_{\mathbf{s},t} & {= f^{- 1}\left( \ldots + \epsilon_{\mathbf{s},t} + \ldots \right),} \\
{\mathbf{ϵ}}_{\mathbf{t}} & {\sim \operatorname{MVNormal}\left( \mathbf{0},\mathbf{\Sigma}_{\epsilon} \right).}
\end{aligned}$$

where $\epsilon_{\mathbf{s},t}$ represent random field deviations at
point $\mathbf{s}$ and time $t$. The random fields are assumed
independent across time steps.

Similarly to the spatial random fields, these spatiotemporal random
fields (including all versions described below) are parameterized
internally with a sparse precision matrix ($\mathbf{Q}_{\epsilon}$)

$${\mathbf{ϵ}}_{\mathbf{t}} \sim {MVNormal}\left( \mathbf{0},\sigma_{\epsilon}^{2}\mathbf{Q}_{\epsilon}^{- 1} \right).$$

#### AR(1) spatiotemporal random fields

First-order auto regressive, AR(1), spatiotemporal random fields
(`spatiotemporal = 'ar1'`) add a parameter defining the correlation
between random field deviations from one time step to the next. They are
defined as

$$\begin{aligned}
\mu_{\mathbf{s},t} & {= f^{- 1}\left( \ldots + \delta_{\mathbf{s},t}\ldots \right),} \\
{\mathbf{δ}}_{t = 1} & {\sim \operatorname{MVNormal}\left( \mathbf{0},\mathbf{\Sigma}_{\epsilon} \right),} \\
{\mathbf{δ}}_{t > 1} & {= \rho{\mathbf{δ}}_{t - 1} + \sqrt{1 - \rho^{2}}{\mathbf{ϵ}}_{\mathbf{t}},\ {\mathbf{ϵ}}_{\mathbf{t}} \sim \operatorname{MVNormal}\left( \mathbf{0},\mathbf{\Sigma}_{\epsilon} \right),}
\end{aligned}$$ where $\rho$ is the correlation between subsequent
spatiotemporal random fields. The
$\rho{\mathbf{δ}}_{t - 1} + \sqrt{1 - \rho^{2}}$ term scales the
spatiotemporal variance by the correlation such that it represents the
steady-state marginal variance. The correlation $\rho$ allows for
mean-reverting spatiotemporal fields, and is constrained to be
$- 1 < \rho < 1$. Internally, the parameter is estimated as `ar1_phi`,
which is unconstrained. The parameter `ar1_phi` is transformed to $\rho$
with
$\rho = 2\left( {logit}^{- 1}\left( \texttt{𝚊𝚛𝟷\_𝚙𝚑𝚒} \right) - 1 \right)$.

#### Random walk spatiotemporal random fields (RW)

Random walk spatiotemporal random fields (`spatiotemporal = 'rw'`)
represent a model where the difference in spatiotemporal deviations from
one time step to the next are IID. They are defined as

$$\begin{aligned}
\mu_{\mathbf{s},t} & {= f^{- 1}\left( \ldots + \delta_{\mathbf{s},t} + \ldots \right),} \\
{\mathbf{δ}}_{t = 1} & {\sim \operatorname{MVNormal}\left( \mathbf{0},\mathbf{\Sigma}_{\epsilon} \right),} \\
{\mathbf{δ}}_{t > 1} & {= {\mathbf{δ}}_{t - 1} + {\mathbf{ϵ}}_{\mathbf{t}},\ {\mathbf{ϵ}}_{\mathbf{t}} \sim \operatorname{MVNormal}\left( \mathbf{0},\mathbf{\Sigma}_{\epsilon} \right),}
\end{aligned}$$

where the distribution of the spatiotemporal field in the initial time
step is the same as for the AR(1) model, but the absence of the $\rho$
parameter allows the spatiotemporal field to be non-stationary in time.
Note that, in contrast to the AR(1) parametrization, the variance is no
longer the steady-state marginal variance.

### Time-varying regression parameters

Parameters can be modelled as time-varying according to a random walk or
first-order autoregressive, AR(1), process. The time-series model is
defined by `time_varying_type`. For all types:

$$\begin{aligned}
\mu_{\mathbf{s},t} & {= f^{- 1}\left( \ldots + \mathbf{X}_{\mathbf{s},t}^{tvc}{\mathbf{γ}}_{\mathbf{t}} + \ldots \right),}
\end{aligned}$$ where ${\mathbf{γ}}_{\mathbf{t}}$ is an optional vector
of time-varying regression parameters and
$\mathbf{X}_{\mathbf{s},t}^{tvc}$ is the corresponding model matrix with
covariate values. This is defined via the `time_varying` argument,
assuming that the `time` argument is also supplied a column name.
`time_varying` takes a *one-sided* formula. `~ 1` implies a time-varying
intercept.

For `time_varying_type = 'rw'`, the first time step is estimated
independently:

$$\begin{aligned}
\gamma_{t = 1} & {\sim \operatorname{Uniform}( - \infty,\infty),} \\
\gamma_{t > 1} & {\sim \operatorname{Normal}\left( \gamma_{t - 1},\sigma_{\gamma}^{2} \right).}
\end{aligned}$$

In this case, the first time-step value is given an implicit uniform
prior. I.e., the same variable should not appear in the fixed effect
formula since the initial value is estimated as part of the time-varying
formula. The formula `time_varying = ~ 1` implicitly represents a
time-varying intercept (assuming the `time` argument has been supplied)
and, this case, the intercept should be omitted from the main effects
(`formula ~ + 0 + ...` or `formula ~ -1 + ...`).

For `time_varying_type = 'rw0'`, the first time step is estimated from a
mean-zero prior:

$$\begin{aligned}
\gamma_{t = 1} & {\sim \operatorname{Normal}\left( 0,\sigma_{\gamma}^{2} \right),} \\
\gamma_{t > 1} & {\sim \operatorname{Normal}\left( \gamma_{t - 1},\sigma_{\gamma}^{2} \right).}
\end{aligned}$$ In this case, the time-varying variable (including the
intercept) *should* be included in the main effects. We suggest using
this formulation, but leave the `'rw'` option so that legacy code works.

For `time_varying_type = 'ar1'`:

$$\begin{aligned}
\gamma_{t = 1} & {\sim \operatorname{Normal}\left( 0,\sigma_{\gamma}^{2} \right),} \\
\gamma_{t > 1} & {\sim \operatorname{Normal}\left( \rho_{\gamma}\gamma_{t - 1},\sqrt{1 - \rho_{\gamma}^{2}}\sigma_{\gamma}^{2} \right),}
\end{aligned}$$ where $\rho_{\gamma}$ is the correlation between
subsequent time steps. The first time step is given a mean-zero prior.

### Spatially varying coefficients (SVC)

Spatially varying coefficient models are defined as

$$\begin{aligned}
\mu_{\mathbf{s},t} & {= f^{- 1}\left( \ldots + \mathbf{X}_{\mathbf{s},t}^{svc}\zeta_{\mathbf{s}} + \ldots \right),} \\
{\mathbf{ζ}} & {\sim \operatorname{MVNormal}\left( \mathbf{0},\mathbf{\Sigma}_{\zeta} \right),}
\end{aligned}$$

where $\mathbf{ζ}$ is a random field representing a spatially varying
coefficient. Usually, $\mathbf{X}_{\mathbf{s},t}^{svc}$ would represent
a prediction matrix that is constant spatially for a given time $t$ as
defined by a one-sided formula supplied to `spatial_varying`. For
example `spatial_varying = ~ 0 + x`, where `0` omits the intercept.

The random fields are parameterized internally with a sparse precision
matrix ($\mathbf{Q}_{\zeta}$)

$${\mathbf{ζ}} \sim {MVNormal}\left( \mathbf{0},\sigma_{\zeta}^{2}\mathbf{Q}_{\zeta}^{- 1} \right).$$

### IID random or multi-level intercepts

Multilevel/hierarchical intercepts are defined as

$$\begin{aligned}
\mu_{\mathbf{s},t} & {= f^{- 1}\left( \ldots + \alpha_{g} + \ldots \right),} \\
\alpha_{g} & {\sim \operatorname{Normal}\left( 0,\sigma_{\alpha}^{2} \right),} \\
 & 
\end{aligned}$$

where $\alpha_{g}$ is an example optional “random” intercept—an
intercept with mean zero that varies by level $g$ and is constrained by
$\sigma_{\alpha}$. This is defined by the `formula` argument via the
`(1 | g)` syntax as in lme4 or glmmTMB. There can be multiple random
intercepts, despite only showing one above. E.g., `(1 | g1) + (1 | g2)`,
in which case they are assumed independent and uncorrelated from each
other.

### Offset terms

Offset terms can be included through the `offset` argument in
[`sdmTMB()`](https://pbs-assess.github.io/sdmTMB/reference/sdmTMB.md).
These are included in the linear predictor as

$$\begin{aligned}
\mu_{\mathbf{s},t} & {= f^{- 1}\left( \ldots + O_{\mathbf{s},t} + \ldots \right),}
\end{aligned}$$

where $O_{\mathbf{s},t}$ is an offset term—a **log transformed**
variable without a coefficient (assuming a log link). The offset is
*not* included in the prediction. Therefore, if `offset` represents a
measure of effort, for example, the prediction is for one unit of effort
(`log(1) = 0`).

## Observation model families

Here we describe the main observation families that are available in
sdmTMB and comment on their parametrization, statistical properties,
utility, and code representation in sdmTMB.

### Binomial

$$\operatorname{Binomial}(N,\mu)$$ where $N$ is the size or number of
trials, and $\mu$ is the probability of success for each trial. If
$N = 1$, the distribution becomes the Bernoulli distribution.
Internally, the distribution is parameterized as the [robust
version](https://kaskr.github.io/adcomp/group__R__style__distribution.html#gaecb5a18095a320b42e2d20c4b120f5f5)
in TMB, which is numerically stable when probabilities approach 0 or 1.
Following the structure of
[`stats::glm()`](https://rdrr.io/r/stats/glm.html), lme4, and glmmTMB, a
binomial family can be specified in one of 4 ways:

1.  the response may be a factor (and the model classifies the first
    level versus all others)
2.  the response may be binomial (0/1)
3.  the response can be a matrix of form `cbind(success, failure)`, or
4.  the response may be the observed proportions, and the `weights`
    argument is used to specify the Binomial size ($N$) parameter
    (`probabilty ~ ..., weights = N`).

Code defined [within
TMB](https://kaskr.github.io/adcomp/group__R__style__distribution.html#gaee11f805f02bc1febc6d7bf0487671be).

Example: `family = binomial(link = "logit")`

### Beta

$$\operatorname{Beta}\left( \mu\phi,(1 - \mu)\phi \right)$$ where $\mu$
is the mean and $\phi$ is a precision parameter. This parametrization
follows Ferrari & Cribari-Neto (2004) and the betareg R package
(Cribari-Neto & Zeileis 2010). The variance is
$\mu(1 - \mu)/(\phi + 1)$.

Code defined [within
TMB](https://kaskr.github.io/adcomp/group__R__style__distribution.html#ga5324c83759d5211c7c2fbbad37fa8e59).

Example: `family = Beta(link = "logit")`

### Gamma

$$\operatorname{Gamma}\left( \phi,\frac{\mu}{\phi} \right)$$ where
$\phi$ represents the Gamma shape and $\mu/\phi$ represents the scale.
The mean is $\mu$ and variance is $\mu \cdot \phi^{2}$.

Code defined [within
TMB](https://kaskr.github.io/adcomp/group__R__style__distribution.html#gab0e2205710a698ad6a0ed39e0652c9a3).

Example: `family = Gamma(link = "log")`

### Gaussian

$$\operatorname{Normal}\left( \mu,\phi^{2} \right)$$ where $\mu$ is the
mean and $\phi$ is the standard deviation. The variance is $\phi^{2}$.

Example: `family = Gaussian(link = "identity")`

Code defined [within
TMB](https://kaskr.github.io/adcomp/dnorm_8hpp.html).

### Lognormal

sdmTMB uses the bias-corrected lognormal distribution where $\phi$
represents the standard deviation in log-space:

$$\operatorname{Lognormal}\left( \log\mu - \frac{\phi^{2}}{2},\phi^{2} \right).$$
Because of the bias correction, ${\mathbb{E}}\lbrack y\rbrack = \mu$ and
${Var}\left\lbrack \log y \right\rbrack = \phi^{2}$.

Code defined [within
sdmTMB](https://github.com/pbs-assess/sdmTMB/blob/18a39eabc111e2179fa589f942c8820d87ad10df/src/utils.h#L47-L54)
based on the TMB [`dnorm()`](https://rdrr.io/r/stats/Normal.html) normal
density.

Example: `family = lognormal(link = "log")`

### Negative Binomial 1 (NB1)

$$\operatorname{NB1}(\mu,\phi)$$

where $\mu$ is the mean and $\phi$ is the dispersion parameter. The
variance scales linearly with the mean
${Var}\lbrack y\rbrack = \mu + \mu/\phi$(Hilbe 2011). Internally, the
distribution is parameterized as the [robust
version](https://kaskr.github.io/adcomp/group__R__style__distribution.html#gaa23e3ede4669d941b0b54314ed42a75c)
in TMB.

Code defined [within
sdmTMB](https://github.com/pbs-assess/sdmTMB/blob/18a39eabc111e2179fa589f942c8820d87ad10df/src/sdmTMB.cpp#L577-L582)
based on NB2 and borrowed from glmmTMB.

Example: `family = nbinom1(link = "log")`

### Negative Binomial 2 (NB2)

$$\operatorname{NB2}(\mu,\phi)$$

where $\mu$ is the mean and $\phi$ is the dispersion parameter. The
variance scales quadratically with the mean
${Var}\lbrack y\rbrack = \mu + \mu^{2}/\phi$(Hilbe 2011). The NB2
parametrization is more commonly seen in ecology than the NB1.
Internally, the distribution is parameterized as the [robust
version](https://kaskr.github.io/adcomp/group__R__style__distribution.html#gaa23e3ede4669d941b0b54314ed42a75c)
in TMB.

Code defined [within
TMB](https://kaskr.github.io/adcomp/group__R__style__distribution.html#ga76266c19046e04b651fce93aa0810351).

Example: `family = nbinom2(link = "log")`

### Poisson

$$\operatorname{Poisson}(\mu)$$ where $\mu$ represents the mean and
${Var}\lbrack y\rbrack = \mu$.

Code defined [within
TMB](https://kaskr.github.io/adcomp/group__R__style__distribution.html#gaa1ed15503e1441a381102a8c4c9baaf1).

Example: `family = poisson(link = "log")`

### Student-t

$$\operatorname{Student-t}(\mu,\phi,\nu)$$

where $\nu$, the degrees of freedom (`df`), is a user-supplied fixed
parameter. Lower values of $\nu$ result in heavier tails compared to the
Gaussian distribution. Above approximately `df = 20`, the distribution
becomes very similar to the Gaussian. The Student-t distribution with a
low degrees of freedom (e.g., $\nu \leq 7$) can be helpful for modelling
data that would otherwise be suitable for Gaussian but needs an approach
that is robust to outliers (e.g., Anderson *et al.* 2017).

Code defined [within
sdmTMB](https://github.com/pbs-assess/sdmTMB/blob/18a39eabc111e2179fa589f942c8820d87ad10df/src/utils.h#L37-L45)
based on the [`dt()`](https://rdrr.io/r/stats/TDist.html) distribution
in TMB.

Example: `family = student(link = "log", df = 7)`

### Tweedie

$$\operatorname{Tweedie}(\mu,p,\phi),\ 1 < p < 2$$

where $\mu$ is the mean, $p$ is the power parameter constrained between
1 and 2, and $\phi$ is the dispersion parameter. The Tweedie
distribution can be helpful for modelling data that are positive and
continuous but also contain zeros.

Internally, $p$ is transformed from
${logit}^{- 1}\left( \texttt{𝚝𝚑𝚎𝚝𝚊𝚏} \right) + 1$ to constrain it
between 1 and 2 and is estimated as an unconstrained variable.

The [source
code](https://kaskr.github.io/adcomp/tweedie_8cpp_source.html) is
implemented as in the [cplm](https://CRAN.R-project.org/package=cplm)
package (Zhang 2013) and is based on Dunn & Smyth (2005). The TMB
version is defined
[here](https://kaskr.github.io/adcomp/group__R__style__distribution.html#ga262f3c2d1cf36f322a62d902a608aae0).

Example: `family = tweedie(link = "log")`

### Gamma mixture

This is a 2 component mixture that extends the Gamma distribution,

$$(1 - p) \cdot \operatorname{Gamma}\left( \phi,\frac{\mu_{1}}{\phi} \right) + p \cdot \operatorname{Gamma}\left( \phi,\frac{\mu_{2}}{\phi} \right),$$
where $\phi$ represents the Gamma shape, $\mu_{1}/\phi$ represents the
scale for the first (smaller component) of the distribution,
$\mu_{2}/\phi$ represents the scale for the second (larger component) of
the distribution, and $p$ controls the contribution of each component to
the mixture (also interpreted as the probability of larger events).

The mean is $(1 - p) \cdot \mu_{1} + p \cdot \mu_{2}$ and the variance
is
$(1 - p)^{2} \cdot \mu_{1} \cdot \phi^{2} + (p)^{2} \cdot \mu_{2} \cdot \phi^{2}$.

Here, and for the other mixture distributions, the probability of the
larger mean can be obtained from
`plogis(fit$model$par[["logit_p_extreme"]])` and the ratio of the larger
mean to the smaller mean can be obtained from
`1 + exp(fit$model$par[["log_ratio_mix"]])`. The standard errors are
available in the TMB sdreport: `fit$sd_report`.

If you wish to fix the probability of a large (i.e., extreme) mean,
which can be hard to estimate, you can fix this value and pass this to
the family:

``` r
sdmTMB(...,
  family = gamma_mix(link = "log", p_extreme = 0.01)
)
```

See also `family = delta_gamma_mix()` for an extension incorporating
this distribution with delta models.

### Lognormal mixture

This is a 2 component mixture that extends the lognormal distribution,

$$(1 - p) \cdot \operatorname{Lognormal}\left( \log\mu_{1} - \frac{\phi^{2}}{2},\phi^{2} \right) + p \cdot \operatorname{Lognormal}\left( \log\mu_{2} - \frac{\phi^{2}}{2},\phi^{2} \right).$$

Because of the bias correction,
${\mathbb{E}}\lbrack y\rbrack = (1 - p) \cdot \mu_{1} + p \cdot \mu_{2}$
and
${Var}\left\lbrack \log y \right\rbrack = (1 - p)^{2} \cdot \phi^{2} + p^{2} \cdot \phi^{2}$.

As with the Gamma mixture, $p$ controls the contribution of each
component to the mixture (also interpreted as the probability of larger
events).

Example: `family = lognormal_mix(link = "log")`. See also
`family = delta_lognormal_mix()` for an extension incorporating this
distribution with delta models. Like with the gamma mixture, fixed
probabilities of extreme events ($p$ in notation above) can be passed
in, e.g.

``` r
sdmTMB(...,
  family = delta_lognormal_mix(p_extreme = 0.01)
)
```

### Negative binomial 2 mixture

This is a 2 component mixture that extends the NB2 distribution,

$$(1 - p) \cdot \operatorname{NB2}\left( \mu_{1},\phi \right) + p \cdot \operatorname{NB2}\left( \mu_{2},\phi \right)$$

where $\mu_{1}$ is the mean of the first (smaller component) of the
distribution, $\mu_{2}$ is the mean of the larger component, and $p$
controls the contribution of each component to the mixture.

Example: `family = nbinom2_mix(link = "log")`

### Delta models

sdmTMB allows for several different kinds of delta (also known as
hurdle) models. These families are implemented by specifying the family
as a delta distribution. For example:

``` r
sdmTMB(
  ...,
  family = delta_gamma()
)
```

The list of supported families is included in the documentation on
[additional
families](https://pbs-assess.github.io/sdmTMB/reference/families.html),
[delta
models](https://pbs-assess.github.io/sdmTMB/articles/delta-models.html),
and [Poisson-link delta
models](https://pbs-assess.github.io/sdmTMB/articles/poisson-link.html).
By default, the `delta_*` families don’t use the Poisson link, but this
structure can be specified with `delta_gamma(type = "poisson-link")`.

In the “standard” delta model implementation, sdmTMB constructs two
internal models (formally, two “linear predictors”), with the first
model representing presence-absence and the second model represting the
positive component (such as catch rates in fisheries applications).
Default links associated with each family can be inspected with
[`delta_lognormal()`](https://pbs-assess.github.io/sdmTMB/reference/families.md)
and equivalent functions. For the standard delta models, the default
first linear predictor link is logit. For the Poisson-link type, the
first linear predictor link is log.

The `formula`, `spatial`, `spatiotemporal`, and `share_range` arguments
of [`sdmTMB()`](https://pbs-assess.github.io/sdmTMB/reference/sdmTMB.md)
can be specified independently as a 2-element list. For example, the
spatial random field might be estimated for only the first linear
predictor with:

``` r
sdmTMB(
  ...,
  family = delta_gamma(),
  spatial = list("on", "off")
)
```

Or for we may want separate main-effects formulas. For example:

``` r
sdmTMB(
  formula = list(
    y ~ depth + I(depth^2), 
    y ~ 1),
  family = delta_gamma()
)
```

All other arguments are shared between the linear predictors.

Currently if delta models contain smoothers, both components must have
the same main-effects formula.

## Gaussian random fields

### Matérn parameterization

The Matérn defines the covariance $\Phi\left( s_{j},s_{k} \right)$
between spatial locations $s_{j}$ and $s_{k}$ as

$$\Phi\left( s_{j},s_{k} \right) = \tau^{2}/\Gamma(\nu)2^{\nu - 1}\left( \kappa d_{jk} \right)^{\nu}K_{\nu}\left( \kappa d_{jk} \right),$$

where $\tau^{2}$ controls the spatial variance, $\nu$ controls the
smoothness, $\Gamma$ represents the Gamma function, $d_{jk}$ represents
the distance between locations $s_{j}$ and $s_{k}$, $K_{\nu}$ represents
the modified Bessel function of the second kind, and $\kappa$ represents
the decorrelation rate. The parameter $\nu$ is set to 1 to take
advantage of the Stochastic Partial Differential Equation (SPDE)
approximation to the GRF to greatly increase computational efficiency
(Lindgren *et al.* 2011). Internally, the parameters $\kappa$ and $\tau$
are converted to range and marginal standard deviation $\sigma$ as
$\text{range} = \sqrt{8}/\kappa$ and
$\sigma = 1/\sqrt{4\pi\exp\left( 2\log(\tau) + 2\log(\kappa) \right)}$.

In the case of a spatiotemporal model with both spatial and
spatiotemporal fields, if `share_range = TRUE` in
[`sdmTMB()`](https://pbs-assess.github.io/sdmTMB/reference/sdmTMB.md)
(the default), then a single $\kappa$ and range are estimated with
separate $\sigma_{\omega}$ and $\sigma_{\epsilon}$. This often makes
sense since data are often only weakly informative about $\kappa$. If
`share_range = FALSE`, then separate $\kappa_{\omega}$ and
$\kappa_{\epsilon}$ are estimated. The spatially varying coefficient
field always shares $\kappa$ with the spatial random field.

### Projection $\mathbf{A}$ matrix

The values of the spatial variables at the knots are multiplied by a
projection matrix $\mathbf{A}$ that bilinearly interpolates from the
knot locations to the values at the locations of the observed or
predicted data (Lindgren & Rue 2015)

$${\mathbf{ω}}^{*} = \mathbf{A}{\mathbf{ω}},$$ where ${\mathbf{ω}}^{*}$
represents the values of the spatial random fields at the observed
locations or predicted data locations. The matrix $\mathbf{A}$ has a row
for each data point or prediction point and a column for each knot.
Three non-zero elements on each row define the weight of the
neighbouring 3 knot locations for location $\mathbf{s}$. The same
bilinear interpolation happens for any spatiotemporal random fields

$${\mathbf{ϵ}}_{t}^{*} = \mathbf{A}{\mathbf{ϵ}}_{t}.$$

### Anisotropy

TMB allows for anisotropy, where spatial covariance may be asymmetric
with respect to latitude and longitude ([full
details](https://kaskr.github.io/adcomp/namespaceR__inla.html)).
Anisotropy can be turned on or off with the logical `anisotropy`
argument to
[`sdmTMB()`](https://pbs-assess.github.io/sdmTMB/reference/sdmTMB.md).
There are a number of ways to implement anisotropic covariance (Fuglstad
*et al.* 2015), and we adopt a 2-parameter rotation matrix $\textbf{𝐇}$.
The elements of $\textbf{𝐇}$ are defined by the parameter vector
$\mathbf{x}$ so that $H_{1,1} = x_{1}$, $H_{1,2} = H_{2,1} = x_{2}$ and
$H_{2,2} = \left( 1 + x_{2}^{2} \right)/x_{1}$.

Once a model is fitted with
[`sdmTMB()`](https://pbs-assess.github.io/sdmTMB/reference/sdmTMB.md),
the anisotropy relationships may be plotted using the
[`plot_anisotropy()`](https://pbs-assess.github.io/sdmTMB/reference/plot_anisotropy.md)
function, which takes the fitted object as an argument. If a barrier
mesh is used, anisotropy is disabled.

### Incorporating physical barriers into the SPDE

In some cases the spatial domain of interest may be complex and bounded
by some barrier such as by land or water (e.g., coastlines, islands,
lakes). SPDE models allow for physical barriers to be incorporated into
the modelling (Bakka *et al.* 2019). With
[`sdmTMB()`](https://pbs-assess.github.io/sdmTMB/reference/sdmTMB.md)
models, the mesh construction occurs in two steps: the user (1)
constructs a mesh with a call to
[`sdmTMB::make_mesh()`](https://pbs-assess.github.io/sdmTMB/reference/make_mesh.md),
and (2) passes the mesh to
[`sdmTMB::add_barrier_mesh()`](https://pbs-assess.github.io/sdmTMB/reference/add_barrier_mesh.md).
The barriers must be constructed as `sf` objects (Pebesma 2018) with
polygons defining the barriers. See
[`?sdmTMB::add_barrier_mesh`](https://pbs-assess.github.io/sdmTMB/reference/add_barrier_mesh.md)
for an example.

The barrier implementation requires the user to select a fraction value
(`range_fraction` argument) that defines the fraction of the usual
spatial range when crossing the barrier (Bakka *et al.* 2019). For
example, if the range was estimated at 10 km, `range_fraction = 0.2`
would assume that the range was 2 km across the barrier. This would let
the spatial correlation decay 5 times faster with distance. From
experimentation, values around 0.1 or 0.2 seem to work well but values
much lower than 0.1 can result in convergence issues.

[This website](https://haakonbakkagit.github.io/btopic128.html) by
Francesco Serafini and Haakon Bakka provides an illustration with INLA.
The implementation within TMB was borrowed from code written by Olav
Nikolai Breivik and Hans Skaug at the [TMB Case
Studies](https://github.com/skaug/tmb-case-studies) Github site.

## Optimization

### Optimization details

The sdmTMB model is fit by maximum marginal likelihood. Internally, a
TMB (Kristensen *et al.* 2016) model template calculates the marginal
log likelihood and its gradient, and the negative log likelihood is
minimized via the non-linear optimization routine
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) in R (Gay 1990;
R Core Team 2021). Random effects are estimated at values that maximize
the log likelihood conditional on the estimated fixed effects and are
integrated over via the Laplace approximation (Kristensen *et al.*
2016).

Like AD Model Builder (Fournier *et al.* 2012), TMB allows for
parameters to be fit in phases and we include the `multiphase` argument
in
[`sdmTMB::sdmTMBcontrol()`](https://pbs-assess.github.io/sdmTMB/reference/sdmTMBcontrol.md)
to allow this. For high-dimensional models (many fixed and random
effects), phased estimation may be faster and result in more stable
convergence. In sdmTMB, phased estimation proceeds by first estimating
all fixed-effect parameters contributing to the likelihood (holding
random effects constant at initial values). In the second phase, the
random-effect parameters (and their variances) are also estimated.
Fixed-effect parameters are also estimated in the second phase and are
initialized at their estimates from the first phase.

In some cases, a single call to
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) may not be
result in convergence (e.g., the maximum gradient of the marginal
likelihood with respect to fixed-effect parameters is not small enough
yet), and the algorithm may need to be run multiple times. In the
[`sdmTMB::sdmTMBcontrol()`](https://pbs-assess.github.io/sdmTMB/reference/sdmTMBcontrol.md)
function, we include an argument `nlminb_loops` that will restart the
optimization at the previous best values. The number of `nlminb_loops`
should generally be small (e.g., 2 or 3 initially), and defaults to 1.
For some sdmTMB models, the Hessian may also be unstable and need to be
re-evaluated. We do this optionally with the
[`stats::optimHess()`](https://rdrr.io/r/stats/optim.html) routine after
the call to [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).
The [`stats::optimHess()`](https://rdrr.io/r/stats/optim.html) function
implements a Newton optimization routine to find the Hessian, and we
include the argument `newton_loops` in
[`sdmTMB::sdmTMBcontrol()`](https://pbs-assess.github.io/sdmTMB/reference/sdmTMBcontrol.md)
to allow for multiple function evaluations (each starting at the
previous best value). By default, this is not included and
`newton_loops` is set to 0. If a model is already fit, the function
[`sdmTMB::run_extra_optimization()`](https://pbs-assess.github.io/sdmTMB/reference/run_extra_optimization.md)
can run additional optimization loops with either routine to further
reduce the maximum gradient.

### Assessing convergence

Much of the guidance around diagnostics and glmmTMB also applies to
sdmTMB, e.g. [the glmmTMB vignette on
troubleshooting](https://CRAN.R-project.org/package=glmmTMB).
Optimization with
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) involves
specifying the number of iterations and evaluations (`eval.max` and
`iter.max`) and the tolerances (`abs.tol`, `rel.tol`, `x.tol`,
`xf.tol`)—a greater number of iterations and smaller tolerance
thresholds increase the chance that the optimal solution is found, but
more evaluations translates into longer computation time. Warnings of
non-positive-definite Hessian matrices (accompanied by parameters with
`NA`s for standard errors) often mean models are improperly specified
given the data. Standard errors can be observed in the output of
`print.sdmTMB()` or by checking `fit$sd_report`. The maximum gradient of
the marginal likelihood with respect to fixed-effect parameters can be
checked by inspecting (`fit$gradients`). Guidance varies, but the
maximum gradient should likely be at least $< 0.001$ before assuming the
fitting routine is consistent with convergence. If maximum gradients are
already relatively small, they can often be reduced further with
additional optimization calls beginning at the previous best parameter
vector as described above with
[`sdmTMB::run_extra_optimization()`](https://pbs-assess.github.io/sdmTMB/reference/run_extra_optimization.md).

## References

Anderson, S.C., Branch, T.A., Cooper, A.B. & Dulvy, N.K. (2017).
Black-swan events in animal populations. *Proceedings of the National
Academy of Sciences*, **114**, 3252–3257.

Bakka, H., Vanhatalo, J., Illian, J., Simpson, D. & Rue, H. (2019).
Non-stationary Gaussian models with physical barriers. *arXiv:1608.03787
\[stat\]*. Retrieved from <https://arxiv.org/abs/1608.03787>

Bürkner, P.-C. (2017). [brms: An R package for Bayesian multilevel
models using Stan](https://doi.org/10.18637/jss.v080.i01). *Journal of
Statistical Software*, **80**, 1–28.

Cribari-Neto, F. & Zeileis, A. (2010). Beta regression in R. *Journal of
Statistical Software*, **34**, 1–24.

Dunn, P.K. & Smyth, G.K. (2005). [Series evaluation of Tweedie
exponential dispersion model
densities](https://doi.org/10.1007/s11222-005-4070-y). *Statistics and
Computing*, **15**, 267–280.

Edwards, A.M. & Auger-Méthé, M. (2019). Some guidance on using
mathematical notation in ecology. *Methods in Ecology and Evolution*,
**10**, 92–99.

Ferrari, S. & Cribari-Neto, F. (2004). [Beta regression for modelling
rates and proportions](https://doi.org/10.1080/0266476042000214501).
*Journal of Applied Statistics*, **31**, 799–815.

Fournier, D.A., Skaug, H.J., Ancheta, J., Ianelli, J., Magnusson, A.,
Maunder, M.N., Nielsen, A. & Sibert, J. (2012). AD model builder: Using
automatic differentiation for statistical inference of highly
parameterized complex nonlinear models. *Optimization Methods and
Software*, **27**, 233–249.

Fuglstad, G.-A., Lindgren, F., Simpson, D. & Rue, H. (2015). Exploring a
new class of non-stationary spatial Gaussian random fields with varying
local anisotropy. *Statistica Sinica*, **25**, 115–133.

Gay, D.M. (1990). Usage summary for selected optimization routines.
*Computing Science Technical Report*, **153**, 1–21.

Hilbe, J.M. (2011). *Negative Binomial Regression*. Cambridge University
Press.

Kristensen, K., Nielsen, A., Berg, C.W., Skaug, H. & Bell, B.M. (2016).
[TMB: Automatic differentiation and Laplace
approximation](https://doi.org/10.18637/jss.v070.i05). *Journal of
Statistical Software*, **70**, 1–21.

Lindgren, F. & Rue, H. (2015). [Bayesian spatial modelling with
R-INLA](https://doi.org/10.18637/jss.v063.i19). *Journal of Statistical
Software*, **63**, 1–25.

Lindgren, F., Rue, H. & Lindström, J. (2011). An explicit link between
Gaussian fields and Gaussian Markov random fields: The stochastic
partial differential equation approach. *J. R. Stat. Soc. Ser. B Stat.
Methodol.*, **73**, 423–498.

Pebesma, E. (2018). Simple Features for R: Standardized Support for
Spatial Vector Data. *The R Journal*, **10**, 439–446. Retrieved from
<https://doi.org/10.32614/RJ-2018-009>

R Core Team. (2021). *R: A language and environment for statistical
computing*. R Foundation for Statistical Computing, Vienna, Austria.
Retrieved from <https://www.R-project.org/>

Thorson, J.T. (2019). [Guidance for decisions using the Vector
Autoregressive Spatio-Temporal (VAST) package in stock, ecosystem,
habitat and climate
assessments](https://doi.org/10.1016/j.fishres.2018.10.013). *Fisheries
Research*, **210**, 143–161.

Wood, S.N. (2017). *Generalized additive models: An introduction with
R*, 2nd edn. Chapman and Hall/CRC.

Wood, S. & Scheipl, F. (2020). *gamm4: Generalized Additive Mixed Models
using ’mgcv’ and ’lme4’*. Retrieved from
<https://CRAN.R-project.org/package=gamm4>

Zhang, Y. (2013). [Likelihood-based and Bayesian methods for Tweedie
compound Poisson linear mixed
models](https://doi.org/10.1007/s11222-012-9343-7). *Statistics and
Computing*, **23**, 743–757.
