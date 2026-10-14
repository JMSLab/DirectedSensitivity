# Problem specification

Fill out this markdown.
After completing the verification steps below, hand the filled markdown to an LLM instructed to follow the workflow in [`WORKFLOW.md`](../../workflow/WORKFLOW.md).
Provide precise mathematical definitions in each section unless instructed otherwise.

The _object of interest_ $c(\nu) \in \mathbb{R}$ depends on an (unknown) parameter $\nu \in \mathcal{N}$, where the space $\mathcal{N}$ is the _nesting model_ and may be infinite-dimensional.
The data-dependent _criterion_ $\hat{q}(\nu) \in \mathbb{R}$ describes the fit of the nesting model to the data.
The _working model_ is a finite-dimensional subset $\mathcal{A}$ of $\mathcal{N}$.
The working model's fit is $\min_{\alpha\in \mathcal{A}}\hat q(\alpha)$, and the baseline estimate $\hat{\alpha} \in \mathcal{A}$ is the minimizer.

The objectives of this template are:
1. Precisely define the following objects:
   - The nesting model $\mathcal{N}$,
   - The working model $\mathcal{A}$,
   - The object of interest $c(\cdot)$, and
   - The criterion $\hat{q}(\cdot)$.
2. Specify the location of data, code, and estimates.

## 1. The nesting model

Define the nesting model $\mathcal{N}$, which encodes economic restrictions on unknown parameters, including unknown functions.
Define the domain and range of all unknown functions formally.
Specify any restrictions on the function class, such as continuity and shape restrictions.
The nesting model may also include a finite set of other known or unknown parameters.
Define the domain or values of these parameters.

**Example**. The nesting model $\mathcal{N}$ is the set of all demand functions $D: [P_0, P_1] \rightarrow [0, 1]$ that are:
- Continuous, and
- Strictly decreasing.

The bounds $P_0$ and $P_1$ are given in [Location of data, code, and estimates](#6-location-of-data-code-and-estimates).

| Your version: |
| --- |

A naive quasi-hyperbolic filer decides on each day $t = 0, \dots, 64$ of the filing season whether to file.
The last dated day of the season, the deadline, is dropped, since the model takes filing as certain on it.
Each day the filer draws a net benefit of filing $X$, in dollars: the value of having the return filed today rather than not, net of the day's hassle.
The refund is excluded because the one-period discount factor is fixed at $1$, so the refund is worth the same whenever it arrives.
Draws are independent across days and filers.

The primitive $\nu = (\beta, G)$ is the present bias $\beta$ and the distribution function $G$ of $X$, the same on every day.
The filer files on day $t$ if the day's draw exceeds the threshold $w_t$, where

$$w_{64} = \beta \mathbb{E}[X],$$

$$w_t = \beta \mathbb{E}\big[X \mathbf{1}\lbrace X > w_{t+1} / \beta \rbrace\big] + G(w_{t+1} / \beta) w_{t+1},$$

so the filer files on day $t$ with probability $p_t = 1 - G(w_t)$.
The division by $\beta$ is the naivety: the filer expects his next self to apply no present bias.

$\mathcal{N}$ is the set of pairs $(\beta, G)$ with $0 < \beta \le 1$ and $G$ a distribution function on $\mathbb{R}$ such that:

- Its mass lies in $[-25, 25]$: $G(-25) \le 10^{-8}$ and $1 - G(25) \le 10^{-8}$.
- Its mean exists and lies in the same interval: $\lvert \mathbb{E}[X \mathbf{1}\lbrace X > 25 \rbrace] \rvert \le 10^{-6}$ and $\lvert \mathbb{E}[X \mathbf{1}\lbrace X \le -25 \rbrace] \rvert \le 10^{-6}$.
- Filing is costly on average: $\mathbb{E}[X] \le 0$.

<br>

## 2. The working model

Define the working model $\mathcal{A} \subset \mathcal{N}$.
Give its parameters, their domain, and the map from parameters to primitives.

**Example**. The working model is the isoelastic family

$$\mathcal{A} = \left\lbrace D_{a,b}(p)=\exp(a-b\log(p)): b > 0, a \le b \log P_0 \right\rbrace.$$

Parameters $a$ and $b$ map to the primitive $\alpha = D_{a,b} \in \mathcal{A}$.

| Your version: |
| --- |

$\mathcal{A}$ is the logistic family with scale one,

$$\mathcal{A} = \left\lbrace (\beta, G_\kappa) : G_\kappa(x) = \big(1 + \exp(-(x + \kappa))\big)^{-1}, 0 < \beta \le 1, \kappa \ge 0, (\beta, G_\kappa) \in \mathcal{N} \right\rbrace,$$

where $\kappa$ is the average cost of filing, so the mean net benefit is $-\kappa$ and its variance is $\pi^2 / 3$.
Parameters $(\beta, \kappa)$ map to the primitive $(\beta, G_\kappa)$.

<br>

## 3. The object of interest

Define the object of interest $c(\cdot): \mathcal{N} \to \mathbb{R}$.
It must be well-defined and finite for every $\nu \in \mathcal{N}$.

**Example**. Writing $\nu = D$, the object is the loss in consumer surplus from a price increase $P_0 \to P_1$:

$$c(D) = \int_{P_0}^{P_1} D(p) \mathrm{d}p.$$

Since $D \in [0,1]$ on the bounded interval $[P_0, P_1]$, the integral lies in $[0, P_1 - P_0]$ and is finite for every $D \in \mathcal{N}$.

| Your version: |
| --- |

$c(\nu) = \beta$ is the present-bias parameter of the filer.
It is finite for every $\nu \in \mathcal{N}$.

<br>

## 4. The criterion

Define the criterion $\hat{q}(\nu)$.
The criterion is a data-dependent measure of fit, where smaller values indicate a better fit.
The criterion must be finite and well-defined for every $\nu \in \mathcal{N}$.

**Example**. With observed shares $s_i$ at prices $p_i \in [P_0, P_1]$, for $i = 1, \dots, n$,

$$\hat{q}(D) = \frac{1}{n} \sum_{i=1}^{n} \big( s_i - D(p_i) \big)^2.$$

This evaluates $D\in[0,1]$ at finitely many points, so it is finite for every $D \in \mathcal{N}$.
The shares $s_i$ and prices $p_i$ are given under [Location of data, code, and estimates](#6-location-of-data-code-and-estimates).


| Your version: |
| --- |

$\hat{q}(\nu)$ is the negative log-likelihood of the observed filing decisions.
With $n_t$ the filers still at risk at the start of day $t$ and $f_t$ those who filed that day,

$$\hat{q}(\nu) = -\sum_{t=0}^{64} \Big( f_t \log p_t + (n_t - f_t) \log(1 - p_t) \Big).$$

It is finite for every $\nu \in \mathcal{N}$ at which every $p_t$ lies in $(0, 1)$; elsewhere $\hat{q}(\nu) = +\infty$.

The data are given in [Location of data, code, and estimates](#6-location-of-data-code-and-estimates).

<br>

## 5. Computational constraints (optional)

State the largest number of evaluations of $c(\cdot)$ and $\hat{q}(\cdot)$ you are willing to spend, for when either is expensive.
Fill in the blanks below.

At most ___ evaluations of $c(\cdot)$ and at most ___ evaluations of $\hat{q}(\cdot)$.


## 6. Location of data, code, and estimates

Provide the location of data inputs.
Provide the location of code that implements the object of interest and criterion (or numerical approximations of these functions).
Provide the values of any parameters referenced in the prior sections.
Provide the location or values of the parameters estimated under the working model.

**Example**.
- $\lbrace (p_i, s_i) \rbrace_{i=1}^{n}$ — `demand.csv`, columns `price` and `share`.
- $c(D)$ — implemented in `consumer_surplus.py`.
- $\hat{q}(D)$ — implemented in `estimator.py`.
- $P_0, P_1$ — equal 100, 200.
- $\hat{a}, \hat{b}$ parameter estimates - `estimates.csv`, columns `intercept` and `slope`.

| Your version: |
| --- |

- $(n_t, f_t)$ — `Distribution Code/Condensed_Data.xlsx` in `temp/MMS2023/jvac067_martinez_meier_sprenger_distribution_code.zip`, columns `frequency` and `filed`, the first 65 rows (`specification.py:SEASON_DAYS`).
- $c(\nu)$ and $\hat{q}(\nu)$ — `specification.py:{Target,Criterion}`. 
  A candidate is `specification.py:Primitive(present_bias, distribution)`, with `distribution` any object with methods `Cdf`, `PartialMean`, and `Mean`; `specification.py:Logistic` implements $G_\kappa$ with location $-\kappa$.
- The restrictions of $\mathcal{N}$ — `specification.py:ConstraintSet`.
- $\hat\alpha = (\hat\beta, G_{\hat\kappa})$ — `output/estimate.json`, produced by `estimate.py`, the maximum-likelihood estimate from `specification.py:FitMaximumLikelihood`.
- Replication of the published estimate — `output/replicate.txt`, produced by `replicate.py`.

<br>

## 7. Consistency and completeness checklist

After you complete the template, ask your LLM to verify the following prior to launching the workflow:

- [x] [The nesting model](#1-the-nesting-model), [The working model](#2-the-working-model), [The object of interest](#3-the-object-of-interest), and [The criterion](#4-the-criterion) are each filled with precise mathematical definitions.
- [x] The nesting model $\mathcal{N}$ is non-empty.
- [x] The working model is a genuine subset of the nesting model, $\mathcal{A} \subset \mathcal{N}$, and is finite-dimensional.
- [x] The object of interest $c(\nu)$ is well-defined and finite for every $\nu \in \mathcal{N}$.
- [x] The criterion $\hat{q}(\nu)$ is well-defined and finite for every $\nu \in \mathcal{N}$.
- [x] Every data or code symbol used in [The nesting model](#1-the-nesting-model) through [Computational constraints](#5-computational-constraints-optional) is resolved in [Location of data, code, and estimates](#6-location-of-data-code-and-estimates).

## 8. Workflow configuration

Should the LLM attempt a Lean formalization of mathematical claims?
Formalization may add time and token cost, and is not required to verify numerical certificates.

- [x] Yes
- [ ] No
