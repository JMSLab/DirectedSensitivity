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

Throughout, $k$ denotes the charity ECU.
The primitive $\nu$ is the distribution $F$ of the altruism parameter $a$ among households solicited for ECU.

Everything else is fixed at the estimates of DellaVigna, List, and Malmendier (2012): the eleven scalar parameters $\hat\theta$, which include the probability $h_0$ of being at home, the probability $r$ of seeing a flyer, the response $\eta$ of home presence to the utility of a visit, the parameter $g_i$ of the utility of giving, and the social pressure $s_k$; the pressure-discharging gift $G^S = 10$; and the altruism distribution for La Rabida, at its censored normal estimate.

At $\hat\theta$, a household with altruism $a$ gives

$$g^{\ast}_k(a) = 0 \text{ for } a \le a^{--}_k, \quad \frac{a}{1-s_k} - g_i \text{ for } a^{--}_k < a \le a^{-}_k, \quad G^S \text{ for } a^{-}_k < a \le a^{+}, \quad a - g_i \text{ for } a > a^{+},$$


with thresholds $a^{--}_k = (1-s_k) g_i$, $a^{-}_k = (1-s_k)(g_i + G^S)$, and $a^{+} = g_i + G^S$.
The utility of a visit is

$$u_k(a) = a \big[\log\big(g_i + g^{\ast}_k(a)\big) - \log g_i\big] - g^{\ast}_k(a) - s_k \big(G^S - g^{\ast}_k(a)\big) \mathbf{1}\lbrace g^{\ast}_k(a) < G^S \rbrace,$$


which is continuous, equals $-s_k G^S < 0$ on $[0, a^{--}_k]$, is strictly increasing above $a^{--}_k$, and crosses zero at a unique $a^0_k$.
A household that sees a flyer is at home with probability $h^{\ast}_k(a) = \min\lbrace 1, \max\lbrace 0, h_0 + \eta u_k(a) \rbrace \rbrace$.
The thresholds $(a^{--}_k, a^{-}_k, a^{+}, a^0_k)$ are $(9.081, 17.64, 20.61, 14.65)$.

$\mathcal{N}$ is the set of probability distributions $F$ on $[0, \bar a_k]$ that consist of an atom $\pi \in [0,1]$ at $a = 0$ and a density $f \ge 0$ on $(0, \bar a_k]$, with $\pi + \int_0^{\bar a_k} f(a) \mathrm{d}a = 1$.
There are no further restrictions: not unimodality, monotonicity, or continuity of $f$.

The values of $\hat\theta$ and the bound $\bar a_k$ are given in [Location of data, code, and estimates](#6-location-of-data-code-and-estimates).

<br>

## 2. The working model

Define the working model $\mathcal{A} \subset \mathcal{N}$.
Give its parameters, their domain, and the map from parameters to primitives.

**Example**. The working model is the isoelastic family

$$\mathcal{A} = \left\lbrace D_{a,b}(p)=\exp(a-b\log(p)): b > 0, a \le b \log P_0 \right\rbrace.$$

Parameters $a$ and $b$ map to the primitive $\alpha = D_{a,b} \in \mathcal{A}$.

| Your version: |
| --- |

$\mathcal{A}$ is the censored normal family

$$\mathcal{A} = \left\lbrace F_{\mu,\sigma} : \mu \in \mathbb{R}, \sigma > 0 \right\rbrace,$$


where $F_{\mu,\sigma}$ has density $\sigma^{-1}\phi\big((a-\mu)/\sigma\big)$ on $(0, \bar a_k]$ and an atom $\Phi(-\mu/\sigma) + 1 - \Phi\big((\bar a_k - \mu)/\sigma\big)$ at $a = 0$: a normal distribution censored below at zero, with its mass above $\bar a_k$ also assigned to the atom.
Parameters $(\mu, \sigma)$ map to the primitive $F_{\mu,\sigma}$, and $F_{\mu,\sigma} \in \mathcal{N}$.

<br>

## 3. The object of interest

Define the object of interest $c(\cdot): \mathcal{N} \to \mathbb{R}$.
It must be well-defined and finite for every $\nu \in \mathcal{N}$.

**Example**. Writing $\nu = D$, the object is the loss in consumer surplus from a price increase $P_0 \to P_1$:

$$c(D) = \int_{P_0}^{P_1} D(p) \mathrm{d}p.$$

Since $D \in [0,1]$ on the bounded interval $[P_0, P_1]$, the integral lies in $[0, P_1 - P_0]$ and is finite for every $D \in \mathcal{N}$.

| Your version: |
| --- |

$c(F)$ is the average welfare effect, in dollars, of a visit by a fund-raiser for ECU without a flyer, per household approached, households not at home included,

$$c(F) = h_0 \int_{[0, \bar a_k]} u_k(a) \mathrm{d}F(a).$$


On $[0, \bar a_k]$, $u_k$ is bounded between $-s_k G^S$ and $u_k(\bar a_k)$, so $c(F)$ is finite for every $F \in \mathcal{N}$.

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

$\hat{q}(F)$ is the minimum-distance objective of DellaVigna, List, and Malmendier (2012),

$$\hat{q}(F) = \big(\hat m - m(F)\big)^{\prime} W \big(\hat m - m(F)\big),$$


where $\hat m$ is the vector of 70 empirical moments, $W$ is the inverse of the diagonal of their variance-covariance matrix, and $m(F)$ is the vector of model moments at $\hat\theta$.
With $\mathbf{1}_A$ the indicator of a set $A$, the 22 moments of ECU are

$$m_{\mathrm{PH,NF}} = h_0, \quad m_{\mathrm{PH,F}} = (1-r) h_0 + r \int h^{\ast}_k \mathrm{d}F, \quad m_{\mathrm{PH,OO}} = (1-r) h_0 + r \int h^{\ast}_k \mathbf{1}_{(a^0_k, \infty)} \mathrm{d}F,$$



$$m_{\mathrm{PG,NF}} = h_0 S, \quad m_{\mathrm{PG,F}} = (1-r) h_0 S + r \int h^{\ast}_k \mathbf{1}_{(a^{--}_k, \infty)} \mathrm{d}F, \quad m_{\mathrm{PG,OO}} = (1-r) h_0 S + r \int h^{\ast}_k \mathbf{1}_{(a^0_k, \infty)} \mathrm{d}F,$$



$$m_{\mathrm{POO,OO}} = r h_0 F\big([0, a^0_k]\big),$$


with $S = F\big((a^{--}_k, \bar a_k]\big)$, and, for each of five giving bins $b$ with altruism interval $I_b$,

$$m_{\mathrm{PG,NF},b} = h_0 F(I_b), \quad m_{\mathrm{PG,F},b} = (1-r) h_0 F(I_b) + r \int h^{\ast}_k \mathbf{1}_{I_b} \mathrm{d}F, \quad m_{\mathrm{PG,OO},b} = (1-r) h_0 F(I_b) + r \int h^{\ast}_k \mathbf{1}_{I_b \cap (a^0_k, \infty)} \mathrm{d}F.$$


The bins are gifts in $(0, G^S)$, exactly $G^S$, $(10, 20]$, $(20, 50]$, and $(50, \infty)$, with $I_b = (a^{--}_k, a^{-}_k]$, $(a^{-}_k, a^{+}]$, $(g_i + 10, g_i + 20]$, $(g_i + 20, g_i + 50]$, and $(g_i + 50, \infty)$.
The remaining 48 moments, 26 from a survey and 22 for La Rabida, do not depend on $F$.

Each moment is bounded and affine in $F$, so $\hat{q}(F)$ is finite for every $F \in \mathcal{N}$.

The values of $\hat m$, $W$, and $\hat\theta$ are given in [Location of data, code, and estimates](#6-location-of-data-code-and-estimates).

<br>

## 5. Computational constraints (optional)

State the largest number of evaluations of $c(\cdot)$ and $\hat{q}(\cdot)$ you are willing to spend, for when either is expensive.
Fill in the blanks below.

At most 50,000 evaluations of $c(\cdot)$ and at most 20,000 evaluations of $\hat{q}(\cdot)$.


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

- $\hat m$ and the variance-covariance matrix behind $W$ — `charity_code_and_appendix/structural estimation/Moments70.mat` in `temp/DLM2012/charity_code_and_appendix.zip`, variables `Moments` and `VCcontrol`.
- $\hat\theta$, $G^S$, and $\bar a_k$ — `specification.py:{SCALARS,GS,A_MAX}`.
- $c(F)$ and $\hat{q}(F)$ — `specification.py:TargetECU` and `specification.py:Criterion`, evaluated at `Primitive(BENCHMARK, lar=Baseline().lar, ecu=F)`.
  A candidate `F` is any `specification.py:AltruismDistribution`; `CensoredNormal` implements $F_{\mu,\sigma}$.
- $\hat\alpha = F_{\hat\mu, \hat\sigma}$ — `output/estimate.json`, produced by `estimate.py`, and read back by `specification.py:Baseline`.
- Replication of the published fit and welfare effects — `output/replicate.txt`, produced by `replicate.py`.

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
