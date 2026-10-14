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

The primitive is the lobbying technology $\nu = (\gamma, f_f, f_a)$: the exponent $\gamma$ with which money buys lobbying effort, and the impact functions $f_f$ and $f_a$ of lobbying effort for and against a policy.

Fix a policy.
Let $\pi \in [0,1]$ be the probability it is enacted absent lobbying, $L_f$ and $L_a$ the sets of lobbying coalitions that value it positively and negatively, $v_i$ the value to coalition $i$, $d_i \in \lbrace 0,1 \rbrace$ its participation, and $s_i \ge 0$ its spending.
Lobbying effort on each side is

$$T_f = \sum_{i \in L_f} s_i^{\gamma}, \qquad T_a = \sum_{j \in L_a} s_j^{\gamma}.$$

The policy is decided by a contest among three parties: its supporters, with impact $f_f(T_f)$; its opponents, with impact $f_a(T_a)$; and Nature, with constant impact $1$.
Each party prevails with probability proportional to its impact, and if Nature prevails the status quo decides, enacting the policy with probability $\pi$.
The policy is therefore enacted with probability

$$P(\pi, T_f, T_a) = \frac{\pi + f_f(T_f)}{1 + f_f(T_f) + f_a(T_a)}.$$

This is the additive contest success function with a draw, equation (6) of Jia, Skaperdas, and Vaidya (2013), https://doi.org/10.1016/j.ijindorg.2012.06.007, with the draw resolved by $\pi$.
Nature's impact of $1$ is a normalization, since multiplying all three impacts by a constant leaves $P$ unchanged.

$\mathcal{N}$ is the set of triples $(\gamma, f_f, f_a)$ with $\gamma \in (0,1)$ and $`f_f, f_a : \mathbb{R}_+ \to \mathbb{R}_+`$ that are:

- Continuously differentiable, with finite derivatives.
- Zero at zero: $f_f(0) = f_a(0) = 0$, so that $P(\pi, 0, 0) = \pi$.
- Strictly increasing.
- Concave.

Every such $P$ lies in $[0,1]$ and is increasing in $\pi$ and $T_f$ and decreasing in $T_a$.
Write $P_f \equiv \partial P/\partial T_f$ and $P_a \equiv \partial P/\partial T_a$.
Because $P$ is increasing and concave in $f_f$ and decreasing and convex in $f_a$, concavity of the impact functions makes $P$ concave in $T_f$ and convex in $T_a$.

This curvature makes each coalition's payoff $P(\pi,T_f,T_a) v_i - s_i$ concave in its own spending: a supporter composes a concave increasing $P$ with the concave $s_i^{\gamma}$, and an opponent composes a convex decreasing $P$ with the same concave $s_j^{\gamma}$ and then multiplies by $v_j < 0$.
The second-stage equilibrium is therefore characterized by first-order conditions.
A participating $i \in L_f$ spends $s_i = (\gamma v_i P_f)^{1/(1-\gamma)}$ and a participating $j \in L_a$ spends $s_j = (\gamma |v_j| (-P_a))^{1/(1-\gamma)}$, so equilibrium effort $(\bar{T}_f, \bar{T}_a)$ solves

$$T_f = C_f \cdot P_f(\pi,T_f,T_a)^{\frac{\gamma}{1-\gamma}}, \qquad T_a = C_a \cdot \big(-P_a(\pi,T_f,T_a)\big)^{\frac{\gamma}{1-\gamma}},$$

$$C_f = \sum_{i \in L_f} d_i (\gamma v_i)^{\frac{\gamma}{1-\gamma}}, \qquad C_a = \sum_{j \in L_a} d_j (\gamma |v_j|)^{\frac{\gamma}{1-\gamma}}.$$

A solution exists for every $\nu \in \mathcal{N}$ but need not be unique, so $(\bar{T}_f, \bar{T}_a)$ is *defined* as the point reached by Gauss-Seidel iteration on the display above, started at $(0,0)$ and run for the number of steps set in the configuration.
This makes $c$ and $\hat q$ well defined on all of $\mathcal{N}$.
Participation $d$ is the profile that maximizes the sum of payoffs among the pure-strategy Nash equilibria of the entry game over the $2^4$ profiles, each profile evaluated at its own second-stage equilibrium.

Everything outside the lobbying technology is held fixed at the published estimates and is not part of $\nu$: the initial enactment probability $\pi_k = \Phi(Z_k \delta + \sigma_{\xi} \xi_k)$, the values $\log|v_{\ell k}| = X_{\ell k}\alpha + \sigma_{\eta \ell}\eta_{\ell k}$ with sign given by the observed position of coalition $\ell$ on policy $k$, and the entry cost.
Locations and values are in [Location of data, code, and estimates](#6-location-of-data-code-and-estimates).

<br>

## 2. The working model

Define the working model $\mathcal{A} \subset \mathcal{N}$.
Give its parameters, their domain, and the map from parameters to primitives.

**Example**. The working model is the isoelastic family

$$\mathcal{A} = \left\lbrace D_{a,b}(p)=\exp(a-b\log(p)): b > 0, a \le b \log P_0 \right\rbrace.$$

Parameters $a$ and $b$ map to the primitive $\alpha = D_{a,b} \in \mathcal{A}$.

| Your version: |
| --- |

$\mathcal{A}$ is the family with linear impact functions,

$$\mathcal{A} = \left\lbrace (\gamma, f_f, f_a) : f_f(T) = \beta_f T, f_a(T) = \beta_a T, \beta_f > 0, \beta_a > 0, \gamma \in (0,1) \right\rbrace,$$

with parameters $(\beta_f, \beta_a, \gamma)$ mapping to the primitive.
It gives the ratio-form contest

$$P(\pi, T_f, T_a) = \frac{\pi + \beta_f T_f}{1 + \beta_f T_f + \beta_a T_a}.$$

Linear impact functions with positive slopes are continuously differentiable, zero at zero, strictly increasing, and concave, so $\mathcal{A} \subset \mathcal{N}$.

<br>

## 3. The object of interest

Define the object of interest $c(\cdot): \mathcal{N} \to \mathbb{R}$.
It must be well-defined and finite for every $\nu \in \mathcal{N}$.

**Example**. Writing $\nu = D$, the object is the loss in consumer surplus from a price increase $P_0 \to P_1$:

$$c(D) = \int_{P_0}^{P_1} D(p) \mathrm{d}p.$$

Since $D \in [0,1]$ on the bounded interval $[P_0, P_1]$, the integral lies in $[0, P_1 - P_0]$ and is finite for every $D \in \mathcal{N}$.

| Your version: |
| --- |

$c(\nu)$ is the average return to lobbying of the coal coalition, in percent: its expected net benefit from lobbying per dollar of total lobbying cost.
For policy $k$ and a draw $r$ of $(\xi, \eta)$, let $`\bar{d}_{k,r}(\nu)`$ be the selected participation profile and $`\bar{P}_{k,r}(\nu) = P(\pi_{k,r}, \bar{T}_{f}, \bar{T}_{a})`$ the equilibrium enactment probability.
For a coalition $\ell$ that participates in $`\bar{d}_{k,r}(\nu)`$, let $`\bar{P}^{-\ell}_{k,r}(\nu)`$ be the enactment probability at the same profile with $\ell$ removed, the other coalitions' participation unchanged and their efforts re-solved, and let $`\bar{s}_{\ell k,r}(\nu)`$ be the equilibrium spending of $\ell$.
Its net benefit from lobbying is

$$b_{\ell k,r}(\nu) = \big(\bar{P}_{k,r}(\nu) - \bar{P}^{-\ell}_{k,r}(\nu)\big) v_{\ell k,r} - \bar{s}_{\ell k,r}(\nu) - c_{\ell k},$$

where $v_{\ell k,r}$ is its value, signed by its position, and $c_{\ell k}$ the entry cost.
With $d_k$ the observed participation profile and $K$ the set of policies on which the coal coalition is observed to lobby (267 of them) and that have at least one draw matching $d_k$,

$$c(\nu) = \frac{100}{|K|} \sum_{k \in K} \mathbb{E}\Big[\ \frac{b_{\mathrm{coal},k}(\nu)}{\bar{s}_{\mathrm{coal},k}(\nu) + c_{\mathrm{coal},k}} \ \Big|\ \bar{d}_{k}(\nu) = d_k \Big],$$

where the expectation is the average over the simulation draws whose selected participation profile matches the observed one, and $c(\nu) = 0$ if $K$ is empty.

Each ratio is finite, because $|\bar{P} - \bar{P}^{-\ell}| \le 1$, the value is finite, and the denominator is at least the entry cost $c_{\ell k} > 0$; so $c(\nu)$ is finite for every $\nu \in \mathcal{N}$.

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

$\hat q(\nu)$ is the negative of the penalized simulated log-likelihood: the fit of the observed enactment $y_k$ and participation $d_k$, penalized by the gap between observed and model-predicted total lobbying expenditure for each of the four coalitions,

$$\hat{q}(\nu) = -\frac{1}{n} \sum_{k=1}^{n} \log \max\big\lbrace \hat{L}_k(\nu),\ \underline{L} \big\rbrace + \frac{\lambda}{n} \sum_{\ell=1}^{4} \left(1 - \frac{\hat{\varphi}_\ell(\nu)}{s_\ell}\right)^{2},$$

$$\hat{L}_k(\nu) = \frac{1}{S}\sum_{r=1}^{S} \mathbf{1}\big\lbrace \bar{d}_{k,r}(\nu) = d_k \big\rbrace \Big[ y_k \bar{P}_{k,r}(\nu) + (1 - y_k)\big(1 - \bar{P}_{k,r}(\nu)\big) \Big],$$

$$\hat{\varphi}_\ell(\nu) = \sum_{k=1}^{n} \frac{1}{S}\sum_{r=1}^{S} \big(\bar{s}_{\ell k, r}(\nu) + c_{\ell k}\big) \bar{d}_{\ell k, r}(\nu),$$

where $s_\ell$ is the observed total lobbying expenditure of coalition $\ell$, $c_{\ell k}$ the entry cost, $\underline{L}$ a floor on the simulated likelihood, and $S$ the number of simulation draws per policy, held fixed across evaluations so that $\hat q$ is a deterministic function of $\nu$.
A draw for which the entry game has no pure-strategy equilibrium contributes zero to both $\hat L_k$ and $\hat \varphi_\ell$.

The first term is bounded above by $-\log \underline{L}$ and below by $0$.
The second term is a finite sum of squares: equilibrium spending is finite for every $\nu \in \mathcal{N}$, because $P$ is bounded and concave in $T_f$, which forces $P_f \to 0$ as $T_f \to \infty$ and hence a finite fixed point of the effort equations.
So $\hat q$ is finite for every $\nu \in \mathcal{N}$.

<br>

## 5. Computational constraints (optional)

State the largest number of evaluations of $c(\cdot)$ and $\hat{q}(\cdot)$ you are willing to spend, for when either is expensive.
Fill in the blanks below.

At most 2,000 evaluations of $c(\cdot)$ and at most 2,000 evaluations of $\hat{q}(\cdot)$.

<br>

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

- $\lbrace (y_k, d_k, Z_k, X_k) \rbrace_{k=1}^{538}$ and the coalition positions — `Supplementary/DataandCode/Data/` in `temp/K2016/rdv029_Supplementary_Data.zip`, the files `data_outcome.txt` (enactment $y_k$, participation $d_k$), `data_explanatory.txt` and `data_polling.txt` (the $Z_k$ terms: constant, more regulation, less regulation, more spending, salience, opinion), and `data_explanatory2.txt` (the position of each coalition, which fixes the sign of $v_{\ell k}$, and the relevance dummy $X_{\ell k}$).
  Read by `specification.py:LoadData`.
- $\lbrace s_\ell \rbrace_{\ell=1}^{4}$ — `rst_spd` in `Supplementary/DataandCode/FORTRAN/dataread.f90` in `temp/K2016/rdv029_Supplementary_Data.zip`, read by `specification.py:LoadSpending`.
- $\delta, \sigma_\xi, \alpha, \sigma_\eta$, held fixed — `K2016_table7_estimates.csv`, columns `parameter` and `estimate`.
- $c_{\ell k} = 5{,}000$ dollars, $\lambda = 50$, $\underline{L} = 10^{-16}$, $S = 5{,}000$, and the number of Gauss-Seidel steps — `config.toml`.
- The draws $(\xi_{k,r}, \eta_{k,r})$ — Kang's, from `Supplementary/Datafortheweb/RandomShockCode/shock.m` and `Supplementary/DataandCode/FORTRAN/lhood.f90` in `temp/K2016/rdv029_Supplementary_Data.zip`, generated by `specification.py:Draws`.
  Policy $k = 0, \dots, 537$ uses block $b = k \bmod 10$, the points $100 + 5{,}000 b$ through $5{,}099 + 5{,}000 b$ of the unscrambled Halton sequence in bases 2, 3, 5, 7 and 11; the first coordinate $u$ gives $\xi = 2u - 1$, the others give $\eta = \Phi^{-1}(u)$, and both are rounded to six decimals.
- $c(\nu)$ and $\hat q(\nu)$ — `specification.py:{Target,Criterion}`, both fields of the dictionary returned by `specification.py:Simulate`.
  A primitive is any object exposing `Value`, `MarginalFor` (which returns $P_f$), `MarginalAgainst` (which returns $-P_a$) and `spending_exponent`; `specification.py:RatioContest` is the working model.
- $\hat\alpha = (\hat\beta_f, \hat\beta_a, \hat\gamma)$ — the published estimates, rows `effectiveness_for`, `effectiveness_against` and `spending_exponent` of `K2016_table7_estimates.csv`, returned by `specification.py:Baseline`.
- Replication of the published Tables 8, 10, 12 and 13 at the published parameters — `output/replicate.txt`, produced by `replicate.py`.

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
