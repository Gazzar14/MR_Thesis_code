## Bounding ACE via Column Generation

### Dynamic Column Generation via the Exclusion Restriction

 The number of
latent response types can become extremely large as the cardinalities of
the instrument, exposure, and outcome increase. Let $|G|=k$ denote the
number of instrument levels, $|X|=m$ the number of exposure levels, and
$|Y|=l$ the number of outcome levels. A latent response type consists of
an exposure response function $X^G$ and an outcome response function
$Y^X$, giving a total of

$$
|U| = m^k l^m.
$$

For example, with $k=10$, $m=10$, and $l=2$, the number of latent
response types is

$$
10^{10}2^{10}=10{,}240{,}000{,}000{,}000.
$$

Explicitly enumerating all such response types to construct the LP is
therefore computationally impractical. While various papers have proposed
numerical approaches to bounding the ACE, some numerical methods are prone
to stability issues and getting stuck in local optima. We instead
retain the linear programming formulation and use column generation based
on the Dantzig--Wolfe decomposition.

The Restricted Master Problem (RMP) is initialized with a limited set of
columns, and additional latent response types are generated as needed by
solving the associated pricing problem.

The implementation exploits the structure of the latent response functions
under the exclusion restriction. In particular, a latent response type can
be represented as

$$
u=(x^G,y^X),
$$

where

$$
x^G=(x^{g_1},\ldots,x^{g_k})
$$

specifies the exposure that would be observed at each instrument level, and

$$
y^X=(y^{x_1},\ldots,y^{x_m})
$$

specifies the potential outcome associated with each exposure level.

Under the exclusion restriction, the outcome response function does not
directly depend on the instrument. This means that we can fix a candidate
$Y^X$ and then optimize the contribution of the exposure response function
$X^G$ to the reduced cost separately for each instrument level.

We first obtain the dual solution by solving the current RMP. Given this
dual solution, the pricing problem searches for a latent response type with
the smallest reduced cost. The reduced cost indicates whether introducing
the corresponding column into the RMP can improve the current objective.

The algorithm splits this problem into two stages. First, it enumerates the
possible outcome response functions $Y^X$. For each fixed $Y^X$, the
exposure response function $X^G$ can be constructed without enumerating the
$m^k$ possible exposure response functions.

Specifically, for each instrument level $g\in G$, the algorithm considers
each possible exposure $x\in X$. Since $Y^X$ is fixed, choosing $x$
determines the corresponding outcome $y^x$, and hence the relevant observed
data cell $(g,x,y)$. The exposure level associated with the largest relevant
dual value is selected independently for each $g$, thereby constructing the
exposure response function $X^{G*}$ that minimizes the reduced cost
conditional on the current $Y^X$.

Thus, each possible $Y^X$ gives rise to one candidate latent response type,

$$
u(Y^X)=\left(X^{G*},Y^X\right),
$$

whose reduced cost is then evaluated. After all possible $Y^X$ functions
have been considered, the candidate with the smallest reduced cost is
selected. If this reduced cost is below $-\epsilon$, the corresponding
column is added to the RMP; otherwise, the pricing step terminates.

Consequently, the implementation enumerates $l^m$ possible outcome response
functions, while the exposure response function for each candidate is
obtained through $k$ independent searches over the $m$ possible exposure
levels. For binary outcomes, this corresponds to enumerating $2^m$ outcome
response functions rather than the full $m^k2^m$ joint latent response
types.

The resulting column-generation procedure can be summarized as follows:

1. **Initialization:** Initialize the RMP with artificial variables assigned
   a large objective penalty and solve the RMP to obtain the corresponding
   dual variables.

2. **Pricing:** Enumerate the possible outcome response functions $Y^X$.
   For each $Y^X$, construct the corresponding exposure response function
   $X^{G*}$ by independently selecting, for each $g\in G$, the exposure
   $x\in X$ with the largest relevant dual contribution. Evaluate the
   reduced cost of the resulting latent response type.

3. **Column selection:** Select the candidate latent response type with the
   smallest reduced cost across all enumerated $Y^X$ functions.

4. **Column addition and termination:** If the minimum reduced cost is less
   than $-\epsilon$, add the corresponding column to the RMP and resolve it
   to obtain updated dual variables. Repeat the pricing and column-selection
   steps. The procedure terminates when the minimum reduced cost is greater
   than or equal to $-\epsilon$.

This algorithm does not reduce the size of the underlying latent response
space. Instead, it avoids explicitly constructing that space during the LP
solution. Its computational advantage comes from exploiting the separability
of the exposure-response component of the pricing problem. The implementation
replaces explicit enumeration of the $m^k$ possible exposure response
functions with $k$ independent searches over the $m$ exposure levels for
each candidate $Y^X$.

### Convergence Tolerance and Finite Samples

The column-generation procedure uses a numerical tolerance $\epsilon$ to
determine when the pricing problem has sufficiently converged. After
considering all possible $Y^X$ response functions, the algorithm terminates
when the smallest reduced cost is greater than or equal to $-\epsilon$,
where $\epsilon$ is a small positive tolerance, such as $10^{-6}$.

If $\epsilon=0$ and the pricing problem is solved exactly, the absence of a
negative reduced cost establishes optimality of the current RMP for the
corresponding LP. With a positive $\epsilon$, the procedure instead
terminates when no column with a reduced cost below $-\epsilon$ is
identified. Thus, the resulting solution is subject to the specified
reduced-cost tolerance and is not guaranteed to be exactly optimal for
the LP.

In practice, sampling variation can result in an empirical distribution
that does not lie exactly within the population IV model's feasible
polytope. Imposing the empirical probabilities as exact equality
constraints may cause infeasibility even if the corresponding population
distribution satisfies the IV model.

To allow for small discrepancies arising from finite-sample variation, the
implementation permits the empirical constraints to be relaxed by a
tolerance $\delta$. Instead of imposing

$$
Ax=b,
$$

the RMP solves the pair of inequalities

$$
b-\delta\mathbf{1}\leq Ax\leq b+\delta\mathbf{1}.
$$

The value of $\delta$ determines the extent to which the empirical
probabilities are allowed to differ from the probabilities generated by
the latent response types. Consequently, the resulting optimization
problem is a relaxed version of the LP defined using the exact empirical
probabilities.

If $L^*$ and $U^*$ denote the lower and upper bounds under the exact
constraints, and $L_\delta$ and $U_\delta$ denote the corresponding optima
under the relaxed constraints, then

$$
L_\delta \leq L^*
\qquad\text{and}\qquad
U_\delta \geq U^*.
$$

The use of $\delta$ is distinct from the column-generation tolerance
$\epsilon$. Column generation determines how the LP is solved without
explicitly enumerating all latent response types, while $\epsilon$
determines the numerical stopping criterion for the pricing problem and
$\delta$ determines the extent to which the empirical constraints are
relaxed.

Setting $\delta=0$ recovers the original LP constraints, while setting
$\epsilon=0$ removes the reduced-cost tolerance from the column-generation
stopping criterion. Thus, when both $\epsilon$ and $\delta$ are zero, the
procedure targets the sharp bounds of the original LP, subject to the
numerical precision of the LP solver.

### Other Features

The pricing procedure can also incorporate restrictions on the set of
possible outcome response functions. For example, under a Monotone
Treatment Response assumption, the potential outcomes are required to
satisfy

$$
y^{x_1}\leq y^{x_2}\leq\cdots\leq y^{x_m}
$$

for ordered exposure levels. Rather than considering all $l^m$ possible
outcome response functions, only those satisfying the monotonicity
restriction are enumerated.

The number of weakly increasing sequences of length $m$ taking values in
$l$ ordered outcome levels is

$$
\binom{m+l-1}{m}.
$$

For binary outcomes, this reduces to $m+1$ possible response functions.
The shape restriction can therefore reduce the enumeration required by
the pricing step.


# References
Dantzig, George B. andWolfe, Philip. Decomposition Principle for Linear Programs. Operations
Research, 8(1):101–111, 1960. doi: 10.1287/opre.8.1.101.

Kilbertus, Niki, Kusner, Matt J., and Silva, Ricardo. A Class of Algorithms for General Instrumental
Variable Models, 2020. URL https://arxiv.org/abs/2006.06366.
