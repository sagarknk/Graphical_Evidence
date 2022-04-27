# Graphical_Evidence

This repository contains code to implement our proposed technique to compute marginal likelihood in Gaussian graphical models. 

Our proposed technique works for a broad class of priors. Specifically, the requirements are: (a) the priors on the diagonal terms 
on the precision matrix can be written as gamma or scale mixtures of gamma random variables and (b) those on the off-diagonal terms 
can be represented as normal or scale mixtures of normal. 

In this work, we mainly focus on three priors:
    1. Bayesian graphical lasso (BGL) 
    2. Graphical horseshoe (GHS) and 
    3. G-Wishart

We validate our proposed technqiue when the prior on the precision matirx is Wishart. As the marginal likelihood is available in a closed form 
in this case, this provides a useful validation. We also compare the maringal likelihood estimates obtained using our proposed technique against
the estimates obtained by annealed importance sampling (Neal,2001), nested sampling (Skilling, 2006) and harmonic mean estimates (Newton and Raftery, 1994). 

There are four directories in this repository, one for each prior (BGL, GHS, G-Wishart and Wishart). More detailed instructions are within the README.txt files
of the respective directories. 
