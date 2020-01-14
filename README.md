# Linear model selection using DPMP on model space.

This implements code implements variable selection for $p>n$ in linear regression model using the Dirichlet process mixture prior on model space. 

The DPMP acts as a post processing analysis of the usual Gibbs-Sampling model exploration strategy. Few Gibbs steps are enough to have a good idea of models posterior distribution.