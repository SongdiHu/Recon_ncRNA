# Recon_ncRNA - achARNement 2.0

In 2016, Tremblay-Savard et al. proposed an approach (achARNement 1.0) using substitution and basepair costs that would focus on simultaneously reconstructing the ancestor of 2 related ncRNA families, assuming that they were created by a duplication followed by subfunctionalization. This approach aimed to infer ancestors that can conform to both structures but not be as specialized to any of them as the extant sequences. To simplify the algorithm, their approach limited the number of dependencies between the two structures that were considered in the cost calculations. 

We present an improvement of that methodology, which considers more structural information (stacking relationships), a more elaborate free energy calculation (using RNAeval from the ViennaRNA package), and makes use of tree decomposition to partition the cost calculation into subproblems.
