# ARTF-ELAP
This algorithms consists in a a three-phase approach for computing a fundamental cycle basis of minimal length in a graph.

## Long description:
This algorithm addresses the combinatorial optimisation problem (COP) of computing a fundamental cycle basis (FCB) of minimum length in a connected graph, the length being the sum of the number of arcs of the cycles of the FCB. This COP being NP-hard, we propose a three-phase heuristic for computing an approximate solution. The first phase is constructive and consists in computing an initial FCB. The second is an improving iterative procedure permitting to reduce the length of the former FCB and leading, in general, to a non fundamental cycle basis (NFCB), whereas the third phase consists in transforming the NFCB into a final FCB better than the initial one.
