# A Quantum-Inspired Single-cell Analysis Framework

A quantum-inspired graph clustering method that reformulates modularity maximisation as a QUBO problem, maps it to an Ising Hamiltonian, and solves it on a 400-pulse coherent Ising machine (CIM).

This repository enables full reproduction of all figures and tables in the paper:

Repository contents
- Complete, runnable Python code for the two baseline methods (Leiden and spectral clustering)  
- Complete R pipeline for Slingshot trajectory inference and dynbenchmark evaluation  
- Pre-computed QCD cluster labels generated on a real 400-pulse coherent Ising machine for all five benchmark datasets  

The core QCD solver that interfaces directly with the photonic hardware is subject to institutional technology-transfer restrictions and is not publicly released. However, because downstream trajectory inference and benchmarking depend solely on the final cluster assignments, the entire analysis remains fully reproducible using the provided pre-computed QCD labels.

Researchers interested in running QCD on their own coherent Ising machine are welcome to contact the corresponding author for collaboration or licensed access.
