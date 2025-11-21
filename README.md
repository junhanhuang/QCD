# A Quantum-Inspired Single-cell Analysis Framework

A quantum-inspired graph clustering method that reformulates modularity maximisation as a QUBO problem, maps it to an Ising Hamiltonian, and solves it on a 400-pulse coherent Ising machine (CIM).

This repository enables full reproduction of all figures and tables in the paper:

Repository contents
- Complete, runnable Python code for the two baseline methods (Leiden and spectral clustering)  
- Complete R pipeline for Slingshot trajectory inference and dynbenchmark evaluation  
- Pre-computed cluster labels generated on a real 400-pulse coherent Ising machine for all five benchmark datasets using the quantum-inspired single-cell analysis framework

The core solver of the quantum-inspired single-cell analysis framework, which interfaces directly with the photonic hardware, is subject to institutional technology-transfer restrictions and is not publicly released. However, because downstream trajectory inference and benchmarking depend solely on the final cluster assignments, the entire analysis remains fully reproducible using the provided pre-computed labels.

Researchers interested in running the quantum-inspired single-cell analysis framework on their own coherent Ising machine are welcome to contact the corresponding author for collaboration or licensed access.
