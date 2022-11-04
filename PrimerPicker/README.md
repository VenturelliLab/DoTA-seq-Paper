# MplexPrimerPicker
Automatic picking of multiplex PCR primers for DOTA-seq

1. Starting with a list of target sequences to design primers for in fasta format, use 1.Primer3-primergenerator.ipynb to generate candidate primers for each gene.

2. Use 2.Calculate-dimer-free-energy.ipynb to calculate free energies of dimerization between all primerset pairs. The script uses ntthal, which is a subcomponent of Primer3 to do the free energy calculations between all primer combinations between every combination of two primer sets (4 primers in total).

3. Use 3.Simulated-annealing-primer-pooling.ipynb to generate pools of primersets that have low primer dimerization potentials. The algorithm uses parallel computing to run simulated annealing on a large number of randomly generated starting primer pools and reports the 10 best performaning primersets based on average deltaG of dimerization.

