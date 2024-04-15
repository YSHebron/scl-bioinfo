# scl-bioinfo

A bioinformatics project that aims to extend existing computational methods to predict (problematic) Protein Complexes from noisy PPINs and improve their performance. Output will be a prediction pipeline containing an improved algorithm. Now done with benchmarking and starting improvements.

Root folder will contain the entrypoint to our pipeline.

## Working Notes

* Moved long-runner.sh and simple-trainer.py to personal backup.
* Moved randomized yeast data to personal backup.

## TODOs

* Filter integrated dataset before PC2P.
* Apply negatome filtering.
* Produce other pre- and post-processing steps for PC2P or the BSS family.
* Create evaluation that will benchmark main algo based on performance of other classical clustering methods (like MCL).
* Try other BSS family algorithms (GCC and CUBCO/+).
* Try other unsupervised (parameter-free) algorithms not in BSS family, use only latest ones.
* Optimization of PC2P using igraph.

## Usage

We will first detail how to use the intermediate programs that prepare the datasets.
*