# scl-bioinfo

A proteomics project that aims to **augment, extend, and integrate** parameter-free protein complex prediction methods to improve their performance on various protein interaction networks.[^1] See paper at [...].

[^1]: By "parameter-free", we are referring to the lack of tunable hyperparameters (comparable to inflation and expansion coefficients for Markov Clustering, etc.) for the component clustering algorithms. Parameters such as input files and match score thresholds are still available for ease of use and experimentation.

* Augment - add necessary pre- and post-processing steps to clustering algorithms.
* Extend - modify clustering algorithm itself as needed, possibly reducing time complexity.
* Integrate - perform ensemble clustering on multiple clustering algorithms.

We hypothesize that parameter-free and purely topological clustering algorithms can contend with parameter-laden clustering algorithms on various areas, especially in adapting to noisy PPINs and predicting overlapping complexes. General focus is on improving recall with minimal loss of precision, and minimizing or eliminating the need for additional biological information to boost performance.

The current output is a prediction pipeline named **P5COMP** that contains 3 component clustering algorithms, namely P5COMP [[1]](#1), CUBCO+ [[2]](#2), and ClusterOne [[3]](#3). **P5COMP** stands for "Parameter-free Pipeline for Predicting Problematic Complexes". P5COMP's contributions and improvements will be based on performance metrics such as F-scores and AUC-PR with emphasis on the biological significance of results. Benchmarked algorithms will include the unmodified component algorithms when ran independently and classical protein complex prediction methods.

## Requirements

**Important**: All commands that will be provided in this documentation must be run from the repository root.

This software is written in Python 3.11 but may be compatible with other versions. This requires various packages, for which a Conda environment `scl` is provided for easy setup.

To set up the `scl` environment using the bundled `environment.yml` file, do:

```sh
conda env create -f environment.yml
```

Activate the new environment using `conda activate scl`. Verify that the environment was installed correctly using `conda env list`.

The `scl-ray` environment from the `environment_ray.yml` can also be set up if the user wishes to run `code/PC2P/PC2P.py` (a component algorithm) with the `-p ray` parameter. This does PC2P parallelization using Ray instead of Multiprocessing. It is necessary to have a separate environment for Ray because Ray is still experimental for Python 3.11.

## Usage

A command line tool for Linux systems are provided, namely `pipeline1.sh`.[^2] For quick help on its usage, try running `./pipeline1.sh -h`. This will display the following message:

[^2]: The pipeline CLI for Windows `pipeline1.ps1` is still under construction.

```sh
usage: ./pipeline1.sh [-p [ppinfile]] [-r [reffile]] [-o [outputdir]] [-n [negfile]] [-f [filter]] [-h]
    
Runs P5COMP on the given PPIN file (ppinfile) and evaluates against the given gold standard (reffile).
Final predicted clusters will be written in outputdir.
Important: Protein names (PID) should be in gene name (ordered locus) or KEGG format (ex. YLR075W) to match gold standards.

options:
    -p [ppinfile]       path to PPIN file (.txt) where each row is (u v w) (required)
    -r [reffile]        path to gold standard or reference complexes file (.txt) (required)
    -o [outputdir]      path to output directory (required)
    -n [negfile]        path to negatome (.txt) where each row is (u v) (optional)
    -f [filter]         filtering type (perpair or perprotein)
    -a [attribs]        attributes for evaluation file name of format 'algo-goldstd-ppin', ex: P5COMP-CYC-Collins
    -r [resfile]        path to results file for evaluation
    -h                  show this help information

```

### Gold Standards and PPINs

This pipeline comes packaged with Yeast gold standards and PPINs for performance evaluation. You may find them in `data/Yeast`. The gold standards used are CYC2008 and SGD, while the PPINs used are Collins, Gavin, KroganCore, KroganExt, and BIM (union of BioGRID, IntACT, and MINT) listed in increasing number of edges. In addition, Negatome 2.0 (a list of high-confidence non-interacting proteins) is also provided to help filter out false positives from the PPINs.

### Example

Here we provide an example of how to use `pipeline1.sh`, from entry to clustering to performance evaluation. We have placed copies of the Yeast Gold Standards, PPINs, and Negatome in the `eval` directory for more convenient access to them.

Now follows is a sample command to run `pipeline1.sh` with the Collins PPIN and evaluate the results against the CYC2008 gold standard. The pipeline's produced clusters can be found in `results` named as `P5COMP-CYC-Collins_clusters.txt` written in `(<cluster_score>_<cluster_length>): proteins` line format. The evaluation results i.e. the performance metrics can be found in comma-separated tabular format in `results/results.csv`. The `-f` attribute indicates that this run will use `perpair` filtering. The alternative is `perprotein` which yields nearly the same results but with worse runtime.

```sh
./pipeline1.sh -p eval/Collins.txt -r eval/CYC2008.txt -o results \
    -n eval/Negatome.txt -f perpair -a "P5COMP-CYC-Collins" -R results/results.csv
```

If needed, save the Results of an important run to a different directory to avoid them getting overwritten when the pipeline is rerun.

Note that `pipeline2.sh` is used for evaluating PC2P, CUBCO+, and ClusterOne and has the same interface as `pipeline1.sh`.

### Batch Evaluation

A batch file is provided to perform batch evaluation of all Method, Gold Standard, and PPIN combinations. One can run this by simply entering `./batch.sh`, a shell script that wraps `metrics.sh` in terms of number of batches (default set to 1).

## Proposed Pipeline

This section illustrates in brief the proposed P5COMP pipeline. See the SCL article for more details.

```sh
DataPreparation >> { Preprocessing >> {Denoising, Hub Decomposition} | Clustering >> {PC2P, CUBCO+, ClusterOne} | Postprocessing >> {Return Hub Proteins, Ensemble Clustering} } >> Evaluation
```

* `DataPreparation` (with filtering) will involve `data/Yeast (PPIN Name) > data/Yeas (FilteredPPINs)t`.
* `Preprocessing`, `Clustering`, and `Postprocessing` will involve `data/Yeast > data/Results`. Human protein clusters are nice-to-haves.
* `Evaluation` will involve `data/Results > data/Analysis`.

Shown below is the current structure of the `data` directory:

```sh
data/
├─ Yeast/ # raw and prepared Yeast PPINs and goldstds
├─ Human/ # raw and prepared Human PPINs and goldstds
├─ Results/ # predicted clusters
├─ Analysis/ # performance metrics evaluation and plots
```

## TODOs

* Granular evaluation of component algorithms with DECOMP.
* Try other CNP-based algorithms (e.g. GCC) and try other unsupervised (parameter-free) algorithms not in CNP family.
* Verify stochastic elements of CUBCO+.
* Try swapping component algorithms or adding to component algorithms with other clustering algos.
* Benchmark against the results of INTEG [[4]](#4) by Yong and Wong (if it can be run, it's written in Perl).
* Try score-based ensemble clustering (currently we are using voting-based ensembling).
* Regenerate precision-recall plots.
* Optimization of CNP-based algorithms using `igraph`.
* Programs only currently tested on Yeast PPINs. Try to scale to Human PPINs and other algorithms.
* Further parallelization of pipeline.
* Use strategic imports from a `main` python file instead of doing the `argparse` route. `main` should be the new entrypoint instead of the shell or Powershell scripts. Alternatively, try setting up everything in Jupyter Notebooks (note that .ipynb can also run shell commands from the notebook environment).
* Make this a pip-installable software.

---

## References

<a id="1">[1]</a> 
Sara Omranian, Angela Angeleska, and Zoran Nikoloski. “PC2P: Parameter-Free Network-Based Prediction of Protein Complexes”. In: Bioinformatics 37.1 (Apr. 2021), pp. 73–81. issn: 1367-4803. doi: 10.1093/bioinformatics/btaa1089. url: https://doi.org/10.1093/bioinformatics/btaa1089 (visited on 06/08/2024).

<a id="2">[2]</a> 
Sara Omranian and Zoran Nikoloski. “CUBCO+: Prediction of Protein Complexes Based on Min-Cut Network Partitioning into Biclique Spanned Subgraphs”. In: Applied Network Science 7.1 (Dec. 2022), pp. 1–12. issn: 2364-8228. doi: 10.1007/s41109-022-00508-5. url: https://appliednetsci.springeropen.com/articles/10.1007/s41109- 022-00508-5 (visited on 06/08/2024).

<a id="3">[3]</a> 
Tam ́as Nepusz, Haiyuan Yu, and Alberto Paccanaro. “Detecting Overlapping Protein Complexes in Protein-Protein Interaction Networks”. In: Nature Methods 9.5 (May 2012), pp. 471–472. issn: 1548-7105. doi: 10.1038/nmeth.1938. url: https://www.nature.com/articles/nmeth.1938 (visited on 06/08/2024).

<a id="4">[4]</a> 
Chern Han Yong and Limsoon Wong. “Prediction of Problematic Complexes from PPI Networks: Sparse, Embedded, and Small Complexes”. In: Biology Direct 10.1 (Aug. 2015), p. 40. issn: 1745-6150. doi: 10.1186/s13062-015-0067-4. url: https://doi.org/10. 1186/s13062-015-0067-4 (visited on 06/09/2024).

---

Scientific Computing Laboratory\
Department of Computer Science\
University of the Philippines Diliman

2nd Semester A.Y. 2023-2024

© Course Materials and Advising by Dr. John Justine Villar
