# scl-bioinfo

A proteomics project that aims to **augment, extend, and integrate** parameter-free protein complex prediction methods to improve their performance on various protein interaction networks. Focus is on adapting to noisy PPINs and predicting overlapping complexes through purely topological graph clustering. Effort is to be made to not include additional biological information such as ontologies and other annotations, and to improve recall without sacrificing precision.

* Augment - add necessary pre- and post-processing steps to clustering algorithms.
* Extend - modify clustering algorithm itself as needed, possibly targeting time complexity.
* Integrate - perform ensemble clustering on multiple clustering algorithms.

Output will be a prediction pipeline that may integrate more than one component unsupervised clustering algorithm. Improvements will be based on general performance metrics such as F-scores and AUC-PR with emphasis on biological significance of results. Benchmarked algorithms will include the unmodified component algorithms when ran independently and classical protein complex prediction methods.

## Usage

A pipeline CLI for Linux systems are provided, namely `pipeline.sh`.[^1] For quick help on their usage, try running `./pipeline.sh -h` if you're on Linux. This will display the following message:

[^1]: The pipeline CLI for Windows `pipeline.ps1` is still under construction.

```sh
usage: ./pipeline.sh [-p [ppinfile]] [-r [reffile]] [-o [outputdir]] [-n [negfile]] [-h]
    
Runs P5COMP on the given PPIN file (ppinfile) and evaluates against the given gold standard (reffile).
Final predicted clusters will be written in outputdir.
Important: Protein names (PID) should be in gene name (ordered locus) or KEGG format (ex. YLR075W) to match gold standards.

options:
    -p [ppinfile]       path to PPIN file (.txt) where each row is (u v s) (required)
    -r [reffile]        path to gold standard or reference complexes file (.txt) (required)
    -o [outputdir]      path to output directory (required)
    -n [negfile]        path to negatome (.txt) where each row is (u v) (optional)
    -h                  show this help information"
```
No help message is currently available for `pipeline.ps1` on Windows. A proper tutorial is still being written. The rest of the Usage section is outdated.

For `code/PC2P` and other algos pending to be implemented here, the pipeline steps are described below. Note that preprocessing (e.g. filtering), clustering, postprocessing (e.g. ensembling), and evaluation steps are still separate. Please take note of the required file name formats.

* Prepare the datasets using `filter_then_score.py` and other included pre-processing programs.
* In particular for `filter_then_score.py`, it performs perprotein and perpair, and/or direct filtering on the selected dataset. Selected dataset has no required filename. `filter_then_score.py` has parameters `ppinname`, `gldstdname`, . We currently only have two `gldstdname` for Yeast, `CYC` and `SGD`. Use `Corum` for human. The outputfile will have the following format: `<ppinname>_<gldstdname>_<filtering>_weighted.txt`. File name may differ based on other PPIN denoising or manipulation step. Each line here has the format `p1 p2 score`. A tab-delimited format (`.tsv`) are also available particularly for use by FINCH-Clustering.
* Afterwards, run `code/PC2P.py`, or any clustering algorithm (may be scripted or ensembled) on the prepared dataset. The required input format is `<ppinname>_<gldstdname>_<filtering>_weighted.{txt,csv,tsv}`. The output file will be named `<ppinname>_<gldstdname>_<filtering>_predicted.txt` (note the change in stem suffix) and each line denoting a single predicted cluster will have the format `(size_score): p1 p2 ...` where `size` is the number of proteins in the cluster and `score` is the weighted sum of the PPI scores of the proteins in the cluster.
* Insert ensemble clustering instructions.
* Finally, to analyze the (working) results, run `code/PC2P_eval.py`. The required input file format is `<ppinname>_<gldstdname>_<filtering>_predicted.txt`. There are currently two output files produced: the first one is `<ppinname>_<gldstdname>_<filtering>_eval.txt` (note the change in stem suffix) which contains precision-recall datapoints for each score threshold (based on recall levels) and the second one is an overall `auc_only.txt` file that contains a summary of the performance metrics (NOTE: rename to `summary.txt`). The output file line format is still not fixed and may be reformatted to aid in the plotting. Currently, no plots are being produced.

```sh
# Sample run from start to finish, no shell scripting.
# filter_then_score.py is to be run via a shell script.
# Run clustering (mp parallel)
python code/PC2P/PC2P.py data/Yeast/FilteredPPINs/Collins_CYC_direct_weighted.txt data/Results/Dummy -p
# Run evaluation
python code/PC2P/PC2P_eval.py data/Results/Dummy/Dummy_CYC_testonly_predicted.txt data/Yeast/CYC_complexes.txt data/Analysis/Dummy
```

* Note that the program names here are tentative due to the ongoing nature of the research.
* If needed, save the Results and Analysis outputs of an important run in a different directory to avoid them getting rewritten when the pipeline is rerun.
* For further help, run the programs with the `--help` flag.

### Proposed Pipeline

```sh
DataPreparation >> { Preprocessing >> {Hub Decomposition} | Clustering >> {PC2P, CUBCO+, ClusterOne}* | Postprocessing >> {Return Hub Proteins, Ensemble Clustering} } >> Evaluation
```

* `DataPreparation` (with filtering) will involve `data/Yeast (PPIN Name) > data/Yeas (FilteredPPINs)t`.
* `Preprocessing`, `Clustering`, and `Postprocessing` will involve `data/Yeast > data/Results`. Human protein clusters are nice-to-haves.
* `Evaluation` will involve `data/Results > data/Analysis`.
* Chosen ensemble clustering algorithms may involve DECOMP and MCL as these are technically unsupervised and effectively parameter free or at least tuning free (MCL parameters will be set without further tuning).

Shown below is the current structure of the `data` directory:

```sh
data/
├─ Yeast/ # raw and prepared Yeast PPINs and goldstds
├─ Human/ # raw and prepared Human PPINs and goldstds
├─ Results/ # predicted clusters
├─ Analysis/ # performance metrics evaluation and plots
```

## TODOs

* Granular.
* Create evaluation that will benchmark main algo based on performance of other classical clustering methods (like MCL).
* Try other BSS family algorithms (GCC and CUBCO/+).
* Try other unsupervised (parameter-free) algorithms not in BSS family, use only latest ones.
* Optimization of CNP-based algorithms using igraph.
* Programs only currently tested on Yeast PPINs. Try to test on Human PPINs.
* Further parallelization of pipeline.
* Use strategic imports from a `main` python file instead of doing the `argparse` route.

---

Department of Computer Science\
University of the Philippines Diliman

2nd Semester A.Y. 2023-2024

© Course Materials and Advising by Dr. John Justine Villar
