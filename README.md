# scl-bioinfo

A bioinformatics project that aims to extend parameter-free computational protein complex prediction methods. Output will be a prediction pipeline containing an improved algorithm that may integrate more than one unsupervised clustering algorithm. Improvements will be based on general performance metrics, currently AUC-PR and F-measaure, but may be composited to include more metrics. Benchmarked algorithms will include the original unmodified algorithms and classical protein complex prediction methods. The Integrated pipeline by Yong and Wong will remain part of the benchmark.

Root folder will eventually contain the entrypoint to our pipeline.

## Usage

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
DataPreparation >> { Preprocessing >> {Hub Decomposition} | Clustering >> {PC2P, FINCH, ONCQS, MCL}* | Postprocessing >> {Return Hub Proteins, Ensemble Clustering} } >> Evaluation
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

## Working Notes

* Moved long-runner.sh and simple-trainer.py to personal backup.
* Moved randomized yeast data to personal backup.
* In preparation for further experiments on generalized clustering algorithms, we have made separate folders inside the `data/Results/` and `data/Analysis/` labeled by which clustering algorithm produced them.

## TODOs

* Wrap up TODOs in initial Gantt Chart.
* Filter integrated dataset before PC2P.
* Apply negatome filtering. Sources: https://github.com/Zhaonan99/NVDT, https://github.com/mattilalab/petrov-et-al-2020, https://github.com/lafleur1/embeddingsPPI/blob/2743925fb9d6489dbe68fc6d0cc70eea8c955a9f/ppiDB/Negatome/Negatome2_combinedStringent.txt
* Produce other pre- and post-processing steps for PC2P or the BSS family.
* Create evaluation that will benchmark main algo based on performance of other classical clustering methods (like MCL).
* Try other BSS family algorithms (GCC and CUBCO/+).
* Try other unsupervised (parameter-free) algorithms not in BSS family, use only latest ones.
* Optimization of PC2P using igraph.
* Programs only currently tested on Yeast PPINs, if possible try to make room for Human PPINs.

## On Data

**Important:** `data/curated` has been completely git ignored due to the large sizes of its contents. If you need said data, please contact the authors. If you are an author, please see our `Datasets` Teams folder. A ZIP file there is also being constantly updated containing the complete datasets that we use. The directory and ZIP file contains the following:

* BIOGRID_3.0.67_2010
* BIOGRID_4.4.227_2023
* Cayetano
* CYC2008_2008
* DIP_2017
* GO_Yeast_2023
* iRefIndex_2023
* SGD_GSE3431_2010
* STRING_2023
* Yong and Wong

`data/curated` contains raw, preprocessed, enriched, as well as reference Protein Interaction and Complex datasets for the target organism, the yeast *S. Cerevisiae*. A statistical overview of the datasets may be found in `data/overview`. Data collection is guided by the related literature, mainly Srihari (2017), Omranian (2022), and Cayetano (2022).

Data here will be fed to the **prediction pipeline's** entry point at `code`. Intermediate data may be produced and fed back to this directory. Legacy data are obtained for fair comparison with previous papers.

`data/curated` currently contains the yeast datasets below, organized into PPI datasets and protein complex gold standards (large files are gitignored).

Yeast PPI Datasets (and GO annotation):

* BIOGRID 3.0.67 (2010) (large files)
* BIOGRID 4.4.227 (2023) (large files)
* DIP PPIN (2017) (large files)
* iRefIndex (2023)
* STRING (2023)
* GO S. Cerevisiae (2023) - topology and annotations

Yeast Gold Standards:

* CYC2008 (2008) (`data/Yeast/CYC_complexes.txt`)
* SGD GSE3431 (2010) (`data/Yeast/SGD_complexes.txt`)

`data/Yeast/CYC_complexes.txt` and `data/Yeast/CYC_complexes_integrated.txt` are the same. If you need the gitignore'd files, you may find it in our Teams.

### Cayetano (2022) Datasets

`data` also contains the raw DIP and SWC datasets used by Cayetano, procured directly from [https://github.com/avancayetano/cs199-bioinformatics]. You may find these in `data/curated/Cayetano`.

### Omranian (2022) Datasets

Omranian's comparative study used 4 Yeast datasets and 2 Human datasets. The Yeast datasets four are comprised of Collins, Gavin, KroganCore, and KroganExt, which are used for the protein complex prediction. The predicted complexes are then benchmarked against the fifth one, the gold standard CYC2008. You may find the first four in `data/curated/Omranian`, while CYC2008 are already in `data`.

**NOTICE**: The Omranian datasets on Yeast will now be our main datasets. Further modifications can be done to edit the scoring or compositiong of the datasets.

### On Edge Weights

Note that the edge weights between proteins may indicate any of the following (as specified by the dataset):

* Affinity or strength of the interaction
* Reliability of the interaction (obtained experimentally)
* Distance of the interaction, as used with shortest paths computation

Please consult the data source for specifications.

### Dataset Updates

We are expanding our testing from just S. Cerevisiae PPINs to also include H. Sapiens, D. Melanogasater, E. Coli, and C. Elegans datasets.

#### PPIN Sources for S. Cerevisiae

[Insert data sources here from IntAct or some other Interactome database]
