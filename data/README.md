# Data

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

* CYC2008 (2008)
* SGD GSE3431 (2010)

If you need the gitignore'd files, you may find it in our Teams.

`data/results` will contain the results of this project, which are performance evaluations and such.

## Cayetano (2022) Datasets

`data` also contains the raw DIP and SWC datasets used by Cayetano, procured directly from [https://github.com/avancayetano/cs199-bioinformatics]. You may find these in `data/curated/Cayetano`.

## Omranian (2022) Datasets

Omranian's comparative study used five (5) datasets. The first four are comprised of Collins, Gavin, KroganCore, and KroganExt, which are used for the protein complex prediction. The predicted complexes are then benchmarked against the fifth one, the gold standard CYC2008. You may find the first four in `data/curated/Omranian`, while CYC2008 are already in `data`.

## On Edge Weights

Note that the edge weights between proteins may indicate any of the following (as specified by the dataset):

* Affinity or strength of the interaction
* Reliability of the interaction (obtained experimentally)
* Distance of the interaction, as used with shortest paths computation

Please consult the data source for specifications.
