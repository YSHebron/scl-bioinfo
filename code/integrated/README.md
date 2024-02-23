
# Integrated Pipeline (SWC+DECOMP+SSS)

This program is described in Yong CH, Wong L, "Discovery of problematic protein complexes from PPI networks: sparse, embedded, and small complexes" [1].

This program uses four third-party clustering program, CMC [2], Coach [3], ClusterOne [4], and IPCA [5]. In its default setting this program will use these four clustering programs to generate clusters. Thus, if you make use of this program in its default setting, please cite, at the minimum, the above five publications [1-5].

If you modify this program (via the batch file) to change the clustering programs used, please cite [1] as well as the relevant references for the clustering programs used [2-5].

## REQUIREMENTS

This program is provided as a set of Perl scripts. Thus it requires a Perl interpreter to be installed, such as [ActivePerl](www.activestate.com/activeperl).

This program comes bundled with four third-party clustering programs:

1. CMC [2] <http://www.comp.nus.edu.sg/~wongls/projects/complexprediction/CMC-26may09/>
2. Coach [3] <http://www1.i2r.a-star.edu.sg/~xlli/coach.zip>
3. ClusterOne [4] (licensed under GNU General Public License v3, sourcecode available at <https://github.com/ntamas/cl1>)
4. IPCA [5] <http://netlab.csu.edu.cn/bioinfomatics/limin/IPCA/>

Java is required to run the ClusterOne algorithm.

This program was developed and tested in the Windows XP/7/8 environment; it may work in other Windows versions as well, or even other operating systems.

## INSTALLATION

Unzip the software package into its own directory, and run it from the command line.

## EASIEST WAY TO RUN THIS PROGRAM

The easiest way to run this program in Windows is to use the two prepared batch files.

### 1. Complex Prediction Batch File

To generate complexes using all the reference complexes for learning, use `run_complexprediction.bat`, with six (6) string arguments:

1. the input data file,
2. the reference complexes file,
3. the GO scheme file,
4. the GO annotations file,
5. the output directory,
6. and the parameter set (either "yeast" or "human").

For **yeast**, using the supplied input files, and the `yeast_predict_output/` output directory, the command would be:

```shell
.\run_complexprediction.bat data_yeast.txt complexes_yeast.txt go_scheme_all_mine.txt my_go_associations.sgd_propagated.txt yeast_predict_output yeast
```

NOTE: complexes_yeast.txt is equivalent to the CYC2008 reference complex.

For **human**, using the supplied input files, and `human_predict_output/` output directory, the command would be:

```shell
.\run_complexprediction.bat data_human.txt complexes_human.txt go_scheme_all_mine.txt my_go_associations.human_propagated.txt human_predict_output human
```

The output files (scored edges and clusters) are stored in the output directory:

- .\SWC\ : output files from running SWC
- .\DECOMP\ : output files from running DECOMP
- .\SSS\ : output files from running SSS
- .\clusters_integrated.txt : final integrated clusters

### 2. Complex Prediction AND Evaluation Batch File

To run cross-validation experiments to get the performance figures as shown in the paper, use `run_xval.bat`, with seven (7) arguments:

1. the input data file,
2. the reference complexes file,
3. the GO scheme file,
4. the GO annotations file,
5. the cross-validation file,
6. the output directory, and
7. the parameter set (either "yeast" or "human").

For yeast, with a "yeast_xval_output\" output directory, the command would be:

```shell
.\run_xval.bat data_yeast.txt complexes_yeast.txt go_scheme_all_mine.txt my_go_associations.sgd_propagated.txt xval_yeast.txt yeast_xval_output yeast
```

For human, with a "human_xval_output\" output directory, the command would be:

```shell
.\run_xval.bat data_human.txt complexes_human.txt go_scheme_all_mine.txt my_go_associations.human_propagated.txt xval_human.txt human_xval_output human
```

The output files (scored edges, clusters, and results) are stored in the output directory:

- `clusters integrated iterN.txt` : integrated clusters for every cross-validation iteration
- `results clusters integrated.txt` : results for integrated clusters
- `results clusters swc.txt` : results for SWC clusters
- `results clusters decomp.txt` : results for DECOMP clusters
- `results clusters sss.txt` : results for SSS clusters

## ADVANCED USE

- The parameter set option (ie. the sixth input parameter in `run_complexprediction.bat`, and the seventh input parameter in `run_xval.bat`, either "yeast" or "human") simply sets four parameters in the batch file (ngo, nhub, decompweight, and sssweight, described in the paper). To change these values you can simply modify the batch files to do so.
- To run the program with your own data files / reference complexes, simply run the batch files, ensuring that the file formats are the same as those provided.
- To generate your own topological data from the primary data sources (PPI, STRING, literature co-occurrence) is not very straightforward. I have some perl scripts that do this but these are not provided in this package. Please email me at <chernycherny@hotmail.com> and I can send them to you.

## CONTACT

Please contact me at <chernycherny@hotmail.com> for comments or questions.

## References

[1] Yong CH, Wong L: Discovery of problematic protein complexes from PPI networks: Sparse, embedded, and small complexes. Biology Direct 10:40 (2015).

[2] Liu G, Wong L, Chua HN: Complex discovery from weighted PPI networks. Bioinformatics 25(15), 1891-1897 (2009).

[3] Wu M, Li X, Kwoh CK, Ng SK.: A core-attachment based method to detect protein complexes in PPI networks. BMC Bioinformatics 10, 169 (2009).

[4] Nepusz T, Yu H, Paccanaro A: Detecting overlapping protein complexes in protein-protein interaction networks. Nat Methods 9, 471-472 (2012).

[5] Li M, Chen J, Wang J, Hu B, Chen G: Modifying the DPClus algorithm for identifying protein complexes based on new topological structures. BMC Bioinformatics 9, 398 (2008).

[6] Yong CH, Liu G, Chua HN, Wong L: Supervised maximum-likelihood weighting of composite protein networks for complex prediction. BMC Syst Biol 6(Suppl 2), 13 (2012).

[7] Liu G, Yong CH, Chua HN, Wong L.: Decomposing PPI networks for complex discovery. Proteome Sci 9(Suppl 1), 15 (2011).

[8] Yong CH, Maruyama O, Wong L: Discovery of small protein complexes from PPI networks with size-specific supervised weighting. BMC Syst Biol 8(Suppl 5), 3 (2014).
