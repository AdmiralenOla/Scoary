# CHANGELOG
v1.6.16 (Aug 2017)
- Bug fixes to vcf2scoary

v1.6.15 (Jul 2017)
- Apparently 1.6.14 did not fix the pip issue for all users. Deleted the function that seemed to cause the pip crash

v1.6.14 (Jul 2017)
- Fixed a bug where Scoary could not be upgraded using pip

v1.6.13 (Jul 2017)
- Fixed a bug where Scoary would not handle converted VCF files due to different headers
- Changed the parameters of --test
- Fixed the GUI that was broken since 1.6.12 due to introduction of --include_input_columns
- Include_input_columns now correctly handles the ALL keyword
- Added the ExampleVCFTrait.csv file to test VCF functionality

v1.6.12 (Jun 2017)
- Convert VCF files to Roary/Scoary format, allowing analysis on a wide range of variants (SNPs, indels, structural variations etc)
- Grab columns from the Roary input and put in the output (To get strain-specific protein names, for example)
- Scoary now comes with a manual, located under docs/tex/scoary_manual.pdf
- The log now includes the original command line

v1.6.11 (Apr 2017)
- Blank values in trait files will now correctly be read as missing. Fixes (#54)
- Added --no_pairwise option for simple set differences / categorical enrichment analysis without causal hypothesis (As requested among others in (#53)
- Modified GUI with no_pairwise and slightly modified look.
- Added ExampleTree.nwk to exampledata
- Added support for travis. (CI tests will be further developed)
- Added example cases in README.md
- Fixed broken links in README.md

v1.6.10 (Jan 2017)
- Scoary now creates a log file (both in terminal and GUI mode)
- Fixed a bug where empirical p-values would exceed 1.0
- Fixed a bug where Scoary would crash when pruning many isolates from the internally calculated phylogenetic tree

v1.6.9 (Dec 2016)
- Scoary now handles missing data specified in the traits file as "NA", "-" or ".".
- Now also handles missing isolates (rows) in the traits file.
- Stylistic changes, including line lengths formats to conform with PEP8.

v1.6.8 (29th Nov 2016)
- Bugfix: Addition of numpy.float64 and long were causing problems on some systems. Now converted to float beforehand. (Fixes #43)

v1.6.7 (28th Nov 2016)
- Scoary article is now out! Updated Readme and citation.

v1.6.6 (9th Nov 2016)
- Changed citation.
- Scoary is now available through Zenodo.

v1.6.5 (25th Oct 2016)
- Added the README_pypi.md file so that setup.py can be used to build/install locally if pip not installed.

v1.6.4 (18th Oct 2016)
- Made the permutations algorithm faster by implementing a dynamic cutoff

v1.6.3 (14th Oct 2016)
- Removed OR calculation from PhyloTree class
- Now prints collapse option to stdout
- Cleaned up threading method

v1.6.2 (13th Oct 2016)
- Collapsing correlated genes is now optional, controlled with the --collapse flag
- Fixed a bug where the collapsing duplicated some genes
- Collapsing is default to False in GUI mode

v1.6.1 (7th Oct 2016)
- Added a GUI, callable as scoary_GUI
- Added post-analysis label-switching permutations, and the option to filter on empirical p-values
- Added multiprocessing to speed up permutations
- Added a second column in the phenotype exampledata
- Adds an outdir option
- Now available through PyPi

v1.5.1 (22nd Sep 2016)
- Users can now provide a custom tree that will be used for analysis instead of the internally calculated one.
- Changed the way trees are written to file. Leaf nodes will no longer have "double" distance. The tree files are now also rooted, which allows them to be read by other software
- Fixed some typos and added a lot of exception throwing for when the user is trying strange things, such as pointing to files that don't exist.

v1.5.0 (13th Sep 2016)
- Scoary is now installable via pip! (Thanks go to Anders Goncalves da Silva). The scoary.py script will now be deprecated, but is still available for legacy use. See [Installation] (#installation)
- The program now also prints out the filtration options being used for the current analysis.

v1.4.2 (13th Sep 2016)
- Fixed a bug that would cause Scoary to crash if ran without any -c options.

v1.4.1 (8th Sep 2016)
- Fixed a bug where results in the output file did not have quotes around them. This could interfere with annotations that have delimiters (like commas) in them.

v1.4.0 (7th Sep 2016)
- The correction options are now a lot more sophisticated, and allow multiple restrictions with individual p-values to be set.
- Fixed a bug that would sometimes result in a too strict Bonferroni p-value in the results.
- No-time now also applies to the reduced gene presence/absence file.
- The results sorting now works a little differently. When using only pairwise comparisons filtering, results will be sorted by these p-values. When other filters (Individual, Bonferroni, Benjamini-Hochberg) are used, results will be sorted by these instead.
- Fixed a bug where the program would crash if you specified a very large max_hits.
- The citation option has been added, and includes a nifty ASCII logo.
- Some optimization when using pairwise comparison to filter results.

v1.3.7 (2nd Sep 2016)
- The no-time argument can now be used to avoid output files (results and tree file) to come with a timestamp in the name. Should make it easier to implement Scoary in automated pipelines. (Credits: Marco Galardini)
- Comma is now the default delimiter in input and output files. The user can specify another input/output delimiter with the delimiter argument. (Note that the two input files and the output files will all have the same delimiter)
- Bug fixes that caused python3 problems in 1.3.6 (Credits: Marco Galardini)

v1.3.6 (28th Jul 2016)
- Simulated example data is included in the exampledata folder. This is primarily intended as a guide to how input files can look, as well as giving users a quick view at what the program can do. Running Scoary with the --test flag will overrun all other options (except --version) and automatically run the exampledata with default options.

v1.3.5 (PRE-RELEASE) (5th Jul 2016)
- You can now use the -w option with -r to write a reduced gene presence/absence file containing only the subset isolates. This ensures that the program will run much faster if you have a large dataset (1000s of isolates) but only want to analyze a subset. Scoary automatically opens and analyzes the newly written file.
- This is a pre-release version. There might still be bugs in the code, in which case I would be grateful if you report them.

v1.3.4 (16th Jun 2016)
- Scoary no longer crashes when using Scipy 0.16 instead of 0.17.
- More information about what's going on is printed. (Useful for very large datasets that take long to analyze)

v1.3.3 (9th Jun 2016)
- BUG FIX: Tree calculation had been broken since 1.3.2 yesterday. Sorry about that.

v1.3.2 (8th Jun 2016)
- Scoary SHOULD now be compatible with Python 3.x. Please report back if you disagree.

v1.3.1 (8th Jun 2016)
- Scoary is now packaged as a proper python package according to the PEP8 style guide (Thanks to Marco Galardini)
- A percentage counter is now shown when calculating pairs

v1.3.0 (1st Jun 2016)
- Major changes to the pairwise comparisons algorithm. Scoary now calculates the maximum number of contrasting pairs, and given that maximum number tests the maximum number of pairs that SUPPORT A -> B (AB-ab pairs) and the maximum number of pairs that OPPOSE A -> B (Ab-aB pairs). The opposite is true for genes where the odds ratio is < 1 (i.e. that indicate A -> b).
- The p-values reported from the pairwise comparisons is now a range. It reports the best (lowest) p-value, which comes from the maximum number of supporting pairs and the minimum number of opposing (given a set total), as well as the worst (highest) p-value, which comes from the minimum number of supporting pairs and the maximum number of non-supporting, given a set total number of pairings. It does this at each node in the tree.
- Scoary can now print the UPGMA tree that is calculated internally from the Hamming distances of the gene_presence_absence matrix. Do this by using the -u flag.
- The elapsed time will now print when finished.

v1.2.3 (30th May 2016)
- Odds ratios should now be correct again. These were behaving strangely since 1.2.0. Apologies.

v1.2.2 (28th May 2016)
- Another bug fix related to the restrict_to option where Scoary would crash if it was NOT set. 

v1.2.1 (26th May 2016)

- Bug fix. The --restrict_to option had become broken in 1.2.0 because one function was passing the full set of isolate names rather than the restricted set downstream in the analyses. This has been fixed.

v1.2.0 (23rd May 2016)

- Major changes. Scoary now implements the pairwise comparisons algorithm to account for population structure.
- Scoary now imports 4 classes: Matrix, QuadTree, PhyloTree and Tip. They two former are used for storing pairwise distances between isolates, and the two latter are used for the pairwise comparisons algorithm. Scoary_methods contain some new functions and most of the old ones have been altered to allow the pairwise comparisons implementation. However, there are no changes to how Scoary calculates the previously implemented statistics.
- There should now be a significant speed increase, as Scoary now stores Fisher's exact test p-values rather than re-calculating for every gene. Previously, a new calculation was made for every gene, even when the 2x2 table was identical.

v1.1.2 (4th May 2016)

- Fixed another bug related to Benjamini-Hochberg p-value adjustments. (Thanks again to cimendes). The numbers should now correspond to R's "p.adjust" method IF the number of tests are the same. (Scoary runs the correction on ALL genes, not just those with p<.05)
- Results are now printed with "" quotation marks around each cell to avoid weird cell breaks if annotations contain semicolons when opening in spreadsheets.

v1.1.1 (3rd May 2016)

- Fixed a bug where adjusted p-values were not being calculated and output correctly. (Thanks to cimendes).
- Benjamini-Hochberg p-values are now calculated from least to most significant (step-up) instead of step-down.
- Holm-Sidak p-values has been removed.
- The Scoary version is now displayed in the help message and if calling Scoary with the --version flag. (Thanks to EricDeveaud)

v1.1 (29th Mar 2016)

- Scoary now imports all methods from Scoary_methods in order to circumvent errors when trying to run Scoary under Python 3.x.
- Genes that have the same naïve p-value now have adjusted p-values penalized with the same factor rather than an increasing one.
