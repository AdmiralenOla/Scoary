![Scoary - Microbial Pan-GWAS](https://cloud.githubusercontent.com/assets/14874487/15772489/b026105a-2971-11e6-9d1e-da4035502869.png)

Scoary is designed to take the gene_presence_absence.csv file from [Roary] (https://sanger-pathogens.github.io/Roary/) as well as a traits file created by the user and calculate the assocations between all genes in the accessory genome and the traits. It reports a list of genes sorted by strength of association per trait.


## Contents
- [What's new] (#whats-new)
- [Dependencies] (#dependencies)
- [Installation] (#installation)
- [Usage] (#usage)
- [Input] (#input)
- [Output] (#output)
- [Options] (#options)
- [Population structure] (#population-structure)
- [Example data] (#example-data)
- [License] (#license)
- [Etymology] (#etymology)
- [Bugs] (#bugs)
- [FAQ] (#faq)
- [Coming soon] (#coming-soon)
- [Acknowledgements] (#acknowledgements)
- [Feedback] (#feedback)
- [Citation] (#citation)
- [Contact] (#contact)

## What's new?
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


## Dependencies

- Python (Tested with versions 2.7 and 3.5)
- [SciPy] (http://www.scipy.org/install.html) (Tested with versions 0.16, 0.17, 0.18)


## Installation

Scoary is a standalone python script and does not require any installation. Simply download and extract the zip archive or clone the git repository:

    git clone https://github.com/AdmiralenOla/Scoary

and you're ready to use Scoary.
If you want to add it to your $PATH variable:

    export PATH="/Path/to/Scoary:$PATH"


## Usage

    scoary.py -g gene_presence_absence.csv -t traits.csv

## Input
Scoary requires two input files: The gene_presence_absence.csv file from [Roary] (https://sanger-pathogens.github.io/Roary/) and a list of traits to test associations to. 

The **gene_presence_absence.csv** file will look something like this:
![gene_presence_absence.csv output] (http://sanger-pathogens.github.io/Roary/images/gene_presence_and_absence.png)

_Source: http://sanger-pathogens.github.io/Roary/_

Make sure you know the delimiter in the file. (By default this is ','). Scoary needs to know.

The **traits.csv** file needs to be formatted in a specific way. 
- It must use the same delimiter as the gene_presence_absence.csv file
- The names of your isolates need to be identical in the two files
- The rows should correspond to your isolates, the columns to the different traits
- Traits needs to be dichotomized. Use "0" to indicate absence and "1" to indicate presence of the trait
- All isolates and traits should be uniquely named and not contain any weird characters (e.g. %;,/&[]@? etc)
- The top left cell should be left blank

It should look something like this:

|         | Trait1 | Trait2 | ... | TraitM |
| ------- | ------ | ------ | --- | ------ |
| Strain1 | 1      | 0      | ... | 1      |
| Strain2 | 1      | 1      | ... | 0      |
| Strain3 | 0      | 0      | ... | 1      |
| ...     | ...    | ...    | ... | ...    |
| StrainN | 1      | 0      | ... | 0      |

You can see an example of how the input files could look in the exampledata folder.

## Output
Scory outputs a single csv file per trait in the traits file. It uses comma "," as a delimiter. The results consists of genes that were found to be associated with the trait, sorted according to significance. (By default, Scoary reports all genes with a naive p-value < 0.05, but the user can change the cut-off value and use adjusted p-values instead)

The results file contains the following columns:

| Column name | Explanation |
| ----------- | ----------- |
| Gene | The gene name |
| Non-unique gene name | The non-unique gene name |
| Annotation | Annotation |
| Number_pos_present_in | The number of trait-positive isolates this gene was found in |
| Number_neg_present_in | The number of trait-negative isolates this gene was found in |
| Number_pos_not_present_in | The number of trait-positive isolates this gene was not found in |
| Number_neg_not_present_in | The number of trait-negative isolates this gene was not found in |
| Sensitivity | The sensitivity if using the presence of this gene as a diagnostic test to determine trait-positivity |
| Specificity | The specificity if using the non-presence of this gene as a diagnostic test to determine trait-negativity |
| Odds_ratio | [Odds ratio] (https://en.wikipedia.org/wiki/Odds_ratio) |
| p_value | The naïve p-value for the null hypothesis that the presence/absence of this gene is unrelated to the trait status |
| Bonferroni_p | A p-value adjusted with Bonferroni's method for multiple comparisons correction (An [FWER] (https://en.wikipedia.org/wiki/Familywise_error_rate) type correction) |
| Benjamini_H_p | A p-value adjusted with Benjamini-Hochberg's method for multiple comparisons correction (An [FDR] (https://en.wikipedia.org/wiki/False_discovery_rate) type correction) |
| Max_pairwise_comparisons | The maximum number of pairs that contrast in both gene and trait characters that can be drawn on the phylogenetic tree without intersecting lines (Read & Nee, 1995; Maddison, 2000) |
| Max_supporting_pairs | The maximum number of these pairs (Max_pairwise_comparisons) that support A->B or A->b, depending on the odds ratio. |
| Max_opposing_pairs | The maximum number of these pairs (Max_pairwise_comparisons) that oppose A->B or A->b, depending on the odds ratio. |
| Best_pairwise_comp_p | The p-value corresponding to the highest possible number of supporting pairs and the lowest possible number of opposing pairs, e.g. the lowest p-value you could end up with when picking a set of maximum number of pairs. |
| Worst_pairwise_comp_p | The p-value corresponding to the lowest possible number of supporting pairs and the highest possible number of opposing pairs, e.g. the highest p-value you could end up with when picking a set of maximum number of pairs. |



## Options
Scoary can take a number of optional arguments to tweak the output and make sure it performs as intended:
```
usage: scoary.py [-h] [-t TRAITS] [-g GENES]
                 [-p P_VALUE_CUTOFF [P_VALUE_CUTOFF ...]]
                 [-c [{I,B,BH,PW,EPW} [{I,B,BH,PW,EPW} ...]]] [-m MAX_HITS]
                 [-r RESTRICT_TO] [-w] [-s START_COL] [-u]
                 [--delimiter DELIMITER] [--no-time] [--test] [--citation]
                 [--version]

Scoary version 1.4.0 - Screen pan-genome for trait-associated genes

optional arguments:
  -h, --help            show this help message and exit
  -t TRAITS, --traits TRAITS
                        Input trait table (comma-separated-values). Trait
                        presence is indicated by 1, trait absence by 0.
                        Assumes strain names in the first column and trait
                        names in the first row
  -g GENES, --genes GENES
                        Input gene presence/absence table (comma-separated-
                        values) from ROARY. Strain names must be equal to
                        those in the trait table
  -p P_VALUE_CUTOFF [P_VALUE_CUTOFF ...], --p_value_cutoff P_VALUE_CUTOFF [P_VALUE_CUTOFF ...]
                        P-value cut-off. SCOARY will not report genes with
                        higher p-values than this. Set to 1.0 to report all
                        genes. Accepts standard form (e.g. 1E-8). Provide a
                        single value or exactly as many values as correction
                        criteria and in corresponding order. (See example
                        under correction). Default = 0.05
  -c [{I,B,BH,PW,EPW} [{I,B,BH,PW,EPW} ...]], --correction [{I,B,BH,PW,EPW} [{I,B,BH,PW,EPW} ...]]
                        Use the indicated p-value for cut-off. I=Individual
                        (naive) p-value. B=Bonferroni adjusted p-value. BH
                        =Benjamini-Hochberg adjusted p. PW=Best (lowest)
                        pairwise comparison p. EPW=Entire range of pairwise
                        comparison p-values. You can enter as many correction
                        criteria as you would like. These will be associated
                        with the p_value_cutoffs you enter. For example "-c
                        Individual PWB -p 0.1 0.05" will apply a naive p-value
                        cutoff of 0.1 AND additionally require that the entire
                        range of pairwise comparison values are below 0.05 for
                        this gene. Default = Individual p-value. (I)
  -m MAX_HITS, --max_hits MAX_HITS
                        Maximum number of hits to report. SCOARY will only
                        report the top max_hits results per trait
  -r RESTRICT_TO, --restrict_to RESTRICT_TO
                        Use if you only want to analyze a subset of your
                        strains. Scoary will read the provided comma-separated
                        table of strains and restrict analyzes to these.
  -w, --write_reduced   Use with -r if you want Scoary to create a new gene
                        presence absence file from your reduced set of
                        isolates. Note: Columns 1-14 (No. sequences, Avg group
                        size nuc etc) in this file do not reflect the reduced
                        dataset. These are taken from the full dataset.
  -s START_COL, --start_col START_COL
                        On which column in the gene presence/absence file do
                        individual strain info start. Default=15. (1-based
                        indexing)
  -u, --upgma_tree      This flag will cause Scoary to write the calculated
                        UPGMA tree to a newick file
  --delimiter DELIMITER
                        The delimiter between cells in the gene
                        presence/absence and trait files, as well as the
                        output file.
  --no-time             Output file in the form TRAIT.results.csv, instead of
                        TRAIT_TIMESTAMP.csv. When used with the -w argument
                        will output a reduced gene matrix in the form
                        gene_presence_absence_reduced.csv rather than
                        gene_presence_absence_reduced_TIMESTAMP.csv
  --test                Run Scoary on the test set in exampledata, overriding
                        all other parameters.
  --citation            Show citation information, and exit.
  --version             Display Scoary version, and exit.

by Ola Brynildsrud (olbb@fhi.no)
```
#### The -r parameter
The **-r** parameter is particularly useful, as you can use it to restrict your analysis to a subset of your isolates without altering the gene_presence_absence or trait files. Simply provide a single-line csv file (delimited by ",") with the names of the isolates you would like to include in the current analysis.

This can be useful for example if you have multiple clades in your dataset but would like to restrict analysis to just one clade. Maybe the trait determinant is not the same in the two clades? Or perhaps you have missing data for some isolates?

The provided file can look something like this:

```
Strain1,Strain2,Strain4,Strain9
```

This will restrict the current analysis to isolates 1,2,4 and 9, and will omit all others.

#### The -w flag
Using the **-w** flag with **-r** will make Scoary write a reduced gene presence/absence file containing only those isolates specified with **-r**. This makes the program run much faster if you are analyzing a small subset of a large dataset.

#### The -s parameter
The **-s** parameter is used to indicate to Scoary which column in the gene_presence_absence.csv file is the _first_ column representing an isolate. By default it is set to 15 (1-based indexing).

#### The -p, -m and -c parameters
These parameters control your output. **-m** sets a hard cut-off on the number of hits reported. With **-p** you can set that no gene with a higher p-value will be reported. (Tip: Set this to 1.0 to report every single gene). You can mix these parameters with **-c**. If you only wanted genes with a Bonferroni-adjusted p-value < 1E-10 you could use _-p 1E-10 -c B_.

##### Combining filtering options
From version 1.4.0, you can also mix different restrictions together. For example, you may want to specify that the entire range of pairwise comparisons p-values be < 1E-5, but you still doubt some of your results. You could try to filter your results more strictly by also requiring an Individual (naive) p-value of less than 0.01. You would then use _-c EPW I -p 1E-5 0.01_. You need to enter the -c options and the -p options in the corresponding order. 

Alternatively, you can specify a single (one) p-value, and this will be taken as the filter for all the specified -c options. For example _-c EPW BH -p 0.05_ will filter the results to only include genes where the entire range of pairwise comparison as well as the Benjamini-Hochberg p-values are > 0.05

#### The -u flag
Calling Scoary with the **-u** flag will cause it to write a newick file of the UPGMA tree that is calculated internally. The tree is based on pairwise Hamming distances in the gene_presence_absence matrix.

## Population structure
Scoary implements the pairwise comparisons algorithm (Read & Nee, 1995; Maddison, 2000) to identify the maximum number of non-intersecting pairs of isolates that contrast in the state of both gene and trait. It does this by creating an UPGMA tree from the information contained in the gene_presence_absence matrix, annotating tips with gene and trait status, and recursively traversing the tree for each gene that were significant in the initial analysis. (i.e. those with p<0.05 if settings are left at default.)

This tells you something about the **number of times** the gene and trait co-emerged in the evolutionary history of your sample. Consider the two following trees. In both scenarios, the gene perfectly predict the trait status, with 10 positive and 11 negative isolates, corresponding to a naïve p-value of 2.8E-6. However, in the first tree, there is a maximum of two non-intersecting contrasting pairs, which must be considered relatively weak evidence for a causal link between this gene and the trait. There can be many other evolutionary events at the stars that explain the observed distribution equally well as this gene. In the second tree, however, there is a maximum of seven possible non-intersecting contrasting pairs, which implies that this gene and the trait co-emerged seven times. This would be considered far stronger evidence for a causal link between the gene and the trait.

![A not-so-significant link between gene and trait](https://cloud.githubusercontent.com/assets/14874487/15569716/8c66322a-2332-11e6-8500-5d27828417c7.png)

![A very significant link between gene and trait](https://cloud.githubusercontent.com/assets/14874487/15569715/8c61f962-2332-11e6-90e7-5c37976071c8.png)

One must also consider that there might be multiple ways of picking the maximum number of contrasting pairs, and of all these possible sets of pairings, some might provide more support for A->B than others. Consider the following tree:

![A best possible pairing](https://cloud.githubusercontent.com/assets/14874487/15708517/271bd3b2-27ff-11e6-9190-8c655622bbfd.png)

The above tree has a maximum of 6 contrasting pairs, and in this tree the pairs have been chosen so that all pairs support A->B. (The presence of the gene caused the presence of the phenotype). However, in this particular tree we could also have picked 6 contrasting pairs where not all pairs supports this. See for example this pairing:

![A worst possible pairing](https://cloud.githubusercontent.com/assets/14874487/15708516/271bf496-27ff-11e6-81d4-7309f2c274cc.png)

The above tree has the same topology and terminal states, and the same number of contrasting pairs, but now we have chosen pairs so that 5 pairs support A->B while 1 pair oppose it (It suggests that A->b / a->B). This is a worst possible pairing which maintains the maximum number of contrasting pairs.

Scoary reports the best (lowest) and worst (highest) p-values, corresponding to the first and the second scenario, respectively. The p-value corresponds to a binomial test using the number of supporting pairs as successes and p=0.5 for each state. A p<0.05 would thus typically be considered as a rejection of the null hypothesis that the expressed phenotype is not associated with the gene.

## Example data
In the exampledata folder you can see an example of how the input files would typically look. This simulated and completely fictitious data set consists of 100 isolates with around 3000 core genes and a total pan-genome of around 9000 genes. 

Here I have used tetracycline resistance as the phenotype for which we would like to know the genetic basis. In this example, only a single gene controls the expression of tetracycline resistance, although the penetrance of the gene is not 100% (i.e. isolates can have the gene and still be susceptible towards tetracycline), and other, unmeasured factors (for example point mutations) can also induce resistance. 

This particular example very clearly identifies the causal gene (looking at the pairwise comparison p-values), whereas in real experiments the results are sometimes a lot messier.

Running Scoary with the --test flag is equivalent to the following command:
```
python ./scoary.py -t ./exampledata/Tetracycline_resistance.csv -g ./exampledata/Gene_presence_absence.csv -u
```

## License
Scoary is freely available under a GPLv3 license.

## Etymology
Scoary is an anagram of "scoring" and "Roary", the pan-genome pipeline. It was named as an homage to Roary.

## Bugs
- Please report bugs here (Issues) or to me directly at olbb@fhi.no

## FAQ
- **How can you justify p=0.5 in your pairwise comparisons method? Is this species-specific?**

The reasoning is as follows: Scoary first finds the maximum number of independent contrasting pairs in a phylogenetic tree, irrespective of gene-trait status. Thus, AB-ab pairs should be equally likely as Ab-aB pairs if your null hypothesis is true. Your null hypothesis in this case, is that there is no detectable association between A/a and B/b. If AB-ab pairs are much more common than Ab-aB pairs then you can be confident that the true p was not 0.5. And if this is the case then then there seems to be an association between your A/a (your gene) and your B/b (your phenotype). A justification for this way of testing can be found in Read and Nee, 1995.
- **Why is my "Best_pairwise_comp_p" higher than my "Worst_pairwise_comp_p"?**

The "best" and "worst" labels are attached to the odds ratio of the gene in the non-population structure-corrected analysis. For example, you may find an odds ratio of 2.0 for a particular gene, meaning presence of the gene was tied to presence of the phenotype. But when you inspect your pairwise comparisons p-values you see that the "best" p-value was 0.2 and the "worst" was 1.0E-5. This means that in your phylogenetic tree, an enrichment of Ab-aB pairs was more common. In other words, the presence of this gene actually seems associated to a _silencing_ of the phenotype, in spite of your original odds ratio. Note that the odds ratio can be inflated for example by sampling of very closely related isolates.

## Coming soon
Please feel free to suggest improvements, point out bugs or methods that could be better optimized.

## Acknowledgements
- Marco Galardini cleaned my code and made many nifty improvements.
- The QuadTree and UPGMA implementation was heavily based on code by Christian Storm Pedersen.
- Inês Mendes pointed out a number of bugs related adjusted p-values and isolate restriction.
- Eric Deveaud added versioning.

## Feedback
I greatly appreciate any feedback, even negative. If you like (or dislike) Scoary, please feel free to tell your friends and colleagues about it. If you don't have friends or colleagues, please feel free to rant about it on your blog or social media profile.

## Citation
Manuscript not yet published.

## Contact
Ola Brynildsrud (ola.brynildsrud@fhi.no)
