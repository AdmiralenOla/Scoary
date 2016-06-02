![Scoary - Microbial Pan-GWAS](https://cloud.githubusercontent.com/assets/14874487/15738049/ad3ffdf6-28a9-11e6-854a-fbceaf2c2328.png)

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
- [License] (#license)
- [Etymology] (#etymology)
- [Coming soon] (#coming-soon)
- [Acknowledgements] (#acknowledgements)
- [Citation] (#citation)
- [Contact] (#contact)

## What's new?

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

- Python (2.7.x)
- [SciPy] (http://www.scipy.org/install.html)


## Installation

Scoary is a standalone python script and does not require any installation. Simply clone the git repository:

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

## Output
Scory outputs a single csv file per trait in the traits file. It uses semicolon ";" as a delimiter to avoid conflict with gene annotations that include commas. The results consists of genes that were found to be associated with the trait, sorted according to significance. (By default, Scoary reports all genes with a p-value < 0.05, but the user can change the cut-off value and use adjusted p-values instead)

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
usage: Scoary.py [-h] [-t TRAITS] [-g GENES] [-p P_VALUE_CUTOFF]
                 [-c {Individual,Bonferroni,Benjamini-Hochberg}] [-m MAX_HITS]
                 [-r RESTRICT_TO] [-u] [-s START_COL] [--delimiter DELIMITER]
                 [--version]

Scoary version v1.3.0 - Screen pan-genome for trait-associated genes

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
  -p P_VALUE_CUTOFF, --p_value_cutoff P_VALUE_CUTOFF
                        P-value cut-off. SCOARY will not report genes with
                        higher p-values than this. Set to 1.0 to report all
                        genes. Default = 0.05
  -c {Individual,Bonferroni,Benjamini-Hochberg}, --correction {Individual,Bonferroni,Benjamini-Hochberg}
                        Instead of cutting off at the individual test p-value
                        (option -p), use the indicated corrected p-value for
                        cut-off. Default = use individual test p-value.
  -m MAX_HITS, --max_hits MAX_HITS
                        Maximum number of hits to report. SCOARY will only
                        report the top max_hits results per trait
  -r RESTRICT_TO, --restrict_to RESTRICT_TO
                        Use if you only want to analyze a subset of your
                        strains. SCOARY will read the provided comma-separated
                        table of strains and restrict analyzes to these.
  -u, --upgma_tree      This flag will cause Scoary to write the calculated
                        UPGMA tree to a newick file
  -s START_COL, --start_col START_COL
                        On which column in the gene presence/absence file do
                        individual strain info start. Default=15. (1-based
                        indexing)
  --delimiter DELIMITER
                        The delimiter between cells in the gene
                        presence/absence and trait files. NOTE: Even though
                        commas are the default they might mess with the
                        annotation column, and it is therefore recommended to
                        save your files using semicolon or tab (" ") instead.
                        SCOARY will output files delimited by semicolon
  --version             Display Scoary version, and exit.

```
#### The -r parameter
The **-r** parameter is particularly useful, as you can use it to restrict your analysis to a subset of your isolates without altering the gene_presence_absence or trait files. Simply provide a single-line csv file (delimited by ",") with the names of the isolates you would like to include in the current analysis.

This can be useful for example if you have multiple clades in your dataset but would like to restrict analysis to just one clade. Maybe the trait determinant is not the same in the two clades? Or perhaps you have missing data for some isolates?

The provided file can look something like this:

```
Strain1,Strain2,Strain4,Strain9
```

This will restrict the current analysis to isolates 1,2,4 and 9, and will omit all others.

#### The -s parameter
The **-s** parameter is used to indicate to Scoary which column in the gene_presence_absence.csv file is the _first_ column representing an isolate. By default it is set to 15 (1-based indexing).

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

## License
Scoary is freely available under a GPLv3 license.

## Etymology
Scoary is an anagram of "scoring" and "Roary", the pan-genome pipeline. It was named as an homage to Roary.

## Coming soon
Please feel free to suggest improvements, point out bugs or methods that could be better optimized.

## Acknowledgements
- The QuadTree and UPGMA implementation was heavily based on code by Christian Storm Pedersen
- Inês Mendes pointed out a number of bugs related adjusted p-values and isolate restriction.
- Eric Deveaud added versioning.

## Citation
Manuscript not yet published.

## Contact
Ola Brynildsrud (ola.brynildsrud@fhi.no)
