# Scoary
Pan-genome wide association studies

Scoary is designed to take the gene_presence_absence.csv file from [Roary] (https://sanger-pathogens.github.io/Roary/) as well as a traits file created by the user and calculate the assocations between all genes in the accessory genome (all genes that are present in *i* genomes where 1 <= *i* < N) and the traits. Scoary reports a list of genes sorted by strength of association per trait.


## Contents
- [What's new] (#whats-new)
- [Dependencies] (#dependencies)
- [Installation] (#installation)
- [Usage] (#usage)
- [Input] (#input)
- [Outout] (#output)
- [Options] (#options)
- [License] (#license)
- [Etymology] (#etymology)
- [Citation] (#citation)
- [Contact] (#contact)

## What's new?

Current release - v1.0 (5th Feb 2016)


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

The gene_presence_absence.csv file will look something like this:
![gene_presence_absence.csv output] (http://sanger-pathogens.github.io/Roary/images/gene_presence_and_absence.png)

_Source: http://sanger-pathogens.github.io/Roary/_

Make sure you know the delimiter in the file. (By default this is ','). Scoary needs to know.

The traits file needs to be formatted in a specific way. 
- It must use the same delimiter as the gene_presence_absence.csv file
- The names of your isolates need to be identical in the two files
- The rows should correspond to your isolates, the columns to the different traits
- Traits needs to be dichotomized. Use "0" to indicate absence and "1" to indicate presence of the trait
- All isolates and traits should be uniquely named and not contain any weird characters
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
| Holm-Sidak_p | A p-value adjusted with Holm-Sidak's method for multiple comparisons correction (An [FWER] (https://en.wikipedia.org/wiki/Familywise_error_rate) type correction) |
| Benjamini_H_p | A p-value adjusted with Benjamini-Hochberg's method for multiple comparisons correction (An [FDR] (https://en.wikipedia.org/wiki/False_discovery_rate) type correction) |



## Options
Scoary can take a number of optional arguments to tweak the output and make sure it performs as intended:
```
usage: SCOARY.py [-h] -t TRAITS -g GENES [-p P_VALUE_CUTOFF]
                 [-c {Individual,Bonferroni,Holm-Sidak,Benjamini-Hochberg}]
                 [-m MAX_HITS] [-r RESTRICT_TO] [-s START_COL]
                 [--delimiter DELIMITER]

Screen pan-genome for trait-associated genes

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
  -c {Individual,Bonferroni,Holm-Sidak,Benjamini-Hochberg}, --correction {Individual,Bonferroni,Holm-Sidak,Benjamini-Hochberg}
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

```
#### The -r parameter
The **-r** parameter is particularly useful, as you can use it to restrict your analysis to a subset of your isolates without altering the gene_presence_absence or trait files. Simply provide a single-line csv file (delimited by ",") with the names of the isolates you would like to include in the current analysis.

This can be useful for example if you have multiple clades in your dataset but would like to restrict analysis to just one clade. Maybe the trait determinant is not the same in the two clades? Or maybe you have data with multiple clones of the same isolate? This will invalidate the statistical analysis, which assumes phylogenetic independency.

The provided file can look something like this:

```
Strain1,Strain2,Strain4,Strain9
```

This will restrict the current analysis to isolates 1,2,4 and 9, and will omit all others.

#### The -s parameter
The **-s** parameter is used to indicate to Scoary which column in the gene_presence_absence.csv file is the _first_ column representing an isolate. By default it is set to 15 (1-based indexing).

## License
Scoary is freely available under a GPLv3 license.

## Etymology
Scoary is an anagram of "scoring" and "Roary", the pan-genome pipeline. It was named as an homage to Roary.

## Citation
Manuscript not yet published

## Contact
Ola Brynildsrud (ola.brynildsrud@fhi.no)
