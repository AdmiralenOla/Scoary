# Scoary
Pan-genome wide association studies

Scoary is designed to take the gene_presence_absence.csv file from [Roary] (https://sanger-pathogens.github.io/Roary/) as well as a traits file created by the user and calculate the assocations between all genes in the accessory genome (all genes that are present in *i* genomes where 1 <= *i* < N) and the traits. Scoary reports a list of genes sorted by strength of association per trait.


## Contents
- [What's new] (#whats-new)
- [Dependencies] (#dependencies)
- [Installation] (#installation)
- [Usage] (#usage)

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

## Options

## Etymology
Scoary is an anagram of "scoring" and "Roary", the pan-genome pipeline. It was named as an homage to Roary.

## Citation
Manuscript not yet published

## Contact
Ola Brynildsrud (ola.brynildsrud@fhi.no)
