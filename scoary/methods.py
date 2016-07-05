#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Script for associating Roary output with phenotypic trait (+/-)
# Idea and implementation: Ola Brynildsrud (olbb@fhi.no)

# Null hypothesis: Gene is equally distributed in Trait+ and Trait- isolates

import argparse
import sys
import csv
import time
from scipy import stats as ss
from scipy import spatial
from .classes import Matrix
from .classes import QuadTree
from .classes import PhyloTree
import scoary

SCOARY_VERSION = scoary.__version__

# Python 2/3 annoyances
try:
    xrange
except NameError:
    xrange = range

def main():
    """
    The main function of Scoary.
    """
    # Parse arguments.
    parser = argparse.ArgumentParser(
            description='Scoary version %s - Screen pan-genome for '
                        'trait-associated genes' % SCOARY_VERSION,
                        epilog='by Ola Brynildsrud (olbb@fhi.no)')
    parser.add_argument('-t', '--traits',
                        required=True,
                        help='Input trait table (comma-separated-values). '
                        'Trait presence is indicated by 1, trait absence by 0. '
                        'Assumes strain names in the first column and trait '
                        'names in the first row')
    parser.add_argument('-g', '--genes',
                        required=True,
                        help='Input gene presence/absence table '
                        '(comma-separated-values) from ROARY. '
                        'Strain names must be equal to those in the trait '
                        'table')
    parser.add_argument('-p', '--p_value_cutoff',
                        help='P-value cut-off. SCOARY will not report genes '
                        'with higher p-values than this. Set to 1.0 to report '
                        'all genes. Accepts standard form (e.g. 1E-8). '
                        'Default = 0.05',
                        default=0.05,
                        type=float)
    parser.add_argument('-c', '--correction',
                        help='Instead of cutting off at the individual test '
                        'p-value (option -p), use the indicated corrected '
                        'p-value for cut-off. '
                        'Default = use individual test p-value.',
                        choices=['Individual',
                                 'Bonferroni',
                                 'Benjamini-Hochberg'],
                        default='Individual')
    parser.add_argument('-m', '--max_hits',
                        help='Maximum number of hits to report. SCOARY will '
                        'only report the top max_hits results per trait',
                        type=int)
    parser.add_argument('-r', '--restrict_to',
                        help='Use if you only want to analyze a subset of your'
                        ' strains. SCOARY will read the provided '
                        'comma-separated table of strains and restrict '
                        'analyzes to these.')
    parser.add_argument('-w', '--write_reduced',
                        help='Use with -r if you want Scoary to create a new '
                        'gene presence absence file from your reduced set of '
                        'isolates.',
                        default=False,
                        action='store_true')
    parser.add_argument('-s', '--start_col',
                        help='On which column in the gene presence/absence '
                        'file do individual strain info start. Default=15. '
                        '(1-based indexing)',
                        default=15,
                        type=int)
    parser.add_argument('-u', '--upgma_tree',
                        help='This flag will cause Scoary to write the '
                        'calculated UPGMA tree to a newick file',
                        default=False,
                        action='store_true')
    parser.add_argument('--delimiter',
                        help='The delimiter between cells in the gene '
                        'presence/absence and trait files. ',
                        default=',',
                        type=str)
    parser.add_argument('--version', help='Display Scoary version, and exit.',
                        action='version',
                        version=SCOARY_VERSION)

    args = parser.parse_args()
    
    if (args.p_value_cutoff > 1.0) or (args.p_value_cutoff <= 0.0):
        sys.exit("P must be between 0.0 and 1.0 or exactly 1.0")
    if (len(args.delimiter) > 1):
        sys.exit("Delimiter must be a single character string. There is no support for tab.")

    starttime = time.time()
    
    with open(args.genes, "rU") as genes, open(args.traits, "rU") as traits:

        if args.restrict_to is not None:
            allowed_isolates = [isolate
                                for line in
                                open(args.restrict_to,"rU")
                                for isolate in line.rstrip().split(",")]
        else:
            # Despite the confusing name
            # this actually means all isolates are allowed
            # and included in the analysis
            allowed_isolates = None
            if args.write_reduced:
                sys.exit("You cannot use the -w argument without specifying a subset (-r)")
            
        print("Reading gene presence absence file")    
        genedic_and_matrix = Csv_to_dic_Roary(genes,
                                              args.delimiter,
                                              startcol=args.start_col - 1,
                                              allowed_isolates=allowed_isolates,
                                              writereducedset=args.write_reduced)
        genedic = genedic_and_matrix["Roarydic"]
        zeroonesmatrix = genedic_and_matrix["Zero_ones_matrix"]
        strains = genedic_and_matrix["Strains"]
        print("Creating Hamming distance matrix based on gene presence/absence")
        TDM = CreateTriangularDistanceMatrix(zeroonesmatrix, strains)
        QT = PopulateQuadTreeWithDistances(TDM)
        print("Building UPGMA tree from distance matrix")
        upgmatree = upgma(QT)
        print("Reading traits file")
        traitsdic = Csv_to_dic(traits, args.delimiter, allowed_isolates)

        print("Finished loading files into memory.")
        print("Tallying genes and performing statistical analyses")

        RES_and_GTC = Setup_results(genedic, traitsdic)
        RES = RES_and_GTC["Results"]
        GTC = RES_and_GTC["Gene_trait_combinations"]

        if args.upgma_tree:
            StoreUPGMAtreeToFile(upgmatree)

        StoreResults(RES,
                     args.max_hits,
                     args.p_value_cutoff,
                     args.correction, upgmatree, GTC)
        print("\nFinished. Checked a total of %d genes for associations to %d trait(s). "
              "Total time used: %d seconds." % (len(genedic),
                                                len(traitsdic),
                                                int(time.time()-starttime)))

    sys.exit(0)


def CreateTriangularDistanceMatrix(zeroonesmatrix, strainnames):
    """
    Converts a raw matrix of 0s and 1s (indicating gene presence and absence,
    rows correspond to genes and columns to allowed strains) to an (upper)
    triangular matrix of pairwise Hamming distances.
    The distance d(i,i) is set to 1 for all i.
    """
    try:
        hamming_distances = list(spatial.distance.pdist(zeroonesmatrix, 'hamming'))
    except TypeError:
        sys.exit("Could not locate scipy.spatial.distance.pdist. Perhaps you have an old version of SciPy installed?")
    nstrains = int((1 + (1 + 8*len(hamming_distances))**0.5)/2)
    TriangularDistanceMatrix = []
    Strain_names = []
    for i in xrange(nstrains):
        add = [1]  # Adding the maximum relative hamming distance to prevent Quadtree algorithm to pick a pair where i = j
        add += hamming_distances[:(nstrains-i-1)]
        hamming_distances = hamming_distances[(nstrains-i-1):]
        TriangularDistanceMatrix.append(add)
        Strain_names.append(strainnames[i])

    return {"matrix": TriangularDistanceMatrix, "names": Strain_names}


def PopulateQuadTreeWithDistances(TDM):
    """
    Takes a triangular distance matrix, creates a quadtree and populates it
    with the hamming distances between isolates.
    First creates a Quadmatrix, so not really optimized.
    """
    Quadmatrix = Matrix(dim=len(TDM["matrix"]))
    for i in xrange(Quadmatrix.dim):
        for j in xrange(i, Quadmatrix.dim):
            try:
                Quadmatrix[i][j] = Quadmatrix[j][i] = TDM["matrix"][i][(j-i)]
            except IndexError:
                sys.exit("There was an error trying to populate the Quadtree "
                         "with pairwise distances. Please report this bug to olbb@fhi.no")
    PopulatedQuadtree = QuadTree(Quadmatrix.dim, names=TDM["names"])
    for i in xrange(Quadmatrix.dim):
        PopulatedQuadtree.insert_row(i, Quadmatrix[i])
    return PopulatedQuadtree

def ReduceSet(genefile, delimiter, startcol=14, allowed_isolates=None):
    csvfile = csv.reader(genefile, skipinitialspace=True, delimiter=delimiter)
    header = next(csvfile)
    allowed_indexes = range(startcol)
    for c in xrange(len(header)):
        if header[c] in allowed_isolates:
            allowed_indexes.append(c)
    
    print("Writing gene presence absence file for the reduced set of isolates")
    reducedfilename = "gene_presence_absence_reduced_" + time.strftime("_%d_%m_%Y_%H%M") + ".csv"
    with open(reducedfilename, "w") as csvout:
        wtr = csv.writer(csvout, delimiter = delimiter)
        newheader = [header[a] for a in allowed_indexes]
        wtr.writerow(newheader)
        for r in csvfile:
            wtr.writerow( tuple(r[a] for a in allowed_indexes) )
    print("Finished writing reduced gene presence absence list to file " + reducedfilename)
    return reducedfilename

def Csv_to_dic_Roary(genefile, delimiter, startcol=14, allowed_isolates=None, writereducedset=False):
    """
    Converts a gene presence/absence file into dictionaries
    that are readable by Roary
    """
    r = {}
    if writereducedset:
        file = open(ReduceSet(genefile,delimiter,startcol,allowed_isolates),"rU")
        csvfile = csv.reader(file, skipinitialspace=True, delimiter=delimiter)
    else:
        csvfile = csv.reader(genefile, skipinitialspace=True, delimiter=delimiter)
    header = next(csvfile)
            
    roaryfile = True

    strains = header[startcol:]
    strain_names_allowed = [val for val in strains if val in allowed_isolates] if allowed_isolates is not None else strains
    zero_ones_matrix = []

    try:
        genecol = header.index("Gene")
        nugcol = header.index("Non-unique Gene name")
        anncol = header.index("Annotation")
    except ValueError:
        print("Warning: Could not properly detect the correct names for all columns in the ROARY table.")
        genecol = 0
        nugcol = 1
        anncol = 2

    for line in csvfile:
        q = line
        try:
            r[q[genecol]] = {"Non-unique Gene name": q[nugcol], "Annotation": q[anncol]} if roaryfile else {}
        except IndexError:
            sys.exit("ERROR: Could not read gene presence absence file. Verify that this file is a proper Roary file "
            "using the specified delimiter (default is ',').")
        # The zero_ones_line represents the presence (1) or absence (0) of a gene. It is used for calculating distances between strains.
        zero_ones_line = []

        for strain in xrange(len(strains)):
            if (allowed_isolates is not None) and strains[strain] not in allowed_isolates:
                continue
            if q[startcol + strain] in ["", "0", "-"]:
                # If the gene is not present, AND The isolate is allowed
                r[q[genecol]][strains[strain]] = 0
                zero_ones_line.append(0)
                # Add a 0 to indicate non-presence
            else:
                # Gene is present if any other value than "", "0" or "-" is in the cell
                r[q[genecol]][strains[strain]] = 1
                zero_ones_line.append(1)
                # Add a 1 to indicate presence of the current gene in this strain

        # Since we are only interested in the differences between strains, no need to append the zero_ones_line if it is all 1's (core gene) or all 0's (not in collection)
        if 1 in zero_ones_line and 0 in zero_ones_line:
            zero_ones_matrix.append(zero_ones_line)

    # Transpose list for distance calculation purposes
    if writereducedset:
        file.close()
    zero_ones_matrix = list(map(list, zip(*zero_ones_matrix)))
    return {"Roarydic": r,
            "Zero_ones_matrix": zero_ones_matrix,
            "Strains": strain_names_allowed}


def Csv_to_dic(csvfile, delimiter, allowed_isolates):
    """
    Converts an input traits file (csv format) to dictionaries readable by Roary
    """
    tab = list(zip(*csv.reader(csvfile, delimiter=delimiter)))
    r = {}
    if len(tab) < 2:
        sys.exit("Please check that your traits file is formatted properly and contain at least one trait")
    for num_trait in xrange(1, len(tab)):
        p = dict(zip(tab[0], tab[num_trait]))
        if "" in p:
            name_trait = p[""]
            del p[""]
        elif "Name" in p:
            name_trait = p["Name"]
            del p["Name"]
        else:
            sys.exit("Make sure the top-left cell in the traits file is either empty or 'Name'. Do not include empty rows")

        # Filter so that only allowed isolates are included
        if allowed_isolates is not None:
            p = {strain: indicator for (strain, indicator) in list(p.items())
                 if strain in allowed_isolates}
        r[name_trait] = p

    return r


def Setup_results(genedic, traitsdic):
    """
    This is the method that actually does all the counting of genes,
    calculation of p-values and post-test adjustment of p-values.
    The majority of running time is currently spent doing Fisher's exact test.
    The running time is somewhat shortened by storing the p-values of known
    2x2 cell distributions, so that identical Fisher's tests will not have
    to be run.
    """
    # Need to create one results dictionary for each trait

    all_traits = {}
    gene_trait_combinations = {}
    fisher_calculated_values = {}
    for trait in traitsdic:
        all_traits[trait] = {}
        # Also need a list of all p-values to perform stepwise FDR methods
        p_value_list = []
        print("Gene-wise counting and Fisher's exact tests for trait: " + trait)

        # We also need a number of tests variable for each genedic.
        # (The number of tests may not be the same, depending on how many genes are in 0 or all strains.)
        number_of_tests = len(genedic)
        initial_number_of_tests = number_of_tests  # This number needed for status bar

        # Additionally, we need a dictionary to, for each gene
        # hold the gene-trait status (e.g "AB" or "ab" of each strain)
        gene_trait_combinations[trait] = {}

        i = 0  # Progress
        for gene in genedic:
            # Status:
            sys.stdout.write("\r{:.2%}".format(float(i)/initial_number_of_tests))
            sys.stdout.flush()
            i += 1  # Progress

            stats = Perform_statistics(traitsdic[trait], genedic[gene])
            stat_table = stats["statistics"]
            gene_trait_combinations[trait][gene] = stats["gene_trait"]
            num_pos = stat_table["tpgp"] + stat_table["tpgn"]
            num_neg = stat_table["tngp"] + stat_table["tngn"]

            if (stat_table["tpgp"] + stat_table["tngp"]) == 0:
                number_of_tests -= 1  # Remove 1 from the number of tests
                continue  # proceed to the next gene

            if (stat_table["tpgn"] + stat_table["tngn"]) == 0:
                number_of_tests -= 1
                continue

            obs_table = [[stat_table["tpgp"],
                          stat_table["tpgn"]],
                         [stat_table["tngp"],
                          stat_table["tngn"]]]
            obs_tuple = (stat_table["tpgp"],
                         stat_table["tpgn"],
                         stat_table["tngp"],
                         stat_table["tngn"])
            if obs_tuple in fisher_calculated_values:
                odds_ratio = fisher_calculated_values[obs_tuple][0]
                p_value = fisher_calculated_values[obs_tuple][1]
            else:
                fisher = ss.fisher_exact(obs_table)
                p_value = fisher[1]
                odds_ratio = fisher[0]
                fisher_calculated_values[obs_tuple] = fisher
            bonferroni_p = p_value * number_of_tests if (p_value * number_of_tests) < 1.0 else 1.0
            p_value_list.append((gene, p_value))

            all_traits[trait][gene] = {
                "NUGN": genedic[gene]["Non-unique Gene name"],
                "Annotation": genedic[gene]["Annotation"],
                "tpgp": stat_table["tpgp"],
                "tngp": stat_table["tngp"],
                "tpgn": stat_table["tpgn"],
                "tngn": stat_table["tngn"],
                "sens": (float(stat_table["tpgp"]) / num_pos * 100) if num_pos > 0 else 0.0,
                "spes": (float(stat_table["tngn"]) / num_neg * 100) if num_neg > 0 else 0.0,
                "OR": odds_ratio,
                "p_v": p_value,
                "B_p": bonferroni_p,
            }

        print("\nAdding p-values adjusted for testing multiple hypotheses")
        # Now calculate Benjamini-Hochberg p-values
        sorted_p_values = sorted(p_value_list, key=lambda x: x[1])  # Sorted list of tuples: (gene, p-value)

        # Find out which p-values are ties
        # Note: Changed from step-down to step-up (from least significant to most significant)
        tie = [sorted_p_values[i-1][1] == sorted_p_values[i][1] for i in xrange(1, len(sorted_p_values))]
        # Initialize dics of corrected p values
        bh_corrected_p_values = {}
        # The least significant gene is entered into the dic
        bh_corrected_p_values[sorted_p_values[len(sorted_p_values)-1][0]] = last_bh = sorted_p_values[len(sorted_p_values)-1][1]

        for (ind, (gene, p)) in reversed(list(enumerate(sorted_p_values[:-1]))):
            bh_corrected_p_values[gene] = min([last_bh, p*number_of_tests/(ind+1.0)]) if not tie[ind] else last_bh
            last_bh = bh_corrected_p_values[gene]

        # Now add values to dictionaries:
        for gene in genedic:
            if gene in all_traits[trait]:
                all_traits[trait][gene]["BH_p"] = bh_corrected_p_values[gene] if bh_corrected_p_values[gene] < 1.0 else 1.0

    return {"Results": all_traits, "Gene_trait_combinations": gene_trait_combinations}


def Perform_statistics(traits, genes):
    """
    The method that adds the presence/absence status for each trait-status combination
    """
    r = {"tpgp": 0, "tpgn": 0, "tngp": 0, "tngn": 0}  # tpgn = trait positive, gene negative
    gene_trait = {}
    for t in traits:  # For each strain
        try:
            if int(traits[t]) == 1 and genes[t] == 1:
                r["tpgp"] += 1
                gene_trait[t] = "AB"
            elif int(traits[t]) == 1 and genes[t] == 0:
                r["tpgn"] += 1
                gene_trait[t] = "aB"
            elif int(traits[t]) == 0 and genes[t] == 1:
                r["tngp"] += 1
                gene_trait[t] = "Ab"
            elif int(traits[t]) == 0 and genes[t] == 0:
                r["tngn"] += 1
                gene_trait[t] = "ab"
            else:
                sys.exit("There was a problem with comparing your traits and gene presence/absence files."
                "Make sure you have formatted the traits file to specification and only use 1s and 0s, and make sure "
                "the Roary file contains empty cells for non-present genes and non-empty text cells for present genes")
        except KeyError:
            print("\nError occured when trying to find " + str(t) + " in the genes file.")
            sys.exit("Make sure strains are named the same in your traits file as in your gene presence/absence file")
    return {"statistics": r, "gene_trait": gene_trait}


def StoreResults(Results, max_hits, p_cutoff, correctionmethod, upgmatree, GTC):
    """
    A method for storing the results. Calls StoreTraitResult for each trait column in the input file
    """
    for Trait in Results:
        print("\nStoring results: " + Trait)
        StoreTraitResult(Results[Trait], Trait, max_hits, p_cutoff, correctionmethod, upgmatree, GTC)


def StoreTraitResult(Trait, Traitname, max_hits, p_cutoff, correctionmethod, upgmatree, GTC):
    """
    The method that actually stores the results. Only accepts results from a single trait at a time
    """
    with open(Traitname + time.strftime("_%d_%m_%Y_%H%M") + ".csv", "w") as outfile:
        # Sort genes by p-value.
        sort_instructions = SortResultsAndSetKey(Trait)

        num_results = max_hits if max_hits is not None else len(Trait)

        cut_possibilities = {"Individual": "p_v", "Bonferroni": "B_p", "Benjamini-Hochberg": "BH_p"}

        outfile.write("Gene;Non-unique gene name;Annotation;Number_pos_present_in;Number_neg_present_in;Number_pos_not_present_in;"
        "Number_neg_not_present_in;Sensitivity;Specificity;Odds_ratio;Naive_p;Bonferroni_p;Benjamini_H_p;Max_Pairwise_comparisons;"
        "Max_supporting_pairs;Max_opposing_pairs;Best_pairwise_comp_p;Worst_pairwise_comp_p\n")

        print("Calculating max number of contrasting pairs for each significant gene")

        for x in xrange(num_results):
            sys.stdout.write("\r{:.2%}".format(float(x)/num_results))
            sys.stdout.flush()

            # Start with lowest p-value, the one which has key 0 in sort_instructions
            currentgene = sort_instructions[x]
            if (Trait[currentgene][cut_possibilities[correctionmethod]] > p_cutoff):
                sys.stdout.write("\r100.00%")
                sys.stdout.flush()
                break

            Max_pairwise_comparisons = ConvertUPGMAtoPhyloTree(upgmatree,
                                                               GTC[Traitname][currentgene])
            max_total_pairs = Max_pairwise_comparisons["Total"]
            max_propairs = Max_pairwise_comparisons["Pro"]
            max_antipairs = Max_pairwise_comparisons["Anti"]
            try:
                best_pairwise_comparison_p = ss.binom_test(max_propairs,
                                                           max_total_pairs,
                                                           0.5) / 2
                worst_pairwise_comparison_p = ss.binom_test(max_total_pairs-max_antipairs,
                                                            max_total_pairs,
                                                            0.5) / 2
            except TypeError:
                sys.exit("There was a problem using scipy.stats.binom_test. Ensure you have a recent distribution of SciPy installed.")

            outfile.write('"' + currentgene + '";"' + str(Trait[currentgene]["NUGN"]) + '";"' + str(Trait[currentgene]["Annotation"]) +
            '";"' + str(Trait[currentgene]["tpgp"]) + '";"' + str(Trait[currentgene]["tngp"]) + '";"' + str(Trait[currentgene]["tpgn"]) +
            '";"' + str(Trait[currentgene]["tngn"]) + '";"' + str(Trait[currentgene]["sens"]) + '";"' + str(Trait[currentgene]["spes"]) +
            '";"' + str(Trait[currentgene]["OR"]) + '";"' + str(Trait[currentgene]["p_v"]) + '";"' + str(Trait[currentgene]["B_p"]) +
            '";"' + str(Trait[currentgene]["BH_p"]) + '";"' + str(max_total_pairs) + '";"' + str(max_propairs) + '";"' + str(max_antipairs) +
            '";"' + str(best_pairwise_comparison_p) + '";"' + str(worst_pairwise_comparison_p) + '"\n')


def SortResultsAndSetKey(genedic):
    """
    A method for returning a dictionary where genes are sorted by p-value
    """
    return {i: gene for (i, gene) in enumerate(sorted(genedic,
                                                      key=lambda x: genedic[x]["p_v"])) }


def upgma(d):
    """
    Returns a UPGMA tree from a QuadTree distance matrix d. Heavily based on
    original implementation by Christian Storm Pedersen.
    """
    n = d.dim
    cluster = [[x] for x in d.names]
    size = n * [1]

    while n > 1:
        # While there are more than a single cluster left, find clusters i and j with minimum distance
        i, j = d.argmin()

        # Build list of distances from the new cluster to all the other clusters
        new_cluster = [cluster[i], cluster[j]]
        new_size = size[i] + size[j]
        new_dist = []
        for k in xrange(d.dim):
            if cluster[k] is None:
                new_dist.append(1)
            else:
                new_dist.append((d.get_elm(i,k)*size[i] + d.get_elm(j,k)*size[j]) / new_size)

        # Insert new row/col in d
        new_dist[i] = sys.maxsize
        d.insert_row(i, new_dist)
        d.insert_col(i, new_dist)
        d.insert_row(j, d.dim * [sys.maxsize])
        d.insert_col(j, d.dim * [sys.maxsize])

        cluster[i] = new_cluster
        cluster[j] = None
        size[i] = new_size
        size[j] = 0

        n -= 1

    return new_cluster


def ConvertUPGMAtoPhyloTree(tree, GTC):
    """
    A method that converts the upgma tree (in nested list form) to a PhyloTree.
    It also needs the status (AB, Ab, aB or ab) of all strains that are to be
    included in the analysis.
    """

    # TRAVERSING TREE: For each binary division - go to left until hit tip. Then go back
    num_AB = float(list(GTC.values()).count("AB"))
    num_Ab = float(list(GTC.values()).count("Ab"))
    num_aB = float(list(GTC.values()).count("aB"))
    num_ab = float(list(GTC.values()).count("ab"))
    OR = ((num_AB + 1)/(num_Ab + 1)) / ((num_aB + 1)/(num_ab + 1))  # Use pseudocounts to avoid 0 or inf OR.
    MyPhyloTree = PhyloTree(leftnode=tree[0],
                            rightnode=tree[1],
                            GTC=GTC,
                            OR=OR)

    return {"Total": MyPhyloTree.max_contrasting_pairs,
            "Pro": MyPhyloTree.max_contrasting_propairs,
            "Anti": MyPhyloTree.max_contrasting_antipairs}


def StoreUPGMAtreeToFile(upgmatree):
    """
    A method for printing the UPGMA tree that is built internally from the 
    hamming distances in the gene presence/absence matrix
    """
    treefilename = str("Tree" + time.strftime("_%d_%m_%Y_%H%M") + ".nwk")
    with open(treefilename, "w") as treefile:
        Tree = str(upgmatree)
        Tree = Tree.replace("[", "(")
        Tree = Tree.replace("]", ")")
        treefile.write(Tree)
        print("Wrote the UPGMA tree to file: %s" % treefilename)

if __name__ == '__main__':
    pass
