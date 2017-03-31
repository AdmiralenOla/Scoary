#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Script for associating Roary output with phenotypic trait (+/-)
# Idea and implementation: Ola Brynildsrud (olbb@fhi.no)

# Null hypothesis: Gene is equally distributed in Trait+ and Trait- 
# isolates

import argparse
import sys
import csv
import time
import random
import logging
from multiprocessing import Pool, Value
from scipy import stats as ss
from scipy import spatial
from .citation import citation
from .classes import Matrix
from .classes import QuadTree
from .classes import PhyloTree
from .classes import PimpedFileHandler
from .classes import ScoaryLogger
from .__init__ import __version__
#import scoary

import os
from pkg_resources import resource_string, resource_filename

#SCOARY_VERSION = scoary.__version__
SCOARY_VERSION = __version__

# Python 2/3 annoyances
try:
    xrange
except NameError:
    xrange = range

# Set up log and message flow
log = ScoaryLogger(logging.getLogger('scoary'))
log.setLevel(logging.DEBUG)

logformat = '%(asctime)s    %(message)s'
logdatefmt='%m/%d/%Y %I:%M:%S %p'
formatter = logging.Formatter(fmt=logformat,datefmt=logdatefmt)
   
def main(**kwargs):
    """
    The main function of Scoary.
    """
        
    # If main has been ran from the GUI, then args already exists
    if len(kwargs) == 0:
        args, cutoffs = ScoaryArgumentParser()
    else:
        args = kwargs["args"]
        cutoffs = kwargs["cutoffs"]
        sys.stdout = kwargs["statusbar"]

    # If the citation arg has been passed, nothing should be done except
    # a call to citation
    if args.citation:
        sys.exit(citation())
   
    # If the test argument was used, all settings are overrided and 
    # defaulted
    if args.test:
        args.correction = ['I','EPW']
        args.delimiter = ','
        args.genes = os.path.join(
            resource_filename(__name__, 'exampledata'), 
            'Gene_presence_absence.csv')
        args.max_hits = None
        args.newicktree = None
        args.no_pairwise = False
        args.outdir = './'
        args.permute = 0
        args.p_value_cutoff = [0.05,0.05]
        args.restrict_to = None
        args.start_col = 15
        args.threads = 4
        args.traits = os.path.join(
            resource_filename(__name__, 'exampledata'),
            'Tetracycline_resistance.csv')
        args.upgma_tree = True
        args.write_reduced = False
        args.no_time = False
        args.collapse = False
        cutoffs = {"I": 0.05, "EPW": 0.05}
    
    # Start clock
    starttime = time.time()
    # currenttime is ONLY used to append to filenames, therefore, it
    # is set empty if args.no_time is True
    if args.no_time:
        currenttime = ""
    else:
        currenttime = time.strftime("_%d_%m_%Y_%H%M")

    # Outdir should end with slash    
    if not args.outdir.endswith("/"):
        args.outdir += "/"
    
    # Initiate console from sys.stdout (which might be a statusbar)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.INFO)
    log.addHandler(console)
    
    # Create log file
    log_filename = os.path.join(args.outdir,
                                "scoary%s.log" % currenttime)
    log_handler = PimpedFileHandler(log_filename, mode='w')
    log_handler.setFormatter(formatter)
    log.addHandler(log_handler)
    log.info("==== Scoary started ====")
    
    # Catch all system exceptions from here on out
    try:

        # Perform tests to make sure all input is okay                        
        if args.traits is None or args.genes is None:
            sys.exit(("The following arguments are required: -t/--traits, "
            "-g/--genes"))
        if args.threads <= 0:
            sys.exit("Number of threads must be positive")
        if not os.path.isfile(args.traits):
            sys.exit("Could not find the traits file: %s" % args.traits)
        if not os.path.isfile(args.genes):
            sys.exit("Could not find the gene presence absence file: %s" 
            % args.genes)
        if ((args.newicktree is not None) and 
        (not os.path.isfile(args.newicktree))):
            sys.exit("Could not find the custom tree file: %s"
            % args.newicktree)
        if not all( [(p <= 1.0 and p > 0.0) for p in args.p_value_cutoff] ):
            sys.exit("P must be between 0.0 and 1.0 or exactly 1.0")
        if (len(args.delimiter) > 1):
            sys.exit("Delimiter must be a single character string. "
            "There is no support for tab.")
        if ((len(args.p_value_cutoff) != len(args.correction)) and
        (len(args.p_value_cutoff) != 1)):
            sys.exit("You can not use more p-value cutoffs than correction "
            "methods. Either provide a single p-value that will be applied "
            "to all correction methods, or provide exactly as many as the "
            "number of correction methods and in corresponding sequence. "
            "e.g. -c I EPW -p 0.1 0.05 will apply an individual p-value "
            "cutoff of 0.1 AND a pairwise comparisons p-value cutoff of "
            "0.05.")
        if "P" in cutoffs and args.permute == 0:
            sys.exit("Cannot use empirical p-values in filtration without "
            "performing permutations. Use '--permute X' where X is a "
            "number equal to or larger than 10")
        if args.permute < 10 and args.permute != 0:
            sys.exit("The absolute minimum number of permutations is 10 "
            "(or 0 to deactivate)")
        elif args.permute > 0:
            try:
                random.shuffle
            except:
                sys.exit("Unable to find random.shuffle. Can not proceed "
                "with permutations")
        if "P" in cutoffs:
            if cutoffs["P"] < (1.0/args.permute):
                sys.exit("Permutation cutoff too low for this number of "
                "permutations")
        if args.permute > 10000:
            log.info("Note: You have set Scoary to do a high number of "
            "permutations. This may take a while.")
        if args.no_pairwise:
            log.info("Performing no pairwise comparisons. Ignoring all "
            "tree related options (user tree, population aware-correction, "
            "permutations).")
            args.permute = 0
            args.newicktree = None
            for m in ["PW","EPW","P"]:
                cutoffs.pop(m,None)
            
        # Start analysis
        with open(args.genes, "rU") as genes, \
        open(args.traits, "rU") as traits:
    
            if args.restrict_to is not None:
                # Note: Despite being a dictionary, the values of
                # allowed_isolates are not currently used, only the keys
                allowed_isolates = {isolate : "all"
                                    for line in
                                    open(args.restrict_to,"rU")
                                    for isolate in line.rstrip().split(",")}
            else:
                # Despite the confusing name
                # this actually means all isolates are allowed
                # and included in the analysis
                allowed_isolates = None
                if args.write_reduced:
                    sys.exit("You cannot use the -w argument without "
                    "specifying a subset (-r)")

            log.info("Reading gene presence absence file")
            
            genedic_and_matrix = \
                Csv_to_dic_Roary(genes,
                args.delimiter,
                startcol=int(args.start_col) - 1,
                allowed_isolates=allowed_isolates,
                writereducedset=args.write_reduced,
                time=currenttime,
                outdir=args.outdir) 
            genedic = genedic_and_matrix["Roarydic"]
            zeroonesmatrix = genedic_and_matrix["Zero_ones_matrix"]
            strains = genedic_and_matrix["Strains"]
    
            # Create or load tree (No need for tree if --no_pairwise)
            if (args.newicktree) is None and not (args.no_pairwise):
                log.info("Creating Hamming distance matrix based on gene "
                "presence/absence")
                TDM = CreateTriangularDistanceMatrix(zeroonesmatrix,
                                                     strains)
                QT = PopulateQuadTreeWithDistances(TDM)
                log.info("Building UPGMA tree from distance matrix")
                upgmatree = upgma(QT)
            elif (args.no_pairwise):
                # Performing --no_pairwise analysis
                log.info("Ignoring relatedness among input sample and "
                "performing only population structure-naive analysis.")
                upgmatree = None
            else:
                log.info("Reading custom tree file")
                from .nwkhandler import ReadTreeFromFile
                upgmatree, members = ReadTreeFromFile(args.newicktree)
                if (sorted(strains) != sorted(members)):
                    # There could be two reasons for this. Either (a) that
                    # the gene presence absence file and the tree file does
                    # not have the same isolates or (b) that the -r flag has
                    # been used.
                    if args.restrict_to is None:
                        sys.exit("CRITICAL: Please make sure that isolates in "
                        "your custom tree match those in your gene "
                        "presence absence file.")
                    else:    
                        # Verify that the all isolates in strains (from GPA
                        # file) are in members (from tree), and if so prune
                        # tree of extraneous isolates. Else, call an error
                        if (all(i in members for i in strains)):
                            # Prune tree of the extraneous isolates
                            log.info("Pruning phylogenetic tree to correspond "
                            "to set of included isolates")
                            Prune = [i for i in members if 
                                          i not in strains]
                            upgmatree = PruneForMissing(upgmatree, Prune)
                        else:
                            sys.exit("CRITICAL: Your provided tree file did "
                            "not contain all the isolates in your gene "
                            "presence absence file.")
            log.info("Reading traits file")
            traitsdic, Prunedic = Csv_to_dic(traits,
                                             args.delimiter,
                                             allowed_isolates,
                                             strains)
            
            log.info("Finished loading files into memory.\n\n")
            log.info("==== Performing statistics ====")
            filtopts = filtrationoptions(cutoffs, args.collapse)
            for line in filtopts:
                log.info(line)

            log.info("Tallying genes and performing statistical analyses")
    
            RES_and_GTC = Setup_results(genedic, traitsdic, args.collapse)
            RES = RES_and_GTC["Results"]
            GTC = RES_and_GTC["Gene_trait_combinations"]
    
            if args.upgma_tree:
                StoreUPGMAtreeToFile(upgmatree,
                                     args.outdir,
                                     time=currenttime)        
    
            StoreResults(RES,
                         args.max_hits,
                         cutoffs,
                         upgmatree,
                         GTC,
                         Prunedic,
                         args.outdir,
                         args.permute,
                         args.threads,
                         args.no_pairwise,
                         time=currenttime,
                         delimiter=args.delimiter)
            log.info("\n")
            log.info("==== Finished ====")
            log.info("Checked a total of %d genes for "
                  "associations to %d trait(s). Total time used: %d "
                  "seconds." % (len(genedic),
                                len(traitsdic),
                                int(time.time()-starttime)))
    
    except SystemExit as e:
        exc_type, exc_value, _ = sys.exc_info()
        log.exception("CRITICAL:")
        log.removeHandler(log_handler)
        log.removeHandler(console)
        sys.exit(exc_value)
    
    if log.critical.called > 0:
        log.info("Scoary finished successfully, but with CRITICAL ERRORS. "
        "Please check your log file.")
    elif log.error.called > 0:
        log.info("Scoary finished successfully, but with ERRORS. Please check "
        "your log file.")
    elif log.warning.called > 0:
        log.info("Scoary finished successfully, but with WARNINGS. Please "
        "check your log file.")
    else:
        log.info("No warnings were recorded.")
    log.removeHandler(log_handler)
    log.removeHandler(console)
    sys.exit(0)

###############################
# FUNCTIONS FOR READING INPUT #
###############################
def Csv_to_dic_Roary(genefile, delimiter, startcol=14, 
    allowed_isolates=None, writereducedset=False, time="",
    outdir="./"):
    """
    Converts a gene presence/absence file into dictionaries
    that are readable by Scoary.
    """
    r = {}
    if writereducedset:
        file = open(ReduceSet(genefile,delimiter,startcol,
                    allowed_isolates,time,outdir),"rU")
        csvfile = csv.reader(file, skipinitialspace=True, 
                             delimiter=delimiter)
    else:
        csvfile = csv.reader(genefile, skipinitialspace=True, 
                             delimiter=delimiter)
    header = next(csvfile)
            
    roaryfile = True
    try:
        header[startcol]
    except IndexError:
        sys.exit("The startcol (-s) you have specified does not seem to "
        "correspond to any column in your gene presence/absence file.")
    strains = header[startcol:]
    
    # Include alternative column spellings
    Roarycols = ["Gene", "Non-unique Gene name", "Annotation",
                 "No. isolates", "No. sequences",
                 "Avg sequences per isolate", "Genome Fragment",
                 "Order within Fragment", "Accessory Fragment",
                 "Accessory Order with Fragment", "QC",
                 "Min group size nuc", "Max group size nuc",
                 "Avg group size nuc", "Order within fragment",
                 "Genome fragment", "Accessory fragment" ]
    
    # Move forwards from startcol to find correct startcol
    if strains[0] in Roarycols:
        for c in xrange(len(strains)):
            if strains[c] not in Roarycols:
                correctstartcol = startcol + c
                break
        log.error("ERROR: Make sure you have set the -s parameter "
        "correctly. You are running with -s %s. This correponds to the "
        "column %s. If this is not an isolate, Scoary might crash or "
        "produce strange results. Scoary thinks you should have run "
        "with -s %s instead" % (str(startcol+1),
                                strains[0],
                                str(correctstartcol + 1)))
    
    # Move backwards from startcol to find correct startcol
    Firstcols = header[:startcol][::-1]
    minus = 0
    Censored_isolates = []
    for c in xrange(len(Firstcols)):
        if Firstcols[c] not in Roarycols:
            minus += 1
            #print(Firstcols[c])
            Censored_isolates.append(Firstcols[c])
        else:
            if minus > 0:
                correctstartcol = startcol - minus
                log.error("ERROR: Make sure you have set the -s "
                "parameter correctly. You are running with -s %s. "
                "Scoary thinks you should have used %s. This excludes the "
                "following, which Scoary thinks are isolates: %s" % 
                (str(startcol+1),
                str(correctstartcol+1), 
                ", ".join(Censored_isolates)))
                break
    
    if allowed_isolates is not None:
        strain_names_allowed = [val for val in strains
                                if val in allowed_isolates.keys()]
    else:
        strain_names_allowed = strains

    zero_ones_matrix = []

    try:
        genecol = header.index("Gene")
        nugcol = header.index("Non-unique Gene name")
        anncol = header.index("Annotation")
    except ValueError:
        log.error("ERROR: Could not properly detect the correct names "
        "for all columns in the ROARY table.")
        genecol = 0
        nugcol = 1
        anncol = 2

    for line in csvfile:
        q = line
        try:
            r[q[genecol]] = ({"Non-unique Gene name": q[nugcol], 
                              "Annotation": q[anncol]} 
                              if roaryfile
                              else {})
        except IndexError:
            sys.exit("CRITICAL: Could not read gene presence absence "
            "file. Verify that this file is a proper Roary file using "
            "the specified delimiter (default is ',').")
        # The zero_ones_line variable represents the presence (1) or 
        # absence (0) of a gene. It is used for calculating distances  
        # between strains.
        zero_ones_line = []

        for strain in xrange(len(strains)):
            if (allowed_isolates is not None):
                if strains[strain] not in allowed_isolates.keys():
                    continue
            if q[startcol + strain] in ["", "0", "-"]:
                # If the gene is not present, AND The isolate is allowed
                r[q[genecol]][strains[strain]] = 0
                zero_ones_line.append(0)
                # Add a 0 to indicate non-presence
            else:
                # Gene is present if any other value than "", "0" or "-"
                # is in the cell
                r[q[genecol]][strains[strain]] = 1
                zero_ones_line.append(1)
                # Add a 1 to indicate presence of the current gene in 
                # this strain

        # Since we are only interested in the differences between 
        # strains, no need to append the zero_ones_line if it is all 1's
        # (core gene) or all 0's (not in collection)
        if 1 in zero_ones_line and 0 in zero_ones_line:
            zero_ones_matrix.append(zero_ones_line)

    # Transpose list for distance calculation purposes
    if writereducedset:
        file.close()
    zero_ones_matrix = list(map(list, zip(*zero_ones_matrix)))
    return {"Roarydic": r,
            "Zero_ones_matrix": zero_ones_matrix,
            "Strains": strain_names_allowed}

def ReduceSet(genefile, delimiter, startcol=14, allowed_isolates=None, 
              time="",outdir="./"):
    """
    Helper function for csv_to_dic_roary.
    Method for writing a reduced gene presence absence file, based only
    on isolates allowed by the restrict_to flag. This can speed up
    Scoary when analyzing subsets of large (e.g. more than a couple 
    of hundred) datasets.
    """
    csvfile = csv.reader(genefile, skipinitialspace=True,
                         delimiter=delimiter)
    header = next(csvfile)
    allowed_indexes = list(range(startcol))
    for c in xrange(len(header)):
        if header[c] in allowed_isolates.keys():
            allowed_indexes.append(c)
    
    log.info("Writing gene presence absence file for the reduced set of "
          "isolates")

    reducedfilename = \
        "%sgene_presence_absence_reduced%s.csv" % (outdir, time)
    
    with open(reducedfilename, "w") as csvout:
        wtr = csv.writer(csvout, delimiter = delimiter)
        newheader = [header[a] for a in allowed_indexes]
        wtr.writerow(newheader)
        for r in csvfile:
            wtr.writerow( tuple(r[a] for a in allowed_indexes) )
    log.info("Finished writing reduced gene presence absence list to "
    "file %s" % str(reducedfilename))
    return reducedfilename

def Csv_to_dic(csvfile, delimiter, allowed_isolates, strains):
    """
    Converts an input traits file (csv format) to dictionaries readable 
    by Scoary.
    """
    tab = list(zip(*csv.reader(csvfile, delimiter=delimiter)))
    r = {}
    # Create dictionary over which trees need pruning due to missing 
    # data
    Prunedic = {}
    if len(tab) < 2:
        sys.exit("Please check that your traits file is formatted "
        "properly and contains at least one trait")
    for num_trait in xrange(1, len(tab)):
        p = dict(zip(tab[0], tab[num_trait]))
        if "" in p:
            name_trait = p[""]
            del p[""]
        elif "Name" in p:
            name_trait = p["Name"]
            del p["Name"]
        else:
            sys.exit("Make sure the top-left cell in the traits file "
            "is either empty or 'Name'. Do not include empty rows")

        # Filter so that only allowed isolates are included
        if allowed_isolates is not None:
            p = {strain: indicator for (strain, indicator) in 
            list(p.items()) if strain in allowed_isolates.keys()}
        # Stop if unknown character found in traits file
        allowed_values = ["0","1","NA",".","-"," ",""]
        if not all([x in allowed_values for x in p.values()]):
            sys.exit("Unrecognized character found in trait file. Allowed "
            "values (no commas): %s" % str(",".join(allowed_values)))
        # Remove isolates with missing values, but only for the 
        # trait for which they are missing
        if ("NA" in p.values() 
        or "-" in p.values() 
        or "." in p.values()
        or " " in p.values()
        or "" in p.values()):
            log.warning("WARNING: Some isolates have missing values for "
            "trait %s. Missing-value isolates will not be counted in "
            "association analysis towards this trait." 
            % str(name_trait))
            p_filt = {strain: indicator for (strain, indicator) in 
                p.items() if indicator not in ["NA","-","."," ",""]}
            Prunedic[name_trait] = [k for (k,v) in p.items() if
                v in ["NA","-","."," ",""]]
            #p = p_filt
        else:
            Prunedic[name_trait] = []
            p_filt = p
        # Remove isolates that did not have rows in the trait file but
        # that were allowed by the GPA file/restrict_to
        #if not all(s in p.keys() for s in strains):
        if not all(s in p.keys() for s in strains):
            log.error("ERROR: Some isolates in your gene presence "
            "absence file were not represented in your traits file. "
            "These will count as MISSING data and will not be included."
            )

            Prunedic[name_trait] += [s for s in strains if
                                     s not in p.keys() and
                                     s not in Prunedic[name_trait]]
        r[name_trait] = p_filt
        Prunedic[name_trait] += [None]

    return r, Prunedic

################################
# FUNCTIONS FOR HANDLING TREES #
################################
def CreateTriangularDistanceMatrix(zeroonesmatrix, strainnames):
    """
    Converts a raw matrix of 0s and 1s (indicating gene presence and 
    absence, rows correspond to genes and columns to allowed strains) to 
    an (upper) triangular matrix of pairwise Hamming distances.
    The distance d(i,i) is set to 1 for all i.
    """
    try:
        hamming_distances = [float(h) for h in list(
        spatial.distance.pdist(zeroonesmatrix, 'hamming'))]
    except TypeError:
        sys.exit("Could not locate scipy.spatial.distance.pdist. "
        "Perhaps you have an old version of SciPy installed?")
    nstrains = int((1 + (1 + 8*len(hamming_distances))**0.5)/2)
    TriangularDistanceMatrix = []
    Strain_names = []
    for i in xrange(nstrains):
        # Adding the maximum relative hamming distance to prevent 
        # Quadtree algorithm to pick a pair where i = j
        add = [1]  
        add += hamming_distances[:(nstrains-i-1)]
        hamming_distances = hamming_distances[(nstrains-i-1):]
        TriangularDistanceMatrix.append(add)
        Strain_names.append(strainnames[i])

    return {"matrix": TriangularDistanceMatrix, "names": Strain_names}

def PopulateQuadTreeWithDistances(TDM):
    """
    Takes a triangular distance matrix, creates a quadtree and populates
    it with the hamming distances between isolates.
    First creates a Quadmatrix, so not really optimized.
    """
    Quadmatrix = Matrix(dim=len(TDM["matrix"]))
    for i in xrange(Quadmatrix.dim):
        for j in xrange(i, Quadmatrix.dim):
            try:
                Quadmatrix[i][j] = \
                Quadmatrix[j][i] = TDM["matrix"][i][(j-i)]
            except IndexError:
                sys.exit("There was an error trying to populate the "
                         "Quadtree with pairwise distances. Please "
                         "report this bug to olbb@fhi.no")
    PopulatedQuadtree = QuadTree(Quadmatrix.dim, names=TDM["names"])
    for i in xrange(Quadmatrix.dim):
        PopulatedQuadtree.insert_row(i, Quadmatrix[i])
    return PopulatedQuadtree

def upgma(d):
    """
    Returns a UPGMA tree from a QuadTree distance matrix d. Heavily 
    based on original implementation by Christian Storm Pedersen.
    """
    n = d.dim
    cluster = [x for x in d.names]
    size = n * [1]

    while n > 1:
        # While there are more than a single cluster left, find clusters
        # i and j with minimum distance
        i, j = d.argmin()

        # Build list of distances from the new cluster to all the other
        # clusters
        new_cluster = [cluster[i], cluster[j]]
        new_size = size[i] + size[j]
        new_dist = []
        for k in xrange(d.dim):
            if cluster[k] is None:
                new_dist.append(1)
            else:
                new_dist.append((d.get_elm(i,k)*size[i] + 
                                 d.get_elm(j,k)*size[j]) / new_size)

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

def PruneForMissing(tree, Prunedic):
    """
    Method for pruning isolates with missing values from the upgma tree
    and gene-trait combinations objects. This happens for a particular
    trait only, and so the full tree and GTC objects are used for traits
    without missing values
    """
    # Traverse tree and prune missing-data isolates
    # Left node is a subtree, go deeper

    # Make slice copy to prevent tampering with original
    tree = tree[:]
    
    if isinstance(tree[0],list):
        tree[0] = PruneForMissing(tree[0], Prunedic)
    # Right node is a subtree, go deeper
    if isinstance(tree[1], list):
        tree[1] = PruneForMissing(tree[1], Prunedic)
 
    # Both isolates have missing data and should be excised
    if (tree[0] in Prunedic) and (tree[1] in Prunedic):
        return
    # Only left node has missing data and should be excised
    elif tree[0] in Prunedic:
        return tree[1]
    # Only right node has missing data and should be excised
    elif tree[1] in Prunedic:
        return tree[0]
    # Neither the left or right nodes have missing data
    else:
        return [tree[0], tree[1]]

def StoreUPGMAtreeToFile(upgmatree, outdir, time=""):
    """
    A method for printing the UPGMA tree that is built internally from 
    the hamming distances in the gene presence/absence matrix
    """
    treefilename = str(outdir + ("Tree%s.nwk" % time))
    with open(treefilename, "w") as treefile:
        Tree = str(upgmatree)
        Tree = Tree.replace("[", "(")
        Tree = Tree.replace("]", ")")
        treefile.write(Tree + ";")
        log.info("Wrote the UPGMA tree to file: %s" % treefilename)       

#####################################
# FUNCTIONS FOR CALCULATING RESULTS #
#####################################
def Setup_results(genedic, traitsdic, collapse):
    """
    This is the method that actually does all the counting of genes,
    calculation of p-values and post-test adjustment of p-values.
    The majority of running time is currently spent doing Fisher's exact
    test. The running time is somewhat shortened by storing the p-values
    of known 2x2 cell distributions, so that identical Fisher's tests 
    will not have to be run.
    """
    # Need to create one results dictionary for each trait

    all_traits = {}
    gene_trait_combinations = {}
    fisher_calculated_values = {}
    for trait in traitsdic:
        all_traits[trait] = {}
        # Also need a list of all p-values to perform stepwise FDR corr
        p_value_list = []
        log.info("Gene-wise counting and Fisher's exact tests for trait: " 
              "%s" % str(trait))

        # We also need a number of tests variable for each genedic.
        # (The number of tests may not be the same, depending on how 
        # many genes are in 0 or all strains.)
        number_of_tests = len(genedic)
        # Need initial number of tests for status bar
        initial_number_of_tests = number_of_tests

        # Additionally, we need a dictionary to, for each gene
        # hold the gene-trait status (e.g "AB" or "ab" of each strain)
        gene_trait_combinations[trait] = {}
        distributions = {}
                
        i = 0  # Progress
        for gene in genedic:
            # Status:
            sys.stdout.write("\r{:.2%}".format(float(i) / 
                                               initial_number_of_tests))
            sys.stdout.flush()
            i += 1  # Progress

            stats = Perform_statistics(traitsdic[trait], genedic[gene])
            stat_table = stats["statistics"]

            num_pos = stat_table["tpgp"] + stat_table["tpgn"]
            num_neg = stat_table["tngp"] + stat_table["tngn"]

            if (stat_table["tpgp"] + stat_table["tngp"]) == 0:
                # No included isolates have the gene
                # Remove 1 from the number of tests
                number_of_tests -= 1  
                # proceed to the next gene
                continue  

            if (stat_table["tpgn"] + stat_table["tngn"]) == 0:
                # All included isolates have the gene
                number_of_tests -= 1
                continue

            # IF gene_trait_combination exists already, ie the current 
            # gene is perfectly correlated with another gene, collapse
            # the two into a single unit
            
            # Should consider adding a softer correlation required than 
            # 100%, i.e. genes that are highly co-distributed should
            #  perhaps be merged earlier?
            if stats["hash"] in distributions and collapse:
                # Find out which gene is correlated
                corr_gene = distributions[stats["hash"]]
                # Collapse the two genes
                newgenename = corr_gene + "-" + gene
                distributions[stats["hash"]] = newgenename
                # Remove 1 from the number of tests
                number_of_tests -= 1
                mergedgenes = True
                
                del(gene_trait_combinations[trait][corr_gene])
                gene_trait_combinations[trait][newgenename] = \
                                                     stats["gene_trait"]
            else:
                distributions[stats["hash"]] = gene
                mergedgenes = False
                gene_trait_combinations[trait][gene] = \
                                                     stats["gene_trait"]
  
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
                 
            if not mergedgenes:
                all_traits[trait][gene] = {
                "NUGN": genedic[gene]["Non-unique Gene name"],
                "Annotation": genedic[gene]["Annotation"],
                "tpgp": stat_table["tpgp"],
                "tngp": stat_table["tngp"],
                "tpgn": stat_table["tpgn"],
                "tngn": stat_table["tngn"],
                "sens": ((float(stat_table["tpgp"]) / num_pos * 100) if 
                    num_pos > 0 else 0.0),
                "spes": ((float(stat_table["tngn"]) / num_neg * 100) if 
                    num_neg > 0 else 0.0),
                "OR": odds_ratio,
                "p_v": p_value}
                p_value_list.append((gene, p_value))
            else:
                del(all_traits[trait][corr_gene])
                all_traits[trait][newgenename] = {
                "NUGN": "",
                "Annotation": "Merged genes",
                "tpgp": stat_table["tpgp"],
                "tngp": stat_table["tngp"],
                "tpgn": stat_table["tpgn"],
                "tngn": stat_table["tngn"],
                "sens": ((float(stat_table["tpgp"]) / num_pos * 100) if 
                    num_pos > 0 else 0.0),
                "spes": ((float(stat_table["tngn"]) / num_neg * 100) if
                    num_neg > 0 else 0.0),
                "OR": odds_ratio,
                "p_v": p_value}
                p_value_list.append((newgenename, p_value))
        
        sys.stdout.write("\r100.00%")
        sys.stdout.flush()
        sys.stdout.write("\n")
        log.info("Adding p-values adjusted for testing multiple hypotheses")
        
        # Now calculate Benjamini-Hochberg and Bonferroni p-values
        # Sorted list of tuples: (gene, p-value)
        # s_p_v = abbreviation for sorted_p_values
        s_p_v = sorted(p_value_list, key=lambda x: x[1])  

        # Find out which p-values are ties
        # Note: Changed from step-down to step-up
        # (from least significant to most significant)
        tie = [s_p_v[i-1][1] == s_p_v[i][1] for i in
            xrange(1, len(s_p_v))]
        # Initialize dics of corrected p values
        # bh_c_p_v = abbreviation for bh_corrected_p_values
        bh_c_p_v = {}
        # The least significant gene is entered into the dic
        bh_c_p_v[s_p_v[len(s_p_v)-1][0]] = last_bh = s_p_v[len(s_p_v)-1][1]

        for (ind, (gene, p)) in reversed(list(enumerate(s_p_v[:-1]))):
            bh_c_p_v[gene] = min([last_bh, p*number_of_tests/(ind+1.0)]\
                                ) if not tie[ind] else last_bh
            last_bh = bh_c_p_v[gene]

        # Now add values to dictionaries:
        for gene in all_traits[trait]:
            all_traits[trait][gene]["B_p"] = min(
                all_traits[trait][gene]["p_v"] * number_of_tests , 1.0)
            all_traits[trait][gene]["BH_p"] = min(bh_c_p_v[gene], 1.0)

    return {"Results": all_traits,
            "Gene_trait_combinations": gene_trait_combinations}

def Perform_statistics(traits, genes):
    """
    Helper function for Setup_results.
    The method that adds the presence/absence status for each 
    trait-status combination.
    """
    r = {"tpgp": 0, "tpgn": 0, "tngp": 0, "tngn": 0}
    # tpgn = trait positive, gene negative, etc
    gene_trait = {}
    distribution_hash = ""
    for t in traits:  # For each strain
        # If missing data, continue without including isolate
        # Note, this should no longer be invoked, since missing-value
        # isolates are not in the traits dic to begin with
        if traits[t] in ["NA","-","."]:
            gene_trait[t] = "NA"
            distribution_hash += str(genes[t])
            continue

        try:
            if int(traits[t]) == 1 and genes[t] == 1:
                r["tpgp"] += 1
                distribution_hash += "1"
                gene_trait[t] = "AB"
            elif int(traits[t]) == 1 and genes[t] == 0:
                r["tpgn"] += 1
                distribution_hash += "0"
                gene_trait[t] = "aB"
            elif int(traits[t]) == 0 and genes[t] == 1:
                r["tngp"] += 1
                distribution_hash += "1"
                gene_trait[t] = "Ab"
            elif int(traits[t]) == 0 and genes[t] == 0:
                r["tngn"] += 1
                distribution_hash += "0"
                gene_trait[t] = "ab"
            else:
                sys.exit("There was a problem with comparing your "
                "traits and gene presence/absence files. Make sure you "
                "have formatted the traits file to specification and "
                "only use 1s and 0s, as well as NA, - or . for missing "
                "data. Also make sure the Roary file contains empty "
                "cells for non-present genes and non-empty text cells "
                "for present genes.")
        except KeyError:
            sys.stdout.write("\n")
            log.critical("CRITICAL: Could not find %s in the genes "
            "file." % str(t) )
            #log.warning("CRITICAL: Could not find %s in the genes file" % str(t))
            sys.exit("Make sure strains are named the same in your "
            "traits file as in your gene presence/absence file")
    return {"statistics": r, "hash": int(distribution_hash, 2),
            "gene_trait": gene_trait}

#################################
# FUNCTIONS FOR CREATING OUTPUT #
#################################
def StoreResults(Results, max_hits, cutoffs, upgmatree, GTC, Prunedic,
                 outdir, permutations, num_threads, no_pairwise,
                 time="", delimiter=","):
    """
    A method for storing the results. Calls StoreTraitResult for each 
    trait column in the input file.
    """
    for Trait in Results:
        sys.stdout.write("\n")
        log.info("Storing results: " + Trait)
        StoreTraitResult(Results[Trait], Trait, max_hits, cutoffs, 
                         upgmatree, GTC, Prunedic, outdir, permutations, 
                         num_threads, no_pairwise, time, delimiter)

def StoreTraitResult(Trait, Traitname, max_hits, cutoffs, upgmatree, 
                     GTC, Prunedic, outdir, permutations, num_threads, 
                     no_pairwise, time="", delimiter=","):
    """
    The method that actually stores the results. Only accepts results 
    from a single trait at a time
    """
    permutations = int(permutations)
    if num_threads > 1:
        multithreaded = True
    else:
        multithreaded = False
    
    fname = (outdir + Traitname + time + '.results.csv')

    with open(fname, "w") as outfile:
        # Sort genes by p-value.
        sort_instructions = SortResultsAndSetKey(Trait)
        if max_hits is None:
            max_hits = len(Trait)
        num_results = min(max_hits, len(Trait))

        cut_possibilities = {
            "I": "p_v",
            "B": "B_p",
            "BH": "BH_p",
            "PW": "Plowest",
            "EPW": "Pboth",
            "P": "Empirical_p"
        }
        
        columns = ["Gene","Non-unique gene name","Annotation",
        "Number_pos_present_in","Number_neg_present_in",
        "Number_pos_not_present_in","Number_neg_not_present_in",
        "Sensitivity","Specificity","Odds_ratio","Naive_p",
        "Bonferroni_p",
        "Benjamini_H_p"]

        if not no_pairwise:
            columns += ["Max_Pairwise_comparisons",
                        "Max_supporting_pairs","Max_opposing_pairs",
                        "Best_pairwise_comp_p","Worst_pairwise_comp_p"]
        
        if permutations >= 10:
            columns.append("Empirical_p")
        
        outfile.write(delimiter.join('"' + c +'"' for c in columns) + 
                      "\n")
        
        if permutations >= 10:
            log.info("Calculating max number of contrasting pairs for "
            "each significant gene and performing %s permutations" 
            % str(permutations))
        elif no_pairwise:
            log.info("Skipping population structure-aware analyses.")
        else:
            log.info("Calculating max number of contrasting pairs for "
            "each nominally significant gene")
        

        if not no_pairwise:    
            # If some isolates have missing values, prune tree and GTC of 
            # these isolates
            if len(Prunedic[Traitname]) > 0:
                upgmatree = PruneForMissing(upgmatree, Prunedic[Traitname])
    
        # Use threads to resolve pairwise comparisons and permutations
        # Each subprocess will be responsible for a domain of the 
        # results list
        
            if multithreaded:
                domains = [list(xrange(tnum, num_results, num_threads)) 
                           for tnum in xrange(num_threads)]
                # The Progress variable holds the total progress status of 
                # the subprocesses. It is locked, but can be written to by 
                # all subprocs
                Progress = Value("h",0)
                pool = Pool(processes = num_threads,initializer=initProcess,
                            initargs=(Progress,))
                argumentlist = {"si":sort_instructions,"tree":upgmatree,
                                "GTC": GTC[Traitname],"cutoffs": cutoffs,
                                "cp": cut_possibilities,"perm":permutations,
                                "Trait":Trait, "Threaded": True}
                # Need to make a copy of dict to prevent then from pointing 
                # to the same object
                all_args = [ (domains[x], dict(argumentlist)) 
                             for x in xrange(len(domains)) ] 
        
                Threadresults = pool.imap(PairWiseComparisons,all_args)

                pool.close()
                pool.join()
            else:
                domains = list(xrange(num_results))
                argumentlist = {"si":sort_instructions,"tree":upgmatree,
                                "GTC": GTC[Traitname],"cutoffs": cutoffs,
                                "cp": cut_possibilities,"perm":permutations,
                                "Trait":Trait, "Threaded": False}
                all_args = (domains, argumentlist)
                Threadresults = [PairWiseComparisons(all_args)]
            sys.stdout.write("\r100.00%")
            sys.stdout.flush()
            sys.stdout.write("\n")

            # Wait for each thread to finish and weave results from all 
            # threads

            # Finally, create a filteredresults that merges results data 
            # from Trait[currentgene] with Threadresult[currentgene]
            Filteredresults = {}
            Threadresults = list(Threadresults)

            for thread in xrange(num_threads):
                for currentgene in Threadresults[thread]:
                    Filteredresults[currentgene] = \
                        Threadresults[thread][currentgene]
                    Filteredresults[currentgene].update(Trait[currentgene])
       
            # Only sort by pairwise values if they are the only present 
            # correction methods
            if ("I" in cutoffs) or ("B" in cutoffs) or ("BH" in cutoffs):
                sort_instructions = SortResultsAndSetKey(Filteredresults)
            elif ("PW" in cutoffs) or ("EPW" in cutoffs):
                sort_instructions = \
                    SortResultsAndSetKeyPairwise(Filteredresults,cutoffs)
            elif ("P" in cutoffs):
                sort_instructions = \
                    SortResultsAndSetKey(Filteredresults,key="Empirical_p")
            else:
                log.info("No filtration applied")
        
        # If no_pairwise mode was run - merge up after this clause.
        else:
            Filteredresults = Trait
            sort_instructions = SortResultsAndSetKey(Filteredresults)

        log.info("Storing results to file")       
        num_filteredresults = min(num_results, len(Filteredresults))
        for x in xrange(num_filteredresults):
            currentgene = sort_instructions[x]
            
            # Final check that the currentgene passes all filters
            if all( [Filteredresults[currentgene]\
                    [cut_possibilities[method]] <= cutoffs[method]
                    for method in cutoffs] ):
                outrow = [currentgene, 
                str(Filteredresults[currentgene]["NUGN"]),
                str(Filteredresults[currentgene]["Annotation"]), 
                str(Filteredresults[currentgene]["tpgp"]),
                str(Filteredresults[currentgene]["tngp"]),
                str(Filteredresults[currentgene]["tpgn"]),
                str(Filteredresults[currentgene]["tngn"]),
                str(Filteredresults[currentgene]["sens"]),
                str(Filteredresults[currentgene]["spes"]),
                str(Filteredresults[currentgene]["OR"]),
                str(Filteredresults[currentgene]["p_v"]),
                str(Filteredresults[currentgene]["B_p"]),
                str(Filteredresults[currentgene]["BH_p"])]
                if not no_pairwise:
                    outrow += [
                    str(Filteredresults[currentgene]["max_total_pairs"]),
                    str(Filteredresults[currentgene]["max_propairs"]),
                    str(Filteredresults[currentgene]["max_antipairs"]),
                    str(Filteredresults[currentgene]["Pbest"]),
                    str(Filteredresults[currentgene]["Pworst"])]
            
                    # If permutations have been performed, print empirical 
                    # p-values as well
                    if permutations >= 10:
                        outrow.append(
                        str(Filteredresults[currentgene]["Empirical_p"]))
                        
                outfile.write(
                delimiter.join('"' + c + '"' for c in outrow) + "\n"
                )

def initProcess(share):
    """
    Initiates the progressbar as a shared variable between multiple
    processes.
    """
    #scoary.Progress = share
    global Progress
    Progress = share

def PairWiseComparisons(nestedlist):
    """
    A threaded version of the pairwise comparisons correction method. 
    Also calculates permutations.
    """
    domain = nestedlist[0]
    sort_instructions = nestedlist[1]["si"]
    tree = nestedlist[1]["tree"]
    GTC = nestedlist[1]["GTC"]
    cutoffs = nestedlist[1]["cutoffs"]
    cutposs = nestedlist[1]["cp"]
    permutations = nestedlist[1]["perm"]
    Trait = nestedlist[1]["Trait"]
    Threaded = nestedlist[1]["Threaded"]
    
    resultscontainer = {}
    num_tot_res = len(Trait)
    progress = 0.0

    for genenumber in domain:
        if Threaded:
            #with scoary.Progress.get_lock():
            #    scoary.Progress.value += 1
            with Progress.get_lock():
                Progress.value += 1
            sys.stdout.write(
            #"\r{:.2%}".format(float(scoary.Progress.value)/num_tot_res)
            "\r{:.2%}".format(float(Progress.value)/num_tot_res)
            )
            sys.stdout.flush()
        else:
            progress += 1.0
            sys.stdout.write("\r{:.2%}".format(progress/num_tot_res))
            sys.stdout.flush()
        # Start with lowest p-value, the one which has key 0 in 
        # sort_instructions
        currentgene = sort_instructions[genenumber]
        resultscontainer[currentgene] = {}
        
        Max_pairwise_comparisons = ConvertUPGMAtoPhyloTree(
            tree, GTC[currentgene])

        resultscontainer[currentgene]["max_total_pairs"] = \
            Max_pairwise_comparisons["Total"]
        resultscontainer[currentgene]["max_propairs"] = \
            Max_pairwise_comparisons["Pro"]
        resultscontainer[currentgene]["max_antipairs"] = \
            Max_pairwise_comparisons["Anti"]
        
        # Switch definition of "best" and "worst" if antipairs outnumber
        # propairs
        if (Max_pairwise_comparisons["Pro"] >= 
            Max_pairwise_comparisons["Anti"]):
            Pbest = "Pbest"
            Pworst = "Pworst"
        else:
            Pbest = "Pworst"
            Pworst = "Pbest"
        try:
            resultscontainer[currentgene][Pbest] = ss.binom_test(
                resultscontainer[currentgene]["max_propairs"],
                resultscontainer[currentgene]["max_total_pairs"],
                0.5)
            resultscontainer[currentgene][Pworst] = ss.binom_test(
                resultscontainer[currentgene]["max_total_pairs"] - 
                resultscontainer[currentgene]["max_antipairs"],
                resultscontainer[currentgene]["max_total_pairs"],
                0.5)
            resultscontainer[currentgene]["Plowest"] = min(
                resultscontainer[currentgene]["Pbest"], 
                resultscontainer[currentgene]["Pworst"])
            resultscontainer[currentgene]["Pboth"] = max(
                resultscontainer[currentgene]["Pbest"],
                resultscontainer[currentgene]["Pworst"])
        except TypeError:
            sys.exit("There was a problem using " 
            "scipy.stats.binom_test. Ensure you have a recent "
            "distribution of SciPy installed.")

        # Loop can break early if filtration includes non-pairwise 
        # correction measures, because if these are all above the p
        # then we are not interested in the pairwise results
        
        resultscontainer[currentgene]["p_v"] = Trait[currentgene]["p_v"]
        resultscontainer[currentgene]["B_p"] = Trait[currentgene]["B_p"]
        resultscontainer[currentgene]["BH_p"] = Trait[currentgene]["BH_p"]
        
        if decideifbreak(cutoffs, resultscontainer[currentgene]):
            # Genes that are cut from the break
            # should not be in the results
            del(resultscontainer[currentgene])
            break    

        # This is also the place to add permutations - 
        # AFTER it has been verified that the current gene passes all 
        # filtration filters
        if permutations >= 10:
            empirical_p_perm = Permute(tree=tree,
                                       GTC=GTC[currentgene],
                                       permutations=permutations,
                                       cutoffs=cutoffs)
            resultscontainer[currentgene]["Empirical_p"] = \
                empirical_p_perm
    
    return resultscontainer
    
def Permute(tree, GTC, permutations, cutoffs):
    """
    A method used to permute the GTC dataset by label-switching a given 
    number of times and calculate corresponding pairwise comparisons
    p-values for each PhyloTree. Finally, it returns the percentile of 
    the unpermuted dataset p-value in the set of ordered p-values 
    resulting from the permutation.
    """
    proceed = True
    # Guard the original GTC against changes by making a copy
    GTCcopy = dict(GTC)
    # r is the number of replicates with a larger test statistic than 
    # the unpermuted
    r = 0
    Unpermuted_tree = ConvertUPGMAtoPhyloTree(tree,GTCcopy)
    
    # Find out if we're after propairs or antipairs - which is more 
    # common and would thus give the high test statistic?

    if Unpermuted_tree["Pro"] >= Unpermuted_tree["Anti"]:
        Pro = "Pro"
    else:
        Pro = "Anti"
    Unpermuted_estimator = (float(Unpermuted_tree[Pro]) / 
                           Unpermuted_tree["Total"])

    if permutations < 10:
        # The program is never supposed to go here since fewer than 10
        # permutations are globally disallowed
        sys.stdout.write("Number of permutations too few. The absolute "
        "minimum is 10.")
        proceed = False

    if proceed:
        for i in xrange(permutations):
            # Make a permutation using random.shuffle
            PermutedGTC = PermuteGTC(GTCcopy)
            # Send new set to phylotree
            NewPhyloTree = ConvertUPGMAtoPhyloTree(tree, PermutedGTC)
            if (float(NewPhyloTree[Pro]) / NewPhyloTree["Total"] >= 
            Unpermuted_estimator):
                r += 1

            # Check how many estimators are higher than the unpermuted
            # If, after more than 30 iterations, r indicates a p > 0.1,
            # abort
            if i >= 30:
                if (1 - ss.binom.cdf(r,i,0.1)) < 0.05:
                    emp_p = (r+1.0)/(i+2.0)
                    break
        else:
            emp_p = (r+1.0)/(permutations+1.0)
    else:
        return None
    
    return emp_p

def PermuteGTC(GTC):
    """
    Returns a permutation of the gene-trait combination dic where trait 
    labels have been swapped
    """
    # Set of trait labels to distribute
    trait_labels = [s[-1] for s in GTC.values()]
    # Shuffle list of labels
    random.shuffle(trait_labels)
    for isolate in GTC:
        # Assign the latter character of GTC to isolates but keep the 
        # gene status intact
        GTC[isolate] = str(GTC[isolate][0] + trait_labels.pop())
    return GTC

def ConvertUPGMAtoPhyloTree(tree, GTC):
    """
    A method that converts the upgma tree (in nested list form) to a 
    PhyloTree. It also needs the status (AB, Ab, aB or ab) of all
    strains that are to be included in the analysis. As of 1.6.0 this
    function is also called from permutations
    """

    # TRAVERSING TREE: For each binary division - go to left until hit 
    # tip. Then go back
    MyPhyloTree = PhyloTree(leftnode=tree[0],
                            rightnode=tree[1],
                            GTC=GTC)

    return {"Total": MyPhyloTree.max_contrasting_pairs,
            "Pro": MyPhyloTree.max_contrasting_propairs,
            "Anti": MyPhyloTree.max_contrasting_antipairs}

##############################
# CURRENTLY UNUSED FUNCTIONS #
##############################
def FindClosestIndex(sortedlist, value, index=None):
    """
    Takes a sorted list and a value and finds which index position the 
    value would have if entered into the list. Does this through a
    binary search
    """
    # METHOD DEPRECATED WHEN P-VALUE STATISTICS NOT NEEDED
    if index is None:
        index = 0
        
    listlen = len(sortedlist)
    if listlen == 0:
        return None
    if listlen == 1:
        if value <= sortedlist[0]:
            return 0
        else:
            return 1
    # Bisect list
    middleindex = listlen // 2
    middlevalue = sortedlist[middleindex]
    if value == middlevalue:
        # Ties in list. Place value exactly at index
        return middleindex
    elif value < middlevalue:
        # Index in left part of list
        index += FindClosestIndex(sortedlist[:middleindex],
                                  value,
                                  index)
    else:
        # Index in right part of list
        index += (middleindex + 
                  FindClosestIndex(sortedlist[middleindex:],
                                   value,
                                   index))
    
    return index   
   
#########################
# MISC HELPER FUNCTIONS #
#########################            
def SortResultsAndSetKey(genedic,key="p_v"):
    """
    A method for returning a dictionary where genes are sorted by 
    p-value
    """
    return {i: gene for (i, gene) in enumerate(sorted(genedic,
                                       key=lambda x: genedic[x][key])) }
                                                      
def SortResultsAndSetKeyPairwise(genedic,correctionmethod):
    """
    A method for returning a dictionary where genes are sorted by 
    pairwise comparison p-values
    """
    if "EPW" in correctionmethod:
        return {i: gene for (i, gene) in enumerate(sorted(genedic,
                                   key=lambda x: genedic[x]["Pboth"])) }

    elif "PW" in correctionmethod:
        return {i: gene for (i, gene) in enumerate(sorted(genedic,
                                 key=lambda x: genedic[x]["Plowest"])) }
    else:
        sys.exit("Something went wrong when using this set of "
        "correction methods. Please report this bug")
                                                      
def filtrationoptions(cutoffs, collapse):
    """
    Converts between easy keys (I, B, etc) and long names of filtration 
    options, for printing to the runlog
    """
    translation = {"I": "Individual (Naive)", "B": "Bonferroni", 
                   "BH":"Benjamini-Hochberg",
                   "PW":"Pairwise comparison (Best)", 
                   "EPW": "Pairwise comparison (Entire range)",
                   "P": "Empirical p-value (permutation-based)"}
    filters = [str(translation[k]) + ":    " + str(v)
               for k,v in cutoffs.items()] 
    filters.append("Collapse genes:    " + str(collapse) + "\n\n")
    filters = ["-- Filtration options --"] + filters

    return filters
    
def decideifbreak(cutoffs, currentgene):
    """
    Decides if the loop should continue or if all filtration thresholds 
    are violated and therefore should break.
    """
    if (("I" in cutoffs) or ("B" in cutoffs) or ("BH" in cutoffs)):
        if "I" in cutoffs and currentgene["p_v"] > cutoffs["I"]:
            return True
        elif "B" in cutoffs and currentgene["B_p"] > cutoffs["B"]:
            return True
        elif "BH" in cutoffs and currentgene["BH_p"] > cutoffs["BH"]:
            return True
        else:
            # All p-values are too low. Cannot break yet.
            return False
    else:
        # Can't break early if none of the basic filtrations are used, 
        # because there might be other genes downstream that pass
        # pairwise and empirical filtration 
        return False

def ScoaryArgumentParser():
    """
    Function for parsing the arguments separate from the main function.
    Makes it easier for other scripts to run Scoary without artificially
    constructing a command line
    """
    # Parse arguments.
    parser = argparse.ArgumentParser(
            description='Scoary version %s - Screen pan-genome for '
                        'trait-associated genes' % SCOARY_VERSION,
                        epilog='by Ola Brynildsrud (olbb@fhi.no)')
    parser.add_argument('-t', '--traits',
                        help='Input trait table '
                        '(comma-separated-values). Trait presence is '
                        'indicated by 1, trait absence by 0. Assumes '
                        'strain names in the first column and trait '
                        'names in the first row')
    parser.add_argument('-g', '--genes',
                        help='Input gene presence/absence table '
                        '(comma-separated-values) from ROARY. '
                        'Strain names must be equal to those in the '
                        'trait table')
    parser.add_argument('-o', '--outdir',
                        help='Directory to place output files.'
                        ' Default = .',
                        default='./')
    parser.add_argument('-p', '--p_value_cutoff',
                        help='P-value cut-off / alpha level. '
                        'For Fishers, Bonferronis, and '
                        'Benjamini-Hochbergs tests, SCOARY will not '
                        'report genes with higher p-values than this. '
                        'For empirical p-values, this is treated as an '
                        'alpha level instead. I.e. 0.02 will filter '
                        'all genes except the lower and upper '
                        'percentile from this test. Run with "-p 1.0" '
                        'to report all genes. Accepts standard form '
                        '(e.g. 1E-8). Provide a single value (applied '
                        'to all) or exactly as many values as '
                        'correction criteria and in corresponding '
                        'order. (See example under correction). '
                        'Default = 0.05',
                        nargs='+',
                        default=[0.05],
                        type=float)
    parser.add_argument('-c', '--correction',
                        help='Apply the indicated filtration measure. '
                        'I=Individual (naive) p-value. '
                        'B=Bonferroni adjusted p-value. '
                        'BH=Benjamini-Hochberg adjusted p. '
                        'PW=Best (lowest) pairwise comparison. '
                        'EPW=Entire range of pairwise comparison '
                        'p-values. P=Empirical p-value from '
                        'permutations. You can enter as many '
                        'correction criteria as you would like. These '
                        'will be associated with the p_value_cutoffs '
                        'you enter. For example "-c I EPW -p 0.1 0.05" '
                        'will apply the following cutoffs: Naive '
                        'p-value must be lower than 0.1 AND the entire '
                        'range of pairwise comparison values are below '
                        '0.05 for this gene. Note that the empirical '
                        'p-values should be interpreted at both tails. '
                        'Therefore, running "-c P -p 0.05" will apply '
                        'an alpha of 0.05 to the empirical (permuted) '
                        'p-values, i.e. it will filter everything '
                        'except the upper and lower 2.5 percent of the '
                        'distribution. '
                        'Default = Individual p-value. (I)',
                        choices=['I',
                                 'B',
                                 'BH',
                                 'PW',
                                 'EPW',
                                 'P'],
                        nargs='*',
                        default=['I'])
    parser.add_argument('-e','--permute',
                        help='Perform N number of permutations of the '
                        'significant results post-analysis. Each '
                        'permutation will do a label switching of the '
                        'phenotype and a new p-value is calculated '
                        'according to this new dataset. After all N '
                        'permutations are completed, the results are '
                        'ordered in ascending order, and the '
                        'percentile of the original result in the '
                        'permuted p-value distribution is reported.',
                        type=int,
                        default=0)
    parser.add_argument('-m', '--max_hits',
                        help='Maximum number of hits to report. SCOARY '
                        'will only report the top max_hits results per '
                        'trait',
                        type=int)
    parser.add_argument('-r', '--restrict_to',
                        help='Use if you only want to analyze a subset '
                        'of your strains. Scoary will read the '
                        'provided comma-separated table of strains and '
                        'restrict analyzes to these.')
    parser.add_argument('-w', '--write_reduced',
                        help='Use with -r if you want Scoary to create '
                        'a new gene presence absence file from your '
                        'reduced set of isolates. Note: Columns 1-14 '
                        '(No. sequences, Avg group size nuc etc) in '
                        'this file do not reflect the reduced dataset. '
                        'These are taken from the full dataset.',
                        default=False,
                        action='store_true')
    parser.add_argument('-s', '--start_col',
                        help='On which column in the gene '
                        'presence/absence file do individual strain '
                        'info start. Default=15. (1-based indexing)',
                        default=15,
                        type=int)
    parser.add_argument('-u', '--upgma_tree',
                        help='This flag will cause Scoary to write the '
                        'calculated UPGMA tree to a newick file',
                        default=False,
                        action='store_true')
    parser.add_argument('-n', '--newicktree',
                        help='Supply a custom tree (Newick format) for '
                        'phylogenetic analyses instead instead of '
                        'calculating it internally.',
                        default=None)
    parser.add_argument('--no_pairwise',
                        help='Do not perform pairwise comparisons. In'
                        'this mode, Scoary will perform population '
                        'structure-naive calculations only. (Fishers '
                        'test, ORs etc). Useful for summary operations '
                        'and exploring sets. (Genes unique in groups, '
                        'intersections etc) but not causal analyses.',
                        default=False,
                        action='store_true')
    parser.add_argument('--delimiter',
                        help='The delimiter between cells in the gene '
                        'presence/absence and trait files, as well as '
                        'the output file. ',
                        default=',',
                        type=str)
    parser.add_argument('--collapse',
                        help='Add this to collapse correlated genes '
                        '(genes that have identical distribution '
                        'patterns in the sample) into merged units. ',
                        default=False,
                        action='store_true')
    parser.add_argument('--threads',
                        help='Number of threads to use. Default = 1',
                        type=int,
                        default=1)
    parser.add_argument('--no-time',
                        help='Output file in the form '
                        'TRAIT.results.csv, instead of '
                        'TRAIT_TIMESTAMP.csv. When used with the -w '
                        'argument will output a reduced gene matrix in '
                        'the form gene_presence_absence_reduced.csv '
                        'rather than '
                        'gene_presence_absence_reduced_TIMESTAMP.csv ',
                        default=False,
                        action='store_true')
    parser.add_argument('--test',
                        help='Run Scoary on the test set in '
                        'exampledata, overriding all other parameters.',
                        default=False,
                        action='store_true')
    parser.add_argument('--citation',
                        help='Show citation information, and exit. ',
                        default=False,
                        action='store_true')
    parser.add_argument('--version', help='Display Scoary version, and '
                        'exit.',
                        action='version',
                        version=SCOARY_VERSION)

    args = parser.parse_args()
    
    if len(args.p_value_cutoff) == 1:
        cutoffs = {c : args.p_value_cutoff[0] for c in args.correction}
    else:
        cutoffs = dict(list(zip(args.correction, args.p_value_cutoff)))
        
    return args, cutoffs

########
# MAIN #
########
if __name__ == '__main__':
    pass
