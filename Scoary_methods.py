#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Script for associating Roary output with phenotypic trait (+/-)
# Idea and implementation: Ola Brynildsrud (olbb@fhi.no)

# Structure of results dictionary:
# Results = {
# 			trait_1 :	{
# 						gene_1 :	{
# 									Annotation: Something
# 									OtherAnnotation : Something
# 									Score_1: a
# 									Score_2: b
# 									p_value: c
# 									p_corr: d
# 									}
# 						gene_2 :	{...}
# 						}
# 			trait_2 : {...}
# 			}

# Null hypothesis: Gene is equally distributed in Trait+ and Trait- isolates


import argparse, os, sys, csv, time
from scipy import stats as ss

SCOARY_VERSION = 'v1.1.1'

def main():
	# Parse arguments.
	parser = argparse.ArgumentParser(description='Scoary version %s - Screen pan-genome for trait-associated genes' %(SCOARY_VERSION), epilog='by Ola Brynildsrud (olbb@fhi.no)')
	parser.add_argument('-t', '--traits', help='Input trait table (comma-separated-values). Trait presence is indicated by 1, trait absence by 0. Assumes strain names in the first column and trait names in the first row' )
	parser.add_argument('-g', '--genes', help='Input gene presence/absence table (comma-separated-values) from ROARY. Strain names must be equal to those in the trait table')
	parser.add_argument('-p', '--p_value_cutoff', help='P-value cut-off. SCOARY will not report genes with higher p-values than this. Set to 1.0 to report all genes. Default = 0.05', default=0.05, type=float)
	parser.add_argument('-c', '--correction', help='Instead of cutting off at the individual test p-value (option -p), use the indicated corrected p-value for cut-off. Default = use individual test p-value.', choices=['Individual', 'Bonferroni', 'Benjamini-Hochberg'], default='Individual')
	parser.add_argument('-m', '--max_hits', help='Maximum number of hits to report. SCOARY will only report the top max_hits results per trait', type=int)
	parser.add_argument('-r', '--restrict_to', help='Use if you only want to analyze a subset of your strains. SCOARY will read the provided comma-separated table of strains and restrict analyzes to these.')
	parser.add_argument('-s', '--start_col', help='On which column in the gene presence/absence file do individual strain info start. Default=15. (1-based indexing)', default=15, type=int)
	parser.add_argument('--delimiter', help='The delimiter between cells in the gene presence/absence and trait files. NOTE: Even though commas are the default they might mess with the annotation column, and it is therefore recommended to save your files using semicolon or tab ("\t") instead. SCOARY will output files delimited by semicolon', default=',', type=str)
	parser.add_argument('--version', help='Display Scoary version, and exit.', default=False,action='store_true')

	args = parser.parse_args()

	# add version for tracing facilities
	if args.version:
		print sys.stdout, "Scoary version: %s" %(SCOARY_VERSION)
		sys.exit(0)
	# check required arguments.
	if args.traits is None or args.genes is None :
		parser.print_usage()
		if args.traits is None :
			print "error: argument -t/--traits is required"
		if args.genes is None:
			print "error: argument -g/--genes is required"
		sys.exit(1)

	with open(args.genes, "rU") as genes, open(args.traits, "rU") as traits:

		if args.restrict_to is not None:
			with open(args.restrict_to, "rU") as restrict_to_isolates:
				allowed_isolates = [isolate for isolate in restrict_to_isolates.next().rstrip().split(",")]
		else:
			allowed_isolates = None # Despite the confusing name, this actually means all isolates are allowed and included in the analysis

		genedic = Csv_to_dic_Roary(genes, args.delimiter, startcol=args.start_col - 1, allowed_isolates=allowed_isolates) # 0-based indexing

		traitsdic = Csv_to_dic(traits, args.delimiter, allowed_isolates)

		print "Finished loading files into memory."
		print "Tallying genes and performing statistical analyses"

		RES = Setup_results(genedic, traitsdic)


		StoreResults(RES, args.max_hits, args.p_value_cutoff, args.correction)


	sys.exit("Finished")



def Csv_to_dic_Roary(genefile, delimiter, startcol=0, allowed_isolates=None):
	r = {}
	csvfile = csv.reader(genefile, skipinitialspace=True)
	header = csvfile.next()

	# If roaryfile = True, then commence reading from column 15 (where the strains start)
	# Code to find out where strains begin. Use if Roary starts including more columns:
	#for x in header:
	#	if x not in ['Gene', 'Non-unique Gene name', 'Annotation', 'No. isolates', 'No. sequences', 'Avg sequences per isolate', 'Genome Fragment', 'Order within Fragment', 'Accessory Fragment', 'Accessory Order with Fragment', 'QC', 'Min group size nuc', 'Max group size nuc', 'Avg group size nuc']:
	#		print x
	#		break

	# This might be expanded if Roary changes its output format. For now:
	roaryfile = True

	strains = header[startcol:]

	try:
		genecol = header.index("Gene")
		nugcol = header.index("Non-unique Gene name")
		anncol = header.index("Annotation")
	except ValueError:
		print "Warning: Could not properly detect the correct names for all columns in the ROARY table."
		genecol = 0
		nugcol = 1
		anncol = 2

	for line in csvfile:
		q = line
		r[q[genecol]] = {"Non-unique Gene name" : q[nugcol], "Annotation" : q[anncol]} if roaryfile else {}

		for strain in xrange(len(strains)):
			if (allowed_isolates is not None) and strains[strain] not in allowed_isolates:
				continue
			if q[startcol + strain] in ["", "0", "-"]:
				# If the gene is not present, AND The isolate is allowed
				r[q[genecol]][strains[strain]] = 0
				# Add a 0 to indicate non-presence
			else:
				# Gene is present if any other value than "", "0" or "-" is in the cell
				r[q[genecol]][strains[strain]] = 1
				# Add a 1 to indicate presence of the current gene in this strain

	return r

def Csv_to_dic(csvfile, delimiter, allowed_isolates):
	tab = zip(*csv.reader(csvfile, delimiter=delimiter))
	r = {}
	if len(tab) < 2:
		sys.exit("Please check that your traits file is formatted properly and contain at least one trait")
	for num_trait in xrange(1,len(tab)):
		p = dict(zip(tab[0], tab[num_trait]))
		if "" in p:
			name_trait = p[""]
			del p[""]
		elif "Name" in p:
			name_trait = p["name"]
			del p["Name"]
		else:
			sys.exit("Make sure the top-left cell in the traits file is either empty or 'Name'. Do not include empty rows")

		# Filter so that only allowed isolates are included
		if allowed_isolates is not None:
			p = {strain: indicator for (strain,indicator) in p.iteritems() if strain in allowed_isolates} # if allowed_isolates is not None
		r[name_trait] = p

	return r

def Setup_results(genedic, traitsdic):

	# Need to create one results dictionary for each trait

	all_traits = {}
	for trait in traitsdic:
		all_traits[trait] = {}
		# Also need a list of all p-values to perform stepwise FDR methods
		p_value_list = []
		print "Gene-wise counting and Fisher's exact tests for trait: " + trait

		# We also need a number of tests variable for each genedic. (The number of tests may not be the same, depending on how many genes are in 0 or all strains.)
		number_of_tests = len(genedic)
		initial_number_of_tests = number_of_tests # This number needed for status bar

		i = 0 # Progress
		for gene in genedic:
			# Status:
			sys.stdout.write("\r{:.2%}".format(float(i)/initial_number_of_tests))
			sys.stdout.flush()
			i += 1 # Progress

			stat_table = Perform_statistics(traitsdic[trait], genedic[gene])
			num_pos = stat_table["tpgp"] + stat_table["tpgn"]
			num_neg = stat_table["tngp"] + stat_table["tngn"]

			if (stat_table["tpgp"] + stat_table["tngp"]) == 0:
				number_of_tests -= 1 # Remove 1 from the number of tests
				continue # proceed to the next gene

			if (stat_table["tpgn"] + stat_table["tngn"]) == 0:
				number_of_tests -= 1
				continue

			obs_table = [[stat_table["tpgp"], stat_table["tpgn"]],[stat_table["tngp"],stat_table["tngn"]]]
			fisher = ss.fisher_exact(obs_table)
			p_value = fisher[1]
			bonferroni_p = p_value * number_of_tests if (p_value * number_of_tests) < 1.0 else 1.0
			p_value_list.append((gene, p_value))

			all_traits[trait][gene] = { \
			"NUGN": genedic[gene]["Non-unique Gene name"], \
			"Annotation": genedic[gene]["Annotation"], \
			"tpgp": stat_table["tpgp"], \
			"tngp": stat_table["tngp"], \
			"tpgn": stat_table["tpgn"], \
			"tngn": stat_table["tngn"], \
			"sens": (float(stat_table["tpgp"]) / num_pos * 100) if num_pos > 0 else 0.0, \
			"spes": (float(stat_table["tngn"]) / num_neg * 100) if num_neg > 0 else 0.0, \
			"OR": fisher[0], \
			"p_v": p_value, \
			"B_p": bonferroni_p, \
			}

		print "\nAdding p-values adjusted for testing multiple hypotheses"
		# Now calculate Benjamini-Hochberg p-values
		sorted_p_values = sorted(p_value_list, key=lambda x: x[1]) # Sorted list of tuples: (gene, p-value)
		# Find out which p-values are ties
		# Note: Changed from step-down to step-up (from least significant to most significant)
		tie = [ sorted_p_values[i-1][1] == sorted_p_values[i][1] for i in xrange(1,len(sorted_p_values)) ][::-1]
		bh_corrected_p_values = {}
		bh_corrected_p_values[sorted_p_values[len(sorted_p_values)-1][0]] = last_bh = sorted_p_values[len(sorted_p_values)-1][1]
		for (ind, (gene,p)) in enumerate(sorted_p_values[::-1][1:]):
			bh_corrected_p_values[gene] = min([last_bh, p*number_of_tests/(ind+1.0)]) if not tie[ind] else last_bh
			last_bh = bh_corrected_p_values[gene]

		# Now add values to dictionaries:
		for gene in genedic:
			if gene in all_traits[trait]:
				all_traits[trait][gene]["BH_p"] = bh_corrected_p_values[gene] if bh_corrected_p_values[gene] < 1.0 else 1.0

	return all_traits


def Perform_statistics(traits, genes):
	r = {"tpgp" : 0, "tpgn": 0, "tngp" : 0, "tngn" : 0} # tpgn = trait positive, gene negative
	for t in traits: # For each strain
		try:
			if int(traits[t]) == 1 and genes[t] == 1:
				r["tpgp"] += 1
			elif int(traits[t]) == 1 and genes[t] == 0:
				r["tpgn"] += 1
			elif int(traits[t]) == 0 and genes[t] == 1:
				r["tngp"] += 1
			elif int(traits[t]) == 0 and genes[t] == 0:
				r["tngn"] += 1
			else:
				sys.exit("There was a problem with comparing your traits and gene presence/absence files. Make sure you have formatted the traits file to specification and only use 1s and 0s, and make sure the Roary file contains empty cells for non-present genes and non-empty text cells for present genes")
		except KeyError:
			print "\nError occured when trying to find " + str(t) + " in the genes file."
			sys.exit("Make sure strains are named the same in your traits file as in your gene presence/absence file")
	return r

def StoreResults(Results, max_hits, p_cutoff, correctionmethod):
	for Trait in Results:
		print "Storing results: " + Trait
		StoreTraitResult(Results[Trait], Trait, max_hits, p_cutoff, correctionmethod)


def StoreTraitResult(Trait, Traitname, max_hits, p_cutoff, correctionmethod): # This method is passed a single trait
	with open(Traitname + time.strftime("_%d_%m_%Y_%H%M") + ".csv", "w") as outfile:
		# Sort genes by p-value.
		sort_instructions = SortResultsAndSetKey(Trait)

		num_results = max_hits if max_hits is not None else len(Trait)

		cut_possibilities = {"Individual": "p_v", "Bonferroni": "B_p", "Benjamini-Hochberg": "BH_p"}

		outfile.write("Gene;Non-unique gene name;Annotation;Number_pos_present_in;Number_neg_present_in;Number_pos_not_present_in;Number_neg_not_present_in;Sensitivity;Specificity;Odds_ratio;p_value;Bonferroni_p;Benjamini_H_p\n")

		for x in xrange(num_results):
			# Start with lowest p-value, the one which has key 0 in sort_instructions
			currentgene = sort_instructions[x]
			if (Trait[currentgene][cut_possibilities[correctionmethod]] > p_cutoff):
				break
			outfile.write(currentgene + ";" + str(Trait[currentgene]["NUGN"]) + ";" + str(Trait[currentgene]["Annotation"]) + ";" + str(Trait[currentgene]["tpgp"]) + ";" + str(Trait[currentgene]["tngp"]) + ";" + str(Trait[currentgene]["tpgn"]) + ";" + str(Trait[currentgene]["tngn"]) + ";" + str(Trait[currentgene]["sens"]) + ";" + str(Trait[currentgene]["spes"]) + ";" + str(Trait[currentgene]["OR"]) + ";" + str(Trait[currentgene]["p_v"]) + ";" + str(Trait[currentgene]["B_p"]) + ";" + str(Trait[currentgene]["BH_p"]) + "\n")

def SortResultsAndSetKey(genedic): # This returns a dictionary where genes are sorted by p_value.
	return {i:gene for (i,gene) in enumerate(sorted(genedic, key=lambda x: genedic[x]["p_v"])) }

