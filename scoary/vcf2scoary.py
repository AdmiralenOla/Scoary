#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Script to search vcf files for mutations within specific coordinates
# Input: 
# -A vcf file
#
# Output:
# -A Roary-like file with mutations sorted in rows, strains as columns and presence/absence in cells
# -Columns: Chromosome, Position, variant (eg C->T), type (eg missense, synonymous, frameshift etc)


# Reading VCF
# File metainfo starts as ##key=value
# These are always formed and should be caught
# example ##fileformat=VCFv4.3 - give warning if format is off
# Columns 8 MANDATORY
# CHROM POS ID REF ALT QUAL FILTER INFO
# OPTIONAL COLUMNS
# FORMAT SAMPLE1 SAMPLE2 etc
# All data lines are tab-delimited
# CHROM : string, no whitespace
# POS : integer. Can have many lines with same pos. Pos=0 or N+1 for telomere positions
# ID : semicolon-delimited list of strings
# REF : string, ACGTN (can be multiple)
# ALT : comma-separated list, ACGTN* (* = allele is missing due to overlapping deletion)
# (NOTE: Suggest splitting ALT variants into different lines to preserve binarity)
# QUAL : float
# FILTER : PASS or semicolon-delimited list
# INFO : semicolon-delimited list of key=value pairs or flags
# FORMAT (optional) : colon-delimited list.
# Genotype fields - Genotype always first field
# GT encoded as allele values separated by | or /. 0 = reference. 1 = first ALT. 2 = second alt etc
# NOTE: Haploid calls (bacteria) have only 1 value
# NOTE: / means genotype unphased. | means genotype phased
# INFO field SVtypes : DELetion, INSertion, DUPlication, INVersion, CNV

import sys
import argparse
import os
import csv
import re
import traceback

__version__ = '0.1b'
__author__ = 'Ola Brynildsrud'
__credits = ['Ola Brynildsrud']
__email__ = 'olbb@fhi.no'

def main():
    """
    Converts VCF files (version 4.x) to Scoary format
    """
    ##########################################################################
    # Parse command line arguments

    parser = argparse.ArgumentParser(
        description='This script takes in vcf files and creates a '
                    'presence/absence matrix of mutations in the '
                    'Roary/Scoary format',
        epilog='by Ola Brynildsrud (olbb@fhi.no)')
    parser.add_argument(
        '--out',
        action='store',
        default='./mutations_presence_absence.csv',
        help='The path to the output file')
    parser.add_argument(
        '--types',
        action='store',
        default='ALL',
        help='The types of variants to include in the output. NOTE: This '
             'works if TYPE=XX can be found in the INFO column of the vcf '
             'file. The special keyword ALL includes all types. This is '
             'the default setting. Common types are snp, mnp, ins, del '
             'and complex. Give as comma-separated list. '
             'Example: --types snp,ins,del')
    parser.add_argument(
        '--version',
        action='version',
        version=__version__)
    parser.add_argument(
        '--force',
        action='store_true',
        default=False,
        help='Force overwriting of output file. (If it already '
             'exists)')
    parser.add_argument(
        'vcf',
        action='store',
        metavar='<VCF_file>',
        help='The VCF file to convert to Roary/Scoary format')

    args = parser.parse_args()
    if args.types is not "ALL":
        args.types = args.types.split(",")

    if os.path.isfile(args.out) and not args.force:
        sys.exit("Outfile already exists. Change name of outfile or "
                 "run with --force")
    if not os.path.isfile(args.vcf):
        sys.exit("Unable to locate input file %s" % args.vcf)

    with open(args.vcf,'rU') as vcffile, open(args.out,'w') as outfile:
        lines = csv.reader(vcffile, delimiter='\t', quotechar='"')
        metainfo = {"##INFO" : {},
                    "##FILTER" : {},
                    "##FORMAT" : {},
                    "##ALT" : {},
                    "##contig" : {},
                    "##META" : {},
                    "##SAMPLE" : {},
                    "##PEDIGREE" : {}
        }
        #for line in lines:
        while True:
            try:
                line = next(lines)
            except StopIteration:
                print(traceback.print_exc())
                sys.exit("ERROR: There appears to be only metainformation "
                         "(lines starting with ##) in your VCF file.")
            # Get metainfo from file
            if line[0][:2] == '##':
                infoline = re.split('=',line[0], maxsplit=1)
                # Capture list output for complex tags
                if infoline[0] in metainfo:
                    ID=re.search(r'ID=(\w+)',infoline[1]).group(1)
                    infolist = re.split(',(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)',infoline[1].strip("<>"))
                    metainfo[infoline[0]][ID] = {}
                    # Enter all elements in infolist into appropriate dic
                    for e in infolist:
                        esplit = e.split("=")
                        metainfo[infoline[0]][ID][esplit[0]] = esplit[1]

                else:
                    metainfo[infoline[0]] = infoline[1]
            else:
                # Have reached the data section of the file
                data = {"header": line}
                break

        try:
            vcfversion = metainfo["##fileformat"].split("v")[1]
            if int(vcfversion[0]) != 4:
                print("WARNING: A VCF format other than 4.x detected."
                      " File parsing may proceed with errors.")
            else:
                print("VCF version %s detected" % vcfversion)
        except:
            print("WARNING: Could not detect VCF format. Expected "
                  "v4.x. File parsing may proceed with errors.")
            print(traceback.print_exc())

        # Check that genotype fields have a single allele
        if metainfo["##FORMAT"]["GT"]["Number"] != "1":
            sys.exit("ERROR: Expected a single allele per genotype. Scoary "
                     "only works for haploid organisms.")

        # Have now caught all metainformation. Now get column information       
        #header = next(line)
        #print header
        data["header"] = data["header"][:9] + ["DUMMY"] + data["header"][9:]
        outfile.write(','.join('"' + c + '"' for c in data["header"]) + '\n')

        while True:
            try:
                line = next(lines)
            except StopIteration:
                print("Reached the end of the file")
                sys.exit(0)
            # Check if line is allowed:
            if args.types is not "ALL":
                vartype = re.search(r'TYPE=(\w+)',line[7]).group(1)
                if vartype not in args.types:
                    continue
            
            # Split line if ALT contains more than one variant
            if "," in line[4]:
                orgline = line[:]
                alts = line[4].split(",")
                c = 1
                for a in alts:
                     newline = orgline[:]
                     newline[4] = a
                     # Only get GT
                     newline[9:] = \
                         [cell.split(":")[0] for cell in orgline[9:]]
                     # Fix dummy comparisons
                     newline[9:] = fixdummy(newline[9:], c)
                     newline = newline[:9] + ["True"] + newline[9:]
                     c += 1
                     writeLine(newline, outfile)

            # Genotype fields need to be 0 or 1
            # GT is always first in colon-separated list
            else:
                newline = line[:9] + ["False"] + line[9:]
                writeLine(newline, outfile)

def writeLine(line, outfile):
    writeline = line[:9] + [cell.split(":")[0] for cell in line[9:]]
    outfile.write(','.join('"' + c + '"' for c in writeline) + '\n')

def fixdummy(line,c):
    newline = line[:]
    try:
        for x in range(len(line)):
            if line[x] == ".":
                # Missing data get entered as reference / no presence
                newline[x] = "0"
            elif int(line[x]) == c:
                newline[x] = "1"
            else:
                newline[x] = "0"
    except ValueError:
        print(newline, c)
        sys.exit(-1)
    return newline

########
# MAIN #
########
if __name__ == '__main__':
    main()
