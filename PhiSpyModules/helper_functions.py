# -*- coding: utf-8 -*-
import os
import sys
import argparse
import re
import pkg_resources
from argparse import RawTextHelpFormatter
from argparse import ArgumentTypeError as err
from .pathtype import PathType


def print_list():
    f = None
    try:
        # with pip we use resource streams that may be files or from archives
        f = pkg_resources.resource_stream('PhiSpyModules', 'data/trainingGenome_list.txt')
    except:
        sys.stderr.write('cannot find list')
        sys.exit(-1)
    for line in f:
        line = line.decode().strip()
        temp = re.split('\t', line)
        if int(temp[3]) == 1:
            print("{}\t{}".format(temp[2], 'data/' + temp[1]))
    f.close()

def is_valid_file(x):
    if not x or not os.path.exists(x):
        raise argparse.ArgumentTypeError("Checking for validity: {0} does not exist".format(x))
    return x


def get_args():
    usage = 'python3 PhiSpy.py [-opt1, [-opt2, ...]] infile'
    parser = argparse.ArgumentParser(
        description="phiSpy is a program for identifying prophages from among microbial genome sequences",
        epilog="(c) 2008-2018 Sajia Akhter, Katelyn McNair, Przemysław Decewicz, Rob Edwards, San Diego State University, San Diego, CA")
    parser.add_argument('infile', type=is_valid_file, help='Input file in genbank format', nargs='?')
    parser.add_argument('-m', '--make_training_data', type=str,
                             help='Create training data from a set of annotated genome files. Requires is_phage=1 qualifier in prophage\'s CDSs')
    parser.add_argument('-t', '--training_set', action='store', default='data/trainSet_genericAll.txt',
                             help='Choose the most closely related set to your genome.')
    parser.add_argument('-l', '--list', action='store_true', default=False,
                             help='List the available training sets and exit')
    #parser.add_argument('-c', '--choose', type=bool, default=False, const=True, nargs='?',
    #                         help='Choose a training set from a list (overrides -t)')
    parser.add_argument('-e', '--evaluate', type=bool, default=False, const=True, nargs='?',
                             help='Run in evaluation mode -- does not generate new data, but reruns the evaluation')
    parser.add_argument('-n', '--number', default=5, type=int,
                             help='Number of consecutive genes in a region of window size that must be prophage genes to be called. [Default: 5]')
    parser.add_argument('-u', '--min_contig_size', default=5000, type=int,
                        help='Minimum contig size (in bp) to be included in the analysis. Smaller contigs will be dropped. [Default: 5000]')
    parser.add_argument('-w', '--window_size', default=30, type=int,
                             help='Window size of consecutive genes to look through to find phages. [Default: 30]')
    parser.add_argument('-g', '--nonprophage_genegaps', default=10, type=int,
                             help='The number of non phage genes betweeen prophages. [Default: 10]')
    parser.add_argument('-r', '--randomforest_trees', default=500, type=int,
                             help='Number of trees generated by Random Forest classifier. [Default: 500]')
    parser.add_argument('--expand_slope', action='store_true', default=False, 
                             help='Use the product of the slope of the Shannon scores in making test sets')
    parser.add_argument('--kmers_type', default='all', choices=['all', 'codon', 'simple'], type=str,
                             help='Type of kmers used for calculating Shannon scores. [Default: all]')
    #parser.add_argument('-i', '--input_dir', help='The input directory that holds the genome')
    parser.add_argument('-o', '--output_dir', help='The output directory to write the results')
    parser.add_argument('-qt', '--quiet', type=bool, default=False, const=True, nargs='?',
                             help='Run in quiet mode')
    parser.add_argument('-k', '--keep', type=bool, default=False, const=True, nargs='?',
                             help='Do not delete temp files')
    parser.add_argument('-v', '--version', type=bool, default=False, const=True, nargs='?', help="Print the version and exit")
    return parser.parse_args()
