#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = 'Venkatramanan Krishnamani'
__copyright__ = 'Copyright 2014'
__version__ = '0.1'

#  The MIT License (MIT)
# =======================
#
# The source code in this file is copyrighted, but you are
# free to use and copy it as long as you don't change or
# remove any of the copyright notices.
#
# -----------------------------------------------------------------------------------
# Pairwise Protein-Interaction Calculator
# Copyright (C) 2013 by Venkatramanan Krishnamani <venky.krishna@me.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this
# software and associated documentation files (the "Software"), to deal in the Software
# without restriction, including without limitation the rights to use, copy, modify,
# merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to the following
# conditions:
#
# The above copyright notice and this permission notice shall be included in all copies
# or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
# PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.

import sys
import os
import re
try:
    import numpy as np
    import matplotlib as mpl
    mpl.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
except Exception, err:
    print "Please install a required python packages (matplotlib, numpy)"
    sys.stderr.write('*** IMPORT ERROR: %s\n' % str(err))
    sys.exit()

mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['axes.grid'] = True
mpl.rcParams['grid.color'] = 'gray'
mpl.rcParams['axes.labelsize'] = 'x-large'
mpl.rcParams['xtick.labelsize'] = 'x-large'
mpl.rcParams['ytick.labelsize'] = 'x-large'
mpl.rcParams['xtick.major.pad'] = 10
mpl.rcParams['ytick.major.pad'] = 10
mpl.rcParams['axes.titlesize'] = 'x-large'

window_size_i = 18
source_sequence_fi = ''
destination_aa_codon_table_fi = ''
source_aa_codon_table_fi = ''
three_to_one = {'Cys':'C', 'His':'H', 'Ile':'I', 'Met':'M', 'Ser':'S', 'Val':'V', 'Ala':'A', 'Gly':'G', 'Pro':'P', 'Thr':'T', 'Arg':'R', 'Phe':'F', 'Tyr':'Y', 'Trp':'W', 'Asp':'D', 'Asn':'N', 'Glu':'E', 'Gln':'Q', 'Lys':'K', 'End':'X', 'Leu':'L'}

def print_centered(string):
    r = 0
    c = 0
    if sys.stdout.isatty():
        r, c = os.popen('stty size', 'r').read().split()
    else:
        r, c = (100,100)
    for n in range(((int(c)+1)-len(string))/2):
        sys.stdout.write(' ')
    print string

def read_arguments():
    '''
    Reads and Parses command line arguments
    '''
    import argparse
    global source_sequence_fi
    global destination_aa_codon_table_fi
    global source_aa_codon_table_fi
    global window_size_i

    parser = argparse.ArgumentParser(usage="\ncodon-optimize.py -f <sequence-for-optimization> -s <source-organism-codon-table> -o <expression-organism-codon-table> [-w <window-size>]", epilog="\n"+__author__+", "+__copyright__+"\n")
    parser.add_argument("-f", "--sequencefile", action="store", dest="source_sequence", default=False, help="<Mandatory> Input DNA/mRNA sequence (Format: FASTA))", type=str)
    parser.add_argument("-s", "--sourcetable", action="store", dest="source_codon_table", default=False, help="<Mandatory> Filename of codon usage table of source organism (Format: CodonFrequency output in GCG Wisconsin Package)", type=str)
    parser.add_argument("-o", "--expressiontable", action="store", dest="expression_codon_table", default=False, help="<Mandatory> Filename of codon usage table of expression organism (Format: CodonFrequency output in GCG Wisconsin Package)", type=str)
    parser.add_argument("-w", "--windowsize", action="store", dest="window_size", default='18', help='<Optional> Window size for MinMax analysis', type=int)
    args = parser.parse_args()

    mandatories = ['source_sequence', 'source_codon_table', 'expression_codon_table']
    for m in mandatories:
        if not args.__dict__[m]:
            parser.print_help()
            sys.exit("\n*** ERROR: Missing a mandatory argument. ***\n")
    window_size_i = args.window_size
    source_sequence_fi = args.source_sequence
    source_aa_codon_table_fi = args.source_codon_table
    destination_aa_codon_table_fi = args.expression_codon_table

def read_cdn_usage(filename):
    fHandle = open (filename, 'r')
    translation_table = {}
    codon_table = {}
    aminoacid_codon_table = {}
    #read basic information
    for line in fHandle.readlines():
        if not line == '\n':
            split = line.split()
            translation_table[split[1]] = three_to_one[split[0]]
            codon_table[split[1]] = float(split[3])
            try:
                aminoacid_codon_table[three_to_one[split[0]]].append({'codon':split[1], 'frequency':float(split[3])})
            except KeyError:
                aminoacid_codon_table[three_to_one[split[0]]] = [{'codon':split[1], 'frequency':float(split[3])}]
    return (aminoacid_codon_table, translation_table, codon_table)

def sortby_most_common_codon(aminoacid_codon_table):
    for key in aminoacid_codon_table.keys():
        aminoacid_codon_table[key].sort(key=lambda x: x['frequency'], reverse=True)
    return aminoacid_codon_table

def average_codon_usage(aminoacid_codon_table):
    average_codon_usage = {}
    for key in aminoacid_codon_table.keys():
        sum_of_frequency = 0
        count = 0
        for item in aminoacid_codon_table[key]:
            sum_of_frequency = sum_of_frequency + item['frequency']
            count += 1
        average_codon_usage[key] = sum_of_frequency/count
    return average_codon_usage

def read_sequence(sequence_file):
    sequenceFile = open(sequence_file, 'r')
    sequence = ''
    sequence_name = ''
    for line in sequenceFile.readlines():
        if re.match(r'^>', line):
            sequence_name = line.rstrip()[1:]
        else:
            sequence += line.rstrip()
    sequenceFile.close()
    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
    return (codons, sequence_name)

def translate_sequence(sequenceFile, translation_table):
    protein_sequence = []
    for codon in sequenceFile:
        protein_sequence.append(translation_table[codon])
    return protein_sequence

def getminmax(dna_sequence, aminoacid_codon_table, translation_table, codon_table):
    aminoacid_sequence = []
    actual_list = []
    max_list = []
    min_list = []
    average_list = []
    max_percent = []
    min_percent = []
    maxsorted_codon_table = sortby_most_common_codon(aminoacid_codon_table)
    average_codon_table = average_codon_usage(aminoacid_codon_table)
    for i in range(0,len(dna_sequence)-window_size_i):
        sum_actual_freq = 0
        sum_max_freq = 0
        sum_min_freq = 0
        sum_average_freq = 0
        for cdn in dna_sequence[i:i+window_size_i]:
            sum_actual_freq = sum_actual_freq + codon_table[cdn]
            sum_max_freq = sum_max_freq + maxsorted_codon_table[translation_table[cdn]][0]['frequency']
            sum_min_freq = sum_min_freq + maxsorted_codon_table[translation_table[cdn]][-1]['frequency']
            sum_average_freq = sum_average_freq + average_codon_table[translation_table[cdn]]
        actual_list.append(sum_actual_freq/window_size_i)
        max_list.append(sum_max_freq/window_size_i)
        min_list.append(sum_min_freq/window_size_i)
        average_list.append(sum_average_freq/window_size_i)
        if not len(actual_list) == len(max_list) == len(min_list) == len(average_list):
            print "Error the list sizes are different"
    for i in range(0, len(actual_list)):
        max_percent.append((actual_list[i]-average_list[i])*100/(max_list[i]-average_list[i]))
        min_percent.append((average_list[i]-actual_list[i])*100/(average_list[i]-min_list[i]))
    minmax = []
    for i in range(len(min_percent)):
        if min_percent[i] < 0:
            minmax.append(max_percent[i])
        else:
            minmax.append(-1*min_percent[i])
    return minmax

def get_colors (list):
    colors = []
    for l in list:
        if l < 0:
            colors.append('r')
        else:
            colors.append('b')
    return colors

def find_nearest(array,value):
    return array[min(range(len(array)), key=lambda i: abs(array[i]['frequency']-value))]

def optimize_codon_map(sequence, source_codon_table, source_translation_table, destination_aa_codon_table):
    optimized_sequence = []
    for codon in sequence:
        optimized_sequence.append(find_nearest(destination_aa_codon_table[source_translation_table[codon]], source_codon_table[codon])['codon'])
    return optimized_sequence

def write_fasta_file(filename, sequence, name, width):
    outfile = open(filename, 'w')
    outfile.write(">Expression Optimized %s\n" % name)
    previous_index = 0
    for i in range(len(sequence)/width + 1):
        outfile.write("%s\n" % ''.join(sequence[previous_index:previous_index+width]))
        previous_index = previous_index + width
    outfile.close()

# Plotting
def minmax_plot(data, filename, title, xlim, xlabel, ylabel, fontsize):
    plot0 = plt.figure("All", figsize=(12, 8))
    ax = plot0.add_subplot(111)
    ax.plot([1, 2])
    ax.bar(range(len(data)), data, 1, color=get_colors(data), edgecolor = "none")
    ax.locator_params(nbins=3)
    ax.set_xlabel(xlabel, fontsize=fontsize-2)
    ax.set_ylabel(ylabel, fontsize=fontsize-2)
    ax.set_xlim([0,xlim])
    ax.set_ylim([-100,100])
    ax.set_xticks([i for i in range(0,xlim+1,100)])
    ax.set_yticks([i for i in range(-100,100+1,20)])
    ax.set_title(title+"\n", fontsize=fontsize)
    pp = PdfPages(filename)
    pp.savefig(plot0, dpi=300)
    pp.close()
    plt.close(plot0)

if __name__ == '__main__':
    '''
    Main Function Call
    '''
    rows = 0
    columns = 0
    if sys.stdout.isatty():
        rows, columns = os.popen('stty size', 'r').read().split()
    else:
        rows, columns = (100,100)
    # Intro message
    print '\n'
    introMessage = "Codon Optimization - Version " + __version__ +"\n--------------------------------------------------------------\nOptimize sequence for expression by codon usage pattern of source organism.\n\n###\n"
    for s in introMessage.split('\n'):
        for n in range(((int(columns)+1)-len(s))/2):
            sys.stdout.write(' ')
        print s
    read_arguments()

print ">>> Optimizing " + source_sequence_fi + "..."
#Read Codon table of source and destination organisms
(source_aminoacid_codon_table, source_translation_table, source_codon_table) = read_cdn_usage(source_aa_codon_table_fi)
(destination_aminoacid_codon_table, destination_translation_table, destination_codon_table) = read_cdn_usage(destination_aa_codon_table_fi)

(original_sequence, original_sequence_name) = read_sequence(source_sequence_fi)
expression_optimized_sequence = optimize_codon_map(original_sequence, source_codon_table,
                                                     source_translation_table,
                                                     destination_aminoacid_codon_table)
original_sequence_source_minmax = getminmax(original_sequence,
                                 source_aminoacid_codon_table,
                                 source_translation_table,
                                 source_codon_table)
original_sequence_destination_minmax = getminmax(original_sequence,
                             destination_aminoacid_codon_table,
                             destination_translation_table,
                             destination_codon_table)
optimized_sequence_destination_minmax = getminmax(expression_optimized_sequence,
                                       destination_aminoacid_codon_table,
                                       destination_translation_table,
                                       destination_codon_table)

#### Translate Protein Sequence for Verification
original_protein_sequence = translate_sequence(original_sequence, source_translation_table)
optimized_protein_sequence = translate_sequence(expression_optimized_sequence, destination_translation_table)

write_fasta_file(source_sequence_fi[:-4]+'_optimized.seq', expression_optimized_sequence, original_sequence_name, 25)

xlim = max(len(original_sequence_source_minmax), len(original_sequence_destination_minmax), len(optimized_sequence_destination_minmax))
minmax_plot(original_sequence_source_minmax, source_sequence_fi[:-4]+'_source_minmax.pdf', original_sequence_name+' (Source expression)', xlim, '', '%MinMax', 18)
minmax_plot(original_sequence_destination_minmax, source_sequence_fi[:-4]+'_expression_minmax.pdf', original_sequence_name+' (Unoptimized expression)', xlim, '', '%MinMax', 18)
minmax_plot(optimized_sequence_destination_minmax, source_sequence_fi[:-4]+'_optimized_expression_minmax.pdf', original_sequence_name+' (Optimized expression)', xlim, '', '%MinMax', 18)
print ">>> Done"
