import argparse
import csv
import glob
from collections import defaultdict
import random
parser = argparse.ArgumentParser()
parser.add_argument('--genes', help='File with list of genes to plot. Single column')
parser.add_argument('--dir', help='Path to PFF output directory')
parser.add_argument('--gene_sep', help='Separator character for gene IDs. Example GCA_00000001.1_00010, separator is _. Default for Prokka is _')
args = parser.parse_args()

genes_of_interest = []
for row in open(args.genes).readlines():
    genes_of_interest.append(row.strip())

f_list = glob.glob(args.dir+'/*.tsv')

d = defaultdict(dict)
for f in f_list:
   with open(args.dir+'/'+f) as res:
       for row in csv.reader(res, delimiter='\t'):
           for g in genes_of_interest:
               if g in row[2]:
                   genome = row[0].split(args.gene_sep)
                   genome = genome[:-1]
                   genome = '_'.join(genome)
                   d[genome][g] = 1

print('DATASET_BINARY')
print('SEPARATOR COMMA')
print('DATASET_LABEL,PFF Data')
print('COLOR,#ff0000')
print('FIELD_SHAPES'+',1'*len(genes_of_interest))
print('FIELD_LABELS,'+','.join(genes_of_interest))


colors = []
for g in genes_of_interest:
    colors.append("#{:06x}".format(random.randint(0, 0xFFFFFF)))

print('FIELD_COLORS,'+','.join(colors))
print('DATA')

for genome, genes in d.items():
    plot_vals = [str(genes.get(g, -1)) for g in genes_of_interest]
    print(genome+','+','.join(plot_vals))


