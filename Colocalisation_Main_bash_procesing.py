import sys
import os
import numpy as np
import pandas as pd
#import vcf as vcf
#import h5py
#import ast
import argparse
#from Colocalisation_functions import *

'''loading parameters from the parameter file'''
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--p',  metavar='p',help='...')
parser.add_argument('-trait',  metavar='t',help='...')

args = parser.parse_args()
print("reading data from :" +args.p)


parameters=pd.read_table(args.p,sep='\t',index_col=0 )
parameters.loc['sum_stats2_file','value']=args.trait
'''get gene list'''

get_gene_list='awk '+'\'NR > 1 {print $1}\' '+ parameters.loc['genes_file','value']+" > "+ parameters.loc['genes_file','value']+"_gene_list"
os.system(get_gene_list)

'''download data'''

os.system("cd "+parameters.loc['sum_stats2_folder','value']+"\
;wget https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/all_snp_gene_associations/"+parameters.loc['sum_stats2_file','value'])

'''grep genes from gwas'''

file=parameters.loc['sum_stats2_folder','value']+parameters.loc['sum_stats2_file','value']
os.system("gzip -cd "+ file+"|head -1 > "+file+"_filtered_genes")

grep_genes_from_gwas="zgrep -F -f "+ parameters.loc['genes_file','value']+"_gene_list "+ file+" >> "+file+"_filtered_genes"
os.system(grep_genes_from_gwas)


'''remove download data'''
os.system("rm "+parameters.loc['sum_stats2_folder','value']+parameters.loc['sum_stats2_file','value'])

