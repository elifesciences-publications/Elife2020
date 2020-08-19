import sys
import os
import numpy as np
import pandas as pd
#import vcf as vcf
#import h5py
import ast
import argparse
from Colocalisation_functions import *

'''loading parameters from the parameter file'''
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--p',  metavar='p',help='...')
parser.add_argument('-trait',  metavar='t',help='...')
parser.add_argument('-run_type',  metavar='r',help='...')

args = parser.parse_args()
print("reading data from :" +args.p)

parameters=pd.read_table(args.p,sep='\t',index_col=0 )
if args.trait is not None: parameters.loc['sum_stats2_file','value']=args.trait

try:os.mkdir(parameters.loc['data_output_folder','value'])
except: 1
try:os.mkdir(parameters.loc['data_output_figures_folder','value'])
except: 1

gene_positions=pd.read_table(parameters.loc['genes_file','value'],sep='\t').set_index('ensembl_gene_id',drop=0)
genes=gene_positions.query(' not replicatedet and inframe_variant')[['gene_name']].index

if args.run_type=="overlap":
    get_summary_stat_data (parameters,bed_flag=False)
elif args.run_type=="overlap_ld":
    get_summary_stat_data_ld(parameters, genes_traits=df)
elif args.run_type=="figures":
    parameters.loc['sum_stats2_file','value']=qtl1.loc[gene,'colocalization'].split(';')[0]
    parameters.loc['sum_stats2_folder','value']='/Users/mirauta/Results/hipsci/QTL_coloc_gwas/gene_coloc'
    plot_manhattan_scatter_overlap(parameters,gene_positions, genes,threshold_replication=-np.log10(0.01),field_manhattan='logp_value',field_scatter='beta')
 
elif args.run_type=="ecaviar":
    _doc="see eCAVIAR requirements for input; Rename   "

    files=np.array(glob.glob(parameters.loc['data_output_folder','value']+'*'));files=files[['_zscore_gwas' in f for f in files]]
    files=np.array([f.replace('_zscore_gwas','*')  for f in files])
    
    run_ecaviar(path_files="/",files=files[:2],ecaviar_path=parameters.loc['ecaviar_path','value'],name_qtl_file='stats1')        
