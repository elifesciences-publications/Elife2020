import sys
import numpy as np
import pandas as pd
from scripts.Create_figures_pQTL_functions import *

#sys.exit()
args=pd.read_table('../hipsci_proteomics_files.tsv',sep='\t',index_col=0)

traits=['eQTL','pQTL' ]
qtl={}
for trait in traits:
    qtl[trait]=pd.read_table(args.loc['path_data'].iloc[0]+args.loc[trait+'_file'].iloc[0],index_col=0,sep='\t')


'''Pannels for Main'''

'''Figure 2'''
plot_fig_2(qtl,path_out=args.loc['path_data'].iloc[0])

'''Figure 3'''
plot_fig3_panel_a(qtl,path_out=args.loc['path_data'].iloc[0])


'''Pannels for Supp figures ''' 

'''Comparing eQTL and pQTL effect sizes'''
plot_sup_fig_eqtl_pqtl_effect(qtl=qtl,path_out=args.loc['data_output_figures_folder'].iloc[0])

