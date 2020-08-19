import sys
import numpy as np
import pandas as pd
from scripts.Create_figures_pQTL_functions import *

args=pd.read_table('./hipsci_proteomics_parameter_file.tsv',sep='\t',index_col=0)


folder_data_raw=args.loc['path_data_h5_files'].iloc[0]
path_data=args.loc['path_data'].iloc[0]
protein_h5_file=args.loc['protein_h5_file'].iloc[0]
pQTL_file=args.loc['pQTL_file'].iloc[0]
peptide_h5_file=args.loc['peptide_h5_file'].iloc[0]
path_out=args.loc['data_output_figures_folder'].iloc[0]

'''
change write flag if no data should be written
- write data: computes the peptide agreement with the sign of the pQTL
- plot_data: creates figures
'''
write_flag="write_data";
if write_flag=="write_data":
    for it,trait in enumerate(['eQTL']):
        if trait=='pQTL':
            proh5file=folder_data_raw+protein_h5_file
        elif trait=='pepQTL':
            proh5file=folder_data_raw+peptide_h5_file
        qtl=pd.read_table(path_data+pQTL_file,sep="\t")

        qtl1=write_peptide_agreement(qtl,folder_data_raw=folder_data_raw,\
                                     peph5file=folder_data_raw+peptide_h5_file,\
                                     proh5file=proh5file)
        qtl1.to_csv(path_data+pQTL_file+"agreement_peptides_all.txt", sep='\t')
#
if write_flag=="plot_data":
    qtl1=pd.read_table(path_data+pQTL_file,sep="\t")
    qtl1['Number of peptides (log2)']=qtl1['npeptidesqtl']
    qtl1['Fraction peptides in agreement']=qtl1['pepagreeqtl']/qtl1['npeptidesqtl']
    qtl1['agr']=qtl1['pepagreeqtl']/qtl1['npeptidesqtl']
    
    for i, ii in enumerate(['replicatedet and inframe_variant','not replicatedet and inframe_variant']):
        plt.figure(figsize=(4,4))
        ax1=plt.subplot(1,1,1)
        ax1.spines['top'].set_visible(False);    ax1.spines['right'].set_visible(False);a.yaxis.set_ticks_position('left'); a.xaxis.set_ticks_position('bottom');a.set_axis_bgcolor('white');
    
        qtl2=qtl1.query(ii)
        qtl2['abs _eff_size']=abs(qtl2['protein_eff_size'])
        qtl2['logqv']=-np.log10(qtl2['qv'])
        qtl2[['Fraction peptides in agreement','qv','protein_eff_size','abs _eff_size']].corr(method='spearman')
        qtl2[['gene_name','Fraction peptides in agreement','qv','protein_eff_size']].corr(method='spearman')
        qtl2['gene_name']
    
        x=np.log2(qtl1['Number of peptides (log2)']+1);y=qtl1['Fraction peptides in agreement']
        x2= np.log2(qtl2['Number of peptides (log2)']+1);y2=qtl2['Fraction peptides in agreement']
        plt.plot((0.95,max(x.max(),x2.max())),(0.5,0.5),'k--',lw=0.5)
        rrange=np.arange(1,2**max(x.max(),x2.max()))
        ax1.fill_between(np.log2(rrange+1), scst.binom.isf(0.95,p=0.5,n=rrange)/rrange, scst.binom.isf(1-0.95,p=0.5,n=rrange)/rrange, 
                 facecolor='grey', interpolate=True,alpha=0.5)
        plt.scatter(x2, y2, s=20,color='red', alpha=1,label='linked to inframe variants')
        plt.ylim(-0.1,1.1)    
        plt.xticks(np.arange(1,x.max()),2**(np.arange(1,x.max())-1))
        plt.ylabel('Fraction peptides\nin agreement',fontsize=16)
        plt.xlabel('Number of peptides assessed',fontsize=16)
        plt.tight_layout()
         
        plt.savefig(path_out+"QTL_VE_peptide"+str(ii)+"new.png",dpi=600)
        plt.savefig(path_out+"QTL_VE_peptide"+str(ii)+"new.svg",dpi=600)
