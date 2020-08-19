import sys
from ast import *
import pandas as pd
import pandas
import numpy as np
import seaborn as sb
import scipy.stats as scst
import scipy  as sc
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from sklearn.linear_model import LinearRegression,LogisticRegression
import statsmodels.sandbox.stats.multicomp as scst2
import matplotlib.pyplot as plt
import re
from numpy.polynomial.polynomial import polyfit
#
#from scripts.postprocess_functions_genome_results import *
#from scripts.postprocess_functions_plots import *
colors=['red','darkgreen','green','blue','orange','magenta','salmon','yellowgreen','yellow']

def plot_fig_2(qtl,path_out=""):
    for fontsize in [16]:
        fig = plt.figure(figsize=(5, 4))
        fig.patch.set_facecolor('white')
        ax2 = fig.add_subplot(121)
        
        trait='pQTL';temp=qtl[trait]
        x=np.array([temp.shape[0],temp.shape[0],\
                    temp[(temp['mrna_nominal_p_value_bonf']<0.01)&(temp['protein_eff_size']*temp['mrna_eff_size']>0)].shape[0]])
        ax2.bar(1,x[2], color='steelblue',label="replicated")
        ax2.bar(1,x[1]-x[2],bottom=x[2], color='lightblue',label="specific")
        ax2.annotate("eQTL",  xy=(0.9, x[2]/2));temp=x[2]/2
        ax2.annotate("no eQTL",  xy=(0.85, x[1]/2+temp));
        plt.xticks([1],['pQTL discovery'],fontsize=fontsize)
        plt.ylim(0,1000)
        plt.ylabel('Number of genes',fontsize=fontsize)
        colors=['darkblue','lightblue','darkkhaki','olive']
        
        ax1 = fig.add_subplot(122)
        for ax in [ax1,ax2]:
            ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False); ax.yaxis.set_ticks_position('left'); 
        
        trait='eQTL';temp=qtl[trait]
        x={'al':temp.shape[0],\
            'rep1':temp[(temp['protein_nominal_p_value_bonf']<0.01)&(temp['mrna_eff_size']*temp['protein_eff_size']>0)].shape[0]}
        
        ax1.bar(3.5,x['rep1'], color='forestgreen')
        ax1.bar(3.5,x['al']-x['rep1'],bottom=x['rep1'], color='palegreen')
        ax1.annotate("pQTL",  xy=(3.32, x['rep1']/2));temp=x['rep1']/2
        ax1.annotate("no pQTL",  xy=(3.32, x['al']/2+temp));temp=x['al']/2
        plt.xticks([3.5],['eQTL discovery'],fontsize=fontsize)
        
        plt.tight_layout()
        
        plt.savefig(path_out+"Main_Fig2.png",dpi=600)
        plt.savefig(path_out+"Main_Fig2.svg",dpi=600)

def plot_fig3_panel_a(qtl,path_out=""):
    qtl1=qtl['pQTL']
    qtl1['replicated_isoform']=qtl1['transcript_nominal_p_value_bonf']<0.01
    qtl1['pav']=(qtl1['inframe']==qtl1['inframe'])|(qtl1['frameshift']==qtl1['frameshift'])
    qtl1['replicated_any']=qtl1['replicated']+(qtl1['replicated']==0)*qtl1['replicated_isoform']*2
    qtl1['replicated_any']=qtl1['replicated_any']+(qtl1['replicated_any']==0)*qtl1['pav']*3
    qtl1['logpv']=-np.log10(qtl1['empirical_feature_p_value'])
    
    
    plt.figure(figsize=(4,4))
    #a =   plt.subplot(1,1,1);
    #a.spines['top'].set_visible(False);    a.spines['right'].set_visible(False);a.yaxis.set_ticks_position('left'); a.xaxis.set_ticks_position('bottom');a.set_axis_bgcolor('white');
    
    my_vals1 = {"is eQTL":qtl1.query('replicated').shape[0]\
    ,"\nis tQTL\nonly": qtl1.query('not replicated and replicated_isoform').shape[0], \
    "no\nmRNA or PAV\nmechanisms": qtl1.query('not replicated and not replicated_isoform and not pav').shape[0],\
    "     in LD with PAV": qtl1.query('not replicated and not replicated_isoform and pav').shape[0]}
     
    my_pal = {"\nis tQTL\nonly": "royalblue","is eQTL": "steelblue","     in LD with PAV": "cadetblue"," ": "lightblue","no\nmRNA or PAV\nmechanisms":"lightblue"}
    
    def absolute_value(val):    return ""+str(np.array(list(my_vals1.keys()))[np.where(np.array(list(my_vals1.values()))==int(np.round(val/100.*qtl1.shape[0], 0)))[0]][0]) \
                                          +"\n("+str(int(np.round(val/100.*qtl1.shape[0], 0)))+")"
    def absolute_value(val):    return int(np.round(val/100.*qtl1.shape[0], 0))
    
    plt.pie(list(my_vals1.values()) , radius=1.1, colors=[my_pal[k] for k in my_vals1.keys()],    wedgeprops=dict(width=0.75, edgecolor='w'),\
           autopct=absolute_value ,        explode=(0, 0, 0.2,0.2),textprops={'fontsize': 10})#, labels=[k+"\n("+str(my_vals1[k])+")" for k in  my_vals1.keys()])
    
    #ax.set(aspect="equal", title='Pie plot with `ax.pie`')
    #plt.annotate(xy=(2-0.25,3.25),s="protein PV > 0.01",fontsize=16)
    plt.tight_layout()
    plt.savefig(path_out+'MAIN_fig2_pqtl_categories.png',dpi=600)
    plt.savefig(path_out+'MAIN_fig2_pqtl_categories.svg',dpi=600)
    
    
    my_pal = {-2: "royalblue",-3: "steelblue",-1: "cadetblue",0:"lightblue"}
    #
    plt.figure(figsize=(4.75,2.5))
    a =   plt.subplot(1,1,1);
    a.spines['top'].set_visible(False);
    a.spines['right'].set_visible(False);
    a.yaxis.set_ticks_position('left'); a.xaxis.set_ticks_position('bottom');a.set_axis_bgcolor('white');
    sb.boxplot(x=-qtl1['replicated_any'] ,y=qtl1['logpv'],palette=my_pal);plt.yscale('log')
    yy=[2,5,10,20,40]
    plt.yticks(yy,yy)
    plt.xticks([0,1,2,3],["is eQTL","is tQTL\nonly","LD with\nPAV","no eQTL,\nno PAV LD"],fontsize=14)
    plt.tight_layout()
    plt.ylabel("-log10 PV",fontsize=14)
    plt.xlabel("",fontsize=14)
    plt.savefig(path_out+'MAIN_fig2_pqtl_categories_insert.png',dpi=600)   
    plt.savefig(path_out+'MAIN_fig2_pqtl_categories_insert.svg',dpi=600)   
    


def plot_sup_fig_eqtl_pqtl_effect(qtl,path_out=""):
    colors={'eQTL':"forestgreen",'pQTL':"steelblue"}
    '''mRNA_protein_replication'''
    fig = plt.figure(figsize=(8,8));fig.patch.set_facecolor('white')
    plt.rcParams["legend.handletextpad"] = 0
    plt.rcParams["legend.fontsize"]= 11
    plt.rcParams["legend.frameon"] = False
    ax = fig.add_subplot(222)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False); ax.yaxis.set_ticks_position('left'); 
    plt.plot((qtl['pQTL']['protein_eff_size'].min(),qtl['pQTL']['protein_eff_size'].max()),(0,0),'k--',lw=0.5)
    plt.plot((0,0),(qtl['pQTL']['mrna_eff_size'].min(),qtl['pQTL']['mrna_eff_size'].max()),'k--',lw=0.5)
    
    
    x2=qtl['pQTL'].query('mrna_nominal_p_value_bonf > 0.01')['mrna_eff_size']
    x=qtl['pQTL'].query('mrna_nominal_p_value_bonf > 0.01')['protein_eff_size']
    x12=qtl['pQTL'].query('mrna_nominal_p_value_bonf < 0.01')['mrna_eff_size']
    x1=qtl['pQTL'].query('mrna_nominal_p_value_bonf < 0.01')['protein_eff_size']
    
    plt.plot(qtl['pQTL'].query('mrna_nominal_p_value_bonf > 0.01')['protein_eff_size'],qtl['pQTL'].\
             query('mrna_nominal_p_value_bonf > 0.01')['mrna_eff_size'],'o',color='lightblue',label=' PV > 0.01',markersize=3)
    plt.plot(qtl['pQTL'].query('mrna_nominal_p_value_bonf < 0.01')['protein_eff_size'],qtl['pQTL'].\
             query('mrna_nominal_p_value_bonf < 0.01')['mrna_eff_size'],'+',color='steelblue',label=' PV < 0.01')
    plt.legend(loc='best', bbox_to_anchor=(0.5, -0.1, 0.5, 0.5),title="eQTL",fontsize=11,markerfirst=False)
    
    b, m = polyfit((x-np.nanmean(x))/np.nanstd(x),(x2-np.nanmean(x2))/np.nanstd(x2), 1)   
    plt.plot(x, b + m * x, '-',color='lightblue',label=None)
    ax.annotate(r"$r^2$"+" "+str(np.around(m**2,2)),  xy=(-1.5,-0.2));    
    b, m = polyfit((x1-np.nanmean(x1))/np.nanstd(x1),(x12-np.nanmean(x12))/np.nanstd(x12), 1)
    plt.plot(x1, b + m * x1, '-',color='steelblue',label=None)
    ax.annotate(r"$r^2$"+" "+str(np.around(m**2,2)),  xy=(-1.5,-1.7));    
               
    plt.ylabel('pQTL effect size at RNA level',rotation=90,fontsize=14)
    plt.xlabel('pQTL effect size at protein level',fontsize=14)

    ax = fig.add_subplot(221)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False); ax.yaxis.set_ticks_position('left'); 
    plt.plot((qtl['eQTL']['mrna_eff_size'].min(),qtl['eQTL']['mrna_eff_size'].max()),(0,0),'k--',lw=0.5)
    plt.plot((0,0),(qtl['eQTL']['protein_eff_size'].min(),qtl['eQTL']['protein_eff_size'].max()),'k--',lw=0.5)
    
    x=qtl['eQTL'].query('protein_nominal_p_value_bonf > 0.01')['mrna_eff_size']
    x2=qtl['eQTL'].query('protein_nominal_p_value_bonf > 0.01')['protein_eff_size']
    x1=qtl['eQTL'].query('protein_nominal_p_value_bonf < 0.01')['mrna_eff_size']
    x12=qtl['eQTL'].query('protein_nominal_p_value_bonf < 0.01')['protein_eff_size']
    plt.plot(x,x2,'o',color='lightgreen', label=' PV > 0.01',markersize=3)
   
    
    plt.plot(qtl['eQTL'].query('protein_nominal_p_value_bonf < 0.01')['mrna_eff_size'],qtl['eQTL'].\
             query('protein_nominal_p_value_bonf < 0.01')['protein_eff_size'],'+',color='forestgreen',label=' PV < 0.01')
    plt.legend(loc='best', bbox_to_anchor=(0.5, 0.0, 0.5, 0.5),\
               title="pQTL",fontsize=11,markerfirst=False)
    b, m = polyfit((x-np.nanmean(x))/np.nanstd(x),(x2-np.nanmean(x2))/np.nanstd(x2), 1)   
    plt.plot(x, b + m * x, '-',color='lightgreen',label=None)
    ax.annotate(r"$r^2$"+" "+str(np.around(m**2,2)),  xy=(-1.9,-0.55));    
    b, m = polyfit((x1-np.nanmean(x1))/np.nanstd(x1),(x12-np.nanmean(x12))/np.nanstd(x12), 1)
    plt.plot(x1, b + m * x1, '-',color='forestgreen',label=None)
    ax.annotate(r"$r^2$"+" "+str(np.around(m**2,2)),  xy=(-1.4,-1.65));    
               
    plt.ylabel('eQTL effect size at RNA level',rotation=90,fontsize=14)
    plt.xlabel('eQTL effect size at protein level',fontsize=14)  
#    ax.annotate("0.1 qv",  xy=(-np.log10(0.05)+0.3,0.61));
    
    ax = fig.add_subplot(223)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False); ax.yaxis.set_ticks_position('left'); 
    temp=qtl['eQTL']
    thrs=np.linspace(1,20,26)
    plt.plot(thrs ,np.array([ np.corrcoef(temp['mrna_eff_size'][temp['qv']<thr], temp['protein_eff_size'][temp['qv']<=thr])[0,1] for thr in 10**-thrs ]))
#    plt.plot((-np.log10(0.1),-np.log10(0.1)),(0.575,np.corrcoef(temp['mrna_eff_size'][temp['qv']<0.1], temp['protein_eff_size'][temp['qv']<=0.1])[0,1]),'k--',lw=0.5)
    ax.spines['top'].set_visible(False);    ax.spines['right'].set_visible(False);ax.yaxis.set_ticks_position('left'); ax.xaxis.set_ticks_position('bottom');ax.set_axis_bgcolor('white');
    plt.ylabel('Correlation of effect sizes \n(mRNA ~ Protein)',fontsize=16)
    plt.xlabel('eQTL FDR (log 10)',fontsize=14)
    ax.annotate("0.1 qv",  xy=(-np.log10(0.05)+0.3,0.61));
    plt.xticks([1,5,10,15,20],[1,5,10,15,20])
     
#    plt.rcParams["legend.handletextpad"] = 1
    ax = fig.add_subplot(224)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False); ax.yaxis.set_ticks_position('left'); 
    for qtlname in ['eQTL','pQTL']:
        field='protein_nominal_p_value_bonf < 0.01' if qtlname=='eQTL' else  'mrna_nominal_p_value_bonf < 0.01' 
        x= np.log(abs(qtl[qtlname].query(field)['protein_eff_size']))-np.log(abs(qtl[qtlname].query(field)['mrna_eff_size']))
        x=x[np.isfinite(x)]
        sb.kdeplot(x,bw=0.025,label=qtlname,shade=True,color=colors[qtlname],alpha=0.51)
        plt.plot((np.median(x),np.median(x)),(0,1.5),'--',color=colors[qtlname],lw=1.5);
    print (x.shape)
#    plt.plot((0,0),(0,1.5),'k--',lw=0.5);
    plt.xlabel('Effect size ratio\n     Protein / mRNA (log)',fontsize=14)
    plt.ylabel('Density',fontsize=14)
    plt.legend(loc=2,fontsize=11)

 
    plt.tight_layout()
    plt.savefig(path_out+'SUPP_fig_eql_pqtl_eff_size.png',bbox_inches='tight',dpi=600)
    plt.savefig(path_out+'SUPP_fig_eql_pqtl_eff_size.svg',bbox_inches='tight',dpi=600)

def plot_distance_TSS(eqtl,pqtl,path_out=""):
    
    
    
    plt.figure(figsize=(8,6))
    plt.subplot(2,1,1)
    xx=eqtl.query('strand==1')
    strand=(xx['strand']==1).astype(int)
    xx['dist']=(xx['poz']-xx['gene_start']*strand-xx['gene_end']*(1-strand))*(-1)**(1-strand)
    x=xx.query('mrna_eff_size>-111')
    mm=max(x['dist'])#sb.histogram(np.clip(np.sign(x)*(abs(x)+1),-10**5.9,10**5.4),shade=1,color="green",bw=0.3)
    rr=np.array([-250000,-50000,0,50000,250000,500000])
    rr2=np.array(['-250k','-50k',0,'+50k','+250',">500k"])
    #sb.jointplot('dist','mrna_eff_size',data=x,color="green")
    plt.plot(np.clip(x['dist'],-250000,500000),abs(x['mrna_eff_size']),'.' ,color="green")
    plt.plot((0,0),(0,1.75),'k--',lw=0.95);plt.plot((-50000,-50000),(0,1.75),'k--',lw=0.95);plt.plot((50000,50000),(0,1.75),'k--',lw=0.95)
    #sb.kdeplot(np.clip(np.sign(x['dist'])*(abs(x['dist'])+1),-250000,250000),shade=1,color="darkblue",bw=0.3)
    plt.xticks(rr,rr2)
    plt.ylabel("eQTL abs. eff. size",fontsize=16)
    plt.xlabel("Distance from gene start (bp)",fontsize=16)
    
    plt.subplot(2,1,2)
    xx=pqtl
    strand=(xx['strand']==1).astype(int)
    xx['dist']=(xx['poz']-xx['gene_start']*strand-xx['gene_end']*(1-strand))*(-1)**(1-strand)
    x=xx.query('mrna_eff_size>-111').query('mrna_nominal_p_value_bonf<=0.01')
    #sb.histogram(np.clip(np.sign(x)*(abs(x)+1),-10**5.9,10**5.4),shade=1,color="green",bw=0.3)
    #sb.jointplot('dist','mrna_eff_size',data=x,color="green")
    plt.plot(np.clip(x['dist'],-250000,500000),abs(x['protein_eff_size']),'.' ,color="darkblue",label="eQTL PV<0.01")
    plt.plot((0,0),(0,1.75),'k--',lw=0.95);plt.plot((-50000,-50000),(0,1.75),'k--',lw=0.95);plt.plot((50000,50000),(0,1.75),'k--',lw=0.95)
    x=xx.query('mrna_eff_size>-111').query('mrna_nominal_p_value_bonf>0.01')
    #sb.histogram(np.clip(np.sign(x)*(abs(x)+1),-10**5.9,10**5.4),shade=1,color="green",bw=0.3)
    #sb.jointplot('dist','mrna_eff_size',data=x,color="green")
    plt.plot(np.clip(x['dist'],-250000,500000),abs(x['protein_eff_size']),'.' ,color="skyblue",label="eQTL PV>0.01")
    plt.plot((0,0),(0,1.75),'k--',lw=0.95);plt.plot((-50000,-50000),(0,1.75),'k--',lw=0.95);plt.plot((50000,50000),(0,1.75),'k--',lw=0.95)
    plt.legend(loc=1)
    plt.xticks(rr,rr2)

    plt.ylabel("pQTL abs. eff. size",fontsize=16)
    
    plt.xlabel("Distance from gene start (bp)",fontsize=16)
    plt.tight_layout()

    plt.savefig(path_out+'SUPP_fig_eql_pqtl_distance.png',bbox_inches='tight',dpi=600)
    plt.savefig(path_out+'SUPP_fig_eql_pqtl_distance.svg',bbox_inches='tight',dpi=600)



def load_mrna_protein_function(folder_data):
    transcript_df=pd.read_table('/Users/mirauta/Data/MS/hipsci/phenotypes/Paired_Expression_ISR_Stringent.txt.gz',index_col=0)
    mrna_df_ensembl=pd.read_table('/Users/mirauta/Data/MS/hipsci/phenotypes/IPSc.ISR.featureCounts.genes.counts.unique.stranded.tsv_counts.HQ_TMM_TPM.tsv.gz',index_col=0)
    phenotype_df0=pandas.read_table(folder_data+'hipsci.proteomics.maxquant.uniprot.TMT_batch_24.20170919_protein_Reporter intensity_ranknorm_regbatch_valid_lines_all_genes.txt',index_col=0)
    peptide_df0=pandas.read_table(folder_data+'hipsci.proteomics.maxquant.uniprot.TMT_batch_24.20170919_peptide_Reporter intensity_ranknorm_regbatch_valid_lines_all_genes.txt.zip',index_col=0)
    linemeta=pandas.read_table(folder_data+'hipsci.proteomics.maxquant.uniprot.TMT_batch_24.20170919_line_map_valid_lines_noreplicates_mrna_simple.txt',index_col=0,header=None)
    
    peptide_dfall=peptide_df0.copy(deep=1)
    phenotype_dfall=phenotype_df0.copy(deep=1)
    
    peptide_df0=peptide_df0[linemeta[1]]
    phenotype_df0=phenotype_df0[linemeta[1]]
    
    mrna_df_ensembl=mrna_df_ensembl[mrna_df_ensembl.columns[np.in1d(np.array(['-'.join(c.split('.')[:2]) for c in mrna_df_ensembl.columns]),       np.array([c.split('@')[0] for c in phenotype_df0.columns]))]]
    mrna_df_ensembl.columns=['-'.join(c.split('.')[:2]) for c in mrna_df_ensembl.columns]
    mrna_df_ensembl=mrna_df_ensembl[linemeta.index]
    
    transcript_df=transcript_df[transcript_df.columns[np.in1d(np.array(['-'.join(c.split('.')[:2]) for c in transcript_df.columns]),       np.array([c.split('@')[0] for c in phenotype_df0.columns]))]]
    transcript_df.columns=['-'.join(c.split('.')[:2]) for c in transcript_df.columns]
    transcript_df=transcript_df[linemeta.index]
    
    return [mrna_df_ensembl,transcript_df,peptide_df0,phenotype_df0,peptide_dfall,phenotype_dfall]


def regout(A,x,y):
    regressor = LinearRegression()
    index=np.isfinite(np.sum(A[np.hstack([x,y])].values,1))
    regressor.fit(X=A[x][index],y=A[y][index])
    y2=A[y]-np.dot(regressor.coef_,A[x].T)
    return [y2,regressor.coef_]
 
def correlation(x,y,typecorr="Pearson"):
    
    index=np.isfinite(x+y)
    if typecorr=="Pearson" : return np.corrcoef(x[index],y[index])[0,1]

def regxfromy(x,y):
    regressor = LinearRegression()
    index=np.isfinite(np.sum(np.hstack([x,y[:,None]]),1))
    regressor.fit(X=x[index],y=y[index])
    residual=y-np.dot(regressor.coef_,x.T)
    residual=residual-np.nanmean(residual)+np.nanmean(y)
    return [residual,regressor.coef_]

def load_string(prometa):
    string=pd.read_table('/Users/mirauta/Data/Annotation/complexes/string_database.9606.protein.links.v10.5.uniprot_ids.txt').drop_duplicates()
    string.columns=['A','B']
    string['A']=np.array([i.split('-')[0]for i in string['A']])
    string['B']=np.array([i.split('-')[0]for i in string['B']])
    string=string.set_index('A',drop=0); temp=np.intersect1d(np.unique(string['A']),np.unique(prometa['no_isoform_info']))
    string=string.loc[temp]
    string=string.set_index('B',drop=0); temp=np.intersect1d(np.unique(string['B']),np.unique(prometa['no_isoform_info']))
    string=string.loc[temp]
    string['pair']=string['A']+'_'+string['B']
    print (string.shape)
    string[['A','B','pair']].to_csv('/Users/mirauta/Data/Annotation/complexes/string_processed.txt',sep='\t')
    return string


def load_complex_prometa(intact_flag=False):
    proteinmeta=pd.read_table('/Users/mirauta/data/MS/hipsci/TMT_batch_24_uniprot_MQv1.6/hipsci.proteomics.maxquant.uniprot.TMT_batch_24.20170919.proteinGroups.txt' , sep='\t')
    proteinmeta['mean']=np.nanmean(proteinmeta[proteinmeta.columns[[('Reporter intensity HPSI' in c)&('Calcium' not in c) for c in proteinmeta.columns]]],1)
    proteinmeta['nlines']=np.nansum(proteinmeta[proteinmeta.columns[[('Reporter intensity HPSI' in c)&('Calcium' not in c) for c in proteinmeta.columns]]]>0,1)
    proteinmeta=proteinmeta[proteinmeta['mean']>0]
    proteinmeta=proteinmeta[np.array(['REV_' not in a.split(";")[0] for a in proteinmeta['Majority protein IDs']])& np.array(['CON_' not in a.split(";")[0] for a in proteinmeta['Majority protein IDs']])]; print (proteinmeta.shape)
    proteinmeta.index=np.hstack([ p.split(";")[0] for p in proteinmeta['Majority protein IDs']])
    
    temp0=np.vstack([len(p.split(";")) for p in proteinmeta['Majority protein IDs']]);print (temp0.shape)
    proteinmeta['Gene names'][proteinmeta['Gene names']!=proteinmeta['Gene names']]=''
    temp2= np.hstack([np.repeat(p.split(";")[0],temp0[ip]) for ip,p in enumerate(proteinmeta['Gene names'])]);print (temp2.shape)
    temp3= np.hstack([np.repeat(proteinmeta['mean'].iloc[ip],temp0[ip]) for ip,p in enumerate(proteinmeta['Gene names'])]);print (temp3.shape)
    temp4= np.hstack([np.repeat(proteinmeta['nlines'].iloc[ip],temp0[ip]) for ip,p in enumerate(proteinmeta['Gene names'])]);print (temp4.shape)
    temp=np.vstack([np.vstack([np.repeat(p.split(";")[0],len(p.split(";"))),p.split(";")]).T  for p in proteinmeta['Majority protein IDs']])
    prometa=pd.DataFrame(temp,index=temp[:,0],columns=['leading','protein_group'])
    prometa['gene_name']=temp2
    prometa['protein_mean_intensity']=temp3
    prometa['protein_nlines']=temp4
    prometa=prometa.set_index('protein_group',drop=0)
    prometa['no_isoform_info']=np.array([c .split('-')[0] for c in prometa.index])
    
    #'''load complexes'''
    complexes=pd.read_table( '/Users/mirauta/Data/Annotation/complexes/allComplexes_CORUM_12_12_17.txt',sep='\t',header=0).set_index('ComplexID',drop=0)
    complexes=complexes.iloc[complexes['subunits(UniProt IDs)'].values==complexes['subunits(UniProt IDs)'].values]
    complexes=complexes.loc[complexes[['ComplexName','subunits(UniProt IDs)']].drop_duplicates().index]
    complexes['ComplexName'][complexes['ComplexName']!=complexes['ComplexName']]=''
    
    complexes_protein=pd.DataFrame(np.vstack([np.hstack([complexes.loc[c, 'subunits(UniProt IDs)'].split(';') for c in complexes.index]),\
              np.hstack([np.repeat(c,len(complexes.loc[c, 'subunits(UniProt IDs)'].split(';'))) for c in complexes.index]),\
               np.hstack([np.repeat(complexes.loc[c,'ComplexName'],len(complexes.loc[c, 'subunits(UniProt IDs)'].split(';'))) for c in complexes.index])]).T,columns=['protein_uniprot_id','complex_id','ComplexName'])
    complexes_protein=complexes_protein.set_index('protein_uniprot_id',drop=0)
    complexes_protein['no_isoform_info']=np.array([c .split('-')[0] for c in complexes_protein.index])
    complexes_protein=complexes_protein.set_index('no_isoform_info',drop=0)
     
    procomplex=np.unique(np.hstack([complexes['subunits(UniProt IDs)'].loc[c].replace(')','').replace('(','').split(';') for c in complexes.index]))
#    print (np.intersect1d(procomplex,prometa['leading'].values).shape)
    procomplex=np.unique(prometa.set_index('no_isoform_info',drop=0).loc[np.intersect1d(np.array([c .split('-')[0] for c in procomplex]),\
                                                    np.unique(prometa['no_isoform_info']))]['no_isoform_info']);
    print (procomplex.shape)
    print (complexes_protein.shape)
    np.unique(complexes_protein.index)
        
    prometa['in_corum']=np.zeros(prometa.shape[0])
    prometa['in_corum']=np.in1d(prometa['no_isoform_info'],procomplex)
    complexes=complexes.set_index('ComplexName',drop=0)
    complexes_protein=complexes_protein.loc[np.intersect1d(complexes_protein['no_isoform_info'],prometa['no_isoform_info'])]
#    complexes_protein['leading']=prometa.set_index('no_isoform_info',drop=0)[['no_isoform_info','leading']].drop_duplicates().loc[complexes_protein.index,'leading']
    complexes_protein=complexes_protein.set_index('complex_id')
    temp0={}
    for i in np.unique(complexes_protein.index):
        temp=complexes_protein.loc[[i]]['no_isoform_info'].values
        if len(temp)>1:
            temp0[i]=np.hstack([np.vstack([np.repeat(t,len(temp)),temp,np.repeat(complexes_protein.loc[[i]]['ComplexName'].iloc[0],len(temp))]) for it,t in enumerate(temp)]).T
    temp=pd.DataFrame(np.vstack(temp0.values()),columns=['A','B','ComplexName'])        
    complex_pairs=temp[temp['A']!=temp['B']]
    temp=complex_pairs.copy(deep=True); temp['A']=complex_pairs['B'].copy(deep=True);temp['B']=complex_pairs['A'].copy(deep=True)
    complex_pairs=pd.concat([complex_pairs,temp]).drop_duplicates()
    complex_pairs['pair']=np.array([complex_pairs.loc[c,'A'].split('-')[0]+'_'+complex_pairs.loc[c,'B'].split('-')[0] for c in complex_pairs.index])
    
    if intact_flag:
        intact=pd.read_table('/Users/mirauta/Data/Annotation/complexes/intact_human.txt').drop_duplicates()
        intact['A']=[c.replace('uniprotkb:','') for c in intact['#ID(s) interactor A']]
        intact['B']=[c.replace('uniprotkb:','') for c in intact['ID(s) interactor B']]
        temp=intact.copy(deep=True); temp['A']=intact['B'].copy(deep=True);temp['B']=intact['A'].copy(deep=True)
        intact=pd.concat([intact,temp]).drop_duplicates()
        intact['A']=np.array([i.split('-')[0]for i in intact['A']])
        intact['B']=np.array([i.split('-')[0]for i in intact['B']])
        intact['pair']=  intact['A']+'_'+ intact['B']
        
        
        print ('loading_string')
        string=pd.read_table('/Users/mirauta/Data/Annotation/complexes/string_database.9606.protein.links.v10.5.uniprot_ids.txt').drop_duplicates()
        string.columns=['A','B']
        string['A']=np.array([i.split('-')[0]for i in string['A']])
        string['B']=np.array([i.split('-')[0]for i in string['B']])
        string=string.set_index('A',drop=0); temp=np.intersect1d(np.unique(string['A']),np.unique(prometa['no_isoform_info']))
        string=string.loc[temp]
        string=string.set_index('B',drop=0); temp=np.intersect1d(np.unique(string['B']),np.unique(prometa['no_isoform_info']))
        string=string.loc[temp]
        string['pair']=string['A']+'_'+string['B']
        print (string.shape)
#        intact_string=pd.concat([intact,string]).drop_duplicates()
    
#    complex_pairs_prgroup=complex_pairs.copy(deep=True)
#    complex_pairs_prgroup=complex_pairs_prgroup.set_index('A',drop=0).loc\
#                                   [np.intersect1d(prometa['protein_group'],complex_pairs_prgroup['A'])]
##    complex_pairs_prgroup['Alead']=prometa.loc[complex_pairs_prgroup.index,'leading'] 
#    
#    complex_pairs_prgroup=complex_pairs_prgroup.set_index('B',drop=0).loc\
#                                   [np.intersect1d(prometa['protein_group'],complex_pairs_prgroup['B'])]
#    complex_pairs_prgroup['Blead']=prometa.loc[complex_pairs_prgroup.index,'leading'] 
#    complex_pairs_prgroup['pairlead']=complex_pairs_prgroup['Alead']+'_'+complex_pairs_prgroup['Blead']

    if intact_flag:
        return [prometa,complexes,procomplex,procomplex,complexes_protein,complex_pairs,intact, string]    
    return [prometa,complexes,procomplex,procomplex,complexes_protein,complex_pairs]    


def get_all_uniprot_ensembl(ensembl='Ensembl_37.75',exon=False):
    annotation_df0=pandas.read_table('/Users/mirauta/Data/Annotation/'+ensembl+'/mart_export.txt.gz',sep='\t').set_index('Gene stable ID',drop=0)
    annotation_df0=annotation_df0.drop_duplicates()
#    
#    annotation_df0['UniProtKB Gene Name ID']==annotation_df0['UniProtKB Gene Name ID']
#    annotation_df0[(annotation_df0['UniProtKB/Swiss-Prot ID']!=annotation_df0['UniProtKB/Swiss-Prot ID'])&\
#    (annotation_df0['UniProtKB/TrEMBL ID']!=annotation_df0['UniProtKB/TrEMBL ID'])]
#    
#    b=annotation_df0[['Gene stable ID', 'Transcript stable ID','Gene name','UniProtKB Gene Name ID']].loc[ annotation_df0 ['UniProtKB Gene Name ID']==annotation_df0 ['UniProtKB Gene Name ID']]
#    b.columns=['Gene stable ID','Transcript stable ID','Gene name','feature_id']
    a=[annotation_df0[['Gene stable ID', 'Transcript stable ID','Exon stable ID','Gene name',field,'Gene start (bp)', 'Gene end (bp)', 'Strand',\
        'Chromosome/scaffold name','Transcript type']]  for field in [ 'UniProtKB/Swiss-Prot ID', 'UniProtKB/TrEMBL ID']]
    a=[a[ifield][a[ifield][field]==a[ifield][field]] for ifield,field in enumerate([ 'UniProtKB/Swiss-Prot ID', 'UniProtKB/TrEMBL ID'])]
    if exon: 
        a=[a[ifield].iloc[np.unique(a[ifield]['Gene stable ID']+a[ifield]['Transcript stable ID']+a[ifield]['Exon stable ID']+a[ifield][field],return_index=1)[1]] for ifield,field in enumerate([ 'UniProtKB/Swiss-Prot ID', 'UniProtKB/TrEMBL ID'])]
        for aa in a:  aa.columns=['ensembl_gene_id', 'ensembl_transcript_id','ensembl_exon_id','gene_name', 'protein_uniprot_id','gene_start','gene_end','strand','chromosome','transcript_type']
        a=pandas.concat(a);
        a=a.iloc[np.unique(a['ensembl_gene_id']+a['ensembl_transcript_id']+a['ensembl_exon_id']+a['protein_uniprot_id'],return_index=1)[1]]
        
    else: 
        a=[a[ifield].iloc[np.unique(a[ifield]['Gene stable ID']+a[ifield]['Transcript stable ID']+a[ifield][field],return_index=1)[1]] for ifield,field in enumerate([ 'UniProtKB/Swiss-Prot ID', 'UniProtKB/TrEMBL ID'])]
        for aa in a:  aa.columns=['ensembl_gene_id', 'ensembl_transcript_id','ensembl_exon_id','gene_name', 'protein_uniprot_id','gene_start','gene_end','strand','chromosome','transcript_type']
        a=pandas.concat(a);
        a=a.iloc[np.unique(a['ensembl_gene_id']+a['ensembl_transcript_id']+a['protein_uniprot_id'],return_index=1)[1]]
        
    a.to_csv('/Users/mirauta/Data/Annotation/'+ensembl+'/'+ensembl+'_processed.txt',  mode='w', sep='\t') 
#    a1=a[a['protein_uniprot_id']==a['protein_uniprot_id']]
    a1=a
    a1.to_csv('/Users/mirauta/Data/Annotation/'+ensembl+'/'+ensembl+'_protein_id_processed.txt',  mode='w', sep='\t') 
          
    return(a)

def fun_split(x,pos=0):    return x.split('_')[pos]




def manhattan_traits(path_data,folder_destination,traits,trait_labels,gene,plot_name):
    plot_manhatan_alone( gene_ensembl_id= gene,folder_name=path_data,\
                         folder_destination=folder_destination+'manuscript_images/'+'/Manhattan/',\
                         plot_name=plot_name,traits=traits,trait_labels=trait_labels,file_name_results_genome='ensembl_gene_id_qtl_results_genome',\
                         qtl_results_file='qtl_results_',colors=np.array(['black','green','orange','blue']),p_value_field='p_value_raw',\
                         figsize=1.75,pdf=False,savefig=True)

def plot_effect_scatter(path_data,figsize=(5,5),main_trait='Protein6', main_feature_id=None,second_trait='Protein6',folder_destination='',data='',gene='',plot_name='',colors=np.array(['black','green','orange','blue'])):
    x={}
    s={}
    plt.figure(figsize=figsize)
    if main_feature_id is None:
        x[main_trait]=data[main_trait][gene]['data/beta'][data[main_trait][gene]['summary_data/min_p_value_feature_id'][:][0]][:]
        s[main_trait]=data[main_trait][gene]['data/snp_id'][data[main_trait][gene]['summary_data/min_p_value_feature_id'][:][0]][:]
    else:
        x[main_trait]=data[main_trait][gene]['data/beta'][main_feature_id][:]
        s[main_trait]=data[main_trait][gene]['data/snp_id'][main_feature_id][:]
    
    transcripts=data[second_trait][gene]['metadata']['feature_id'][:].astype('U');transcripts=[t.split('_')[0] for t in transcripts]
    for trait in [second_trait,'mRNA']:

        x[trait]= [data[trait][gene]['data/beta'][ff][:] for ff in list(data[trait][gene]['data/beta/'].keys())] 
        s[trait]= [data[trait][gene]['data/snp_id'][ff][:] for ff in list(data[trait][gene]['data/beta/'].keys())] 
        
        x[trait]=[x[trait][ss][np.in1d(s[trait][ss],np.intersect1d(s[main_trait],s[trait][ss]))]for ss in np.arange(len(x[trait]))]
        corr=np.array([np.corrcoef(x[trait][ss],x[main_trait][np.in1d(s[main_trait],np.intersect1d(s[main_trait],s[trait][ss]))])[0,1] for ss in np.arange(len(x[trait]))])
        x[trait]=np.array(x[trait])[np.argsort(abs(corr))[-5:]];s[trait]=np.array(s[trait])[np.argsort(abs(corr))[-5:]]
        for ss in np.arange(x[trait].shape[0]):
            plt.plot(x[main_trait][np.in1d(s[main_trait],np.intersect1d(s[main_trait],s[trait][ss]))],x[trait][ss], 'o', mfc='none',color='grey'if trait=='mRNA'else colors[ss+2],markersize=3,\
                     label='gene level mRNA'if trait=='mRNA'else transcripts[ss])
    plt.legend(loc=1,fontsize=8)
    if main_feature_id is None: plt.xlabel(main_trait+' eff. size')
    else: plt.xlabel(main_trait+' '+main_feature_id+' eff. size')
    plt.ylabel(second_trait+' eff. size')
#    print (data[main_trait][gene]['metadata/gene_name'][:][0])
    plt.savefig(folder_destination+'manuscript_images/Manhattan/'+plot_name+ data[main_trait][gene]['metadata/gene_name'][:].astype('U')[0] +'_scatter.png',dpi=600,bbox_inches='tight')
    plt.savefig(folder_destination+'manuscript_images/Manhattan/'+plot_name+ data[main_trait][gene]['metadata/gene_name'][:].astype('U')[0] +'_scatter.svg',dpi=600,bbox_inches='tight')
    
def distribtion_windows():
    lengths=pd.DataFrame(data=np.zeros([len(data['Protein']),4]),columns=['window','nsnps','minpv','minpv_raw'],index=list(data['Protein'].keys()))
    for ig,gene in enumerate(data['Protein'].keys()):
        if ig%100==0: print (ig)
        positions=data['Protein'][gene]['data/position'][data['Protein'][gene]['summary_data/min_p_value_feature_id'][:][0]][:]  
        lengths.loc[gene,'window']=positions.max()-positions.min()
        lengths.loc[gene,'nsnps']=len(positions)
        lengths.loc[gene,'minpv']=data['Protein'][gene]['summary_data']['min_p_value'][:][0]
        lengths.loc[gene,'minpv_raw']=data['Protein'][gene]['data/p_value_raw'][data['Protein'][gene]['summary_data/min_p_value_feature_id'][:][0]][:].min()  
        
        
    plt.hist((lengths['window']),bins=100)
    plt.xlabel('cis window (last - first snp)')

#    
#def plot_effect_scatter_2genes(path_data,data,qtls,genes,folder_destination,plot_name):
#    x={}
#    s={}
#    plt.figure(figsize=(5,5))
#    
#    for i,trait in enumerate(np.array(qtls)):
#        1
#        x[trait]= [data[trait][genes[i]]['data/beta'][ff][:] for ff in list(data[trait][genes[i]]['data/beta/'].keys())] 
#        s[trait]= [data[trait][genes[i]]['data/snp_id'][ff][:].astype('U') for ff in list(data[trait][genes[i]]['data/beta/'].keys())] 
#    
#    x[trait]=[x[trait][ss][np.in1d(s[trait][ss],np.intersect1d(s['Protein'],s[trait][ss]))]for ss in np.arange(len(x[trait]))]
#    
#    for ss in np.arange(x[trait].shape[0]):
#       plt.plot(x['Protein'][np.in1d(s['Protein'],np.intersect1d(s['Protein'],s[trait][ss]))],x[trait][ss], 'o',color='grey'if trait=='mRNA'else colors[ss],markersize=2,label='mRNA'if trait=='mRNA'else 'Trans '+str(ss))
#    plt.legend(loc=2,fontsize=8)
#    plt.xlabel('Protein eff. size')
#    plt.ylabel('mRNA/Transcript\n eff. size')
##    print (data['Protein'][gene]['metadata/gene_name'][:][0])
#    plt.savefig(folder_destination+'manuscript_images/Manhattan/'+plot_name+ data['Protein'][gene]['metadata/gene_name'][:].astype('U')[0] +'_scatter.png',dpi=600,bbox_inches='tight')
#    

def customaxis(ax, c_left='k', c_bottom='k', c_right='none', c_top='none',
               lw=3, size=12, pad=8):

    for c_spine, spine in zip([c_left, c_bottom, c_right, c_top],
                              ['left', 'bottom', 'right', 'top']):
        if c_spine != 'none':
            ax.spines[spine].set_color(c_spine)
            ax.spines[spine].set_linewidth(lw)
        else:
            ax.spines[spine].set_color('none')
    if (c_bottom == 'none') & (c_top == 'none'): # no bottom and no top
        ax.xaxis.set_ticks_position('none')
    elif (c_bottom != 'none') & (c_top != 'none'): # bottom and top
        ax.tick_params(axis='x', direction='out', width=lw, length=7,
                      color=c_bottom, labelsize=size, pad=pad)
    elif (c_bottom != 'none') & (c_top == 'none'): # bottom but not top
        ax.xaxis.set_ticks_position('bottom')
        ax.tick_params(axis='x', direction='out', width=lw, length=7,
                       color=c_bottom, labelsize=size, pad=pad)
    elif (c_bottom == 'none') & (c_top != 'none'): # no bottom but top
        ax.xaxis.set_ticks_position('top')
        ax.tick_params(axis='x', direction='out', width=lw, length=7,
                       color=c_top, labelsize=size, pad=pad)
    if (c_left == 'none') & (c_right == 'none'): # no left and no right
        ax.yaxis.set_ticks_position('none')
    elif (c_left != 'none') & (c_right != 'none'): # left and right
        ax.tick_params(axis='y', direction='out', width=lw, length=7,
                       color=c_left, labelsize=size, pad=pad)
    elif (c_left != 'none') & (c_right == 'none'): # left but not right
        ax.yaxis.set_ticks_position('left')
        ax.tick_params(axis='y', direction='out', width=lw, length=7,
                       color=c_left, labelsize=size, pad=pad)
    elif (c_left == 'none') & (c_right != 'none'): # no left but right
        ax.yaxis.set_ticks_position('right')
        ax.tick_params(axis='y', direction='out', width=lw, length=7,
                       color=c_right, labelsize=size, pad=pad)
        
def reduce_snp_draft(snp_df,threshold=0.8):
    allsnps=snp_df.columns
    cor=abs(np.corrcoef(snp_df.T))>=threshold
#    
##    cor=pd.DataFrame(np.zeros([4,4]),index=np.arange(4),columns=np.arange(4))                    
##    allsnps=cor.columns
##    for i in np.arange(4):cor.loc[i,i]=1
#    cor[1,2]    =1;cor[2,1]    =1
#    cor[2,3]    =1;cor[3,2]    =1
#    cor[1,3]    =1;cor[3,1]    =1
#    cor[2,5]    =1;cor[5,2]    =1
#    sn=allsnps[1:6]
#    cor[1,3]    =1;cor.loc[3,1]    =1
#    print (cor)
    cor=pd.DataFrame(data=cor,index=allsnps,columns=allsnps)
#    cor=cor.loc[sn,sn];allsnps=sn
    cor1=pd.DataFrame(data=np.triu(cor,k=1),index=allsnps,columns=allsnps).sum(1)
    cor2=pd.DataFrame(data=np.triu(cor,k=1),index=allsnps,columns=allsnps)
#    cor3=pd.DataFrame(data=np.tril(cor,k=-1),index=allsnps,columns=allsnps)
    uniquesnps=cor1[cor1==0].index
    duplicatedsnps=cor1[cor1>0].index
    '''the unique ones'''
    rez=pd.DataFrame(data=uniquesnps,index=uniquesnps,columns=['snp_id'])
    ''' the ones that are in LD.. if a,c and a,b but not b,c returns a,b)'''
    rez2=pd.DataFrame(data=duplicatedsnps,index=allsnps[np.argmax(cor2.loc[duplicatedsnps].values*cor2.loc[duplicatedsnps].values.sum(0),1)],columns=['snp_id'])
#    rez2=pd.DataFrame(data=duplicatedsnps,index=allsnps[np.argmax(cor2.loc[duplicatedsnps].values,1)],columns=['snp'])
    rez=pd.concat([rez,rez2])
    rez['lead_snp_id']=rez.index
#    print (rez)
    rez.set_index('snp_id').loc[sn]
    rez.loc[sn]
    sb.heatmap(cor.loc[sn,sn][::-1]+0)
    return(rez)
#snp_df=pd.read_table("/Users/mirauta/Downloads/- (2).txt",sep=',',index_col=0)
#snp_df.columns=['_'.join(s.split('_')[1:]) for s in snp_df.columns]
def reduce_snp(snp_df,threshold=0.8):
    ''' input a snp df  samples(rows) x snps( columns)'''
    ''' returns a df with columns: 'lead_snp_id' the pruned snp names and snp_id'''
    allsnps=snp_df.columns
    cor=abs(np.corrcoef(snp_df.T))>=threshold
    cor1=pd.DataFrame(data=np.triu(cor,k=1),index=allsnps,columns=allsnps).sum(1)
    cor2=pd.DataFrame(data=np.triu(cor,k=1),index=allsnps,columns=allsnps)

    uniquesnps=cor1[cor1==0].index.values
    duplicatedsnps=cor1[cor1>0].index.values
    '''the unique ones'''
    rez=pd.DataFrame(data=uniquesnps,index=uniquesnps,columns=['snp_id'])
    ''' the ones that are in LD.. if a,c and a,b but not b,c returns a,b)'''
    rez2=pd.DataFrame(data=duplicatedsnps,index=allsnps[np.argmax(cor2.loc[duplicatedsnps].values*cor2.loc[duplicatedsnps].values.sum(0),1)],columns=['snp_id'])
    rez=pd.concat([rez,rez2])
    rez['lead_snp_id']=rez.index
    return(rez)

#rez=reduce_snp(snp_df,threshold=0.5)
#rez.shape
#np.unique(rez.index).shape
#values=pd.DataFrame(data=np.arange(5),index=np.unique(rez.index),columns=['pv'])
#extended=pd.DataFrame(index=rez['snp_id'],columns=['pv'])
#extended.loc[rez['snp_id'],'pv']=values.loc[rez.set_index('snp_id').loc[rez['snp_id'],'lead_snp_id'],'pv'].values
#
#         
def trans_bring_cis_pv(trans0,df,prometa):
    for s in trans0.keys():
        trans0[s]=trans0[s].set_index('snp_id',drop=0)
        trans0[s]['cis_protein_uniprot_id']=df[s].set_index('snp_id').iloc[np.unique(df[s]['snp_id'],return_index=1)[1]].loc[trans0[s].index]['feature_id']
        trans0[s]['cis_gene_name']=df[s].set_index('snp_id').iloc[np.unique(df[s]['snp_id'],return_index=1)[1]].loc[trans0[s].index]['gene_name']
        trans0[s]['cis_protein_p_value']=df[s].set_index('snp_id').iloc[np.unique(df[s]['snp_id'],return_index=1)[1]].loc[ trans0[s].index]['p_value_raw']
        trans0[s]['cis_protein_beta']=df[s].set_index('snp_id').iloc[np.unique(df[s]['snp_id'],return_index=1)[1]].loc[ trans0[s].index]['beta']
        trans0[s]['pair']=trans0[s]['feature_id']+'_'+trans0[s]['cis_protein_uniprot_id']
        trans0[s]['index']=trans0[s].index+'_'+trans0[s]['ensembl_gene_id']+'_'+trans0[s]['feature_id']
  
    return trans0
def return_cistranslead_pv(trans0,df,prometa):
    translead={}
    for s in trans0.keys():
        dd=trans0[s].copy(deep=True)
        dd=dd.set_index('index',drop=0)

        a=dd.groupby('cis_protein_uniprot_id').apply(lambda dd:dd.empirical_feature_p_value.argmin())
        translead[s]=dd.loc[a.values]
    return translead


def return_significant_field_qv(trans0,df,prometa,thr,field):
    rez={}
    
    for name in trans0.keys():
        if field=='snp_id':nfeatures=np.unique(trans0[name]['feature_id']).shape[0]
        elif field=='feature_id':nfeatures=np.unique(trans0[name]['snp_id']).shape[0]
        
        dd=trans0[name].copy(deep=True)
        dd[field+'_qv']=np.zeros(dd.shape[0])+np.nan
        dd=dd.set_index('pair',drop=0)
        dd0=dd.set_index(field,drop=0)
        cisfeatures=np.unique(dd[field])
        
        for cis in cisfeatures:
            x=dd0.loc[[cis]]
            extended=np.hstack([x['p_value'].values, np.random.uniform(x['p_value'].values.max(),1,nfeatures-x.shape[0])])
            dd.loc[x['pair'],field+'_qv']= scst2.fdrcorrection0(extended)[1][:x.shape[0]]
        rez[name]=dd[dd[field+'_qv']<thr]
#        a=dd.groupby('cis_protein_uniprot_id').apply(lambda dd:np.where(dd.snp_qv<0.05)[0])
#        translead[name]=dd.iloc[np.hstack(a.values)]
    return rez


def return_transcislead_pv(trans0,df,prometa):
    translead={}
    for s in trans0.keys():
        dd=trans0[s].copy(deep=True)
        dd['index']=dd['feature_id']+'_'+dd['snp_id']
        dd=dd.set_index('index',drop=0)

        a=dd.groupby('feature_id').apply(lambda dd:dd.empirical_feature_p_value.argmin())
        translead[s]=dd.loc[a.values]
    return translead

def load_gwas_old():
    annotation_snp=pd.read_table('/Users/mirauta/Data/Genotypes/gwas_catalog_v1.0.1-associations_e88_r2017-04-24_2.txt',sep='\t',encoding='latin_1').set_index('SNP_ID_Current_Data',drop=0)
    annotation_snp=annotation_snp.iloc[np.where(annotation_snp.index.values==annotation_snp.index.values)[0]]; print (annotation_snp.shape)
    annotation_snp['gene_name']=np.array([np.array(re.split(',| - ',g)) for g in annotation_snp['MAPPED_GENE'].values]); print (annotation_snp.shape)
    annotation_snp['gene_name']=np.array([g[np.array(['LOC' not in gg for gg in g])] for g in annotation_snp['gene_name']]); print (annotation_snp.shape)
    annotation_snp['number_genes']=np.array([len(g) for g in annotation_snp['gene_name']]); print (annotation_snp.shape)
    annotation_snp['index']=annotation_snp['SNP_ID_Current_Data']+'_'+annotation_snp['MAPPED_TRAITs']
    annotation_snp=annotation_snp[annotation_snp['index']==annotation_snp['index']]
    annotation_snp=annotation_snp.iloc[np.unique(annotation_snp['index'],return_index=1)[1]]; print (annotation_snp.shape)
    annotation_snp['SNP_pos']=[i.split('_')[1] for i in annotation_snp['SNP_ID_Current_Data']]
    annotation_snp['SNP_chrom']=[i.split('_')[0] for i in annotation_snp['SNP_ID_Current_Data']]
    
    annotation_snp[[]].to_csv('/Users/mirauta/Results/hipsci/QTL_june2018/gwas_snps.txt')
    return annotation_snp
def sum_ints(text):
    a=np.array(text.replace(',','') .split(' '))
    return a[np.array([aa.isnumeric() for aa in a])& np.array([a[iaa-1]!="from" for iaa,aa in enumerate(a)])].astype(int).sum()
    
def load_gwas2(gwas_file='gwas_catalog_v1.0.1-associations_e92_r2018-04-10.txt'):
    
    annotation_snp=pd.read_table('/Users/mirauta/Data/Genotypes/'+gwas_file,sep='\t',encoding='latin_1').set_index('SNP_ID_CURRENT',drop=0)
    annotation_snp=annotation_snp[annotation_snp['CHR_POS_HG37']==annotation_snp['CHR_POS_HG37']]
    annotation_snp=annotation_snp[annotation_snp['CHR_ID_HG37']==annotation_snp['CHR_ID_HG37']]
    
    annotation_snp['SNP_ID_CURRENT']=annotation_snp['CHR_ID_HG37'].astype('U')+'_'+annotation_snp['CHR_POS_HG37'].astype('U')
    annotation_snp['SNP_ID_CURRENT']=[s.split('.')[0] for s in annotation_snp['SNP_ID_CURRENT']]; print (annotation_snp.shape)
    annotation_snp[['SNP_ID_CURRENT','SNPS']].head()
    
    annotation_snp1=annotation_snp.set_index('SNP_ID_CURRENT',drop=0)
    annotation_snp=annotation_snp.set_index('SNP_ID_CURRENT',drop=0)

    annotation_snp['index']=annotation_snp['SNP_ID_CURRENT'].astype('U')+'_'+annotation_snp['MAPPED_TRAIT'].astype('U')

    print (annotation_snp.shape)
    annotation_snp=annotation_snp[annotation_snp['OR or BETA']==annotation_snp['OR or BETA']]
    print (annotation_snp.shape)
    annotation_snp['Study_size']=[sum_ints(text=text) for text in annotation_snp['INITIAL SAMPLE SIZE'].values]

#    
    annotation_snp=annotation_snp[annotation_snp['Study_size']>1000]
    
    
    annotation_snp=annotation_snp.iloc[np.unique(annotation_snp['index'],return_index=1)[1]]; print (annotation_snp.shape)
    annotation_snp['SNP_pos']=[i.split('_')[1] for i in annotation_snp['SNP_ID_CURRENT']]
    annotation_snp['SNP_chrom']=[i.split('_')[0] for i in annotation_snp['SNP_ID_CURRENT']]
    annotation_snp=annotation_snp[['SNP_ID_CURRENT','SNP_pos','SNP_chrom','MAPPED_TRAIT','P-VALUE','OR or BETA','INITIAL SAMPLE SIZE','SNPS','Study_size']]


    annotation_snp[[]].to_csv('/Users/mirauta/Results/hipsci/QTL_june2018/gwas_snps.txt')
    return [annotation_snp,annotation_snp1]
   
def load_gwas(gwas_file='gwas_catalog_v1.0.1-associations_e92_r2018-04-10.txt',snp_id='SNP_ID_CURRENT',mapped_trait='MAPPED_TRAIT'):
    
    annotation_snp=pd.read_table('/Users/mirauta/Data/Genotypes/'+gwas_file,sep='\t',encoding='latin_1')
    annotation_snp=annotation_snp.set_index(snp_id,drop=0)
    annotation_snp=annotation_snp[annotation_snp['CHR_POS_HG37']==annotation_snp['CHR_POS_HG37']]
    annotation_snp=annotation_snp[annotation_snp['CHR_ID_HG37']==annotation_snp['CHR_ID_HG37']]
    
    annotation_snp['SNP_ID_CURRENT']=annotation_snp['CHR_ID_HG37'].astype('U')+'_'+annotation_snp['CHR_POS_HG37'].astype('U')
    annotation_snp['SNP_ID_CURRENT']=[s.split('.')[0] for s in annotation_snp['SNP_ID_CURRENT']]; print (annotation_snp.shape)
#    annotation_snp[['SNP_ID_CURRENT','SNPS']].head()
    
    annotation_snp=annotation_snp.set_index('SNP_ID_CURRENT',drop=0)
    
#    annotation_snp=annotation_snp.iloc[np.where(annotation_snp.index.values==annotation_snp.index.values)[0]]; print (annotation_snp.shape)
#    annotation_snp['gene_name']=np.array([np.array(re.split(',| - ',g)) for g in annotation_snp['MAPPED_GENE'].values.astype('U')]); print (annotation_snp.shape)
#    annotation_snp['gene_name']=np.array([g[np.array(['LOC' not in gg for gg in g])] for g in annotation_snp['gene_name']]); print (annotation_snp.shape)
#    annotation_snp['number_genes']=np.array([len(g) for g in annotation_snp['gene_name']]); print (annotation_snp.shape)
    annotation_snp['index']=annotation_snp['SNP_ID_CURRENT'].astype('U')+'_'+annotation_snp[mapped_trait].astype('U')
    #    plt.hist(np.clip(-np.log10(annotation_snp['P-VALUE'].astype(float)+10**-100),0,10))
    
    
#    print (np.unique(annotation_snp.index[annotation_snp.index==annotation_snp.index].astype('U')).shape)
#    for c in annotation_snp.columns: 
#        print ("\n>>>")
#        print (annotation_snp[[c]].head())
    try:
        print (annotation_snp.shape)
    #    annotation_snp=annotation_snp[annotation_snp['OR or BETA']==annotation_snp['OR or BETA']]
        print (annotation_snp.shape)
        annotation_snp['Study_size']=[sum_ints(text=text) for text in annotation_snp['INITIAL SAMPLE SIZE'].values]
        annotation_snp=annotation_snp[annotation_snp['Study_size']>1000]
    except:1
#    annotation_snp[['Study_size','INITIAL SAMPLE SIZE']]
#    annotation_snp[annotation_snp['Study_size']<1000][['Study_size','INITIAL SAMPLE SIZE']]
#    annotation_snp.loc['11_116663707',['Study_size','INITIAL SAMPLE SIZE']].iloc[0].loc['INITIAL SAMPLE SIZE']
#    
    
    
    annotation_snp['MAPPED_TRAIT']=annotation_snp[mapped_trait].astype('U')
    annotation_snp=annotation_snp.iloc[np.unique(annotation_snp['index'],return_index=1)[1]]; print (annotation_snp.shape)
    annotation_snp['SNP_pos']=[i.split('_')[1] for i in annotation_snp['SNP_ID_CURRENT']]
    annotation_snp['SNP_chrom']=[i.split('_')[0] for i in annotation_snp['SNP_ID_CURRENT']]
    try:
        annotation_snp=annotation_snp[['SNP_ID_CURRENT','SNP_pos','SNP_chrom','MAPPED_TRAIT','P-VALUE','OR or BETA','INITIAL SAMPLE SIZE','SNPS','Study_size']]
    except:annotation_snp=annotation_snp[['SNP_ID_CURRENT','SNP_pos','SNP_chrom','MAPPED_TRAIT']]


#    annotation_snp[[]].to_csv('/Users/mirauta/Results/hipsci/QTL_june2018/gwas_snps.txt')
    return annotation_snp

#
def update_gwas_38_37():
    ann=pd.read_table('/Users/mirauta/Data/Genotypes/'+'gwas_catalog_v1.0-associations_e95_r2019-03-01.tsv',sep='\t',encoding='latin_1')
    
    file=pd.DataFrame(columns=['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO'])
    file['#CHROM']=ann[ 'CHR_ID']
    file['POS']=ann['CHR_POS']
    file['ID']=['i'+i for i in np.arange(file.shape[0]).astype('U')]
    file['REF']='N'
    file['ALT']='N'
    file['QUAL']=1
    file['FILTER']=1
    file['INFO']=1
    file=file[file['POS']==file['POS']]
    file=file[["x" not in x for x in file['POS']]]
    file['#CHROM']=[j.split(';')[0] for j in file['#CHROM']]
    file['POS']=[j.split(';')[0] for j in file['POS']]
    file.to_csv('/Users/mirauta/Data/Genotypes/gwas_catalog_v1.0-associations_e95_r2019-03-01_snps_to_ensembl.vcf',index=0,sep='\t')
    update_pos=pd.read_table('/Users/mirauta/Data/Genotypes/gwas_38_37.vcf',skiprows=4,sep='\t',encoding='latin_1').set_index('ID')
    old_pos=pd.read_table('/Users/mirauta/Data/Genotypes/gwas_catalog_v1.0-associations_e95_r2019-03-01_snps_to_ensembl.vcf',sep='\t',encoding='latin_1').set_index('ID')
    old_pos['CHR_ID_HG37']=update_pos.loc[old_pos.index,'#CHROM']
    old_pos['CHR_POS_HG37']=update_pos.loc[old_pos.index,'POS']
    old_pos= old_pos[np.isfinite(old_pos['CHR_POS_HG37'])]
    old_pos['Combined38']=old_pos[ '#CHROM']+'_'+old_pos['POS'].astype('U')
    old_pos['Combined']=old_pos['CHR_ID_HG37']+'_'+old_pos['CHR_POS_HG37'].astype('int').astype('U')
    old_pos=old_pos.set_index('Combined38',drop=0)
    old_pos=old_pos.iloc[np.unique(old_pos.index,return_index=1)[1]]
    ann['Combined38']=ann[ 'CHR_ID']+'_'+ann['CHR_POS']
    ann=ann.set_index('Combined38',drop=0)
    for c in ['CHR_ID_HG37','Combined','CHR_POS_HG37']: ann[c]=old_pos.loc[ann.index,c]
    ann.to_csv('/Users/mirauta/Data/Genotypes/'+'gwas_catalog_v1.0-associations_e95_r2019-03-01_37_coordinates.tsv',sep='\t',encoding='latin_1')

def process_ld_chromomes():
    for chrom in np.arange(1,23):
        print (chrom)
        
        ld=pandas.read_table('/Users/mirauta/Data/Annotation/Ensembl_37.75/r08/chr'+str(chrom),sep=' ',header=None)      
        ld.columns=['pos','ld_pos'];         ld=ld.set_index('pos',drop=0)
        ld['chromosome']=chrom
        ld['ld_pos']=ld['pos'].astype('U')+','+ld['ld_pos']
        ld.index=ld['chromosome'].astype('U')+':'+ld['pos'].astype('U')
        ld.to_csv('/Users/mirauta/Data/Annotation/Ensembl_37.75/r08/chr'+str(chrom)+'_processed.txt',sep='\t')    
        
        
        
def get_pep_info(data,gene,df,trait):
        temp=data['Peptide'][gene]['data']
        qtlpeps=np.array(list(temp['snp_id'].keys()))
        temp2=np.array([temp['p_value_raw'][pep][:][temp['snp_id'][pep][:].astype('U')==df[trait].loc[gene,'snp_id']]<0.01 for pep in qtlpeps])
        qtlpeps2=qtlpeps[[False if len(t)==0 else t[0] for t in temp2]]
        return qtlpeps2

def get_subfeature_qtl(data,gene,df11,ann,peptides):
        snp=df11.loc[gene,'snp_id']
        temp=data['Transcript'][gene]['data']
      
        qtlpeps=np.array(list(temp['snp_id'].keys()))
        rez=pd.DataFrame(index=qtlpeps,columns=['beta','p_value','ensembl_transcript_id'])
#        temp2=np.array([temp['p_value_raw'][pep][:][]
        for pep in qtlpeps:
            index=temp['snp_id'][pep][:].astype('U')==snp
            rez.loc[pep,'ensembl_transcript_id']=pep.split('_')[0]
            try:
                rez.loc[pep,'gene_beta']=df11.loc[gene,'beta']
                rez.loc[pep,'gene_replicated_beta']=df11.loc[gene,'replicated_beta']
                rez.loc[pep,'gene_p_value_raw']=df11.loc[gene,'p_value_raw']
                rez.loc[pep,'gene_replicated_p_value_raw']=df11.loc[gene,'replicated_p_value_raw']
                rez.loc[pep,'beta']=temp['beta'][pep][:][index][0]
                rez.loc[pep,'p_value_raw']=temp['p_value_raw'][pep][:][index]
                rez.loc[pep,'empirical_feature_p_value']=temp['empirical_feature_p_value'][pep][:][index]
                rez.loc[pep,'transcript_type']=ann.loc[[rez.loc[pep,'ensembl_transcript_id']],'transcript_type'].iloc[0]
            except:
                rez.loc[pep,'transcript_type']='not_protein'
        rez.index=rez['ensembl_transcript_id']
        try:
            peptidesrez=peptides.loc[rez.index];#peptidesrez=peptidesrez.set_index('peptide',drop=0);
#        peptidesrez=peptidesrez.set_index('ensembl_isoform_id',drop=0);
#        peptidesrez=peptidesrez.loc[rez.index];   
#        peptidesrez=peptidesrez[peptidesrez.index==peptidesrez.index]
#        temp=data['Peptide'][gene]['data']
#        peptidesrez['peptide_beta']=np.array([temp['beta'][pep][:][temp['snp_id'][pep][:].astype('U')==snp][0] for pep in peptidesrez.index]).astype('U')
#        peptidesrez=peptidesrez.set_index('ensembl_isoform_id',drop=0);peptidesrez=peptidesrez.loc[rez.index];        
            rez['peptide']=[sum(peptidesrez.loc[[tr],'peptide']==peptidesrez.loc[[tr],'peptide']) for tr in rez.index]
            rez['peptide2']=[';'.join(peptidesrez.loc[[tr],'peptide'][peptidesrez.loc[[tr],'peptide']==peptidesrez.loc[[tr],'peptide']]) for tr in rez.index]
#        rez['peptide_beta']=[';'.join(peptidesrez.loc[[tr],'peptide_beta'][peptidesrez.loc[[tr],'peptide']==peptidesrez.loc[[tr],'peptide']]) for tr in rez.index]
        except:1
        try: 
            rez1=rez[rez['gene_beta']*rez['beta']>0]        
            
#        print(rez[['peptide2','beta','transcript_type','peptide','gene_beta']])
            return [rez,rez1.loc[np.argmin(rez1['p_value_raw'])]]
        except: return [rez,rez.iloc[0]]

def replicate_isoform(path_data,data,ann0,df,proteinlabel,prometa):
    ind=0
    genes=df['Transcript_'+proteinlabel].index
    rez=df['Transcript_'+proteinlabel].copy(deep=True)
    rez['replicated_p_value_raw_protein']=np.zeros(rez.shape[0])+np.nan
    rez['replicated_beta_protein']=np.zeros(rez.shape[0])+np.nan
    rez['feature_id_protein']=np.zeros(rez.shape[0]).astype('U')
    for gene in genes:
        tr=df['Transcript_'+proteinlabel].loc[gene,'feature_id']
        snp=df['Transcript_'+proteinlabel].loc[gene,'snp_id']
        try: 
            leadpr=prometa.loc[ann0.loc[tr.split('_')[0],'protein_uniprot_id'],'leading']

            x=data['Transcript'][gene]['data']
            print (ann0.loc[[xx.split('_')[0]for xx in list(x['snp_id'].keys())],'protein_uniprot_id'])
            x['p_value_raw'][leadpr][:][x['snp_id'][leadpr][:].astype('U')==snp]
             
            x=data[proteinlabel][gene]['data']
            rez.loc[gene,'replicated_p_value_raw_protein']=\
               x['p_value_raw'][leadpr][:][x['snp_id'][leadpr][:].astype('U')==snp]
            rez.loc[gene,'replicated_beta_protein']=\
               x['beta'][leadpr][:][x['snp_id'][leadpr][:].astype('U')==snp]
            df['Transcript_'+proteinlabel].loc[gene,'p_value_raw']
            rez.loc[gene,'protein_feature_id']=leadpr
            ind=ind+1
        except:1
#        if ind==1111:sys.exit()
    rez.to_csv(path_data+'Transcript_'+proteinlabel+'_isoform_specific_qtl_results.txt', sep='\t',index=False)

#for proteinlabel in ['Protein6','Protein7M']:replicate_isoform(path_data,data,ann0,df,proteinlabel,prometa)
#rez1=pandas.read_table(path_data+'Transcript_'+proteinlabel+'_isoform_specific_qtl_results.txt', sep='\t').set_index('ensembl_gene_id')
##rez1.head()
##plt.plot(rez['replicated_beta'],rez['replicated_beta_protein'],'o')
#
#rez1[rez1['protein_feature_id']!=rez1['protein_feature_id']]

def load_cis_for_trans(trans,df,anntrall,anntrcoding,prometa,proteinlabel):
    print ('load cis')
    ''' get the cis type'''

    trans=trans.set_index('feature_id',drop=0)
    

    trans['feature_position']=prometa.loc[trans['feature_id'],'start']
    trans['feature_chromosome']=prometa.loc[trans['feature_id'],'chromosome']
    trans['feature_chromosome'][trans['feature_chromosome']=='HSCHR6_MHC_DBB']='6'
    
    trans['feature_chromosome'][trans['feature_chromosome']=='HSCHR6_MHC_SSTO']='6'         
    trans['feature_chromosome'][trans['feature_chromosome']=='HSCHR6_MHC_MCF']='6'
    trans['feature_chromosome'][trans['feature_chromosome']=='HG1322_PATCH']='6'
    trans['feature_chromosome'][trans['feature_chromosome']=='HG1436_HG1432_PATCH']='X'
         
    
    trans['snp_chromosome']=np.array([s.split('_')[0]for s in trans['snp_id']])
#trans['snp_chromosome'][trans['snp_chromosome']=='HSCHR6_MHC_DBB']='6'
    trans['snp_position']=np.array([s.split('_')[1]for s in trans['snp_id']])
    trans['trans_is_cis']=(trans['snp_chromosome']==trans['feature_chromosome'])&(abs(trans['snp_position'].astype(float)-trans['feature_position'])<10**6)
    trans['trans_is_cis'][(~np.in1d(trans['feature_chromosome'],np.arange(23).astype('U')))&((abs(trans['snp_position'].astype(float)-trans['feature_position'])<10**6))]=True
    trans=trans.set_index('snp_id',drop=0)
    

    trans['cis_protein']=np.in1d(trans.index,df[proteinlabel+'_'+proteinlabel]['snp_id'])
    trans['cis_rna']=np.in1d(trans.index,df['mRNA_mRNA']['snp_id'])
    
    trans['cis_protein_rep']=np.in1d(trans.index,df[proteinlabel+'_mRNA'][df[proteinlabel+'_mRNA']['replicated_bonf']<0.01]['snp_id'])
    trans['cis_protein_only']=trans['cis_protein']&(~trans['cis_protein_rep'])
    
    trans['cis_rna_rep']=np.in1d(trans.index,df['mRNA_'+proteinlabel][df['mRNA_'+proteinlabel]['replicated_bonf']<0.01]['snp_id'])
    trans['cis_rna_notassesed']=np.in1d(trans.index, np.setdiff1d(df['mRNA_mRNA'] ['snp_id'],df['mRNA_'+proteinlabel] ['snp_id']))
    trans['cis_rna_only']=trans['cis_rna']&(~trans['cis_rna_rep'])&(~trans['cis_rna_notassesed'])
         
    
    trans['cis_protein_ensembl_id']=df[proteinlabel+'_'+proteinlabel].set_index('snp_id').iloc[np.unique(df[proteinlabel+'_'+proteinlabel]['snp_id'],return_index=1)[1]].loc[trans.index]['ensembl_gene_id']
#    df[proteinlabel+'_'+proteinlabel].set_index('snp_id').iloc[np.unique(df[proteinlabel+'_'+proteinlabel]['snp_id'],return_index=1)[1]].loc[trans.index]['ensembl_gene_id']
    trans['cis_rna_ensembl_id']=df['mRNA_mRNA'].set_index('snp_id').iloc[np.unique(df['mRNA_mRNA']['snp_id'],return_index=1)[1]].loc[trans.index]['ensembl_gene_id']
    trans['cis_rna_ensembl_id'][trans['cis_rna_ensembl_id']!=trans['cis_rna_ensembl_id']]=''
    trans['cis_protein_ensembl_id'][trans['cis_protein_ensembl_id']!=trans['cis_protein_ensembl_id']]=''
    
    for f in ['cis_protein','cis_rna']:
        trans=trans.set_index(f+'_ensembl_id',drop=0)
        trans[f+'_gene_name']=anntrall.loc[trans[f+'_ensembl_id'],'gene_name']
        trans[f+'_gene_name'][trans[f+'_gene_name']!=trans[f+'_gene_name']]=trans[f+'_ensembl_id']
        
    trans=trans.set_index('ensembl_gene_id',drop=0)
    trans['gene_name']=anntrall.loc[trans['ensembl_gene_id'],'gene_name']    
    trans['gene_name'][trans['gene_name']!=trans['gene_name']]=trans['ensembl_gene_id']
    trans=trans.set_index('feature_id',drop=0)
    idx=trans['ensembl_gene_id']=='0.0'
    trans['gene_name'][idx]=prometa.loc[trans[idx].index]['superior_feature_id']
    trans=trans.set_index('ensembl_gene_id',drop=0)
    trans['cis_crna_gene_name']=trans['cis_rna_gene_name'].copy(deep=True)
    trans['cis_ncrna_gene_name']=trans['cis_rna_gene_name'].copy(deep=True)
 
    
    trans=trans.set_index('snp_id',drop=0)

         
    trans['cis_protein_snp_id']=df[proteinlabel+'_'+proteinlabel].set_index('snp_id',drop=0).iloc[np.unique(df[proteinlabel+'_'+proteinlabel]['snp_id'],return_index=1)[1]].loc[trans.index]['snp_id']
    trans['cis_rna_snp_id']=df['mRNA_mRNA'].set_index('snp_id',drop=0).iloc[np.unique(df['mRNA_mRNA']['snp_id'],return_index=1)[1]].loc[trans.index]['snp_id']
    
    
    trans['cis_crna']=np.in1d(trans['cis_rna_ensembl_id'],anntrcoding.index)
    trans['cis_ncrna']=~np.in1d(trans['cis_rna_ensembl_id'],np.hstack([[''],anntrcoding.index]))
    #trans['cis_protein_beta']=df[proteinlabel+'_'+proteinlabel].set_index('snp_id').iloc[np.unique(df[proteinlabel+'_'+proteinlabel]['snp_id'],return_index=1)[1]].loc[trans.index]['beta']
    
    
    trans['cis_ncrna_ensembl_id']=trans['cis_rna_ensembl_id'];trans['cis_crna_ensembl_id']=trans['cis_rna_ensembl_id']
    
    trans=trans.set_index('snp_id',drop=0)
    trans['cis_protein_uniprot_id']=df[proteinlabel+'_'+proteinlabel].set_index('snp_id').iloc[np.unique(df[proteinlabel+'_'+proteinlabel]['snp_id'],return_index=1)[1]].loc[trans.index]['feature_id']
    trans['cis_protein_uniprot_id'][trans['cis_protein_uniprot_id']!=trans['cis_protein_uniprot_id']]=''
    trans['cis_protein_gene_name']=df[proteinlabel+'_'+proteinlabel].set_index('snp_id').iloc[np.unique(df[proteinlabel+'_'+proteinlabel]['snp_id'],return_index=1)[1]].loc[trans.index]['gene_name']
    trans['cis_ensembl_gene_id']=trans['cis_protein_ensembl_id'].copy(deep=True)
    trans['cis_ensembl_gene_id'][trans['cis_protein_ensembl_id']=='']=trans['cis_rna_ensembl_id'][trans['cis_protein_ensembl_id']=='']
    trans['cis_gene_name']=trans['cis_protein_gene_name'].copy(deep=True)
    trans['cis_gene_name'][(trans['cis_protein_gene_name']=='')|(trans['cis_protein_gene_name']!=trans['cis_protein_gene_name'])]=trans['cis_rna_gene_name'][(trans['cis_protein_gene_name']=='')|(trans['cis_protein_gene_name']!=trans['cis_protein_gene_name'])]
 
    #trans['cis_rna_uniprot_id']=np.nan
    #trans['cis_rna_uniprot_id'][trans['cis_rna_ensembl_id'].values==trans['cis_rna_ensembl_id'].values]=ann.loc[trans[trans['cis_rna_ensembl_id']==trans['cis_rna_ensembl_id']]['cis_rna_ensembl_id']]['protein_uniprot_id']
    
    
    trans['cis_protein_p_value']=df[proteinlabel+'_'+proteinlabel].set_index('snp_id').iloc[np.unique(df[proteinlabel+'_'+proteinlabel]['snp_id'],return_index=1)[1]].loc[trans.index]['p_value_raw']
    trans['cis_rna_p_value']=df['mRNA_mRNA'].set_index('snp_id').iloc[np.unique(df['mRNA_mRNA']['snp_id'],return_index=1)[1]].loc[trans.index]['p_value_raw']
    trans['cis_rna_c_p_value']=trans['cis_rna_p_value'];trans['cis_ncrna_p_value']=trans['cis_rna_p_value']
    
    trans['cis_protein_beta']=df[proteinlabel+'_'+proteinlabel].set_index('snp_id').iloc[np.unique(df[proteinlabel+'_'+proteinlabel]['snp_id'],return_index=1)[1]].loc[trans.index]['beta']
    trans['cis_rna_beta']=df['mRNA_mRNA'].set_index('snp_id').iloc[np.unique(df['mRNA_mRNA']['snp_id'],return_index=1)[1]].loc[trans.index]['beta']
    trans['cis_rna_c_beta']=trans['cis_rna_beta'];trans['cis_ncrna_beta']=trans['cis_rna_beta']; 
    
    trans['cis_protein_replicated_beta']=df[proteinlabel+'_mRNA'].set_index('snp_id').iloc[np.unique(df[proteinlabel+'_mRNA']['snp_id'],return_index=1)[1]].loc[trans.index]['replicated_beta']
    trans['cis_rna_replicated_beta']=df['mRNA_'+proteinlabel].set_index('snp_id').iloc[np.unique(df['mRNA_'+proteinlabel]['snp_id'],return_index=1)[1]].loc[trans.index]['replicated_beta']
    trans['cis_rna_c_replicated_beta']=trans['cis_rna_replicated_beta'];trans['cis_ncrna_replicated_beta']=trans['cis_rna_replicated_beta']
    trans['index']=trans['snp_id']+'_'+trans['ensembl_gene_id']
    
    idx = trans.groupby('ensembl_gene_id').apply(lambda df:df.empirical_feature_p_value.argmin());idx=idx.values+'_'+idx.index
    leadtrans=trans.set_index('index',drop=0).loc[idx]
    leadtrans['qv']=scst2.fdrcorrection0(leadtrans['empirical_feature_p_value'])[1]
 
    idx = trans[trans['cis_protein']].groupby('ensembl_gene_id').apply(lambda df:df.empirical_feature_p_value.argmin());idx=idx.values+'_'+idx.index
    leadtransprotein=trans[trans['cis_protein']].set_index('index',drop=0).loc[idx]
    leadtransprotein['qv']=scst2.fdrcorrection0(leadtransprotein['empirical_feature_p_value'])[1]   
    
    idx = trans[trans['cis_rna']].groupby('ensembl_gene_id').apply(lambda df:df.empirical_feature_p_value.argmin());idx=idx.values+'_'+idx.index
    leadtransrna=trans[trans['cis_rna']].set_index('index',drop=0).loc[idx]
    leadtransrna['qv']=scst2.fdrcorrection0(leadtransrna['empirical_feature_p_value'])[1]
    
#    trans['qv']=np.nan
#    x=trans
#    extended=np.hstack([x['p_value'].values, np.random.uniform(0.1,1,np.unique(x['snp_id']).shape[0]*   np.unique(trans['feature_id']).shape[0]-x.shape[0])])
#    x['qv']=scst2.fdrcorrection0(extended)[1][:x.shape[0]]
#    temp=x[x['qv']<0.5]
#    trans.loc[temp.index, 'qv']=x['qv'].loc[temp.index]
    
#    transhigh=transhigh.set_index('feature_id',drop=0)
#    transhigh=pd.concat([transhigh.loc[[gene]].iloc[np.argmin(transhigh.loc[[gene],'empirical_feature_p_value'].values)].T for gene in transhigh['feature_id']],1).T
#    
        
    return [trans,leadtrans,leadtransprotein,leadtransrna]

def nanunique(x):
    return np.unique(x[x==x])


def boxplot(snp_df,snp_id,y,color,fig = plt.figure(figsize=(4,4)),ii0=1,i=1,ii=1 ):

    ##    sb.swarmplot(y=y[snp==0],color=colors[itr*2],label=tr.split('_')[0]);  
    y=y.loc[snp_df.index]
    snp=snp_df[snp_id].values
    alleles=snp_id.split('_')[2:];alleles=[alleles[0]+alleles[0],alleles[0]+alleles[1],alleles[1]+alleles[1]] 
    
    #    plt.title(gene_name)
    fig.patch.set_facecolor('white')
    ax2 = fig.add_subplot(ii0,ii,i);    ax2.spines['right'].set_visible(False); ax2.spines['top'].set_visible(False); #ax1.yaxis.set_ticks_position('right'); 
    ax2.plot([y[snp==s].mean() for s in np.unique( snp )],'_',color=color,markersize=20) 
    sb.swarmplot(y=y,x= snp,color=color )
#    sb.boxplot(y=y,x= snp,color=color)

    plt.xticks(np.arange(np.unique(snp).shape[0]),alleles[:np.unique(snp).shape[0]])
    ax2.legend(loc=1,fontsize=8)   


def plots_trans(x,cis, mrna_df_ensembl,peerrna,protein_df,peer,snp_df):
    for ix in np.arange(x.shape[0])[:]:
        print (ix)
        try:
            x1=x.iloc[ix]
            snp_id=x1['snp_id'] 
            y=np.log(protein_df.loc[x1['feature_id']])
            y=regxfromy(peer,y)[0]
            fig=plt.figure(figsize=(16,4))
            boxplot(snp_df,snp_id,y=y,color='steelblue',fig=fig,i=1)
            plt.ylabel(x1['gene_name']+' (trans Protein)')
            try:plt.title(x1['cis_rna_annotation'].split('[')[0][:55]+'\n'+x1['gene_name']+' ('+x1['feature_chromosome']+') <<-- '+x1['cis_gene_name']+\
                      ' ('+x1['snp_chromosome']+')')
            except:1
            
            y2= mrna_df_ensembl.loc[x1['cis_ensembl_gene_id']]          
            y2=regxfromy(peerrna,y2)[0]
            boxplot(snp_df,snp_id,y=y2,color='grey',fig=fig,i=2)
            plt.ylabel(x1['cis_gene_name'] )
            ax2 = fig.add_subplot(1,4,3); 
            plt.plot(y2,y,'o')
            plt.xlabel(x1['cis_gene_name']+'cis RNA')
            plt.ylabel(x1['gene_name'])
            plt.plot(y2[snp_df[snp_id]==0].values,y[snp_df[snp_id]==0].values,'go')
            plt.plot(y2[snp_df[snp_id]>1].values,y[snp_df[snp_id]>1].values,'ro')
            plt.plot(y2[snp_df[snp_id]==1].values,y[snp_df[snp_id]==1].values,'mo')
            plt.tight_layout()
                
            y22= mrna_df_ensembl.loc[x1['ensembl_gene_id']]
            y22=regxfromy(peerrna,y22)[0]
            ax4 = fig.add_subplot(1,4,4); 
            boxplot(snp_df,snp_id,y=y22,color='forestgreen',fig=fig,i=4)
            plt.ylabel(x1['gene_name']+'(trans mRNA)')
#            ax4 = fig.add_subplot(1,4,4);
#            plt.hist([y2[y!=y],y2[y==y]],normed=0,label=['missing protein','detected'])
#            plt.xlabel(x1['cis_rna_gene_name'] +'cis mRNA')
#            plt.legend(loc=1)
            plt.tight_layout()    
        except: 
            print (ix)
            print(x1['gene_name']+' regulated by '+x1['cis_rna_gene_name'])
        
        plt.savefig('/Users/mirauta/Results/hipsci/manuscript_images/Manhattan/Cis_trans_boxplots/'+cis+'_'+x1['gene_name']+'_'+x1['cis_gene_name']+'.png',dpi=600)
        plt.savefig('/Users/mirauta/Results/hipsci/manuscript_images/Manhattan/Cis_trans_boxplots/'+cis+'_'+x1['gene_name']+'_'+x1['cis_gene_name']+'.svg',dpi=600)
        plt.gcf().clear()
        
        

def correlation(x,y):
    return np.corrcoef(x[(x+y)==(x+y)],y[(x+y)==(x+y)])[0,1]


def retrieve_extracellular_domains(membrane,pepmeta):
    print ('retrieve teh exteracellular infor for peptides: distance from extracellular domains margins; positive means inside extracellular')
    pepmeta['extracellular']=np.nan
    for pr in membrane.index:
         temp=np.array(membrane.loc[pr,'Topological domain'].split(';'))
         temp=temp[['Extracellular' in t for t in temp]]
         if len(temp)==0: continue
         extra=[np.array(t.replace(' TOPO_DOM ','').replace('TOPO_DOM ','').split(' Extracellular')[0].split(' ')).astype(int) for t in temp]
         pep=pepmeta[pepmeta['superior_feature_id']==pr]
         pep['len']=[len(p) for p in pep.index]
         pep['end_in_protein']=pep['start_in_protein']+pep['len']
         for reg in extra:
             pepmeta.loc[pep.index,'extracellular']=np.nanmax(np.vstack([pepmeta.loc[pep.index,'extracellular'].values,\
                                                              np.vstack([pep['start_in_protein']-reg[0],     reg[1]-pep['end_in_protein']]).min(0)]),0)
         
    return pepmeta     

def bring_data_trans(x,df,df2,ann1):
    ann1=ann1.set_index('ensembl_gene_id',drop=0)
    x['index']=x.index;
    x=x.set_index('snp_id',drop=0) 
    df1=df.iloc[np.unique(df.index,return_index=1)[1]]        
    x[['cis_p_value','cis_ensembl_gene_id','cis_protein','cis_protein_beta','cis_gene_name']]=df1.loc[x['snp_id'],['empirical_feature_p_value','ensembl_gene_id','feature_id','beta', 'gene_name']]
    df21=df2.iloc[np.unique(df2.index,return_index=1)[1]]        
    x[['cis_replicated_p_value_mRNA']]=df1.loc[x['snp_id'],[ 'replicated_bonf' ]]
    x=x.set_index('ensembl_gene_id',drop=0)
    x['trans_gene_name']=ann1.loc[x['ensembl_gene_id'],'gene_name']
 
    x=x.set_index('snp_id',drop=0)
#    x['pair']=np.array([xx.split('-')[0]+"_" for xx in x['feature_id']]).astype('U')
#    x['pair']=x['pair']+np.array([xx.split('-')[0]for xx in x['cis_protein']])
#    ''' retreive PPI info'''
#    x['PPI_complex']=np.in1d(x['pair'],complex_pairs['pair'])
    x=x.set_index('index',drop=0) 
    return(x)


def plot_cofactors(cofactors, cofactors2,proteinlabel,qtl,field='sign',field2='logpv',log_field2=True,field2_threshold=2,thrs=10.0**-np.linspace(1,6,num=6),\
                  method='logistic',title='eQTL_replication_cofactors',xlim=None):
    import warnings
    warnings.filterwarnings('ignore')
 

    allprotein=pd.read_table('/Users/mirauta/Data/MS/hipsci/phenotypes/hipsci.proteomics.maxquant.uniprot.TMT_batch_24.20170919_protein_Reporter intensity_ranknorm_regbatch_valid_lines_recurrent_genes.txt',index_col=0)
    replicated=allprotein[allprotein.columns[np.array(['bubh_3' in b for  b in allprotein.columns])]]

    try:
        qtl['coef. of error']=qtl['tech_noise2'].copy(deep=True)# (np.array([ np.nanstd(replicated.loc[pr])/np.nanmean(replicated.loc[pr]) for pr in qtl['replicated_feature_id'].values]))
        qtl['protein CV']=qtl['cv_intensity'].copy(deep=True)
        qtl['protein mean']=qtl['mean_intensity'].copy(deep=True)
        qtl['MPC']=qtl['CORUM complex'].copy(deep=True)
        qtl['number of peptides']=qtl['countpep2'].copy(deep=True)
        
        qtl['coef. of error'],regcoef=regout(A=qtl,x=['protein mean','number of peptides'], y='coef. of error')
        qtl['protein CV'],regcoef=regout(A=qtl,x=['protein mean','coef. of error'], y='protein CV')  
        qtl['protein mean'],regcoef=regout(A=qtl,x=['number of peptides'], y='protein mean')
        qtl['sign']=(qtl['replicated_bonf']<10**-2)&(qtl['beta']*qtl['replicated_beta']>0)
        qtl=qtl[np.isfinite(np.sum(qtl[cofactors+['sign']].values.astype(float),1))]
    except:1

    plt.figure(figsize=(8,4))
    for i,indexqtl in enumerate([np.isfinite(qtl['beta'])]):        
        a=plt.subplot(1,2,1);    a.spines['top'].set_visible(False);a.spines['right'].set_visible(False);
        for ithr,thr in enumerate(thrs):
            print (thr)
            index=np.isfinite(np.sum(qtl[cofactors+[field]].values.astype(float),1))&(qtl['qv']<thr); 
            index=index&indexqtl                   
            X=qtl[cofactors][index].values.astype(float)
            y=qtl[field][index].values
            if method =='logistic':
                X=np.array([x>np.nanmean(x) for ix,x in enumerate(X.T)]).T
                #reg = stsmdiscrete.Logit(y>=y.mean(), X); result = reg.fit_regularized() coef=result.params
                reg = LogisticRegression(random_state=0,solver='lbfgs',penalty='l2');
                reg = LogisticRegression(random_state=0,solver='saga',penalty='l1');
                
                result = reg.fit(X=X, y=y>=y.mean());                 coef=result.coef_[0]
            print (coef)
            plt.barh(0.75+ithr*0.075+np.arange(len(coef)+1),np.hstack([coef,0]),  align='center', alpha=1,height=0.5/len(thrs),color='tan')

        plt.plot((0,0),(0.5,len(coef)),'k-',lw=0.5)
        plt.yticks(0.75+np.arange(len(coef)+1), np.hstack([cofactors,'all genes']),rotation=00,fontsize=9)
        if xlim is not None : plt.xlim((-xlim ,xlim ))
       # plt.yticks(np.linspace(0.8,1.4,num=7)-1,np.linspace(0.8,1.4,num=7))
        plt.xlabel("log odds ratio",fontsize=16)
        plt.legend(loc=1,fontsize=8)
    
    a=plt.subplot(1,2,2);    a.spines['top'].set_visible(False);a.spines['right'].set_visible(False);a.spines['left'].set_visible(False)
    plt.yticks([],[])
    
    
    for ithr,thr in enumerate(thrs):
        index=np.isfinite(np.sum(qtl[cofactors+[field]].values.astype(float),1))&(qtl['qv']<thr); 
        index=index&indexqtl
        y=qtl[field][index].values
        X=qtl[cofactors][index].values.astype(float);X=np.array([x>np.nanmean(x) for ix,x in enumerate(X.T)]).T
        plt.barh(0.75+ithr*0.075+np.arange(len(coef)+1),np.hstack([X.sum(0),y.shape[0]]),  align='center', alpha=1,height=0.5/len(thrs),color='palegreen')
        plt.barh(0.75+ithr*0.075+np.arange(len(coef)+1),np.hstack([((X*y[:,None]).sum(0)),y.sum()]),  align='center', alpha=1,height=0.5/len(thrs),color='forestgreen')
        #plt.barh(0.75+ithr*0.075+np.arange(len(coef)+1),np.hstack([-(1-X).sum(0),-y.shape[0]]),  align='center', alpha=1,height=0.5/len(thrs),color='palegreen', hatch='////')
       # plt.barh(0.75+ithr*0.075+np.arange(len(coef)+1),np.hstack([-(((~X)*y[:,None]).sum(0)),-y.sum()]),  align='center', alpha=1,height=0.5/len(thrs),color='forestgreen', hatch='////')
        tempx=np.hstack([((X).sum(0)),y.shape[0]]);
        tempxx=np.hstack([((X*y[:,None]).sum(0)),y.sum()])
        tempxm=np.hstack([((1-X).sum(0)),-y.shape[0]]);
        tempxxm=np.hstack([((~X*y[:,None]).sum(0)),y.sum()])
    #for it, t in enumerate(tempxxm): plt.annotate(s=str(np.around(tempxxm[it]/tempxm[it]*100,1))+'%' ,xy=(-tempxxm[it]-500,it+0.65)) 
    for it, t in enumerate(tempxx): plt.annotate(s=str(np.around(tempxx[it]/tempx[it]*100,1))+'%' ,xy=(tempxx[it]-0.5,it+0.65)) 
        
    plt.plot((0,0),(0.5,len(coef)+1),'k-',lw=0.5)
    plt.yticks([],[])
    plt.xlabel('number of genes',fontsize=16)
    plt.savefig("/Users/mirauta/Results/hipsci/manuscript_images/"+title+".png",bbox_inches='tight',dpi=600)
    plt.savefig("/Users/mirauta/Results/hipsci/manuscript_images/"+title+".svg",bbox_inches='tight',dpi=600)
    plt.show()

    return [X,y,coef,qtl]

#from scipy import interpolate
def estimateqv(pv, pi0=None):
 
    m=pv.shape[0]
    if pi0 is None:pi0=0.99
    psort = sc.argsort(pv)
    pv = pv[psort]
    
    qv =pv+np.nan
    
    
    for i in np.arange(m-2,-1,-1):
        qv[i] =  pi0*m*pv[i]/(i+1.0)
    qv_temp = qv.copy()
    qv = sc.zeros_like(qv)
    qv[psort] = qv_temp

    return qv 
#rr=np.random.uniform(0,1,10**5)
#a0=np.hstack([ 10**-8,10**-7 , rr  ])
#pv=np.hstack([ 10**-7,10**-8,np.random.uniform(10**-7.98,10**-6.995,100), rr  ])
#qv=estimateqv(pv, pi0=None)
#scst.spearmanr(pv[:11115],qv[:11115])
#
#print(estimateqv(a0)[:3])
#print(estimateqv(a)[:3])
#
#print((estimateqv(a0)<0.1).sum())
#print((estimateqv(a)<0.1).sum())
#print((a<0.001).sum())

def force_normal_distribution(phenotype, method='gaussnorm', reference=None):
    _doc='rank transform x into ref/ gaussian;keep the range; keep ties'
    
    if method=='log':
        return np.log(1+phenotype)
    
    if method=='log_standardize':
        temp=np.log(1+phenotype)
        return (temp-np.nanmean(temp))/np.nanstd(temp)

    if method=='arcsin':
        return np.arcsin(np.sqrt(phenotype))
    
    if method=='arcsin_standardize':
        temp=np.arcsin(np.sqrt(phenotype))
        return (temp-np.nanmean(temp))/np.nanstd(temp)
    
    if method=='standardize':
        return (phenotype-np.nanmean(phenotype))/np.nanstd(phenotype)
                        
    indextoupdate = np.isfinite(phenotype)
    y1 = phenotype[indextoupdate]
    yuni,yindex=np.unique(y1, return_inverse=True)
    phenotypenorm=phenotype.copy()
    
    if method =='gaussnorm':

        sref = scst.norm.isf(np.linspace(0.001, 0.999,num=yuni.shape[0])[::-1])
        phenotypenorm[indextoupdate]=sref[yindex]
        return phenotypenorm
    
    elif method=='ranknorm':
        try:
            xref1=np.unique(reference[np.isfinite(reference)])
            sref=np.sort(xref1)[np.linspace(0,xref1.shape[0]-0.001, num=y1.shape[0]).astype(int)]
        except:
            print ('reference missing. provide reference to force_normal_distribution or choose gaussnorm')
            return 1
        phenotypenorm[indextoupdate]=sref[np.argsort(np.argsort(y1))]
        return phenotypenorm
    
    elif method=='ranknorm_duplicates':
        try:
            xref1=np.unique(reference[np.isfinite(reference)])### unique values from reference
            sref=np.sort(xref1)[np.linspace(0,xref1.shape[0]-0.001, num=yuni.shape[0]).astype(int)]
        except: 
            print ('reference missing. provide reference to force_normal_distribution or choose gaussnorm')
            return 1
        phenotypenorm[indextoupdate]=sref[yindex]
        return phenotypenorm  
    
    else:
        print ('methods are: log, log_standardize, standardize, gaussnorm, ranknorm, ranknorm_duplicates, arcsin, arcsin_standardize')





def get_ld(qtl):
    qtl['chromosome_snp']=np.array([ s.split('_')[0] for s in  qtl['snp_id']]).astype(int)
    qtl['position_snp']=np.array([ s.split('_')[1] for s in  qtl ['snp_id']]).astype(int)
    qtl=qtl.set_index('snp_id',drop=0)
    qtl['ld_snps']=np.zeros(qtl.shape[0],dtype='object')
#    print (np.unique(qtl['chromosome_snp']))
    for chrom in np.unique(qtl['chromosome_snp']):
        print (chrom)
        ld=pandas.read_table('/Users/mirauta/Data/Annotation/Ensembl_37.75/r08/chr'+str(chrom),sep=' ',header=None)
        ld.columns=['pos','ld_pos']
        ld=ld.set_index('pos',drop=0)
    
        temp=qtl[qtl['chromosome_snp']==chrom]
        temp['position_snp']=temp['position_snp'].astype(int)
#        temp=temp[['position_snp','ld_snps','ensembl_gene_id','snp_id']].set_index('position_snp',drop=0)
        temp=temp[['position_snp','ld_snps','snp_id']].set_index('position_snp',drop=0)
        temp=temp.iloc[np.unique(temp.index,return_index=1)[1]]
        temp['ld_pos']=ld.loc[temp.index,'ld_pos']
        temp['ld_pos'][temp['ld_pos']==temp['ld_pos']]=\
            temp['position_snp'][temp['ld_pos']==temp['ld_pos']].astype('U')+','+temp['ld_pos'][temp['ld_pos']==temp['ld_pos']]
        temp['ld_pos'][temp['ld_pos']!=temp['ld_pos']]=temp['position_snp'][temp['ld_pos']!=temp['ld_pos']].astype('U')
        temp=temp.set_index('snp_id',drop=0)['ld_pos']
#        print (temp.shape)
#        print (temp.loc[np.intersect1d(temp.index,qtl.index)])
        qtl.loc[temp.index,'ld_snps']=temp.loc[np.intersect1d(temp.index,qtl.index)]
#        print ("finished")
    return qtl

def get_gwas(qtl,annotation_snp,pos_field='CHR_POS_HG37',chrom_field='CHR_ID_HG37',name=''):
    
    chroms=np.unique(annotation_snp[chrom_field])
    annotation_snp[pos_field]=annotation_snp[pos_field].astype('U')
    qtl['chromosome_snp']=qtl['chromosome_snp'].astype('U')
    for chrom in chroms:
        annchrom=annotation_snp[annotation_snp[chrom_field]==chrom]
        qtlchrom=qtl.query('chromosome_snp==@chrom')
        for ig,gene in enumerate(qtlchrom.index):
            qtlgene=qtlchrom.loc[gene]
            snps=qtlgene['ld_snps'].split(',')
            
        #    snps['rna']=qtl.loc[gene,'rna_snp_id_ldblock'].split(';')
        #    qtl['protein_gwas_trait']np.intersect1d(snps['protein'],snps['rna'])
            temp=np.intersect1d(snps,annchrom[pos_field].values)
            if len(temp)==0: continue
            
            if ig%100==0: print (ig)
            temp_ann=annchrom[np.in1d(annchrom[pos_field].values,temp)]
            qtl.loc[gene,'gwas_snp'+name]=';'.join(temp_ann.index)
            qtl.loc[gene,'gwas_snp_rsid'+name]=';'.join(temp_ann['SNPS'])
            try:
#                print (gene)
                try:qtl.loc[gene,'gwas_trait'+name]= ';'.join(temp_ann['MAPPED_TRAITS'])
                except: qtl.loc[gene,'gwas_trait'+name]= ';'.join(temp_ann['MAPPED_TRAIT'])
                qtl.loc[gene,'gwas_p_value'+name]= ';'.join(temp_ann['P-VALUE'])
               
                try:qtl.loc[gene,'PUBMEDID'+name]= ';'.join( temp_ann['PUBMEDID'].astype('U')) 
                except:1
            except: 1
    return qtl 

def write_snp_snap(qtl1,path_data,name,snp_field='snp_id'):

    qtl1['snp_id_dots']=[':'.join(s.split('_')[:2]) for s in qtl1 [snp_field]]
    
    pd.DataFrame(qtl1 ['snp_id_dots']). to_csv(path_data+name+"_snp_lists_SNPsnap.txt",header=False,index=False)
#sys.exit()   
def plot_peptides_linear2(fig,L,peps,snppos,pavpos,qtl1,g):
 

    ax = plt.subplot(1,1,1)
    ax.spines['top'].set_visible(True); ax.spines['right'].set_visible(False);    ax.spines['bottom'].set_visible(False);    ax.xaxis.set_ticks_position('none');  ax.xaxis.set_ticklabels(''); 
    for pe in peps.query('beta<0').index: ax.add_patch(Rectangle((peps.loc[pe,'start_in_protein'], 0), peps.loc[pe,'peplen'], (peps.loc[pe,'beta']), facecolor='darkgoldenrod',    alpha=1))             
    for pe in peps.query('beta>0').index: ax.add_patch(Rectangle((peps.loc[pe,'start_in_protein'], 0), peps.loc[pe,'peplen'], (peps.loc[pe,'beta']), facecolor='skyblue',    alpha=1))             
    for p in pavpos:plt.plot((p,p),(0,(qtl1.loc[g,'protein_eff_size'])) ,'--',color='k',lw=1)#    l, b, w, h = ax.get_position().bounds
    ax.spines['bottom'].set_position('zero')
#    ax.set_position([l*2, b*2, w*2, h])
    ax.plot(np.arange(L),np.repeat((qtl1.loc[g,'protein_eff_size']),L),'--',color="saddlebrown" if qtl1.loc[g,'protein_eff_size']<0 else "steelblue",lw=1)
    
    
    

    plt.plot((0,L),(0,0),'k-',lw=0.5)
#    ax.spines['left'].set_position(('axes', 0.6))
#    ax.spines['right'].set_color('none')
    ax.spines['bottom'].set_position(('axes', min(0,min(peps['beta'].min(),qtl1.loc[g,'protein_eff_size']))))
#    ax.spines['top'].set_color('none')
#    ax.spines['left'].set_smart_bounds(True)
#    ax.spines['bottom'].set_smart_bounds(True)
##    ax.xaxis.set_ticks_position('bottom')
#    ax.yaxis.set_ticks_position('left')
#    ax.invert_yaxis()
#    ax.xaxis.set_ticks_position('top')
#    ax.xaxis.tick_top()
#    ax.invert_yaxis()
    
#    pos1 = ax.get_position() # get the original position 
#    print (pos1)
#    pos2 = [pos1.x0 + 0.3, pos1.y0 +0.3,  pos1.width / 2.0, pos1.height / 5.0] 
#    ax.set_position(pos2)
#    pos=[0, peps['beta'].min(), L,  peps['beta'].max()] 
#    ax.set_position(pos, which='both')
    
    plt.ylabel('Effect size',fontsize=16)
    plt.title('Position in protein '+qtl1.loc[g,'feature_id'],fontsize=16)
    plt.tight_layout()

#fig=plt.figure(figsize=(5,2))
##        plt.title(gene_name)
##        plot_peptides(fig=fig,L=L, peps=peps,snppos=None,pavpos=pozs)
#
#plot_peptides_linear2(fig=fig,L=L,peps=peps,snppos=None,pavpos=pozs,qtl1=qtl1,g=g)
#
#numdata = 100
#t = np.linspace(0, 100, numdata)
#y = 1/t**(1/2.0)
#
#fig = plt.figure()
#ax = fig.add_subplot(1, 1, 1)
#ax.xaxis.set_ticks_position('top')
#ax.yaxis.grid(linestyle = '-', color = 'gray')
#
#ax.plot(t, y, 'g-', linewidth = 1.5)
#
#plt.show()
def plot_peptides_linear3(fig,L,peps,snppos,pavpos,qtl1,g):
 

    ax = plt.subplot(1,1,1)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False);    ax.spines['bottom'].set_visible(True);    ax.xaxis.set_ticks_position('none');  ax.xaxis.set_ticklabels(''); 
    shift=min(peps['beta'].min(),qtl1.loc[g,'protein_eff_size'])-0.05
    for pe in peps.query('beta<0').index: 
        ax.add_patch(Rectangle((peps.loc[pe,'start_in_protein'], -shift), peps.loc[pe,'peplen'], (peps.loc[pe,'beta']), facecolor='darkgoldenrod',    alpha=1))             
    for pe in peps.query('beta>0').index: 
        ax.add_patch(Rectangle((peps.loc[pe,'start_in_protein'], -shift), peps.loc[pe,'peplen'], (peps.loc[pe,'beta']), facecolor='skyblue',    alpha=1))             
    for p in pavpos:plt.plot((p,p),(-shift,(qtl1.loc[g,'protein_eff_size'])-shift) ,'--',color='k',lw=1)#    l, b, w, h = ax.get_position().bounds
    ax.spines['bottom'].set_position('zero')
#    ax.set_position([l*2, b*2, w*2, h])
    plt.plot(np.arange(L),np.repeat((qtl1.loc[g,'protein_eff_size']),L)-shift,'--',color="saddlebrown" if qtl1.loc[g,'protein_eff_size']<0 else "steelblue",lw=1)
#    plt.plot((0,L),(-shift,-shift),'k-',lw=1)
    rr=np.around(np.arange(min(peps['beta'].min(),qtl1.loc[g,'protein_eff_size'])-0.1*(qtl1.loc[g,'protein_eff_size']>0),max(peps['beta'].max(),qtl1.loc[g,'protein_eff_size'])+0.1*(qtl1.loc[g,'protein_eff_size']<0),0.1),1)    
    
    plt.yticks(rr-shift,rr)


#    pos1 = ax.get_position() # get the original position 
#    print (pos1)
#    pos2 = [pos1.x0 + 0.3, pos1.y0 +0.3,  pos1.width / 2.0, pos1.height / 5.0] 
#    ax.set_position(pos2)
#    pos=[0, peps['beta'].min(), L,  peps['beta'].max()] 
#    ax.set_position(pos, which='both')
    
    plt.ylabel('Effect size',fontsize=16)
    plt.xlabel('Position in protein '+qtl1.loc[g,'feature_id'],fontsize=16)
    plt.tight_layout()
#fig=plt.figure(figsize=(5,2))
#plot_peptides_linear2(fig=fig,L=L,peps=peps,snppos=None,pavpos=pozs,qtl1=qtl1,g=g)

def plot_peptides_linear(fig,L,peps,snppos,pavpos,qtl1,g):
 

    ax = plt.subplot(1,1,1)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False);    ax.spines['bottom'].set_visible(True);    ax.xaxis.set_ticks_position('none');  ax.xaxis.set_ticklabels(''); 
    for pe in peps.query('beta<0').index: ax.add_patch(Rectangle((peps.loc[pe,'start_in_protein'], 0), peps.loc[pe,'peplen'], abs(peps.loc[pe,'beta']), facecolor='darkgoldenrod',    alpha=1))             
    for pe in peps.query('beta>0').index: ax.add_patch(Rectangle((peps.loc[pe,'start_in_protein'], 0), peps.loc[pe,'peplen'], abs(peps.loc[pe,'beta']), facecolor='skyblue',    alpha=1))             
    for p in pavpos:plt.plot((p,p),(0,abs(qtl1.loc[g,'protein_eff_size'])) ,'--',color='k',lw=1)#    l, b, w, h = ax.get_position().bounds
    ax.spines['bottom'].set_position('zero')
#    ax.set_position([l*2, b*2, w*2, h])
    plt.plot(np.arange(L),np.repeat(abs(qtl1.loc[g,'protein_eff_size']),L),'--',color="saddlebrown" if qtl1.loc[g,'protein_eff_size']<0 else "steelblue",lw=1)
    
    
    plt.ylabel('Effect size',fontsize=16)
    plt.xlabel('Position in protein '+qtl1.loc[g,'feature_id'],fontsize=16)
    plt.tight_layout()

    
def plot_peptides(fig,L,peps,snppos,pavpos):
    
    CL=330/L
#    patches2 = [  Wedge((.0, .0), 1.05, 180+CL*(i-0.15),180+CL*(i+0.15) ,width=0.05) for i in np.linspace(0,L,10)]
    patches = [  Wedge((.0, .0), 1.025, 180+CL*0,180+CL*L  ,width=0.07)]
    patches = patches+[  Wedge((.0, .0), 0.975, 180+CL*0,180+CL*L  ,width=0.01)]
    if snppos is not None:patches = patches+[  Wedge((.0, .0), 0.95, 180+CL*(snppos-0.25),180+CL*(snppos+0.25) ,width=0.15)]
    patches = patches+[  Wedge((.0, .0), 0.95, 180+CL*(v-0.25),180+CL*( v+0.25) ,width=0.15)for v in pavpos]
    
    patches2 = [  Wedge((0,0), .9, 180+CL*peps.loc[pe,'start_in_protein'],180+CL*(peps.loc[pe,'stop_in_protein'] ), width=(peps.loc[pe,'beta'])) for pe in peps.query('beta>0').index]              # Full sector]
    patches3 = [  Wedge((0,0), .9, 180+CL*peps.loc[pe,'start_in_protein'],180+CL*(peps.loc[pe,'stop_in_protein'] ), width=-(peps.loc[pe,'beta'])) for pe in peps.query('beta<0').index]              # Full sector]
    

    
    if snppos is None: 
        colors=np.array(['grey','grey', 'k'])
        colorsi = np.repeat(2,len(patches))
        colorsi[:2]=np.arange(2)
        
    else:
        colors=np.array(['grey','grey','darkred','k'])
        colorsi = np.repeat(3,len(patches))
        colorsi[:3]=np.arange(3)
        
    
           
     
    #    plt.plot((-1,0),(0,0),'k-',lw=0.5)
    #    plt.annotate('0',xy=(-1,0.0510),fontsize=8);plt.annotate('-1',xy=(0,0.051),fontsize=8)
#    fig=plt.figure(figsize=(4,4))
#    plt.title(gene_name)
    ax = plt.subplot(1,1,1)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False);    ax.spines['bottom'].set_visible(False); ax.spines['left'].set_visible(False); 
    ax.yaxis.set_ticks_position('none');     ax.xaxis.set_ticks_position('none');  ax.yaxis.set_ticklabels('');  ax.xaxis.set_ticklabels(''); 
    p = PatchCollection(patches, alpha=1)
    #p.set_array(np.array(colors))
    p.set_color(colors[colorsi])
    ax.add_collection(p)
    p2 = PatchCollection(patches2, alpha=1)
    p2.set_color('skyblue')
    ax.add_collection(p2)
    
    p3 = PatchCollection(patches3, alpha=1)
    p3.set_color('lightcoral')
    ax.add_collection(p3)   

    #plt.plot(sc.sin(z),sc.cos(z),'-',color='grey');
    #fig.colorbar(p, ax=ax)
    #plt.xlim
    plt.xlim(-1.2,1.2)
    plt.ylim(-1.2,1.2)
#    plt.show()
#
#g='ENSG00000227184'
#g='ENSG00000213965'
#gene=g
#gene_name=qtl['pQTL'].loc[g,'gene_name']
#protein=qtl['pQTL'].loc[g,'feature_id']
#peps=pepmeta.query('superior_feature_id==@protein')[['start_in_protein']]
#peps=peps.loc[np.intersect1d(list(peph5[g]['data/p_value_raw'].keys()),peps.index)]
#peps.sort_values('start_in_protein',inplace=True)
#peps=peps[[len(np.where(peph5[g]['data/snp_id'][pe][:].astype('U')==qtl1.loc[g,'snp_id'])[0])>0 for pe in peps.index]]
#peps['pv']=[peph5[g]['data/p_value_raw'][pe][np.where(peph5[g]['data/snp_id'][pe][:].astype('U')==qtl1.loc[g,'snp_id'])[0][0]]for pe in peps.index]
#peps['beta']=[peph5[g]['data/beta'][pe][np.where(peph5[g]['data/snp_id'][pe][:].astype('U')==qtl1.loc[g,'snp_id'])[0][0]]for pe in peps. index]
#peps['stop_in_protein']=[peps.loc[pe,'start_in_protein']+len(pe)  for pe in peps.index] 
#
#print (ig)
#print (g)
#print (gene_name)
#print(peps['beta'])
#
#L=len(records[np.where(refseqids==qtl1[qtl1['gene_name']==gene_name]['feature_id'].values[0])[0][0]].seq)
#snppos=1;pavpos=[1]
#
#CL=330/L
##    patches2 = [  Wedge((.0, .0), 1.05, 180+CL*(i-0.15),180+CL*(i+0.15) ,width=0.05) for i in np.linspace(0,L,10)]
#patches = [  Wedge((.0, .0), 1.025, 180+CL*0,180+CL*L  ,width=0.07)]
#patches = patches+[  Wedge((.0, .0), 0.975, 180+CL*0,180+CL*L  ,width=0.01)]
#patches = patches+[  Wedge((.0, .0), 0.95, 180+CL*(snppos-0.25),180+CL*(snppos+0.25) ,width=0.15)]
#patches = patches+[  Wedge((.0, .0), 0.95, 180+CL*(v-0.25),180+CL*( v+0.25) ,width=0.15)for v in pavpos]
#
#patches2 = [  Wedge((0,0), .9, 180+CL*peps.loc[pe,'start_in_protein'],180+CL*(peps.loc[pe,'stop_in_protein'] ), width=(peps.loc[pe,'beta'])) for pe in peps.query('beta>0').index]              # Full sector]
#patches3 = [  Wedge((0,0), .9, 180+CL*peps.loc[pe,'start_in_protein'],180+CL*(peps.loc[pe,'stop_in_protein'] ), width=-(peps.loc[pe,'beta'])) for pe in peps.query('beta<0').index]              # Full sector]
#
#colorsi = np.repeat(4,len(patches))
#colors=np.array(['grey','grey','darkred','k'])
#colorsi[:4]=np.arange(4)
#       
# 
##    plt.plot((-1,0),(0,0),'k-',lw=0.5)
##    plt.annotate('0',xy=(-1,0.0510),fontsize=8);plt.annotate('-1',xy=(0,0.051),fontsize=8)
#fig=plt.figure(figsize=(4,4))
#plt.title(gene_name)
#ax = plt.subplot(1,1,1)
#ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False);    ax.spines['bottom'].set_visible(False); ax.spines['left'].set_visible(False); 
#ax.yaxis.set_ticks_position('none');     ax.xaxis.set_ticks_position('none');  ax.yaxis.set_ticklabels('');  ax.xaxis.set_ticklabels(''); 
#p = PatchCollection(patches, alpha=1)
##p.set_array(np.array(colors))
#p.set_color(colors[colorsi])
#ax.add_collection(p)
#p2 = PatchCollection(patches2, alpha=1)
#p2.set_color('skyblue')
#ax.add_collection(p2)
#
#p3 = PatchCollection(patches3, alpha=1)
#p3.set_color('lightcoral')
#ax.add_collection(p3)
#
#plt.xlim(-1.2,1.2)
#plt.ylim(-1.2,1.2)

 

inframe_variant=np.array(['inframe_deletion','inframe_insertion','incomplete_terminal_codon_variant','stop_lost','stop_gained','missense_variant']) 
frameshift_variant=np.array(['frameshift_variant,feature_elongation','feature_truncation']) 
UTR_variant=np.array(['UTR_variant','3_prime_UTR_variant','5_prime_UTR_variant']) 
intron_variant=np.array(['intron_variant']) 
splicing_variant=np.array(['splice_acceptor_variant','splice_donor_variant','splice_region_variant','splicing_variant']) 

variants={'inframe_variant': np.array(['inframe_deletion','inframe_insertion','incomplete_terminal_codon_variant','stop_lost','stop_gained','missense_variant']),
'frameshift_variant':np.array(['frameshift_variant,feature_elongation','feature_truncation']), 
'UTR_variant':np.array(['UTR_variant','3_prime_UTR_variant','5_prime_UTR_variant']),
'intron_variant':np.array(['intron_variant']),'synonymous_variant':np.array(['synonymous_variant']),
'splicing_variant':np.array(['splice_acceptor_variant','splice_donor_variant','splice_region_variant','splicing_variant'])}

plot_variants=np.array(['intron_variant','UTR_variant','splicing_variant',   'synonymous_variant', 'frameshift_variant','inframe_variant'])


def plot_fisher_exact_results(fig,qtl1,it,its=2,thr=0.1,trait='pQTL',xlim=(-0.20,6),ylim=(-0.20,3.5)):
#    try:
    qtl1['VE_LD1']=['no proxy gene variants' not in p for p in qtl1['VE_LD']]
    al0=np.hstack(qtl1['VE_LD'][qtl1['VE_LD1']].values)
    rep0=np.hstack(qtl1['VE_LD'][(qtl1['VE_LD1'])&(qtl1['replicatedet'])].values)
    nonrep0=np.hstack(qtl1['VE_LD'][(qtl1['VE_LD1'])&(~qtl1['replicatedet'])].values)
#    except:
#        al0=np.hstack(qtl1['VE_LD'][(qtl1['VE_LD']!=0)].values)
#        rep0=np.hstack(qtl1['VE_LD'][(qtl1['VE_LD']!=0)&(qtl1['replicatedet'])].values)
#        nonrep0=np.hstack(qtl1['VE_LD'][(qtl1['VE_LD']!=0)&(~qtl1['replicatedet'])].values)
#    print (qtl1['replicatedet'].shape)
#    print (qtl1['replicatedet'].sum())
#    print (qtl1['replicatedet'].shape[0]-qtl1['replicatedet'].sum())
    for x in [nonrep0,rep0,al0]:       
        x[np.in1d(x,inframe_variant)]='inframe_variant';
        x[np.in1d(x,['feature_elongation','feature_truncation'])]='frameshift_variant';
        x[np.in1d(x,['splice_acceptor_variant','splice_donor_variant','splice_region_variant'])]='splicing_variant';
        x[np.in1d(x,['3_prime_UTR_variant','5_prime_UTR_variant'])]='UTR_variant';
        x[x=='stop_retained_variant']='synonymous_variant'
        x[x=='coding_sequence_variant']='0'
#        x[x=='splice_variant']='0'
        x[x=='intergenic_variant']='0'
        x[x=='initiator_codon_variant']='0'
    for iix,x in enumerate([nonrep0,rep0,al0]):
        xx=np.array([nonrep0,rep0,al0])[iix]
        x=xx[xx!='0']
        x=x[x!=0]
   
    x=np.unique(rep0,return_counts=1);xrep={};
    for ix in np.arange(x[1].shape[0]): xrep[x[0][ix]]=x[1][ix]
    try:del xrep['0']
    except:1
    x=np.unique(nonrep0,return_counts=1);xnonrep={};
    for ix in np.arange(x[1].shape[0]):  xnonrep[x[0][ix]]=x[1][ix]
    try:del xnonrep['0']
    except:1
    x=np.unique(al0,return_counts=1);xal={};
    for ix in np.arange(x[1].shape[0]): xal[x[0][ix]]=x[1][ix]
    try:del xal['0']
    except:1
    for k in np.setdiff1d(list(xrep.keys()),list(xnonrep.keys())):    xnonrep[k]=0   
    for k in np.setdiff1d(list(xnonrep.keys()),list(xrep.keys())):   xrep[k]=0   
                         
    al=np.array([xal[k] for k  in plot_variants])
    rep=np.array([xrep[k] for k  in plot_variants])
    nonrep=np.array([xnonrep[k] for k  in plot_variants])
 
    ax = fig.add_subplot(1,its,it+1);ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False); ax.yaxis.set_ticks_position('left'); 
    temp=np.array([scst.fisher_exact(np.array([[xnonrep[k],nonrep.sum()-xnonrep[k]],[xal[k],al.sum()-xal[k]]]),alternative='two-sided')[1] for k  in plot_variants])
    tempval=np.array([scst.fisher_exact(np.array([[xnonrep[k],nonrep.sum()-xnonrep[k]],[xal[k],al.sum()-xal[k]]]),alternative='two-sided')[0] for k  in plot_variants])
    plt.plot(-np.log10(temp),tempval,'o')
    plt.xlim((-np.log10(temp)).min()-(-np.log10(temp)).max()*0.1,(-np.log10(temp)).max()*1.3)
    for iis, s in enumerate(temp):
        if s<thr:
            st=plot_variants[iis].replace('_','\n').replace('variant','')
            st=st[:1].capitalize()+st[1:]
            plt.annotate(xy=(-np.log10(temp)[iis]+0.2 ,tempval[iis]-0.4),s=st,fontsize=14)
    plt.ylabel('Fisher test: oddsratio',fontsize=16)
    plt.xlabel('Fisher test: -log10 PV',fontsize=16)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.plot(xlim,(1,1),'k--',lw=0.5)

def plot_fisher_exact_results_text_number_genes(fig,qtl1,it,its=2,thr=0.1,trait='pQTL',xlim=(-0.20,6),ylim=(-0.20,3.5)):
#    try:
    qtl1['VE_LD1']=['no proxy gene variants' not in p for p in qtl1['VE_LD']]
    qtl1=qtl1[qtl1['VE_LD1']]
    
  
             
#        x[np.in1d(x,['feature_elongation','feature_truncation'])]='frameshift_variant';
#        x[np.in1d(x,['splice_acceptor_variant','splice_donor_variant','splice_region_variant'])]='splicing_variant';
#        x[np.in1d(x,['3_prime_UTR_variant','5_prime_UTR_variant'])]='UTR_variant';
#        x[x=='stop_retained_variant']='synonymous_variant'
#        x[x=='coding_sequence_variant']='0'
##        x[x=='splice_variant']='0'
#        x[x=='intergenic_variant']='0'
#        x[x=='initiator_codon_variant']='0'
#    for var in variants.keys():
#        variants[var]
#        qtl1[var]=[False if qtl1['VE_LD'].iloc[n]==0 else len(np.intersect1d(variants[var],  qtl1['VE_LD'].iloc[n]))>0 for n in np.arange(qtl1.shape[0])]
#    
#    
                        
    test={}
    for k  in plot_variants:
        test[k]=scst.fisher_exact(np.array([[qtl1.query('not replicatedet')[k].sum(),(~qtl1.query('not replicatedet')[k]).sum()],\
                                               [qtl1.query('replicatedet')[k].sum(),(~qtl1.query('replicatedet')[k]).sum()]]),alternative='two-sided')    
    temp=np.array([test[k][1] for k  in plot_variants])
    tempval=np.array([test[k][0] for k  in plot_variants])

    ax = fig.add_subplot(1,its,it+1);ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False); ax.yaxis.set_ticks_position('left'); 
    
    plt.plot(-np.log10(temp),tempval,'o')
    plt.xlim((-np.log10(temp)).min()-(-np.log10(temp)).max()*0.1,(-np.log10(temp)).max()*1.3)
    for iis, s in enumerate(temp):
        if s<thr:
            st=plot_variants[iis].replace('_','\n').replace('variant','')
            st=st[:1].capitalize()+st[1:]
            plt.annotate(xy=(-np.log10(temp)[iis]+0.2 ,tempval[iis]-0.4),s=st,fontsize=14)
    plt.ylabel('Fisher test: oddsratio',fontsize=16)
    plt.xlabel('Fisher test: -log10 PV',fontsize=16)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.plot(xlim,(1,1),'k--',lw=0.5)





def plot_fisher_exact_bar(fig,qtl1,it,thr=0.3,trait='pQTL'):

    try:
        al0=np.hstack(qtl1['VE_LD'][(qtl1['VE_LD']!='0')].values)
        rep0=np.hstack(qtl1['VE_LD'][(qtl1['VE_LD']!='0')&(qtl1['replicatedet'])].values)
        nonrep0=np.hstack(qtl1['VE_LD'][(qtl1['VE_LD']!='0')&(~qtl1['replicatedet'])].values)
    except:
        al0=np.hstack(qtl1['VE_LD'][(qtl1['VE_LD']!=0)].values)
        rep0=np.hstack(qtl1['VE_LD'][(qtl1['VE_LD']!=0)&(qtl1['replicatedet'])].values)
        nonrep0=np.hstack(qtl1['VE_LD'][(qtl1['VE_LD']!=0)&(~qtl1['replicatedet'])].values)

    for x in [nonrep0,rep0,al0]:
        x[np.in1d(x,['inframe_deletion','inframe_insertion','incomplete_terminal_codon_variant','stop_lost','stop_gained','missense_variant'])]='inframe_variant';
        x[np.in1d(x,['feature_elongation','feature_truncation'])]='frameshift_variant';
        x[np.in1d(x,['splice_acceptor_variant','splice_donor_variant','splice_region_variant'])]='splicing_variant';
        x[np.in1d(x,['downstream_gene_variant','upstream_gene_variant'])]='intergenic_variant';
#        x[np.in1d(x,['3_prime_UTR_variant','5_prime_UTR_variant'])]='UTR_variant';
        x[np.in1d(x,['3_prime_UTR_variant'])]='UTR_variant';
        x[x=='stop_retained_variant']='synonymous_variant'
        x[x=='coding_sequence_variant']='0'
        x[x=='intergenic_variant']='0'
        x[x=='splice_variant']='0'
        x[x=='initiator_codon_variant']='0'
        x=x[np.in1d(x,np.unique(x,return_counts=1)[0][np.unique(x,return_counts=1)[1]>20])]

 
    x=np.unique(rep0,return_counts=1);xrep={};  
    for ix in np.arange(x[1].shape[0]): xrep[x[0][ix]]=x[1][ix]
    try:del xrep['0']
    except:1
    x=np.unique(al0,return_counts=1);xal={};
    for ix in np.arange(x[1].shape[0]): xal[x[0][ix]]=x[1][ix]
    try:del xal['0']
    except:1
    x=np.unique(nonrep0,return_counts=1);xnonrep={};
    for ix in np.arange(x[1].shape[0]):  xnonrep[x[0][ix]]=x[1][ix]
    try:del xnonrep['0']
    except:1
    for k in np.setdiff1d(list(xrep.keys()),list(xnonrep.keys())):    xnonrep[k]=0   
    for k in np.setdiff1d(list(xnonrep.keys()),list(xrep.keys())):   xrep[k]=0   

    al=np.array([xal[k] for k  in plot_variants])
    rep=np.array([xrep[k] for k  in plot_variants])
    nonrep=np.array([xnonrep[k] for k  in plot_variants])
    ax = fig.add_subplot(1,2,it+1);ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False); ax.yaxis.set_ticks_position('left'); 
    plt.bar(height=rep,left=np.arange(len(plot_variants)),color='forestgreen' if trait=='eQTL' else "steelblue",label='is eQTL or tQTL' if trait=='pQTL' else "is pQTL")
    plt.bar(height=nonrep,left=np.arange(len(plot_variants)),bottom=rep,color='palegreen' if trait=='eQTL' else "skyblue",\
            label='is not eQTL or tQTL' if trait=='pQTL' else "is not pQTL")
    plt.ylabel('Number of QTLs',fontsize=16)
    plt.xticks(np.arange(len(plot_variants))-0.0005,[p.replace('_variant','')[:1].capitalize()+p.replace('_variant','')[1:].replace('_','\n') for p in plot_variants],rotation=45) 
    plt.xlabel('Variant effect',fontsize=16)
    plt.legend(loc=1,fontsize=11)

#fig = plt.figure(figsize=(8, 4))
#fig.patch.set_facecolor('white')
#for it,trait in enumerate(['pQTL','eQTL']):
#    qtl1=qtl[trait]
#    plot_fisher_exact_bar(fig,qtl1,it,thr=0.3,trait=trait)
#
#plt.tight_layout()

def get_ld(qtl):
    qtl['chromosome_snp']=np.array([ s.split('_')[0] for s in  qtl['snp_id']]).astype(int)
    qtl['position_snp']=np.array([ s.split('_')[1] for s in  qtl ['snp_id']]).astype(int)
    qtl=qtl.set_index('snp_id',drop=0)
    qtl['ld_snps']=np.zeros(qtl.shape[0],dtype='object')
#    print (np.unique(qtl['chromosome_snp']))
    for chrom in np.unique(qtl['chromosome_snp']):
        print (chrom)
        ld=pandas.read_table('/Users/mirauta/Data/Annotation/Ensembl_37.75/r08/chr'+str(chrom),sep=' ',header=None)
        ld.columns=['pos','ld_pos']
        ld=ld.set_index('pos',drop=0)
    
        temp=qtl[qtl['chromosome_snp']==chrom]
        temp['position_snp']=temp['position_snp'].astype(int)
        temp=temp[['position_snp','ld_snps','ensembl_gene_id','snp_id']].set_index('position_snp',drop=0)
        temp=temp.iloc[np.unique(temp.index,return_index=1)[1]]
        temp['ld_pos']=ld.loc[temp.index,'ld_pos']
        temp['ld_pos'][temp['ld_pos']==temp['ld_pos']]=\
            temp['position_snp'][temp['ld_pos']==temp['ld_pos']].astype('U')+','+temp['ld_pos'][temp['ld_pos']==temp['ld_pos']]
        temp['ld_pos'][temp['ld_pos']!=temp['ld_pos']]=temp['position_snp'][temp['ld_pos']!=temp['ld_pos']].astype('U')
        temp=temp.set_index('snp_id',drop=0)['ld_pos']
#        print (temp.shape)
#        print (temp.loc[np.intersect1d(temp.index,qtl.index)])
        qtl.loc[temp.index,'ld_snps']=temp.loc[np.intersect1d(temp.index,qtl.index)]
#        print ("finished")
    return qtl



def get_VE_type(ve):
  
    VE={'inframe':np.array(['inframe_deletion','inframe_insertion','incomplete_terminal_codon_variant','stop_lost','stop_gained','missense_variant']) }
    VE['frameshift']=np.array(['feature_elongation','feature_truncation','frameshift_variant'])
    VE['splicing']=np.array(['splice_acceptor_variant','splice_donor_variant','splice_region_variant'])
    VE['UTR']=np.array(['3_prime_UTR_variant','5_prime_UTR_variant'])
    VE['intron']=np.array(['intron_variant'])
    VE['NMD_transcript']=np.array(['NMD_transcript_variant'])
    VE['synonymous']=np.array(['synonymous_variant','start_retained_variant','stop_retained_variant'])
    vet=ve['MAPPED_TRAITS']
    for f in VE.keys():
        print (f)
        ve[f]=[any([(vet.iloc[i].find(ff)!=-1) for ff in VE[f]]) for i in np.arange(vet.shape[0])]
    return (ve)    

def get_VE(qtl,snp_VE,peptide_transcript1):
  
    try:qtl['chromosome_gene']
    except:qtl['chromosome_gene']=qtl['chromosome_snp']
        
    VE={'inframe':np.array(['inframe_deletion','inframe_insertion','incomplete_terminal_codon_variant','stop_lost','stop_gained','missense_variant']) }
    VE['frameshift']=np.array(['feature_elongation','feature_truncation','frameshift_variant'])
    VE['splicing']=np.array(['splice_acceptor_variant','splice_donor_variant','splice_region_variant'])
    VE['UTR']=np.array(['3_prime_UTR_variant','5_prime_UTR_variant'])
    VE['intron']=np.array(['intron_variant'])
    VE['NMD_transcript']=np.array(['NMD_transcript_variant'])
    VE['synonymous']=np.array(['synonymous_variant','start_retained_variant','stop_retained_variant'])

    qtl['VE_LD']=np.zeros(qtl.shape[0],dtype='object')
    for f in VE.keys():
        qtl[f]=np.repeat('',qtl.shape[0])


    
        
    for n in np.arange(qtl.shape[0]):
        if n%10==0: print (n)
#        print (n)
        pozld=np.array(qtl['ld_snps'].iloc[n].split(',')).astype(int)
    
        pozld=pozld[((pozld>qtl['gene_start'].iloc[n])*(pozld<qtl['gene_end'].iloc[n]))|((pozld<qtl['gene_start'].iloc[n])*(pozld>qtl['gene_end'].iloc[n]))]
        if len(pozld)==0: continue
    
        snp_VE_gene=snp_VE.fetch(qtl['chromosome_gene'].iloc[n], pozld.min(),pozld.max())
        temp=np.array([])
        for i in pozld:
            for reader in snp_VE_gene.fetch(qtl['chromosome_gene'].iloc[n], i-1,i):
                try: reader.INFO['VE']
                except: continue
                variants_ve=np.array([r.split('|')[0] for r in reader.INFO['VE']])
                try:variants_poz=np.in1d([r.split('|')[3] for r in reader.INFO['VE']], peptide_transcript1.loc[[qtl['ensembl_gene_id'].iloc[n]],'ensembl_isoform_id'].values)
                except:variants_poz=np.arange(len(variants_ve))
                
                temp0=variants_ve[variants_poz]
                variant_types=np.hstack([f if len(np.intersect1d(VE[f],temp0))>0 else '' for f in VE.keys()])
                for f in np.unique(variant_types[variant_types!='']):
                    qtl[f].iloc[n]=str(reader.CHROM)+'_'+str(reader.POS)+'_'+str(reader.REF)+'_'+str(reader.ALT[0])+';'+qtl[f].iloc[n]
                temp=np.hstack([temp,temp0])
#                  except:1
        qtl['VE_LD'].iloc[n]=';'.join(np.unique(temp))
        qtl['VE_LD'][qtl['VE_LD']==0]='no proxy gene variants'
        
    return qtl

def write_peptide_agreement(qtl1, folder_data_raw, peph5file='/Users/mirauta/Results/hipsci/QTL_may2019/param_peptide_maxquant_peer_log_7peer_divided_new_selected_ensembl_gene_id_qtl_results_genome.h5',
    proh5file='/Users/mirauta/Results/hipsci/QTL_may2019/param_protein_maxquant_peer_log_newpqtl_ensembl_gene_id_qtl_results_genome.h5'):
    
    proh5=h5py.File(proh5file,'r')
    peph5=h5py.File(peph5file,'r')
    pepmeta=pandas.read_table(folder_data_raw + 'hipsci.proteomics.maxquant.uniprot.TMT_batch_24.20170919_peptide_meta_recurrent_genes.txt',sep='\t',index_col=0)
#    pepmetaall=pandas.read_table(folder_data_raw + 'hipsci.proteomics.maxquant.uniprot.TMT_batch_24.20170919_peptide_meta_all_genes.txt',sep='\t',index_col=0)
    
    overlap=pd.read_table('/Users/mirauta/Data/MS/hipsci/phenotypes/may_2019/notmapped_overlap_peptide_snp_idSNP_multiple_junction_phase_30.01.txt',sep='\t',index_col=0)


    pepmeta=pepmeta.loc[np.setdiff1d(pepmeta.index,overlap['feature_id'])]
    pepquant=pandas.read_table(folder_data_raw + 'hipsci.proteomics.maxquant.uniprot.TMT_batch_24.20170919_peptide_Reporter intensity_ranknorm_regbatch_valid_lines_recurrent_genes.txt',sep='\t',index_col=0)
    
    qtl1['pepagreeqtl']=np.nan
    qtl1['pepagreeqtlsign']=np.nan
    qtl1['npeptidesqtl']=np.nan
     
    pepmeta1=pepmeta.copy(deep=1)
    pepquant1=pepquant[(pepquant>0).sum(1)>=30]
    pepmeta1=pepmeta.loc[np.intersect1d(pepmeta.index,pepquant1.index)]
    pepmeta1['nlines']=(pepquant>0).sum(1).loc[pepmeta1.index]
    pepmeta1['quant']=np.log(pepquant).mean(1).loc[pepmeta1.index]
    
    for igen,gen in enumerate(np.intersect1d(list(peph5.keys()),qtl1.index)):
        if igen%10==0: print (igen)
        
        pro= qtl1.loc[gen,'feature_id'] 
        pronames=np.array(list(proh5[gen]['data/p_value_raw'].keys()))
        try:pro=pronames[ np.in1d(pronames.astype('S32').astype('U'),pro)][0]
        except:pro=pronames[[pro == p[:len(pro)] for p in pronames]][0]
        gene_name=qtl1.loc[gen,'gene_name']
        snp_id=qtl1.loc[gen,'snp_id']
        snp_id='_'.join(snp_id.split('_')[:2])
    #        snp_id=qtl1.loc[gen,'VE_LD']
        temp=np.array(list(proh5[gen]['data/snp_id'][pro][:])).astype('U')
        
        temp=np.array(['_'.join(t.split('_')[:2])for t in temp])
        ind=np.where(temp==snp_id)[0]
        if len(ind)==0:continue
    
        
        peps=pd.DataFrame(index=np.intersect1d(list(peph5[gen]['data/beta'].keys()),pepmeta1.index),columns=['p_value_raw','beta'])
        for f in ['start_in_protein','nlines','quant']:peps[f]= pepmeta1.loc[peps.index,f]
#        peps=peps[peps['nlines']>=np.nanmedian(peps['nlines'])]
        try: peps=peps.loc[np.intersect1d(peps.index,pepmeta1.set_index('superior_feature_id').loc[pro,'feature_id'])]
        except: 
            try:   peps=peps.loc[np.intersect1d(peps.index,pepmeta1.set_index('superior_feature_id').loc[pro.split('-')[0],'feature_id'])]
            except:   
                try:peps=peps.loc[np.intersect1d(peps.index,pepmeta1.set_index('gene_name').loc[gene_name,'feature_id'])]
                except: continue
 
        for f in ['beta','p_value_raw']:peps.loc[gen,f]= proh5[gen]['data'][f][pro][:] [ind][0]
        
        for pep in np.setdiff1d(peps.index,[gen+"lead",gen]):
            temp=np.array(list(peph5[gen]['data/snp_id'][pep][:])).astype('U')
            temp=np.array(['_'.join(t.split('_')[:2])for t in temp])
            ind=np.where(temp==snp_id)[0]
            if len(ind)==0:continue
            for f in ['beta','p_value_raw']:
                peps.loc[pep,f]=peph5[gen]['data'][f][pep][:][ind][0]

#        pepmeta.loc[peps.loc[np.intersect1d(peps.index,pepmeta.index)].index,'beta']=peps.loc[np.intersect1d(peps.index,pepmeta.index),'beta']
        qtl1.loc[gen,'pepagreeqtl']= ((peps['p_value_raw']<1)&(peps['beta']*peps.loc[gen,'beta']>0)).sum()-1
        qtl1.loc[gen,'pepagreeqtlsign']= ((peps['p_value_raw']<0.01)&(peps['beta']*peps.loc[gen,'beta']>0)).sum()-1
        qtl1.loc[gen,'npeptidesqtl']= ((peps['p_value_raw']<1)).sum()-1
    qtl1['fraction']=qtl1['pepagreeqtlsign']/qtl1['npeptidesqtl']
    return (qtl1)
def write_snp_snap(qtl1,path_data,name,snp_field='snp_id'):

    qtl1['snp_id_dots']=[':'.join(s.split('_')[:2]) for s in qtl1 [snp_field]]
    
    pd.DataFrame(qtl1 ['snp_id_dots']). to_csv(path_data+name+"_snp_lists_SNPsnap.txt",header=False,index=False)
#sys.exit()

def get_ld_snp_snap(file  ,snps): 
    try: 
        ld=pd.read_table(file+'.tsv',sep='\t',index_col=0)
#        mapping=pd.read_table(file+'_buddy map_input.tsv',sep='\t',index_col=0)
    except:
        ld={}
        for chrom in np.arange(1,23):
            print (chrom)
            ld[chrom]=pandas.read_table('/Users/mirauta/Data/Annotation/Ensembl_37.75/r08/chr'+str(chrom)+'_processed.txt',sep='\t',index_col=0)        
            ld[chrom]=ld[chrom].loc[np.intersect1d(ld[chrom].index,snps)]
        for chrom in np.arange(1,23):
            print (chrom)
            ld[chrom]['ld_pos_chrom']=ld[chrom]['ld_pos'].replace(regex=r',', value=','+str(chrom)+':')
            ld[chrom]['ld_pos_chrom']=str(chrom)+':'+ld[chrom]['ld_pos_chrom']    
        #    pd.concat(ld.values(),0).to_csv(path_data+'ld_expansion_SNP_snap_distance100_ldb08_new',sep='\t')    
        ld=pd.concat(ld.values(),0)
        ld['n_ld']=ld['ld_pos_chrom'].str.count(pat=':')
        ld=ld[np.isfinite(ld['n_ld'])]
        ld['n_ld']=ld['n_ld'].astype(int)
        ld['ld_pos_chrom_index']=[','.join(np.repeat(l,ld['n_ld'].loc[l])) for il,l in enumerate(ld['ld_pos_chrom'].index)]
        ld.to_csv(file+'.tsv',sep='\t')  

    mapping=pd.DataFrame(np.hstack([e.split(',') for e in ld['ld_pos_chrom_index']]),index=np.hstack([e.split(',') for e in ld['ld_pos_chrom']]),columns=['ld'])
#        mapping=mapping[mapping==mapping]
#        mapping=pd.DataFrame(index=','.join(mapping.index),data=','.join(np.squeeze(mapping.values)))
#        
#        mapping.to_csv(file+'_buddy map_input.tsv',sep='\t') 
    return [ld,mapping]

def get_enrichement_mapping(trait_labels,mapping,back,qtl):
    mappingtrait={}
    for i,trait in enumerate(trait_labels):
        mappingtrait[trait]=mapping[np.in1d(mapping['ld'],qtl[trait].index)]
        print (mappingtrait[trait].shape)

    return mappingtrait

def get_enrichement(trait_labels,mapping,back,qtl,mappingtrait, annotation_snp,ld):          
    backcolumns=list(back[trait_labels[0]].columns[:-1].values)
    rez={}
    for i,trait in enumerate(trait_labels):
        print (trait)
        qtl[trait]['chromosome']=qtl[trait]['chromosome'].astype('U')
        back[trait]['chromosome']=back[trait]['chromosome'].astype('U')        
        rez[trait]=pd.DataFrame(index=np.unique(back[trait]['chromosome']),columns=['ld_exp','ld_exp_back','ld_exp_gwas','snp_gwas_back','snp_gwas','ld_exp_gwas_back']+backcolumns)

        chrs=np.unique(qtl[trait]['chromosome'])
        print(chrs)
        for chrom in chrs:
            
            try:
                temp=ld.loc[ qtl[trait][qtl[trait]['chromosome']==chrom].index,'ld_pos_chrom']
                print (chrom)
            except:continue
            rez[trait]['ld_exp'][chrom]=np.unique(np.hstack([s.split(',') for s in  temp[temp==temp]]))
            rez[trait]['ld_exp_gwas'][chrom]=np.intersect1d(annotation_snp.index,rez[trait]['ld_exp'][chrom])
            rez[trait]['snp_gwas'][chrom]=np.intersect1d(mappingtrait[trait].loc[rez[trait]['ld_exp_gwas'][chrom]].values,temp.index)
            
            temp1=back[trait][back[trait]['chromosome']==chrom].astype('U')
            temp= ld.loc[ np.hstack(temp1.values) ,'ld_pos_chrom']; #temp index contains back random snps
            rez[trait]['ld_exp_back'][chrom]=np.unique(np.hstack([s.split(',') for s in  temp[temp==temp]]))
            rez[trait]['ld_exp_gwas_back'][chrom]=np.intersect1d(annotation_snp.index,rez[trait]['ld_exp_back'][chrom])
            rez[trait]['snp_gwas_back'][chrom]=np.intersect1d(mapping.loc[rez[trait]['ld_exp_gwas_back'][chrom],'ld'],temp.index.astype('U'))
            for column in backcolumns:
                rez[trait].loc[chrom,column]=np.intersect1d(rez[trait]['snp_gwas_back'][chrom],temp1[column].astype('U'))
#    for i,trait in enumerate(trait_labels):
#        rez[trait].to_csv(path_data+'GWAS_rez_'+trait,sep='\t')    
    final={}
    final2={}
    for i,trait in enumerate(trait_labels):
        chrs2=np.array(list(rez[trait]['snp_gwas'].keys()))
        chrs2=chrs2[[isinstance(rez[trait]['snp_gwas'][chrom],np.ndarray) for chrom in chrs2]]
                
        final[trait]=pd.DataFrame(index=chrs2,columns=['back','qtl','back_shape','qtl_shape']+backcolumns)
        final2[trait]={}

        final[trait]['qtl'] = [rez[trait]['snp_gwas'][chrom].shape[0] for chrom in chrs2]
        final2[trait]['qtl']= np.hstack([rez[trait]['snp_gwas'][chrom] for chrom in chrs2])
        final[trait]['qtl_shape'] = [qtl[trait][qtl[trait]['chromosome']==chrom].shape[0] for chrom in chrs2]
        final[trait]['back']= [rez[trait]['snp_gwas_back'][chrom].shape[0] for chrom in chrs2]
        final[trait]['back_shape'] = [back[trait][back[trait]['chromosome']==chrom].shape[0]*100 for chrom in chrs2]
        for column in backcolumns:
            final[trait][column]=[rez[trait][column][chrom].shape[0] for chrom in chrs2]
            final2[trait][column]= np.hstack([rez[trait][column][chrom] for chrom in chrs2])
    
    
    final3=pd.DataFrame(np.hstack([np.hstack([ (((final2[trait]['qtl'].shape[0]+1)/(qtl[trait].shape[0]+1)))/\
                                       (((1+final2[trait][column].shape[0])/(1+back[trait].shape[0])))  for column in backcolumns]) for i,trait in enumerate(trait_labels)]),columns=['value'])
    final3['trait']=np.hstack([ np.repeat(trait,100) for i,trait in enumerate(trait_labels)])

    return [final2,final3]

#temp=get_enrichement(trait_labels,mapping,back,qtl,mappingtrait,annotation_snp,ld)




def load_QTL_vep(path_data='/Users/mirauta/Results/hipsci/QTL_may2019/',name="qtl_results_LD_VEP_gwas",traits=['pQTL','eQTL' ]):
    def changesep(x):  
        if '[' in x: 
            return literal_eval(x.replace("\' \'",'\',\'' ).replace("...,\n",'\',\'' ).replace("\n ...,",'\',\'' ))                
        else: return [x]
    def mergeback(x):return ';'.join(x)
    
    qtl={}
    for trait in traits:
        if trait=='pQTL': 
            try:qtl[trait]=pd.read_table(path_data+'protein_mrna'+'_'+name+'.txt', sep='\t')
            except:qtl[trait]=pd.read_table(path_data+'protein_mrna'+'_'+name+'.tsv', sep='\t')
        elif trait=='eQTL': qtl[trait]=pd.read_table(path_data+'mrna_protein'+'_'+name+'.txt', sep='\t')
        elif trait=='tQTL': qtl[trait]=pd.read_table(path_data+'transcript_protein'+'_'+name+'.txt', sep='\t')
        elif trait=='pepQTL': qtl[trait]=pd.read_table(path_data+'peptide_protein'+'_'+name+'.txt', sep='\t')
        qtl1=qtl[trait]  
#        qtl1['VE_LD']=qtl1['VE_LD'].apply(changesep)
#        qtl1[qtl1['VE_LD']=="None"]
        qtl1['VE_LD2']=qtl1['VE_LD']#.apply(mergeback)
    #    qtl1['VE_LD']
        qtl1=qtl1.set_index('ensembl_gene_id',drop=0)
        qtl1['VE_LD']=[0 if (s==0)|(s!=s) else s.split(';') for s in qtl1['VE_LD2']]
        if trait=="eQTL":qtl1['sign']=qtl1['mrna_eff_size']>0
        elif trait=="pQTL":qtl1['sign']=qtl1['protein_eff_size']>0
        qtl1['variant_consequences']=['' if (s=='0')|((s==0)|(s!=s)) else s.replace('@',';') for s in qtl1['VE_LD2']]
        qtl1['len_ld_snps']=[np.array(len(qtl1['ld_snps'].iloc[n].split(','))) for n in np.arange(qtl1.shape[0])]
#        qtl1['3utr']=[False if qtl1['VE_LD'].iloc[n]==0 else len(np.intersect1d('3_prime_UTR_variant',  qtl1['VE_LD'].iloc[n]))>0 for n in np.arange(qtl1.shape[0])]
#        qtl1['inframe_variant']=[False if qtl1['VE_LD'].iloc[n]==0 else len(np.intersect1d(inframe_variant,  qtl1['VE_LD'].iloc[n]))>0 for n in np.arange(qtl1.shape[0])]
        if trait=='pQTL':
            qtl1['replicated']=((qtl1['mrna_nominal_p_value_bonf']<0.01)&(qtl1['mrna_eff_size']*qtl1['protein_eff_size']>0))
            qtl1['replicatedet']=qtl1['replicated'] |   ((qtl1['transcript_nominal_p_value_bonf']<0.01)&(qtl1['transcript_eff_size']*qtl1['protein_eff_size']>0))
#            except:qtl1['replicatedet']=qtl1['replicated']
        if trait=='pepQTL':
            qtl1['replicated']=((qtl1['mrna_nominal_p_value_bonf']<0.01)&(qtl1['mrna_eff_size']*qtl1['peptide_eff_size']>0))
            qtl1['replicatedet']=qtl1['replicated'] |   ((qtl1['transcript_nominal_p_value_bonf']<0.01)&(qtl1['transcript_eff_size']*qtl1['peptide_eff_size']>0))
        elif trait=='eQTL':
            qtl1['replicated']=((qtl1['protein_nominal_p_value_bonf']<0.01)&(qtl1['mrna_eff_size']*qtl1['protein_eff_size']>0))
            qtl1['replicatedet']=qtl1['replicated'] |   ((qtl1['peptide_nominal_p_value_bonf']<0.01)&(qtl1['peptide_eff_size']*qtl1['mrna_eff_size']>0))     
        elif trait=='tQTL':
            qtl1['replicated']=((qtl1['protein_nominal_p_value_bonf']<0.01)&(qtl1['transcript_eff_size']*qtl1['protein_eff_size']>0))
            qtl1['replicatedet']=qtl1['replicated'] |   ((qtl1['peptide_nominal_p_value_bonf']<0.01)&(qtl1['peptide_eff_size']*qtl1['mrna_eff_size']>0))     
        
        qtl[trait]=qtl1
    variants2={'UTR_variant':  ['UTR'],'frameshift_variant':['frameshift'],'inframe_variant':  ['inframe'],
     'intron_variant':  ['intron'], 'splicing_variant':  ['splicing'], 'synonymous_variant':  ['synonymous']}
    for trait in qtl.keys():
        for var in variants2.keys():
            qtl[trait][var]=qtl[trait][ variants2[var]]==qtl[trait][ variants2[var]]
         
    return qtl

def plot_QTL_vep(qtl,dest_folder):
    
#    fig = plt.figure(figsize=(8, 3.5))
#    fig.patch.set_facecolor('white')
#    for it,trait in enumerate(['pQTL','eQTL']):
#        qtl1=qtl[trait]
#        plot_fisher_exact_results_text_number_genes(fig,qtl1,it,its=2,thr=0.05, xlim=(-0.20,3),ylim=(-0.20,3))
#    plt.tight_layout()
#    plt.savefig(dest_folder+"QTL_VE.png",dpi=600)
    
    fig = plt.figure(figsize=(4, 3.5))
    fig.patch.set_facecolor('white')
    for it,trait in enumerate(['pQTL' ]):
        qtl1=qtl[trait]
        plot_fisher_exact_results_text_number_genes(fig,qtl1,it,its=1,thr=0.5, xlim=(-0.20,6),ylim=(-0.20,6))
    plt.tight_layout()
    plt.savefig(dest_folder+"QTL_VE_protein.png",dpi=600)
    plt.savefig(dest_folder+"QTL_VE_protein.svg",dpi=600)
    
    
    fig = plt.figure(figsize=(8, 4))
    fig.patch.set_facecolor('white')
    for it,trait in enumerate(['pQTL','eQTL']):
        qtl1=qtl[trait]
        plot_fisher_exact_bar(fig,qtl1,it,thr=0.3,trait=trait)
    
    plt.tight_layout()
    plt.savefig(dest_folder+"QTL_VE_hist.png",dpi=600)
    plt.savefig(dest_folder+"QTL_VE_hist.svg",dpi=600)
    
#plot_QTL_vep(qtl,dest_folder="/Users/mirauta/Results/hipsci/manuscript_images/")
#


def sample_equilibrated_simple(fieldeq,samplesi,bins=10):
    thrs=np.nanpercentile(np.hstack([samplesi[fieldeq] for f in samplesi.keys()]),np.arange(bins+1)*(100/bins))
    thrs[-1]=thrs[-1]+0.0001
    thrs=np.sort(np.unique(thrs))
    proteins_bin=np.vstack([[samplesi[(samplesi[fieldeq]>=thr)&(samplesi[fieldeq]<thrs[ithr+1])].shape[0] for ithr,thr in enumerate(thrs[:-1])] for f in samplesi.keys()]).min(0)
    proteins_bin.sum()
    sample={}
    for f in samplesi.keys():
        sample[f]=np.hstack([np.random.choice(samplesi[(samplesi[fieldeq]>=thr)&(samplesi[fieldeq]<thrs[ithr+1])].index,\
                            proteins_bin[ithr],replace=False) if proteins_bin[ithr]>0  \
              else[]for ithr,thr in enumerate(thrs[:-1])])

    samples={}
    for f in samplesi.keys():
        samples[f]=samplesi.loc[sample[f]]     
    
    for f in samplesi.keys():
        print(sample[f].shape)
        print(samples[f].shape)
#        print (samplesi[fieldeq].median())
#        sb.kdeplot( samplesi[fieldeq] ,label=f,bw=0.2)
    return samples

def sample_equilibrated(fieldeq,samplesi,bins=10):
    thrs=np.nanpercentile(np.hstack([samplesi[f][fieldeq] for f in samplesi.keys()]),np.arange(bins+1)*(100/bins))
    thrs[-1]=thrs[-1]+0.0001
    thrs=np.sort(np.unique(thrs))
    proteins_bin=np.vstack([[samplesi[f][(samplesi[f][fieldeq]>=thr)&(samplesi[f][fieldeq]<thrs[ithr+1])].shape[0] for ithr,thr in enumerate(thrs[:-1])] for f in samplesi.keys()]).min(0)
    proteins_bin.sum()
    sample={}
    for f in samplesi.keys():
        sample[f]=np.hstack([np.random.choice(samplesi[f][(samplesi[f][fieldeq]>=thr)&(samplesi[f][fieldeq]<thrs[ithr+1])].index,\
                            proteins_bin[ithr],replace=False) if proteins_bin[ithr]>0  \
              else[]for ithr,thr in enumerate(thrs[:-1])])

    samples={}
    for f in samplesi.keys():
        samples[f]=samplesi[f].loc[sample[f]]     
    
    for f in samplesi.keys():
        print(sample[f].shape)
        print(samples[f].shape)
#        print (samplesi[f][fieldeq].median())
#        sb.kdeplot( samplesi[f][fieldeq] ,label=f,bw=0.2)
    return samples



def plot_manhatan_genes_feature_snp_file( folder_destination='/Users/mirauta/Results/hipsci/manuscript_images/Manhattan/',\
                      plot_name='manhattan',\
                      trans=None,qtl=None,metapro=None, \
                      colors=['darkblue','lightcoral','orange','slateblue','darkseagreen']*10,\
                      figsize=4, cis_gene_name =None,\
                      p_value_field='p_value',log_flag=True,ylim=None,  fplot=None,ax=None):
    
    
    snps=qtl['snp_id']
    snp=qtl.query('gene_name==@cis_gene_name')['snp_id'].values[0]
    trans1=trans.set_index('trans_gene_name').loc[trans.query('qv_nominal<0.2' ).loc[snp,'trans_gene_name']]
    
    
    genes=np.unique(trans1.index)
    rez=pd.DataFrame(np.random.uniform(0.01,1,len(snps)*len(genes)).reshape(-1,len(genes)),columns=genes,index=snps)
    rez['chromosome']=np.array([s.split('_')[0] for s in rez.index]).astype(int)
    rez['position']=np.array([s.split('_')[1] for s in rez.index]).astype(float)
    for g in genes:
        1
        rez.loc[trans1.loc[g,'snp_id'],g]=trans1.loc[g].set_index('snp_id').loc[trans1.loc[g,'snp_id'],p_value_field]
    rez=rez.iloc[np.argsort(rez['chromosome'].astype(int))]
    
    lengths=pd.DataFrame(index=np.unique(rez['chromosome']),columns=['len']); 
    for chr in np.unique(rez['chromosome']):              
        lengths.loc[chr,'len']= np.max(rez['position'] [rez['chromosome'] ==chr])
    lengths['cum']=np.hstack([0,np.cumsum(lengths['len'])[:-1]])
    lengths['cum2']=lengths['cum']+lengths['len']/2;
#    lengths.loc[1,'cum2']=lengths.loc[1,'len']/2
    rez['genome_position']=lengths.loc[rez['chromosome'],'cum'].values+rez['position'].values
    rez=rez.iloc[np.argsort(rez['genome_position'].astype(int))]
    generez=metapro.set_index('gene_name',drop=0).loc[np.hstack([genes,cis_gene_name])][['gene_name','chromosome','start']].drop_duplicates()
    
    generez['start']=generez['start'].values+lengths.loc[generez['chromosome'].values.astype(int),'cum'].values.astype(int)

    fig=plt.figure(figsize=(figsize*2,figsize*len(genes)))
    fig.set_facecolor('white')
    fplot = gridspec.GridSpec(len(genes)*6,8)
    fplot.update(hspace=10, wspace=10)    
    ax = [plt.subplot(fplot[(i*3):(i*3+3),:10]) for i in np.arange(len(genes))]
    
    
  
    for i , a in enumerate(ax):
        a.spines['top'].set_visible(False);    a.spines['right'].set_visible(False);a.yaxis.set_ticks_position('left'); a.xaxis.set_ticks_position('bottom');a.set_axis_bgcolor('white');
        a.plot(rez['genome_position'] ,-np.log10(rez[genes[i]].values),'o',color=colors[i],markersize=4)
        a.plot(rez['genome_position'][rez[genes[i]].values<10**-4] ,-np.log10(rez[genes[i]].values)[rez[genes[i]].values<10**-4],'o',color=colors[i],markersize=4)
#        a.set_ylim(0,17)
        a.set_xticks(lengths['cum2'])
        a.set_yticks([0,5,10])
#        a.set_yticks(np.arange(4)*5)

#    plt.show()
        a.set_xticklabels(lengths.index,fontsize=9)
        a.set_xticklabels(np.hstack([np.arange(1,19),'','20','',22]),fontsize=10)
        
        
        a.add_patch(Rectangle((generez.loc[genes[i],'start'], 0), 2*10**7, 5, facecolor="grey",    alpha=0.99935))
        a.annotate(generez.iloc[i].loc['gene_name'], xy=(generez.iloc[i]['start'], 10) ,fontsize=10)
#        plt.title(generez.iloc[i].loc['gene_name'])
#        a.set_ylabel(generez.iloc[i].loc['gene_name']+"\n -log10PV", labelpad=15,fontsize=16,rotation=90,horizontalalignment= 'center' ,verticalalignment= 'center')
        a.set_ylabel("-log10 PV", labelpad=15,fontsize=16,rotation=90,horizontalalignment= 'center' ,verticalalignment= 'center')
        if i==(len(ax)-1):a.set_xlabel("Genome position",labelpad=10,fontsize=16)
#        a.set_ylabel(" -log10PV", labelpad=10,fontsize=10,rotation=90,horizontalalignment= 'center' ,verticalalignment= 'center')
        
        if i==0: a.annotate(qtl.query('gene_name==@cis_gene_name')['snp_rsid'].values[0], xy=(generez.loc[cis_gene_name,'start'], 8) ,fontsize=10)
         
#    plt.tight_layout()
#    plt.show()   
#    plt.savefig(folder_destination+plot_name+'.svg',dpi=600,bbox='tight')
    plt.savefig(folder_destination+cis_gene_name+'_'+plot_name+'.png',dpi=600,bbox='tight')
    plt.savefig(folder_destination+cis_gene_name+'_'+plot_name+'.svg',dpi=600,bbox='tight')


def fun_get_cis_data_trans(path_data='/Users/mirauta/Results/hipsci/QTL_may2019/', folder_data='/Users/mirauta/Data/MS/hipsci/phenotypes/may_2019/', \
                           trans_file='param_cistrans_protein_maxquant_peer_log_extended_pqtl_qtl_results_feature_snp_1.0'):

    p_value_field='p_value'
    metapro=pd.read_table(folder_data+'hipsci.proteomics.maxquant.uniprot.TMT_batch_24.20170919_protein_meta_recurrent_genes.txt',index_col=0)
    #metapro=pd.read_table(folder_data+'hipsci.proteomics.maxquant.uniprot.TMT_batch_24.20170919_protein_meta_recurrent_genes_recurrent_peptide.txt',sep='\t',index_col=0)
    metapro['protein_id']=[p.split('-')[0] for p in metapro.index]
    #overlap=pd.read_table('/Users/mirauta/Data/MS/hipsci/phenotypes/may_2019/notmapped_overlap_peptide_snp_idSNP_multiple_junction_phase_30.01.txt',sep='\t',index_col=0).set_index('feature_id',drop=0)
    #metapep=pd.read_table('/Users/mirauta/Data/MS/hipsci/phenotypes/jan_including_nonvalidpeptides/hipsci.proteomics.maxquant.uniprot.TMT_batch_24.20170919_peptide_meta_all_genes.txt',index_col=0)
    #metapep=metapep.loc[np.intersect1d(metapep.index,overlap['feature_id'])]
    
    [prometa,complexes,corumallproteins,corumproteins,complexes2,complex_pairs,intact,string] = load_complex_prometa(intact_flag=True)
    complex_pairs=complex_pairs.set_index('A',drop=0); complex_pairs=complex_pairs.loc[np.intersect1d(complex_pairs.index,metapro['protein_id'])]
    complex_pairs=complex_pairs.set_index('B',drop=0); complex_pairs=complex_pairs.loc[np.intersect1d(complex_pairs.index,metapro['protein_id'])]
    complex_pairs=complex_pairs.iloc[np.unique(complex_pairs['pair'],return_index=1)[1]]
    
    intact=intact.set_index('A',drop=0); intact=intact.loc[np.intersect1d(intact.index,metapro['protein_id'])]
    intact=intact.set_index('B',drop=0); intact=intact.loc[np.intersect1d(intact.index,metapro['protein_id'])]
    intact=intact.iloc[np.unique(intact['pair'],return_index=1)[1]]
    
    string=string.set_index('A',drop=0); string=string.loc[np.intersect1d(string.index,metapro['protein_id'])]
    string=string.set_index('B',drop=0); string=string.loc[np.intersect1d(string.index,metapro['protein_id'])]
    string=string.iloc[np.unique(string['pair'],return_index=1)[1]]
    
    qtl=load_QTL_vep(path_data='/Users/mirauta/Results/hipsci/QTL_may2019/',name="qtl_results_LD_VEP_gwas_rsid",traits=['pQTL'])['pQTL'].set_index('snp_id',drop=0)
    
    trans=pd.read_table(path_data+trans_file+'.txt',sep='\t').set_index('snp_id',drop=0)
    trans['qv_nominal']=scst2.fdrcorrection0(trans[p_value_field] )[1]
    print(trans.query('qv_nominal<0.1').shape)
    print(trans.query('qv_nominal<0.2').shape)
    
    trans=trans.set_index('feature_id',drop=0)
    trans['trans_gene_name']=metapro.loc[trans.index,'gene_name']
    trans['trans_protein_id']=metapro.loc[trans.index].index
    trans['trans_ensembl_gene_id']=metapro.loc[trans.index,'ensembl_gene_id']
    
    trans=trans.set_index('snp_id',drop=0)
    trans['cis_gene_name']=qtl.loc[trans.index,'gene_name']
    trans['cis_protein_id']=qtl.loc[trans.index,'protein_id']
    trans['cis_ensembl_gene_id']=qtl.loc[trans.index,'ensembl_gene_id']
    
    def fun_isoform (p): return  p.split('-')[0]
    trans['trans_protein_id_no_isoform']=trans['trans_protein_id'].apply(fun_isoform)
    trans['cis_protein_id_no_isoform']=trans['cis_protein_id'].apply(fun_isoform)
    trans['pair_no_isoform']=trans['cis_protein_id_no_isoform']+'_'+trans['trans_protein_id_no_isoform']
    trans['pair']=trans['cis_protein_id']+'_'+trans['trans_protein_id']
    trans['absbeta']=np.abs(trans['beta'])
    qtl['cis_protein_id_no_isoform']=qtl['feature_id'].apply(fun_isoform)
    complex_pairs1=complex_pairs[np.in1d(complex_pairs['A'],qtl['cis_protein_id_no_isoform'])|np.in1d(complex_pairs['B'],qtl['cis_protein_id_no_isoform'])]
    pairs=np.intersect1d(trans['pair_no_isoform'],complex_pairs1['pair'])
    trans=trans.set_index('pair_no_isoform',drop=0)
    trans['pair_in_corum']=0;trans['pair_in_corum'].loc[pairs]=1
    
    string1=string[np.in1d(string['A'],qtl['cis_protein_id_no_isoform'])|np.in1d(string['B'],qtl['cis_protein_id_no_isoform'])]
    pairs=np.intersect1d(trans['pair_no_isoform'],string1['pair'])
    trans=trans.set_index('pair_no_isoform',drop=0)
    trans['string']=0;trans['string'].loc[pairs]=1
    
    intact1=intact[np.in1d(intact['A'],qtl['cis_protein_id_no_isoform'])|np.in1d(intact['B'],qtl['cis_protein_id_no_isoform'])]
    pairs=np.intersect1d(trans['pair_no_isoform'],intact1['pair'])
    trans=trans.set_index('pair_no_isoform',drop=0)
    trans['intact']=0;trans['intact'].loc[pairs]=1     
         
    trans['pair_annotation']=(trans['string']==1)|(trans['pair_in_corum']==1)
    trans.to_csv(path_data+trans_file+'_cisdata_all.txt',sep='\t')

    signtrans=trans.query('qv_nominal<0.1' );print (signtrans.shape)
    signtrans.to_csv(path_data+trans_file+'_significant_cisdata_qv01.txt',sep='\t')
    
    signtrans=trans.query('qv_nominal<0.2' );print (signtrans.shape)
    qtl=qtl.set_index('snp_id',drop=0)
    signtrans['cis_beta']=qtl.loc[signtrans.index,'beta']
    signtrans.to_csv(path_data+trans_file+'_significant_cisdata_qv02.txt',sep='\t')

    signtrans=trans.query('qv_nominal<0.4' );print (signtrans.shape)
    qtl=qtl.set_index('snp_id',drop=0)
    signtrans['cis_beta']=qtl.loc[signtrans.index,'beta']
    signtrans.to_csv(path_data+trans_file+'_significant_cisdata_qv04.txt',sep='\t')



def plot_gene_transcript_rezolution(qtl,data,genes,gene_names,trait,traitsnp=None,plot_peptides=False,manhattan_flag=False,remove_transcript_flag=False,keep_peptide_flag=False,snp_df=None,peptide_transcript=None,\
              phenotype_df=None,mrna_df_ensembl=None,transcript_df=None,peptide_df=None,anntr=None,args=None):
    
    for ig,gene in enumerate(genes): 
        colors=np.array(['sienna', 'darkkhaki', 'seagreen']*100)
        if traitsnp is None: traitsnp=trait
        gene_name=gene_names[ig]
        try:protein=qtl['pQTL'].loc[gene,'feature_id']
        except:
            try:protein=qtl['eQTL'].loc[gene,'replicated_feature_id']
            except: protein=qtl['tQTL'].loc[gene,'replicated_feature_id']
        
        if gene_name!=gene_name:gene_name=''

        snp_id=qtl[traitsnp].loc[gene,'snp_id']
        snp=snp_df[snp_id].values

        alleles=qtl[trait].loc[gene,'snp_id'].split('_')[2:];alleles=[alleles[0]+alleles[0],alleles[0]+alleles[1],alleles[1]+alleles[1]]  
        ypr=np.log(phenotype_df.loc[protein].values); 
        transcripts=np.array(list(data['tQTL'][gene]['data/features'])).astype('U')

        genepep=peptide_transcript.loc[[gene]].set_index('peptide',drop=0)
        genepep=genepep.loc[np.intersect1d(peptide_df.index,genepep.index)]

        anntrgene=anntr.set_index('ensembl_transcript_id').loc[[t.split('_')[0]for t in transcripts]]
        anntrgene['cod']=anntrgene['transcript_type']=='protein_coding'
        anntrgene.index=anntrgene.index+'_'+anntrgene['gene_name']
        anntrgene['mean']=transcript_df.loc[transcripts].mean(1).loc[anntrgene.index]
        anntrgene=pd.concat([anntrgene.query('cod == False').sort_values(['mean'], ascending=[0])[:4],
                   anntrgene.query('cod == True').sort_values(['mean'], ascending=[0])[:4]],0)
        transcripts=anntrgene.index
#        x=pd.concat([transcript_df.loc[transcripts],np.log(phenotype_df.loc[protein]),\
#                     mrna_df_ensembl.loc[[gene]],snp_df[qtl[traitsnp].loc[[gene],'snp_id']].T],0);print (x.shape)
#        x=x.T.corr()
#        print (qtl1.loc[gene])
#        sb.heatmap(x[::-1])
##        plt.xticks(np.arange(len(transcripts))+0.5,np.around(transcript_df.loc[transcripts].mean(1),2))
#        plt.xticks(np.arange(len(transcripts))+0.5,anntrgene.loc[x.index]['cod'][::-1])
#        plt.show()
#        transcripts=abs(x.loc[qtl[traitsnp].loc[[gene],'snp_id'].values,transcripts]).T.sort_values(qtl[traitsnp].loc[[gene],'snp_id'].values[0])[-3:].index
#        anntrgene.loc[transcripts]
        
    #    transcripts=transcripts[np.in1d([t.split('_')[0]for t in transcripts],genepep.loc[genepepqtl,'ensembl_isoform_id'])]
    #    if len(transcripts)==0: continue
    #    genepep.set_index('ensembl_isoform_id').loc['ENST00000495526']
        if len(transcripts)>10:transcripts=transcripts[transcript_df.loc[transcripts].sum(1)>np.nanpercentile(transcript_df.loc[transcripts].sum(1),00)]
        transcripts=transcripts[transcript_df.loc[transcripts].mean(1)>0.0001]
        if gene_name=='MMAB': transcripts=np.array(['ENST00000540016_MMAB','ENST00000537496_MMAB'])[:2]
#        if gene_name=='CCDC114': transcripts=np.array(['ENST00000315396_CCDC114', 'ENST00000474199_CCDC114'])
        if (gene_name=='PEX6'): transcripts=np.array([])
#        if (gene_name=='NTPCR'): transcripts=np.array(['ENST00000487953','ENST00000366628'])
#        if (gene_name=='RUVBL2'): transcripts=np.array(['ENST00000596247','ENST00000601968'])
        if (gene_name=='CTTN'): transcripts=np.array(['ENST00000301843_CTTN','ENST00000346329_CTTN'])
        if (gene_name=='METTL21A'): transcripts=np.array(['ENST00000425132_METTL21A','ENST00000477919_METTL21A'])
        

                
        anntrgene=anntr.set_index('ensembl_transcript_id').loc[[t.split('_')[0]for t in transcripts]]
        anntrgene['index']=anntrgene.index+'_'+gene_name
        anntrgene['meanTPM']=transcript_df.loc[anntrgene['index']].mean(1).values
        if gene_name=='CTTN':anntrgene=anntrgene[anntrgene['transcript_type']=='protein_coding'];
        transcripts=np.unique(anntrgene['index'])
        ''' the  transcript_peptide table''' 
        transpep=pd.DataFrame(data=np.vstack([np.in1d(np.unique(genepep.index),genepep[genepep['ensembl_isoform_id']==t].index) for t in  anntrgene.index]).T,\
                       index=np.unique(genepep.index),columns=[tr.split('_')[0].replace('ENST00000','Tr-') for tr in anntrgene.index]).T
        corpep=np.array([[correlation(x=peptide_df.loc[pep].values, y=transcript_df.loc[tr].values)  for pep in np.unique(genepep['peptide'])] \
                              for tr in anntrgene['index']])
        corpep[~transpep.values]=np.nan;
        corpep=pd.DataFrame(corpep,  index=anntrgene.index,columns=np.unique(genepep['peptide']));
        
        p_value_field='p_value_raw'#    p_value_field='empirical_feature_p_value'        
        temppeptide=peptide_transcript.loc[gene,['peptide','mean_expression']].drop_duplicates()
    
        if (gene_name=='PEX6'): fig =plt.figure(figsize=(2.75,2.5*2))
        else: fig = plt.figure(figsize=(8,4))

        fig.patch.set_facecolor('white')
        ax1 = fig.add_subplot(122);

                             
        ax1.spines['right'].set_visible(False); ax1.spines['top'].set_visible(False);# ax1.yaxis.set_ticks_position('right'); 
        
        sb.boxplot(y=ypr,x=snp,color='steelblue',width=0.8)
        plt.plot([np.nanmedian(ypr[snp==s]) for s in np.arange(1)],'-',color='steelblue',label='protein ')
    #    for pep in np.unique(genepep['peptide']):
        if gene_name=='CTTN':
            plt.ylim(12,19)
            for pe in ['QDSAAVGFDYK']:
                plt.plot([np.nanmedian(np.log(peptide_df.loc[pe])[snp==s]) for s in np.arange(1)],'o',color='skyblue',label='unique peptide ')
                sb.boxplot(y=np.log(peptide_df.loc[pe]),x=snp,color='skyblue',boxprops=dict(alpha=0.9),width=0.5)
#        if (gene_name!='PEX6'): plt.legend(loc=3 if   gene_name!='CTTN' else 4,fontsize=10)  
        plt.xticks(np.arange(np.unique(snp).shape[0]),alleles[:np.unique(snp).shape[0]],fontsize=16)
        ax1.set_ylabel('log abundance', color='k',fontsize=16)
        
        if gene_name=='CTTN': plt.yticks(np.around(np.linspace(np.nanmin([np.min(ypr),np.log(np.nanmin(peptide_df.loc[pe].values))]),np.max( [np.max(ypr),np.log(np.max(peptide_df.loc[pe].values))]),num=4),1))
        else: 
            plt.yticks(np.around(np.linspace(np.min([np.nanmin(ypr)]),np.max( [np.nanmax(ypr) ]),num=4),1))
#            plt.yticks(np.around(np.linspace(np.min(ypr),np.max(ypr),num=4),1))
        ax2 = fig.add_subplot(121);
                                                          
        ax2.spines['right'].set_visible(False); ax2.spines['top'].set_visible(False);  
        
    #    fig = plt.figure(figsize=(8,8))
    
    
        y=(mrna_df_ensembl.loc[gene].values);
#        y= regxfromy(peer_mrna,y)[0];
        ylist={}
        ylist['gene']=pd.DataFrame(data=np.vstack([y,snp]).T,columns=['value','snp'])
        ylist['gene']['type']='gene'
#        qtl['eQTL'].loc[gene][:14]
        qtl['pQTL'].loc[gene][:24]
        
    #    #y=y-np.nanmean(y)
#        plt.plot([np.nanmedian(y[snp==s]) for s in np.arange(1)],'-',color='forestgreen',label='mRNA' if gene_name=='PEX6' else 'mRNA gene level')
#        sb.boxplot(y=y,x= snp,color='green',width=0.3)
        if gene_name=='ELP5':
            transcripts=transcripts.loc[['ENST00000574255','ENST00000354429']]
    
        noisoform=qtl['pQTL'].loc[gene,'feature_id'].split('-')[0]
        
    
        flatui=['green']
        for itr,tr in enumerate(transcripts):
#            cod=anntrgene.loc[tr.split('_')[0],'transcript_type']=='protein_coding'

            y=(transcript_df.loc[tr].values); #y= regxfromy(peer_mrna,y)[0];
            y=pd.DataFrame(np.vstack([y,snp,np.repeat(tr.split('_')[0],y.shape[0])]).T,columns=['value','snp','isoforms']);y['value']=y.value.astype(float)#y=y-np.nanmean(y)
            ylist[tr]=y[['value','snp']]
            ylist[tr]['snp']=ylist[tr]['snp'].astype(float).astype(int)
            
            if (anntrgene.loc[[tr.split('_')[0]],'protein_id'][0]== qtl['pQTL'].loc[gene,'feature_id'].split('-')[0]):  
                if anntrgene.query('protein_id==@noisoform')['meanTPM'].argmax()==tr.split('_')[0]:
                    flatui=flatui+['c']
                    ylist[tr]['type']= tr.split('_')[0]
            elif (anntrgene.loc[[tr.split('_')[0]],'transcript_type'][0]== "protein_coding"):
                flatui=flatui+['indianred']
            else:
                flatui=flatui+['grey']
#            else: ylist[tr]['type']="Noncoding"
        
        if gene_name=='CTTN': 
            colors=np.array(['c', 'indianred'])
            flatui = ['green','c','indianred']
        elif gene_name=='MMAB': 
            colors=np.array([ 'c','grey','m'])
            flatui = ['green','c','grey']
        elif gene_name=='METTL21A': 
            colors=np.array([ 'c','grey'])
            flatui = ['green','c','grey']
         
        sb.set_palette(flatui)
            
#            y['snp1']=y['snp'].astype(float)+0.5
#            plt.plot([y['value'][y['snp']==s].mean() for s in [0]],'-',color=colors[itr],label=tr.split('_')[0] ) 
#            if (gene_name=='MMAB')|(gene_name=='CCDC114'):
#                plt.plot([y['value'][y['snp']==s].mean() for s in [0]],'o',color=colors[itr] if cod else 'grey',label='isoform '+ str(itr+1) if itr==0 else 'isoform '+ str(itr+1) +' (non-coding)')
#            else:
#                plt.plot([y['value'][y['snp']==s].mean() for s in [0]],'o',color=colors[itr] if cod else 'grey',label='isoform '+ str(itr+1))
#            if cod:
##                print (y['snp'])
##                print (y['snp1'])
#
#                sb.violinplot(y='value',x='snp1',color=colors[itr] ,data=y,\
#                           boxprops=dict(alpha=0.5+0.4* (colors[itr]=='grey')),width=0.25) 
##                sb.boxplot(y='value',x='snp',color=colors[itr] ,data=y,\
##                           boxprops=dict(alpha=0.5+0.4* (colors[itr]=='grey')),width=0.25) 
#            else: sb.violinplot(y='value',x='snp1',color='grey' ,data=y,\
#                             boxprops=dict(alpha=0.5+0.4* (colors[itr]=='grey')),width=0.25) 
        
        yy=pd.concat([ylist[k] for k in ylist.keys()],0)
        


        sb.boxplot(y='value',x='snp',data=yy,hue='type',width=0.75) 
        plt.legend(loc=None,fontsize=0)
             
        yg=(mrna_df_ensembl.loc[gene].values);# y= regxfromy(peer_mrna,y)[0];
        plt.yticks(np.around(np.linspace(np.min([np.min(yg),np.min(transcript_df.loc[transcripts].values)]),np.max( [np.max(yg),np.max(transcript_df.loc[transcripts].values)]),num=4),1))
    #    if gene_name!='CCDC114': sb.boxplot(y=y,x= snp,color='green',width=0.8)
    
        ax2.spines['top'].set_visible(False); ax2.spines['right'].set_visible(False); ax2.yaxis.set_ticks_position('left'); 
        ax2.set_ylabel('log TPM', color='k',fontsize=16)
        plt.xticks(np.arange(np.unique(snp).shape[0]),alleles[:np.unique(snp).shape[0]],fontsize=16)
        plt.xlabel('')
#        if (gene_name!='PEX6'):plt.legend(loc=3,fontsize=10)  

        plt.tight_layout()
    
        plt.savefig('/Users/mirauta/Results/hipsci/manuscript_images/Transcripts/'+'fig3_peptide_'+gene_name+'_boxplot.png',dpi=900)
        plt.savefig('/Users/mirauta/Results/hipsci/manuscript_images/Transcripts/'+'fig3_peptide_'+gene_name+'_boxplot.svg',dpi=900,bbox='tight')
    #    plt.gcf().clear()
    #
    #    fig = plt.figure(figsize=(3,4))
    #    sb.heatmap(corpep .T ,cmap="YlGnBu");
    #    plt.yticks(np.arange(corpep.shape[1])+0.5,[pp[:17]if len(pp)<=17 else pp[:12]+'...' for pp in corpep.columns],fontsize=7)          
    #    plt.xticks(np.arange(2)+0.5,corpep.index,fontsize=9,rotation=0)          
    #    plt.xlabel('Ensembl transcript ID')
    #    
    #    plt.savefig('/Users/mirauta/Results/hipsci/manuscript_images/Transcripts/'+'fig3_peptide_'+gene_name+'_heatmap.png',dpi=900)
    #    plt.gcf().clear()
          
        if plot_peptides:
            fig = plt.figure(figsize=(4,3))
        #    sb.heatmap(corpep.T[:15].iloc[np.argsort(corpep.loc['ENST00000301843'][:15])].T  ,cmap="Blues",cbar_kws={"orientation": "horizontal"});
            sb.heatmap(corpep.T[:25].T  , cbar_kws={"orientation": "horizontal"}, vmin=-.5, vmax=0.5,center=0,cmap="RdBu_r")

            plt.xticks([5,corpep.T[5:15].T.shape[1]-0.5],['Shared peptides','Unique'],fontsize=9,rotation=0)          
            plt.yticks([ .5,1.5],np.hstack(['Isoform 1', 'Isoform2']),fontsize=9,rotation=0)          
            plt.ylabel('')
            plt.savefig('/Users/mirauta/Results/hipsci/manuscript_images/Transcripts/'+'fig3_peptide_'+gene_name+'_heatmap2.png',dpi=900,bbox_inches='tight')
            plt.savefig('/Users/mirauta/Results/hipsci/manuscript_images/Transcripts/'+'fig3_peptide_'+gene_name+'_heatmap2.svg',dpi=900,bbox_inches='tight')
            
            fig = plt.figure(figsize=(4,2))
            ax2 = fig.add_subplot(111);                                                   
            ax2.spines['right'].set_visible(False); ax2.spines['top'].set_visible(False);  ax2.spines['left'].set_visible(False); 
        #    sb.heatmap(corpep.T[:15].iloc[np.argsort(corpep.loc['ENST00000301843'][:15])].T  ,cmap="Blues",cbar_kws={"orientation": "horizontal"});
            
            betas=[data['Peptide'][gene]['data/beta'][p][data['Peptide'][gene]['data/snp_id'][p].astype('U')==snp_id]for p in list(data['Peptide'][gene]['data/beta'].keys())]
            
            sb.swarmplot(betas, size=6)
            plt.plot((-.3,1),(0.0,0.0),'k--',lw=.5)
            plt.plot((0,0),(-0.5,0.5),'k--',lw=.5)
            plt.ylim(-0.2,0.2)
            plt.xlim(-0.3,1)
            plt.xlabel('pepQTL effect size', fontsize=16)
            ax2.annotate("shared by\nisoforms 1 & 2",  xy=(0.1, 0.1));
            ax2.annotate("unique to\nisoform 1",  xy=(0.7, 0.051));
            
#            plt.xticks([5,corpep.T[5:15].T.shape[1]-0.5],['Shared peptides','Unique'],fontsize=9,rotation=0)          
#            plt.yticks([ .5,1.5],np.hstack(['Isoform 1', 'Isoform2']),fontsize=9,rotation=0)          
            plt.ylabel('Peptides', fontsize=16)
            plt.tight_layout()
            
            plt.savefig('/Users/mirauta/Results/hipsci/manuscript_images/Transcripts/'+'fig3_peptide_'+gene_name+'_heatmap3.png',dpi=900,bbox_inches='tight')
            plt.savefig('/Users/mirauta/Results/hipsci/manuscript_images/Transcripts/'+'fig3_peptide_'+gene_name+'_heatmap3.svg',dpi=900,bbox_inches='tight')
        
        plt.show()
