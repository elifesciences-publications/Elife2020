import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge, Polygon,Rectangle
import seaborn as sb
#import vcf as vcf
#import h5py
#import glob as glob
#import subprocess
import pickle
import ast
import pexpect
 
def split0(x,poz): return x.split('_')[poz]
def split1(x,split,poz,join,replace="",poz0=0): return (join).join(x.replace(replace,"").split(split)[poz0:poz])
def inlist(x,searchlist): return any(np.in1d(x,searchlist))

def get_chromosome(chr,gwas,df ):
       
    dfchrom={}
    dfchrom['gwas']=gwas.query('chr==@chr')
    dfchrom['gwas']['snp_pos']=dfchrom['gwas']['snp_pos'].astype(float)
    for i,trait in enumerate(df.keys()):
        dfchrom[trait]=df[trait].query('chromosome==@chr')
        dfchrom[trait]['snp_pos']=dfchrom[trait]['snp_pos'].astype(float)
#    bimchrom=bim.query('chrom==@chr')
#    gwasqtlchrom=qtlinldfile.query('chrom==@chr')
    dfchrom['gwas'].index= dfchrom['gwas']['snp_chr_pos'].astype('U')
    return dfchrom
def get_chromosome_region(snp_gene,dfchrom,region):
    poz=int(snp_gene.split('_')[1])
    chrompoz={}
    chrompoz['gwas']=dfchrom['gwas'].query('snp_pos< (@poz+@region*10**5) and snp_pos> (@poz-@region*10**5)')
    
    for i,trait in enumerate(dfchrom.keys()):
        
        chrompoz[trait]=dfchrom[trait].query('snp_pos< (@poz+@region*10**5) and snp_pos> (@poz-@region*10**5)')
#        print (trait)
#        print (chrompoz[trait].shape)
        
    return chrompoz
 
def write_snp(file_sum_stat_name,snp_gene,chrompoz,gene,path_destination,traits, snp_df,snp_flag=True):
    
    snps=chrompoz['gwas'].index
#    for trait in traits:
##        print (chrompoz[trait].query('gene_ensembl_id==@gene').shape)
#        snps=np.intersect1d(snps,chrompoz[trait].query('ensembl_gene_id==@gene').index)

    if snp_flag: 
        
        temp=snp_df.T[(snp_df==2).sum(0)>(snp_df==0).sum(0)]
        temp[temp==2]=3;        temp[temp==0]=2;temp[temp==3]=0;
        snp_df[temp.index]=temp.T    
        ld=snp_df.corr();print (ld.shape)
        ld.to_csv(path_destination+file_sum_stat_name+'_ld',sep='\t',index=0,header=0)
    else: snp_df=None 
    chrompozgene={}
    for trait in chrompoz.keys():
        chrompozgene[trait]=chrompoz[trait].loc[snps]


    for trait in chrompozgene.keys():
        chrompozgene[trait]['zscore'].to_csv(path_destination+file_sum_stat_name+'_zscore_'+trait,sep='\t',header=0)
    for trait in chrompozgene.keys():
        try:chrompozgene[trait]['logp_value']
        except:chrompozgene[trait]['logp_value']=   -np.log10(chrompozgene[trait]['pvalue'])

    for trait in chrompozgene.keys():
        chrompozgene[trait][['beta','zscore','logp_value']].to_csv(path_destination+file_sum_stat_name+'_data_'+trait,sep='\t',header=1)
 
    return [file_sum_stat_name,snps,snp_df,chrompozgene]

def get_snp_df_from_bed(bed,fam,bim,snps):
    
    bim2=bim.loc[snps]
    bim2=bim2.query('i>0')
    snp_idxs =bim2['i'].values
    snp_names = bim2['snp'].values
    snp_df=pd.DataFrame(data=bed[snp_idxs,:].compute().transpose() ,index=fam.index,columns=snp_names)
    return  snp_df 

#def get_qtl_data_for_colocalisation(parameters):
#
#    
#    return [qtl,df,sample2individual_df,bed,fam,bim,qtl_snps,qtl_snp_genes]

def get_summary_stat_data(parameters,                          bed_flag=False):
                          
    '''load gwas and intersect with qtl variants and genes if available'''
    def split0(x,poz): return x.split('_')[poz]
    def split01(x,poz): return x.split('.')[poz]
    def split1(x,split,poz,join,replace="",poz0=0): return (join).join(x.replace(replace,"").split(split)[poz0:poz])
    def inlist(x,searchlist): return any(np.in1d(x,searchlist))

    fields_sum_stats1=ast.literal_eval(parameters.loc['fields_sum_stats1','value'])
    fields_sum_stats2=ast.literal_eval(parameters.loc['fields_sum_stats2','value'])
    '''get variant data'''
    if bed_flag:
        bim=pd.read_table(parameters.loc['sum_stats1_folder','value']+'hipsci_202_bim.txt.gz',index_col=0,sep=',')
        with open(parameters.loc['sum_stats1_folder','value']+'hipsci_202_bed.txt.pickle', 'rb') as handle:
            bed = pickle.load(handle)
        fam=pd.read_table(parameters.loc['sum_stats1_folder','value']+'hipsci_all_fam.txt.gz',index_col=0,sep=',')

    '''get QTL data'''
    gene_positions=pd.read_table(parameters.loc['genes_file','value'],sep='\t').set_index('ensembl_gene_id',drop=0)
    
    df={}
    for i,q in enumerate(['stats1' ]): 
        df[q]=pd.read_table(parameters.loc['sum_stats1_folder','value']+parameters.loc['sum_stats1_file','value'], sep='\t',compression='gzip',index_col=0)
        
        for f in fields_sum_stats1.keys(): 
            try:df[q][f]=df[q][fields_sum_stats1[f]]
            except:1
#            
#        df1=df[q].set_index('feature_id',drop=0)
#        df1=df1.loc[gene_positions['feature_id']]
#        df1.to_csv(parameters.loc['sum_stats1_folder','value']+"significant_lead_"+parameters.loc['sum_stats1_file','value'], sep='\t',compression='gzip' )        
    
        if fields_sum_stats1['snp_id'] is not None:
            df[q]['ref']=df[q]['snp_id'].apply(split0,args={2})
            df[q]['alt']=df[q]['snp_id'].apply(split0,args={3})
            df[q]['snp_pos']=df[q]['snp_id'].apply(split0,args={1}).astype(int)
            df[q]['snp_chr_pos']=df[q]["chr"].astype('U')+"_"+df[q][ "snp_pos"].astype('U')
            df[q]['snp_chr_pos_gene']=df[q][ "snp_chr_pos"] +"@"+    df[q]['ensembl_gene_id']

        

        df[q]['log_pvalue']=-np.log10(df[q]['pvalue'])
        df[q]['zscore']=df[q]['beta']/df[q]['beta_se']
        print (np.unique(df[q].index).shape)
            
    
    for i,q in enumerate(df.keys()):
        if fields_sum_stats1['ensembl_gene_id'] is not None:  
            df[q]=df[q].set_index('gene_ensembl_id',drop=0)
            df[q]=df[q].loc[np.intersect1d(df[q]['ensembl_gene_id'],gene_positions['ensembl_gene_id'])]
    
    
        
    gwas=pd.read_table(parameters.loc['sum_stats2_folder','value']+parameters.loc['sum_stats2_file','value'])

    for k in fields_sum_stats2.keys():
        try:gwas[k]=gwas[fields_sum_stats2[k]]
        except:1
   
    if fields_sum_stats2['chr'] is None:
        try:
            def fun (x):return x.split('_')[0]
            gwas['chr']=gwas['snp_id'].apply(fun)
            def fun (x):return x.split('_')[1]
            gwas['snp_pos']=gwas['snp_id'].apply(fun)
            def fun (x):return x.split('_')[2]
            gwas['ref']=gwas['snp_id'].apply(fun)
            def fun (x):return x.split('_')[3]
            gwas['alt']=gwas['snp_id'].apply(fun)
        except:print ('fail _snp_id')
    def fun (x):
        try:return float(x)
        except:return  np.nan
    
    gwas['beta_se']= gwas['beta_se'].apply(fun)
    gwas['zscore']=gwas['beta']/gwas['beta_se']
    gwas=gwas[np.isfinite(gwas['zscore'])]
    gwas['snp_chr_pos']=gwas["chr"].astype('U')+"_"+gwas[ "snp_pos"].astype('U')
    gwas['snp_id']=gwas["chr"].astype('U')+"_"+gwas["snp_pos"].astype('U')+"_"+gwas["ref"].astype('U')+"_"+gwas["alt"].astype('U')
    
    gwas=gwas.set_index('snp_chr_pos',drop=0)
    gwas=gwas.loc[np.intersect1d(np.unique(gwas['snp_chr_pos']),np.unique(df[q]['snp_chr_pos']))];
    df[q]=df[q].set_index('snp_chr_pos',drop=0)
    df[q]=df[q].loc[np.intersect1d(np.unique(gwas['snp_chr_pos']),np.unique(df[q]['snp_chr_pos']))];
    print (gwas.shape)
    print (df[q].shape)

    if fields_sum_stats2['ensembl_gene_id'] is not None: gwas['ensembl_gene_id']=gwas['ensembl_gene_id'].apply(split01,args={0})
         
    gene_flag=False
    if fields_sum_stats2['ensembl_gene_id'] is not None:
        gene_flag=True
        
        gwas['snp_chr_pos_gene']=gwas[ "snp_chr_pos"] +"@"+    gwas['ensembl_gene_id']
        gwas=gwas.set_index('snp_chr_pos_gene',drop=0)
        df[q]=df[q].set_index('snp_chr_pos_gene',drop=0)
        gwas=gwas.loc[np.intersect1d(df[q].index,gwas.index)]        
        df[q]=df[q].loc[np.intersect1d(df[q].index,gwas.index)]
        print (gwas.shape)
        print (df[q].shape)
    gwas=gwas.set_index('snp_id',drop=0)
    df[q]=df[q].set_index('snp_id',drop=0)
    gene_positions=gene_positions.loc[np.intersect1d(gene_positions.index,np.unique(gwas['ensembl_gene_id']))] 
        

    chrs=np.unique([snp.split('_')[0] for snp in gene_positions['snp_id']]    )
    dfchromlist={}
    for chr in chrs:         
        
        try:
            dfchromlist[chr]=get_chromosome(chr=chr,gwas=gwas,df=df);           
#            print ("succeeded 1 :  "+ gene+" -"+file)
        except:continue   
    if bed_flag:
        bimchrom={}
        bim=bim.set_index('chrom',drop=0)
        for chr in chrs.astype(int):  
            print (chr)
            bimchrom[str(chr)]=bim.loc[chr].set_index('i',drop=0)  
            bimchrom[str(chr)]['snp']=bimchrom[str(chr)]['chrom_pos']
                    
    for gene in gene_positions.index:
#        sys.exit()
        chr=gene_positions.loc[gene,'snp_id'].split('_')[0]
        snp_gene=gene+'@'+gene_positions.loc[gene,'snp_id']
        
    
#        try:
        window=10
        chrompoz=get_chromosome_region(snp_gene,dfchrom=dfchromlist[chr],region=window)
#            if chrompoz['gwas'].shape[0]!=chrompoz['stats1'].shape[0]:
#                sys.exit()
        if fields_sum_stats2['ensembl_gene_id'] is not None:
            chrompoz['gwas']=chrompoz['gwas'].query('ensembl_gene_id==@gene') 
        if fields_sum_stats1['ensembl_gene_id'] is not None:
            chrompoz['stats1']=chrompoz['stats1'].query('ensembl_gene_id==@gene')            
        
        snps=chrompoz['gwas']['snp_chr_pos']
        for trait in np.setdiff1d(list(chrompoz.keys()),'gwas'):
            snps=np.intersect1d(snps,chrompoz[trait]['snp_chr_pos'])
        if snps.shape[0]<1:
            print ("no intersection data for gene"+gene)
            continue
        
        while snps.shape[0]>1000:
#            print ("restricting the analysis on a smaller window")
            window=window-0.2
            chrompoz=get_chromosome_region(snp_gene,dfchromlist[chr],region=window)
            if fields_sum_stats2['ensembl_gene_id'] is not None:
                chrompoz['gwas']=chrompoz['gwas'].query('ensembl_gene_id==@gene') 
            if fields_sum_stats1['ensembl_gene_id'] is not None:
                chrompoz['stats1']=chrompoz['stats1'].query('ensembl_gene_id==@gene') 
            snps=chrompoz['gwas']['snp_chr_pos']
            for trait in np.setdiff1d(list(chrompoz.keys()),'gwas'):
                snps=np.intersect1d(snps,chrompoz[trait].query('gene_ensembl_id==@gene')['snp_chr_pos'])
        window=window+0.2
        chrompoz=get_chromosome_region(snp_gene,dfchromlist[chr],region=window)
        if fields_sum_stats2['ensembl_gene_id'] is not None:
            chrompoz['gwas']=chrompoz['gwas'].query('ensembl_gene_id==@gene') 
        if fields_sum_stats1['ensembl_gene_id'] is not None:
            chrompoz['stats1']=chrompoz['stats1'].query('ensembl_gene_id==@gene') 
        snps=chrompoz['gwas']['snp_chr_pos']
        for trait in np.setdiff1d(list(chrompoz.keys()),'gwas'):
            snps=np.intersect1d(snps,chrompoz[trait].query('gene_ensembl_id==@gene')['snp_chr_pos'])
        for k in chrompoz.keys():
            chrompoz[k]=chrompoz[k].set_index('snp_chr_pos',drop=0)#            print ("succeeded 2 : "+ gene+" -"+file)
#        sys.exit()
#        except:
#            continue
 
#        try:
        if bed_flag: snp_df=get_snp_df_from_bed(bed,fam,bimchrom[str(chr)],snps=snps)
        [name,snpsgene,snp_df,chrompozgene]=write_snp(file_sum_stat_name=gene+"_"+parameters.loc['sum_stats1_file','value']+'__'+parameters.loc['sum_stats2_file','value'],\
        snp_gene=snp_gene,chrompoz=chrompoz,gene=gene,path_destination=parameters.loc['data_output_folder','value'],  traits=np.setdiff1d(list(chrompoz.keys()),'gwas'),\
                                                                                     snp_df=None,snp_flag = bed_flag)
        print ('wrote gene: '+ gene)
#        except:
#            continue



def get_summary_stat_data_ld(parameters, genes_traits):

    def split1(x): return x.split('_')[1]   
    bim=pd.read_table(parameters.loc['sum_stats1_folder','value']+'hipsci_202_bim.txt.gz',index_col=0,sep=',')
    with open(parameters.loc['sum_stats1_folder','value']+'hipsci_202_bed.txt.pickle', 'rb') as handle:
        bed = pickle.load(handle)
    fam=pd.read_table(parameters.loc['sum_stats1_folder','value']+'hipsci_all_fam.txt.gz',index_col=0,sep=',')
    fam202=pd.read_table(parameters.loc['sum_stats1_folder','value']+'hipsci_202_fam.txt.gz',index_col=0,sep=',')

    gene_positions=pd.read_table(parameters.loc['genes_file','value'],sep='\t').set_index('ensembl_gene_id',drop=0)
    chrs=np.unique([snp.split('_')[0] for snp in gene_positions['snp_id']]    )
    bimchrom={}
    bim=bim.set_index('chrom',drop=0)
    for chr in chrs.astype(int):  
        print (chr)
        bimchrom[str(chr)]=bim.loc[chr].set_index('i',drop=0)  
        bimchrom[str(chr)]['snp']=bimchrom[str(chr)]['chrom_pos']
                    
 
    for ig,file in  enumerate(genes_traits.index) :

        print (ig)
        
        chrompozgene={}
        try:  
            for k in ['gwas','stats1']:
                chrompozgene[k] = pd.read_table(parameters.loc['data_output_folder','value'] + file.replace("*","_data_")+k)
#                chrompozgene[k]['snp_id2']=str(gene_positions.loc[genes_traits.loc[file,'genes'],'chromosome'])+"_"+chrompozgene[k]['snp_chr_pos'].astype('U')
                chrompozgene[k]['snp_id']=gene_positions.loc[genes_traits.loc[file,'genes'],'snp_id']                                
#                chrompozgene[k]['snp_chr_pos']=chrompozgene[k].index;chrompozgene[k].index=chrompozgene[k]['snp_chr_pos'].apply(split1).astype(int)
        except: continue    
        
        chr=gene_positions.loc[genes_traits.loc[file,'genes'],'chromosome']
#        snp_gene=gene+'@'+gene_positions.loc[gene,'snp_id']
        
        snps=chrompozgene[k]['snp_chr_pos']
        try:
            snp_df=get_snp_df_from_bed(bed,fam,bimchrom[str(chr)].set_index('snp',drop=0),snps=snps) .loc[fam202.index]
        except: 
            chrompozgene[k]['snp_chr_pos']=str(gene_positions.loc[genes_traits.loc[file,'genes'],'chromosome'])+"_"+chrompozgene[k]['snp_chr_pos'].astype('U')
            snps=chrompozgene[k]['snp_chr_pos']
            snp_df=get_snp_df_from_bed(bed,fam,bimchrom[str(chr)].set_index('snp',drop=0),snps=snps) .loc[fam202.index]
        
        temp=snp_df.T[(snp_df==2).sum(0)>(snp_df==0).sum(0)]
        temp[temp==2]=3;        temp[temp==0]=2;temp[temp==3]=0;
        snp_df[temp.index]=temp.T    
        ld=snp_df.corr();print (ld.shape)
        ld.to_csv(parameters.loc['data_output_folder','value'] + file.replace("*","_ld"),sep='\t',index=0,header=0)
        
        print ('wrote gene: '+ file)
#        except:
#            continue



def plot_manhattan_scatter_overlap(parameters,gene_positions, genes, \
                                   threshold_replication=2,field_manhattan='logp_value',\
                                   field_scatter='beta'):
    traitcolor={'stats1':"darkblue",'eQTL':"green","gwas":"salmon"}
    plotsnpcolor={'Coloc.':"--k",'Missense':":k","QTL":"-k"}
    plotsnpcolor2={'Coloc.':"*",'Missense':"^","QTL":"o"};
    def split1(x): return x.split('_')[1]   
    
 
    for gene in  genes :
        plotsnp={}
        plotsnp['QTL']=np.array(gene_positions.loc[gene,'snp_id'].split('_')[:2]).astype(int)[1]
         
        try:
            pavsnp=np.array([int(np.array(s.split('_'))[1]) if len(s)>0 else np.nan for s in gene_positions.loc[gene,'inframe'].split(';')])
            plotsnp['Missense']=pavsnp[pavsnp==pavsnp] .astype(int) 
        except:1
    
        file=gene+"_"+parameters.loc['sum_stats1_file','value']+'__'+parameters.loc['sum_stats2_file','value']
        chrompozgene={}
        try:  
            for k in ['gwas','stats1']:
                chrompozgene[k]=pd.read_table(parameters.loc['data_output_folder','value'] + file  +"_data_"+k,index_col=0)
                chrompozgene[k]['index']=chrompozgene[k].index;chrompozgene[k].index=chrompozgene[k]['index'].apply(split1).astype(int)
        except:
            try:
                for k in ['gwas','stats1']:
                    file=str(gene_positions.loc[gene,'chromosome'])+"_"+str(int(gene_positions.loc[gene,'position']))+"---"+\
                            parameters.loc['sum_stats2_file','value']
                    chrompozgene[k]=pd.read_table(parameters.loc['sum_stats2_folder','value'] + "/"+file  +"__data_"+k+"_"+gene,index_col=0)
                    chrompozgene[k]['index']=chrompozgene[k].index;chrompozgene[k].index=chrompozgene[k]['index'].apply(split1).astype(int)
            except: continue    
            
        '''contine only if nominal significnat in gwas'''
        try:  
            if chrompozgene['gwas'].loc[plotsnp['QTL'],'logp_value']<threshold_replication: continue
        except: print (gene)
        
        field=field_manhattan   
    #        field='zscore'
        plt.figure(figsize=(7,4))
        a=plt.subplot(1,1,1);    a.spines['top'].set_visible(False);a.spines['right'].set_visible(False);
        x={}
        for k in chrompozgene.keys():            x[k]=chrompozgene[k][[field]]
    #        for k in chrompozgene.keys():            x[k]=chrompozgene[k][[field]]/chrompozgene[k][field].values.std()
        maxval=max([max(x[k][field]) for k in chrompozgene.keys()])
        a.add_patch(Rectangle((gene_positions.loc[gene,'gene_start'], 0), gene_positions.loc[gene,'gene_end']-gene_positions.loc[gene,'gene_start'], maxval, facecolor="grey",    alpha=0.35))
        
        start=True
        for s in plotsnp.keys():
            i=plotsnp[s]
            if s=="Missense":
                print (i)
                if start: 
                    plt.plot((i[0],i[0]),(0,maxval),plotsnpcolor[s],label=s,lw=0.5)
                    for j in np.arange(1,len(i)):
                        plt.plot((i[j],i[j]),(0,maxval),plotsnpcolor[s],lw=0.5)
                    start=False
                else: 
                    for j in np.arange(0,len(i)):
                        plt.plot((i[j],i[j]),(0,maxval),plotsnpcolor[s],lw=0.5)
           
                
            else:plt.plot((i,i),(0,maxval),plotsnpcolor[s],label=s,lw=0.5)
     
        for k in chrompozgene.keys():
            a.plot(chrompozgene[k].index,x[k],'o',color=traitcolor[k],label=k if k!='gwas' else "GWAS", markersize=4)
            
        shift=(chrompozgene[k].index.max()-chrompozgene[k].index.min())
        plt.xlim(chrompozgene[k].index.min()-shift*0.05,chrompozgene[k].index.max()+shift*0.15)
        plt.xticks(chrompozgene[k].index[[0,int(chrompozgene[k].shape[0]/2),chrompozgene[k].shape[0]-1]],\
                   [j+' Mb' for j in np.around(chrompozgene[k].index[[0,int(chrompozgene[k].shape[0]/2),chrompozgene[k].shape[0]-1]]/10**6,1).values.astype('U')])
        plt.ylabel("-log10 PV",fontsize=16)
        plt.xlabel('Genome position',fontsize=16)
    #        plt.title(str(chrompozgene[k].shape[0])+"_"+str(np.around(gene_positions.loc[gene,file],2))        +"_"+str(gene_positions.loc[gene,'gwas_trait_all_gwas'])[:10]+"\n"+gene_positions.loc[gene,'gene_name']+"_"+file.replace('_colocpv_pQTL',''))
        plt.title(gene_positions.loc[gene,'gene_name']+"\n"+parameters.loc['sum_stats2_file','value'].replace('.txt.gz_filtered_genes',''),fontsize=16)
        plt.tight_layout()
        plt.legend(loc=2,fontsize=12)
        plt.savefig(parameters.loc['data_output_figures_folder','value']+gene_positions.loc[gene,'gene_name']\
                     +"_"+parameters.loc['sum_stats2_file','value'].replace('.txt.gz_filtered_genes','')+'_manhattan.png',dpi=600)   
    
        plt.figure(figsize=(3.25,3.25))
        a=plt.subplot(1,1,1); a.spines['top'].set_visible(False);a.spines['right'].set_visible(False);
        field=field_scatter 
        x={}
    #        for k in chrompozgene.keys():            x[k]=chrompozgene[k][[field]]/chrompozgene[k][field].values.std()
        for k in chrompozgene.keys():            x[k]=chrompozgene[k][[field]]
        plt.plot((x['gwas'].values.min(),x['gwas'].values.max()),(0,0),'k-',lw=0.425)
        try:plt.plot( (0,0),(min(x['stats1'].values.min(),x['eQTL'].values.min()),max(x['stats1'].values.max(),x['eQTL'].values.max())),'k-',lw=0.425)
        except:plt.plot( (0,0),(min(x['stats1'].values.min(),x['stats1'].values.min()),max(x['stats1'].values.max(),x['stats1'].values.max())),'k-',lw=0.425)
        
        for k in np.setdiff1d(np.array(list(chrompozgene.keys())),'gwas'):
            plt.plot(x['gwas'].iloc[np.unique(x['gwas'].index,return_index=1)[1]],\
                     x[k].iloc[np.unique(x[k].index,return_index=1)[1]],'o',color=traitcolor[k],markersize=3,label=k)
    #            plt.plot(x['gwas'].loc[colocsnp],x[k].loc[colocsnp],'o',color=traitcolor[k],markersize=10 )
    #            for s in plotsnp.keys():
    #                try:plt.plot(x['gwas'].loc[plotsnp[s]],x[k].loc[plotsnp[s]],plotsnpcolor2[s],color="red",markersize=10,mfc=traitcolor[k],label='_nolegend_')
    #                except:1
    #                plt.show()
            for s in plotsnp.keys():
                if any(np.in1d('Missese',list(plotsnp.keys()))):
                    try:plt.plot(x['gwas'].loc[plotsnp[s]],x[k].loc[plotsnp[s]],plotsnpcolor2[s],color="red",markersize=10,mfc=traitcolor[k],label=s )
                    except:1
                else: 
                    try:plt.plot(x['gwas'].loc[plotsnp[s]],x[k].loc[plotsnp[s]],plotsnpcolor2[s],color="red",markersize=10,mfc=traitcolor[k],label=s,lw=3 )
                    except:1
                    
    #        plt.plot()
    #                plt.show()
            
    #        plt.xticks(chrompozgene[k].index[[0,int(chrompozgene[k].shape[0]/2),chrompozgene[k].shape[0]-1]],\
    #                   [j+' kb' for j in np.around(chrompozgene[k].index[[0,int(chrompozgene[k].shape[0]/2),chrompozgene[k].shape[0]-1]]/10**6,1).values.astype('U')])
        plt.ylabel("Effect size QTL",fontsize=14)
        plt.xlabel("Effect size GWAS",fontsize=14)
        
    
    #        box = a.get_position()
    #        a.set_position([box.x0, box.y0, box.width, box.height])
    #        a.legend(loc='upper center', bbox_to_anchor=(0.2, 1.05),ncol=1,fontsize=9.3)
        
    #        plt.legend(loc='upper center',ncol=3,fontsize=10)
    #        plt.title(str(chrompozgene[k].shape[0])+"_"+str(np.around(qtlpr.loc[gene,file],2))+"\n"+qtlpr.loc[gene,'gene_name']+"_"+file.replace('_colocpv_pQTL',''))
    #        print (qtlpr.loc[gene][:80])
    #        print (qtlpr.loc[gene,['snp_id','inframe']])
        plt.tight_layout()
        plt.savefig(parameters.loc['data_output_figures_folder','value']+gene_positions.loc[gene,'gene_name']\
                     +"_"+parameters.loc['sum_stats2_file','value'].replace('.txt.gz_filtered_genes','')+'_scatter.png',dpi=600)   
    
#        plt.show()

def scp(user2,host2,tgt,pwd,local,opts='',timeout=30):
    ''' Performs the scp command. Transfers file(s) from local host to remote host '''
    cmd = f'scp {user2}@{host2}:{tgt} {local}'
    print("Executing the following cmd:",cmd,sep='\n')

    tmpFl = '/tmp/scp.log'
    fp = open(tmpFl,'wb')
    childP = pexpect.spawn(cmd,timeout=timeout)
    try:
        childP.sendline(cmd)
        childP.expect([f"{user2}@{host2}'s password:"])
        childP.sendline(pwd)
        childP.logfile = fp
        childP.expect(pexpect.EOF)
        childP.close()
        fp.close()

        fp = open(tmpFl,'r')
        stdout = fp.read()
        fp.close()

        if childP.exitstatus != 0:
            raise Exception(stdout)
    except KeyboardInterrupt:
        childP.close()
        fp.close()
        return

    print(stdout)    

def run_ecaviar(path_files,files,ecaviar_path="/Users/mirauta/Applications_additional/caviar/CAVIAR-C++/eCAVIAR"):
    for file in files:
#        print(file)
        name=file.replace('*','')
            
        try:
            print ("Run on:  "+ path_files + name + "_zscore_gwas"+"\n and\n"+path_files + name + "_zscore_"+name_qtl_file)        
            command_name=[ecaviar_path,\
                          "-o", path_files + name + "_out",\
                          "-l",path_files + name + "_ld",\
                          "-l",path_files + name + "_ld" ,\
                          "-z",path_files + name + "_zscore_gwas",\
                          "-z",path_files + name + "_zscore_"+name_qtl_file,"-c","2"]
#            print (' '.join(command_name))
            subprocess.run(command_name, stdout=subprocess.PIPE)
           
        except: 
            print ("Failed eCAVIAR for:  "+ name)
            continue
