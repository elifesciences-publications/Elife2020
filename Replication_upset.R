#new <- list of genes 1
#old2 <- list of genes affected in old paper
#newFull <- list of genes affected at any level

library(UpSetR)
path_data='/Users/mirauta/Results/hipsci/QTL_may2019/'
path_destination='/Users/mirauta/Results/hipsci/manuscript_images/upset/'
trait_labels=c('Protein7M', 'mRNA','Transcript','Peptide')
peptide<- read.delim(paste(path_data, 'peptide_protein_qtl_results_LD_VEP_gwas_agreement_peptides_all.txt',sep=''),sep='\t')
protein<- read.delim(paste(path_data, 'protein_mrna_qtl_results_LD_VEP_gwas_agreement_peptides_all.txt',sep=''),sep='\t')
mrna<- read.delim(paste(path_data, 'mrna_protein_qtl_results_LD_VEP_gwas_agreement_peptides_all.txt',sep=''),sep='\t')
transcript<- read.delim(paste(path_data, 'transcript_protein_qtl_results_LD_VEP_MAF_rsid_rsid_inframe.txt',sep=''),sep='\t')
# dim(mrna)
# 
x=peptide;
x=x[!is.na(x$peptide_eff_size+x$transcript_eff_size+x$protein_eff_size+x$mrna_eff_size),]
listInput <- list(Protein= x$ensembl_gene_id[(x$protein_nominal_p_value_bonf<0.01)&(x$peptide_eff_size*x$protein_eff_size>0)&(is.finite(x$protein_nominal_p_value_bonf))], 
                  mRNA =x$ensembl_gene_id[(x$mrna_nominal_p_value_bonf<0.01)&(x$peptide_eff_size*x$mrna_eff_size>0)&(is.finite(x$mrna_nominal_p_value_bonf))],
                  Transcript_isoform=x$ensembl_gene_id[(x$transcript_nominal_p_value_bonf<0.01)&(x$peptide_eff_size*x$transcript_eff_size>0)&(is.finite(x$transcript_nominal_p_value_bonf))],
                  Peptide= x$ensembl_gene_id)

png(paste(path_destination,'peptide_centric.png',sep=''), width = 320*1.7, height = 320*0.7)
a=upset(data=fromList(listInput), keep.order = T,set.metadata = NULL,
      matrix.color = "darkblue",
      main.bar.color = "darkblue", mainbar.y.label = "Replicated QTLs",
      mainbar.y.max = NULL, sets.bar.color = "gray23",
      sets.x.label = " # QTLs", point.size = 2.9, line.size = 0.7,
      mb.ratio = c(0.65, 0.35), expression = NULL, att.pos = NULL,
      att.color = "darkblue", order.by = c("freq", "degree"),
      decreasing = c(T, T),  text.scale=c(1.7, 1.7, 1.7, 1.25, 1.7,1.7))
a
dev.off()
x=peptide;
a0=(x$peptide_eff_size*x$mrna_eff_size>0)&(x$mrna_nominal_p_value_bonf<0.01)
a=((x$peptide_eff_size*x$mrna_eff_size<0)|(x$mrna_nominal_p_value_bonf>0.01))&(x$transcript_nominal_p_value_bonf>0)
b=(x$peptide_eff_size*x$transcript_eff_size>0)&(x$transcript_nominal_p_value_bonf<0.01)&(x$mrna_nominal_p_value_bonf>0)
c=(x$peptide_eff_size*x$transcript_eff_size>0)&(x$transcript_nominal_p_value_bonf<0.01)
sum(a0|c,na.rm=1)

sum(a0,na.rm=1)
sum(a&b,na.rm=1)
sum(a,na.rm=1)


x=protein;
x=x[!is.na(x$peptide_eff_size+x$transcript_eff_size+x$protein_eff_size+x$mrna_eff_size),]
listInput <- list(Protein= x$ensembl_gene_id, 
                  mRNA =x$ensembl_gene_id[(x$mrna_nominal_p_value_bonf<0.01)&(x$protein_eff_size*x$mrna_eff_size>0)&(is.finite(x$mrna_nominal_p_value_bonf))],
                  Transcript_isoform=x$ensembl_gene_id[(x$transcript_nominal_p_value_bonf<0.01)&(x$protein_eff_size*x$transcript_eff_size>0)&(is.finite(x$transcript_nominal_p_value_bonf))],
                  Peptide= x$ensembl_gene_id[(x$peptide_nominal_p_value_bonf<0.01)&(x$protein_eff_size*x$peptide_eff_size>0)&(is.finite(x$peptide_nominal_p_value_bonf))])
#sum(apply(fromList(listInput[c('mRNA','Transcript_isoform')]),FUN=sum,MARGIN=1)>0)

png(paste(path_destination,'protein_centric.png',sep=''), width = 320*1.7, height = 320*0.7)
a=upset(data=fromList(listInput), keep.order = T,set.metadata = NULL,
        matrix.color = "darkblue",
        main.bar.color = "darkblue", mainbar.y.label = "Replicated QTLs",
        mainbar.y.max = NULL, sets.bar.color = "gray23",
        sets.x.label = " # QTLs", point.size = 2.9, line.size = 0.7,
        mb.ratio = c(0.65, 0.35), expression = NULL, att.pos = NULL,
       order.by = c("freq", "degree"),
        decreasing = c(T, T),  text.scale=c(1.7, 1.7, 1.7, 1.25, 1.7,1.7))
a
dev.off()

x=protein;
a0=(x$protein_eff_size*x$mrna_eff_size>0)&(x$mrna_nominal_p_value_bonf<0.01)
a=((x$protein_eff_size*x$mrna_eff_size<0)|(x$mrna_nominal_p_value_bonf>0.01))&(x$transcript_nominal_p_value_bonf>0)
b=(x$protein_eff_size*x$transcript_eff_size>0)&(x$transcript_nominal_p_value_bonf<0.01)&(x$mrna_nominal_p_value_bonf>0)
c=(x$protein_eff_size*x$transcript_eff_size>0)&(x$transcript_nominal_p_value_bonf<0.01)
sum(c,na.rm=1)
sum(a0|c,na.rm=1)
sum(a0,na.rm=1)
sum(a&b,na.rm=1)
sum(a,na.rm=1)

x=mrna;
x=x[!is.na(x$peptide_eff_size+x$transcript_eff_size+x$protein_eff_size+x$mrna_eff_size),]
listInput <- list(mRNA= x$ensembl_gene_id, 
                  Protein =x$ensembl_gene_id[(x$protein_nominal_p_value_bonf<0.01)&(x$mrna_eff_size*x$protein_eff_size>0)&(is.finite(x$protein_nominal_p_value_bonf))],
                  Transcript_isoform=x$ensembl_gene_id[(x$transcript_nominal_p_value_bonf<0.01)&(x$transcript_eff_size*x$mrna_eff_size>0)&(is.finite(x$transcript_nominal_p_value_bonf))],
                  Peptide= x$ensembl_gene_id[(x$peptide_nominal_p_value_bonf<0.01)&(x$mrna_eff_size*x$peptide_eff_size>0)&(is.finite(x$peptide_nominal_p_value_bonf))])
png(paste(path_destination,'mrna_centric.png',sep=''), width = 320*1.7, height = 320*0.7)
a=upset(data=fromList(listInput), keep.order = T,set.metadata = NULL,
        matrix.color = "darkblue",
        main.bar.color = "darkblue", mainbar.y.label = "Replicated QTLs",
        mainbar.y.max = NULL, sets.bar.color = "gray23",
        sets.x.label = " # QTLs", point.size = 2.9, line.size = 0.7,
        mb.ratio = c(0.65, 0.35), expression = NULL, att.pos = NULL,
         order.by = c("freq", "degree"),
        decreasing = c(T, T),  text.scale=c(1.7, 1.7, 1.7, 1.25, 1.7,1.7))
a
dev.off()


x=transcript;

listInput <- list(Transcript_isoform= x$ensembl_gene_id, 
                  Protein =x$ensembl_gene_id[(x$protein_nominal_p_value_bonf<0.01)&(x$transcript_eff_size*x$protein_eff_size>0)&(is.finite(x$protein_nominal_p_value_bonf))],
                  mRNA=x$ensembl_gene_id[(x$mrna_nominal_p_value_bonf<0.01)&(x$transcript_eff_size*x$mrna_eff_size>0)&(is.finite(x$mrna_nominal_p_value_bonf))],
                  Peptide= x$ensembl_gene_id[(x$peptide_nominal_p_value_bonf<0.01)&(x$transcript_eff_size*x$peptide_eff_size>0)&(is.finite(x$peptide_nominal_p_value_bonf))])
png(paste(path_destination,'transcript_isoform_centric.png',sep=''), width = 320*1.7, height = 320*0.7)
a=upset(data=fromList(listInput), keep.order = T,set.metadata = NULL,
        matrix.color = "darkblue",
        main.bar.color = "darkblue", mainbar.y.label = "Replicated QTLs",
        mainbar.y.max = NULL, sets.bar.color = "gray23",
        sets.x.label = " # QTLs", point.size = 2.9, line.size = 0.7,
        mb.ratio = c(0.65, 0.35), expression = NULL, att.pos = NULL,
         order.by = c("freq", "degree"),
        decreasing = c(T, T),  text.scale=c(1.7, 1.7, 1.7, 1.25, 1.7,1.7))
#c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
a        
dev.off()

x=transcript;
x=x[!is.na(x$peptide_eff_size+x$transcript_eff_size+x$protein_eff_size+x$mrna_eff_size),]
b=(x$protein_eff_size*x$transcript_eff_size>0)&(x$protein_nominal_p_value_bonf<0.01)
c=(x$peptide_eff_size*x$transcript_eff_size>0)&(x$peptide_nominal_p_value_bonf<0.01)
d=(x$mrna_eff_size*x$transcript_eff_size>0)&(x$mrna_nominal_p_value_bonf<0.01)
sum(b,na.rm=1)
sum(c,na.rm=1)
sum(b|c,na.rm=1)
sum(c&(1-b)&(1-d),na.rm=1)
sum(a0|c,na.rm=1)
sum(a0,na.rm=1)
sum(a&b,na.rm=1)
sum(a,na.rm=1)