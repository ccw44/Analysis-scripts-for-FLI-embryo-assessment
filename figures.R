#!/usr/bin/env Rscript
# Every proteomics figure, produced from the outputs of analysis.R and
# enrichment.R.
#
#   Rscript figures.R results/proteomics
#
# Sections, in output order:
#   1  Fig. 2B, GSEA panels, per-protein supplements   - matrix reanalysis
#   2  Fig. 2C                                        - provider result review
#   3  supplementary ORA and GSEA panels              - provider result review
#   4  exploratory ORA, raw P < 0.05 and |log2FC| >= 1
#   5  exploratory volcano for the same foreground
#
# Sections 2 and 3 read the provider candidate branch; sections 1, 4 and 5 read the
# matrix reanalysis. Each keeps its own candidate set and enrichment background.
args<-commandArgs(TRUE);stopifnot(length(args)==1);b<-normalizePath(args[1])
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})

s<-readRDS(file.path(b,'statistics/analysis_state.rds'))
e<-readRDS(file.path(b,'enrichment/enrichment_state.rds'))

put<-function(x,n)fwrite(x,file.path(b,'source_data',paste0(n,'.tsv')),sep='\t',na='NA')
put_rel<-function(x,p)fwrite(x,file.path(b,p),sep='\t',na='NA')
savefig<-function(p,n,w,h){
 grDevices::cairo_pdf(file.path(b,'figures',paste0(n,'.pdf')),width=w/25.4,height=h/25.4,family='Arial');print(p);dev.off()
}


# ===========================================================================
# 1  Fig. 2B, GSEA panels and per-protein supplements (matrix reanalysis)
# ===========================================================================
theme_set(theme_classic(base_size=8,base_family='Arial')+theme(plot.title=element_text(size=9,face='bold'),plot.subtitle=element_text(size=7),legend.position='bottom',legend.text=element_text(size=7)))
d<-copy(s$main);d[,status_plot:='Not significant'];d[significant_FDR==TRUE,status_plot:=ifelse(log2FC>0,'HS higher','LS higher')]
p<-ggplot(d,aes(log2FC,-log10(pmax(padj,1e-300)),color=status_plot))+geom_point(size=.6,alpha=.6)+geom_hline(yintercept=-log10(.05),linetype=2,linewidth=.3)+geom_vline(xintercept=c(-1,1),linetype=3,linewidth=.3)+ggrepel::geom_text_repel(data=d[significant_FDR==TRUE],aes(label=gene_symbol),size=2.7,seed=20261031,max.overlaps=Inf,min.segment.length=0,show.legend=FALSE)+scale_color_manual(values=c('Not significant'='#9C9C9C','HS higher'='#B46B33','LS higher'='#39799B'))+labs(x='log2 fold change (HS / LS)',y='-log10 BH-adjusted P',color='FDR < 0.05',title='Protein abundance differences',subtitle=sprintf('%d tested | Differential proteins: %d | HS higher: %d | LS higher: %d',nrow(d),sum(d$significant_FDR),sum(d$significant_FDR & d$log2FC>0),sum(d$significant_FDR & d$log2FC<0)))
savefig(p,'Fig2B_protein_volcano',140,105);put(d,'Fig2B_volcano')
k<-copy(e$GSEA[['KEGG']]);k<-k[padj<.05];setorder(k,NES)
k[,label:=vapply(Description,function(x)paste(strwrap(sub(' - Mus musculus.*','',x),width=32),collapse='\n'),character(1))];k[,label:=factor(label,levels=unique(label))]
p<-ggplot(k,aes(NES,label,size=size,color=padj))+geom_vline(xintercept=0,linetype=3,linewidth=.3)+geom_point()+scale_size_continuous(range=c(3,5),name='Gene-set size')+scale_color_gradient(low='#B46B33',high='#39799B',name='BH adjusted P')+labs(x='NES (positive: higher in HS)',y=NULL,title='KEGG GSEA',subtitle='All FDR-significant sets; measured-gene ranking')+theme(legend.position='bottom')
savefig(p,'Fig2C_KEGG_GSEA_primary',170,95);put(k,'Fig2C_KEGG_GSEA_primary')
# The old candidate-list ORA is supplementary, even after background correction.
src<-file.path(b,'figures');for(ext in c('pdf','svg','png')){
 from<-file.path(src,paste0('Fig2C_exploratory_KEGG.',ext));to<-file.path(src,paste0('Supplementary_legacy_candidate_KEGG_ORA.',ext));if(file.exists(from))stopifnot(file.rename(from,to))
}
rows<-d[significant_FDR==TRUE,protein_group];mat<-s$log2[rows,,drop=FALSE]
long<-as.data.table(mat,keep.rownames='protein_group');long<-melt(long,id.vars='protein_group',variable.name='sample_id',value.name='log2_intensity');long<-merge(long,d[,.(protein_group,gene_symbol)],by='protein_group');long[,group:=ifelse(substr(sample_id,1,1)=='H','HS','LS')]
put(long,'significant_protein_sample_values')
p<-ggplot(long,aes(group,log2_intensity,color=group,shape=group))+geom_point(position=position_jitter(width=.07,height=0,seed=20261031),size=1.7)+facet_wrap(~gene_symbol,nrow=1,scales='free_y')+scale_color_manual(values=c(HS='#B46B33',LS='#39799B'))+labs(x=NULL,y='Observed log2 intensity',title='Individual values for FDR-significant protein groups',subtitle='Each point represents one independent embryo; n = 3 per group')+theme(legend.position='none',strip.text=element_text(size=8))
savefig(p,'Supplementary_significant_protein_values',183,75)
ci<-d[significant_FDR==TRUE];ci[,gene_symbol:=factor(gene_symbol,levels=gene_symbol[order(log2FC)])];put(ci,'significant_protein_CIs')
p<-ggplot(ci,aes(log2FC,gene_symbol))+geom_vline(xintercept=0,linetype=3,linewidth=.3)+geom_segment(aes(x=CI95_low,xend=CI95_high,yend=gene_symbol),linewidth=.6,color='#39799B')+geom_point(size=2,color='#39799B')+labs(x='log2 fold change (HS / LS), moderated 95% CI',y=NULL,title='Effect sizes for FDR-significant protein groups',subtitle='Marginal model intervals; not simultaneous selection-adjusted intervals')
savefig(p,'Supplementary_significant_protein_CI',150,95)
for(fam in c('GO_CC','GO_MF')){
 g<-copy(e$GSEA[[fam]]);g<-head(g[padj<.05][order(padj,pval)],15)
 g[,label:=vapply(Description,function(x)paste(strwrap(x,width=50),collapse='\n'),character(1))];g[,label:=factor(label,levels=rev(unique(label)))];put(g,paste0(fam,'_GSEA_display'))
 p<-ggplot(g,aes(NES,label,color=padj,size=size))+geom_vline(xintercept=0,linetype=3,linewidth=.3)+geom_point()+scale_color_gradient(low='#B46B33',high='#39799B',name='BH adjusted P')+labs(x='NES (positive: higher in HS)',y=NULL,title=paste(fam,'GSEA'),subtitle='Top 15 FDR-significant sets; complete results supplied')
 savefig(p,paste0('Supplementary_',fam,'_GSEA'),183,180)
}
writeLines(c('Figure2B: protein-group BH FDR, all 5359 points; labels mark all five significant groups.',
 'Figure2C: three KEGG GSEA sets significant in the new prespecified reanalysis, not old ORA results.',
 'Protein-value points are raw independent sample summaries, with no aggregate uncertainty implied.',
 'Effect intervals are marginal moderated 95% CIs, not simultaneous confidence intervals after selection.',
 'QC: PCA/correlation use the 5028 complete groups without imputation; inclusion counts are in audit/eligibility.tsv.',
 'Heatmaps: all 85 legacy candidates on three pages; row z-scores; gray missing; columns fixed H1-H3/L1-L3.',
 'Final exports are PDF only with editable text; temporary previews are not deliverables.',
 'Static preflight missing-count/uncertainty warnings are addressed by explicit source tables and point/CI definitions.'),file.path(b,'audit/figure_contract_and_QA_notes.txt'))
cat('FINAL FIGURES COMPLETED\n')


# ===========================================================================
# 2  Fig. 2C as published (provider candidate set)
# ===========================================================================
z<-copy(e$ORA[['legacy85_all KEGG']]);z<-head(z[Count>0][order(Pvalue,term)],15)
gt<-e$gene_statistics
counts<-rbindlist(lapply(seq_len(nrow(z)),function(i){
 ids<-strsplit(z$genes[i],';',fixed=TRUE)[[1]];x<-gt[ENTREZID%in%ids]
 stopifnot(nrow(x)==z$Count[i],all(is.finite(x$logFC)))
 data.table(term=z$term[i],HS_higher=sum(x$logFC>0),LS_higher=sum(x$logFC<0),equal=sum(x$logFC==0))
}))
stopifnot(all(counts$equal==0))
z<-merge(z,counts,by='term',sort=FALSE);setorder(z,Pvalue,term)
z[,label:=vapply(Description,function(t)paste(strwrap(sub(' - Mus musculus.*','',t),width=44),collapse='\n'),character(1))]
z[,label:=factor(label,levels=rev(unique(label)))];z[,minus_log10_P:=-log10(pmax(Pvalue,1e-300))]
upperP<-ceiling(max(z$minus_log10_P)*2)/2;upperCount<-max(z$Count)+1
axis_scale<-upperCount/upperP;z[,line_x:=minus_log10_P*axis_scale]
bars<-melt(z,id.vars=c('term','label'),measure.vars=c('HS_higher','LS_higher'),variable.name='direction',value.name='count')
bars[,direction:=factor(direction,levels=c('HS_higher','LS_higher'))]
fwrite(z,file.path(b,'source_data/Fig2C_ORA_reference_style.tsv'),sep='\t')
fwrite(bars,file.path(b,'source_data/Fig2C_ORA_reference_stacked_counts.tsv'),sep='\t')
theme_set(theme_classic(base_size=8,base_family='Arial')+theme(plot.title=element_text(size=10,face='bold'),plot.subtitle=element_text(size=7),axis.text=element_text(size=7),legend.position='bottom',legend.title=element_blank(),legend.text=element_text(size=8),plot.margin=margin(8,12,8,8)))
p<-ggplot()+geom_col(data=bars,aes(count,label,fill=direction),orientation='y',width=.65,position=position_stack(reverse=TRUE))+geom_path(data=z,aes(line_x,label,group=1),color='#222222',linewidth=.45)+geom_point(data=z,aes(line_x,label),shape=22,size=2.1,fill='#E8BD43',color='#222222',stroke=.4)+scale_fill_manual(values=c(HS_higher='#D94A38',LS_higher='#4B8DB9'),labels=c(HS_higher='Higher in HS',LS_higher='Higher in LS'))+scale_x_continuous(name='Mapped candidate genes (count)',limits=c(0,upperCount),breaks=seq(0,upperCount,by=1),expand=expansion(mult=c(0,.015)),sec.axis=sec_axis(~./axis_scale,name='-log10(P-value)',breaks=pretty(c(0,upperP),n=5)))+labs(y=NULL,title='Exploratory ORA–KEGG',subtitle='Bars: abundance direction of contributing candidate genes. Line/squares: unadjusted P-value.\nNo ORA term passed BH FDR < 0.05; top 15 nominal-P terms with hits are shown.')
grDevices::cairo_pdf(file.path(b,'figures/Fig2C_ORA_KEGG_reference_style.pdf'),width=183/25.4,height=165/25.4,family='Arial');print(p);dev.off()
writeLines(c('Final requested format: PDF only.',
 'Active Figure2C: Fig2C_ORA_KEGG_reference_style.pdf.',
 'Red/blue stacked counts denote contributing legacy candidate genes with positive/negative HS-minus-LS gene-level log2FC, not necessarily FDR-significant individual genes.',
 'Bottom x axis is count; upper axis is -log10 nominal ORA P; linear coordinate conversion is explicit in source data and code.',
 'Total count equals the sum of HS-higher and LS-higher counts; assertions passed.',
 'The original multi-panel ORA alternative is retained as supplementary PDF only.'),file.path(b,'audit/final_display_specification.txt'))
cat('REFERENCE-STYLE ORA PDF COMPLETED\n')


# ===========================================================================
# 3  Supplementary ORA and GSEA panels (provider candidate set)
# ===========================================================================
z<-copy(e$ORA[['legacy85_all KEGG']]);z<-head(z[Count>0][order(Pvalue,term)],15)
genes<-e$gene_statistics
z[,label:=vapply(Description,function(t)paste(strwrap(sub(' - Mus musculus.*','',t),width=38),collapse='\n'),character(1))]
z[,label:=factor(label,levels=rev(unique(label)))];z[,display_P:=pmax(Pvalue,1e-300)]
fwrite(z,file.path(b,'source_data/Fig2C_ORA_overview.tsv'),sep='\t')
edges<-rbindlist(lapply(seq_len(nrow(z)),function(i){data.table(term=z$term[i],label=z$label[i],ENTREZID=strsplit(z$genes[i],';',fixed=TRUE)[[1]])}))
edges<-merge(edges,genes[,.(ENTREZID,symbols,logFC,padj)],by='ENTREZID',all.x=TRUE,sort=FALSE)
colorder<-unique(edges[order(logFC,symbols),ENTREZID]);symbols<-setNames(genes$symbols,genes$ENTREZID)
edges[,ENTREZID:=factor(ENTREZID,levels=colorder)];edges[,label:=factor(label,levels=levels(z$label))]
fwrite(edges,file.path(b,'source_data/Fig2C_ORA_gene_contributions.tsv'),sep='\t')
theme_set(theme_classic(base_size=7,base_family='Arial')+theme(plot.title=element_text(size=8,face='bold'),plot.subtitle=element_text(size=7),legend.position='bottom',legend.text=element_text(size=6),axis.text=element_text(size=6)))
overview<-ggplot(z,aes(fold_enrichment,label,color=display_P,size=Count))+geom_point()+scale_color_gradient(low='#B46B33',high='#93B6C7',trans='log10',name='P-value')+scale_size_continuous(range=c(1.4,4),breaks=c(2,4,6),name='Hit genes')+guides(color=guide_colorbar(barwidth=grid::unit(18,'mm'),title.position='top'),size=guide_legend(nrow=1,title.position='top',override.aes=list(color='#4D4D4D')))+labs(x='Fold enrichment',y=NULL,title='a  ORA–KEGG overview')
matrix<-ggplot(edges,aes(ENTREZID,label,fill=logFC))+geom_tile(color='white',linewidth=.15)+scale_x_discrete(labels=symbols)+scale_fill_gradient2(low='#39799B',mid='white',high='#B46B33',midpoint=0,name='HS−LS log2FC')+guides(fill=guide_colorbar(barwidth=grid::unit(27,'mm'),title.position='top'))+labs(x='Contributing candidate genes',y=NULL,title='b  Gene contributions and abundance direction')+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_text(angle=60,hjust=1,vjust=1,size=6))
p<-(overview|matrix)+plot_layout(widths=c(1.3,1.2))+plot_annotation(title='Exploratory ORA–KEGG of legacy protein candidates',subtitle='Color in a shows unadjusted P-value. No ORA term passed BH FDR < 0.05.',theme=theme(plot.title=element_text(size=9,face='bold'),plot.subtitle=element_text(size=7)))
savefig(p,'Fig2C_ORA_KEGG_multidimensional',183,175)
# GSEA stays supplementary; FDR selection is retained but plot colors use nominal P.
k<-copy(e$GSEA[['KEGG']]);show<-head(k[order(padj,pval)],15)
show[,label:=vapply(Description,function(t)paste(strwrap(sub(' - Mus musculus.*','',t),width=45),collapse='\n'),character(1))];show[,label:=factor(label,levels=rev(unique(label)))];show[,FDR_supported:=padj<.05]
fwrite(show,file.path(b,'source_data/Supplementary_KEGG_GSEA_Pvalue_display.tsv'),sep='\t')
p<-ggplot(show,aes(NES,label,color=pval,size=size,shape=FDR_supported))+geom_vline(xintercept=0,linetype=3,linewidth=.3)+geom_point()+scale_color_gradient(low='#B46B33',high='#93B6C7',trans='log10',name='P-value')+scale_shape_manual(values=c('TRUE'=16,'FALSE'=1),name='BH FDR < 0.05')+labs(x='NES (positive: higher in HS)',y=NULL,title='Supplementary KEGG GSEA',subtitle='Top 15 ranked by adjusted P; color shows unadjusted P-value')
savefig(p,'Supplementary_KEGG_GSEA_Pvalue',183,150)
for(fam in c('GO_CC','GO_MF')){
 g<-copy(e$GSEA[[fam]]);g<-head(g[padj<.05][order(padj,pval)],15)
 g[,label:=vapply(Description,function(x)paste(strwrap(x,width=50),collapse='\n'),character(1))];g[,label:=factor(label,levels=rev(unique(label)))]
 p<-ggplot(g,aes(NES,label,color=pval,size=size))+geom_vline(xintercept=0,linetype=3,linewidth=.3)+geom_point()+scale_color_gradient(low='#B46B33',high='#93B6C7',trans='log10',name='P-value')+labs(x='NES (positive: higher in HS)',y=NULL,title=paste(fam,'GSEA'),subtitle='Top 15 FDR-significant sets; color shows unadjusted P-value')
 savefig(p,paste0('Supplementary_',fam,'_GSEA_Pvalue'),183,180)
}
writeLines(c('Fig2C is exploratory ORA KEGG, not GSEA.',
 'Panel a: top15 nominal-P ordered terms with hits; fold enrichment, Count, nominal P-value.',
 'Panel b: same terms and rows; all hit candidate genes shown; uniquely mapped gene-level abundance direction.',
 'No adjusted-P values were deleted or replaced; no ORA term met BH FDR<0.05.',
 'GSEA is supplementary; enrichment color is nominal P-value; FDR remains the selection/annotation criterion.'),file.path(b,'audit/ORA_display_user_update.txt'))
cat('USER REQUESTED ORA DISPLAY COMPLETED\n')


# ===========================================================================
# 4  Exploratory ORA, raw P < 0.05 and |log2FC| >= 1
# ===========================================================================
d<-fread(file.path(b,'statistics/main_all.tsv'));cand<-d[Pvalue<.05 & abs(log2FC)>=1]
stopifnot(nrow(cand)==200,sum(cand$log2FC>0)==149,sum(cand$log2FC<0)==51)
put_rel(cand,'statistics/nominalP05_FC2_candidates.tsv')
mp<-fread(file.path(b,'audit/gene_mapping_entrez.tsv'));mp<-unique(mp[entrez_status=='unique_entrez' & protein_group%in%d$protein_group,.(protein_group,gene_symbol,ENTREZID)])
mem<-fread(file.path(b,'reference/all_gene_set_memberships.tsv'))[family=='KEGG'];mem[,ENTREZID:=as.character(ENTREZID)]
terms<-unique(fread(file.path(b,'reference/all_term_names.tsv'))[,.(term,Description)])
bg<-sort(unique(mem$ENTREZID));stopifnot(all(bg%in%mp$ENTREZID))
sets<-lapply(split(mem$ENTREZID,mem$term),unique);sets<-sets[lengths(sets)>=2 & lengths(sets)<=500]
put_rel(data.table(ENTREZID=bg),'source_data/ORA_nominalP_FC2_KEGG_background.tsv')
tracking<-merge(cand,mp,by='protein_group',all.x=TRUE,allow.cartesian=TRUE);tracking[,KEGG_annotated:=!is.na(ENTREZID)&ENTREZID%in%bg]
put_rel(tracking,'audit/nominalP_FC2_mapping_tracking.tsv')
lists<-list(all=cand$protein_group,up=cand[log2FC>0,protein_group],down=cand[log2FC<0,protein_group]);out<-list();status<-list()
for(nm in names(lists)){
 ids<-sort(unique(mp[protein_group%in%lists[[nm]],ENTREZID]));fg<-intersect(ids,bg)
 put_rel(data.table(ENTREZID=fg),paste0('source_data/ORA_nominalP_FC2_',nm,'_KEGG_foreground.tsv'))
 stopifnot(length(fg)>0)
 z<-rbindlist(lapply(names(sets),function(id){a<-sets[[id]];hit<-intersect(fg,a);data.table(term=id,Count=length(hit),foreground_n=length(fg),background_n=length(bg),set_size=length(a),Pvalue=phyper(length(hit)-1,length(a),length(bg)-length(a),length(fg),lower.tail=FALSE),genes=paste(hit,collapse=';'))}))
 z[,`:=`(padj=p.adjust(Pvalue,'BH'),GeneRatio=paste0(Count,'/',foreground_n),BgRatio=paste0(set_size,'/',background_n),fold_enrichment=(Count/foreground_n)/(set_size/background_n))]
 z<-merge(z,terms,by='term',all.x=TRUE);setorder(z,Pvalue,term)
 stopifnot(all(z$padj>=0 & z$padj<=1),max(abs(z$padj-p.adjust(z$Pvalue,'BH')))<1e-12)
 put_rel(z,paste0('enrichment/ORA_nominalP_FC2_',nm,'_KEGG_all.tsv'));put_rel(z[padj<.05],paste0('enrichment/ORA_nominalP_FC2_',nm,'_KEGG_FDR05.tsv'));out[[nm]]<-z
 status[[nm]]<-data.table(list=nm,input_protein_groups=length(lists[[nm]]),unique_mapped_genes=length(ids),KEGG_foreground=length(fg),background=length(bg),tested_terms=nrow(z),nominal_P05=sum(z$Pvalue<.05),FDR05=sum(z$padj<.05),min_padj=min(z$padj))
}
put_rel(rbindlist(status),'enrichment/ORA_nominalP_FC2_summary.tsv')
gt<-fread(file.path(b,'statistics/gene_level_all.tsv'));gt[,ENTREZID:=as.character(ENTREZID)]
z<-head(out$all[Count>0],15)
counts<-rbindlist(lapply(seq_len(nrow(z)),function(i){ids<-strsplit(z$genes[i],';',fixed=TRUE)[[1]];x<-gt[ENTREZID%in%ids];stopifnot(nrow(x)==z$Count[i],all(is.finite(x$logFC)),all(x$logFC!=0));data.table(term=z$term[i],HS_higher=sum(x$logFC>0),LS_higher=sum(x$logFC<0))}))
z<-merge(z,counts,by='term',sort=FALSE);setorder(z,Pvalue,term);stopifnot(all(z$Count==z$HS_higher+z$LS_higher))
z[,label:=vapply(Description,function(t)paste(strwrap(sub(' - Mus musculus.*','',t),width=44),collapse='\n'),character(1))];z[,label:=factor(label,levels=rev(unique(label)))];z[,minus_log10_P:=-log10(pmax(Pvalue,1e-300))]
upperP<-ceiling(max(z$minus_log10_P)*2)/2;upperCount<-max(z$Count)+1;axis_scale<-upperCount/upperP;z[,line_x:=minus_log10_P*axis_scale]
bars<-melt(z,id.vars=c('term','label'),measure.vars=c('HS_higher','LS_higher'),variable.name='direction',value.name='count');bars[,direction:=factor(direction,levels=c('HS_higher','LS_higher'))]
put_rel(z,'source_data/Fig2C_ORA_KEGG_nominalP05_FC2.tsv');put_rel(bars,'source_data/Fig2C_ORA_KEGG_nominalP05_FC2_stacked_counts.tsv')
theme_set(theme_classic(base_size=8,base_family='Arial')+theme(plot.title=element_text(size=10,face='bold'),plot.subtitle=element_text(size=7),axis.text=element_text(size=7),legend.position='bottom',legend.title=element_blank(),legend.text=element_text(size=8),plot.margin=margin(8,12,8,8)))
subtitle<-sprintf('200 candidate protein groups: raw P < 0.05 and |log2FC| >= 1; %d mapped KEGG genes.\nTop 15 terms by unadjusted P-value; %d of %d tested terms pass BH FDR < 0.05.',status$all$KEGG_foreground,status$all$FDR05,status$all$tested_terms)
p<-ggplot()+geom_col(data=bars,aes(count,label,fill=direction),orientation='y',width=.65,position=position_stack(reverse=TRUE))+geom_path(data=z,aes(line_x,label,group=1),color='#222222',linewidth=.45)+geom_point(data=z,aes(line_x,label),shape=22,size=2.1,fill='#E8BD43',color='#222222',stroke=.4)+scale_fill_manual(values=c(HS_higher='#D94A38',LS_higher='#4B8DB9'),labels=c(HS_higher='Higher in HS',LS_higher='Higher in LS'))+scale_x_continuous(name='Mapped candidate genes (count)',limits=c(0,upperCount),breaks=pretty(c(0,upperCount),n=6),expand=expansion(mult=c(0,.015)),sec.axis=sec_axis(~./axis_scale,name='-log10(P-value)',breaks=pretty(c(0,upperP),n=5)))+labs(y=NULL,title='Exploratory ORA-KEGG: nominal-P / twofold candidates',subtitle=subtitle)
grDevices::cairo_pdf(file.path(b,'figures/Fig2C_ORA_KEGG_nominalP05_FC2.pdf'),width=183/25.4,height=165/25.4,family='Arial');print(p);dev.off()
writeLines(c('New exploratory foreground: current limma protein-group raw P < 0.05 and abs(log2FC) >= 1; 200 groups, 149 HS higher and 51 LS higher.',
 'Protein groups map through existing unique symbol-to-Entrez audit; ambiguous/unmapped groups excluded from gene-level ORA. Duplicate Entrez IDs counted once.',
 'Measured, tested and KEGG-mapped background retained; existing KEGG snapshot reused, eligible sets have 2-500 measured members. Hypergeometric upper-tail test; BH separately across all eligible KEGG terms for all/up/down lists, including zero-hit terms.',
 'Plot red/blue counts use signs of existing gene-level median-collapsed limma refits; they count unique contributing genes, not protein groups or proven activation. All hit counts must equal the sum of both directions.',
 'Plot ordered by raw P, first 15 nonzero-hit terms, no biological hand-selection. Complete results and BH values retained. Foreground selection remains exploratory even if a pathway passes BH.',
 'Original Fig2C_ORA_KEGG_reference_style.pdf and all earlier results preserved unchanged.'),file.path(b,'audit/nominalP_FC2_ORA_specification.txt'))
writeLines(c(capture.output(sessionInfo())),file.path(b,'logs/nominalP_FC2_ORA_sessionInfo.txt'));print(rbindlist(status));print(z[,.(Description,Count,Pvalue,padj)])


# ===========================================================================
# 5  Exploratory volcano for the same foreground
# ===========================================================================
d<-fread(file.path(b,'statistics/main_all.tsv'));d[,candidate:=!is.na(Pvalue)&Pvalue<.05&abs(log2FC)>=1]
n<-sum(d$candidate);up<-sum(d$candidate&d$log2FC>0);down<-sum(d$candidate&d$log2FC<0);stopifnot(n==200,up==149,down==51)
d[,plot_status:=ifelse(candidate,ifelse(log2FC>0,'HS higher','LS higher'),'Not selected')]
# Retain the five FDR-supported protein labels from the original main figure.
labelled<-d[significant_FDR==TRUE]
p<-ggplot(d,aes(log2FC,-log10(pmax(Pvalue,1e-300)),color=plot_status))+geom_point(size=.65,alpha=.65)+geom_hline(yintercept=-log10(.05),linetype=2,linewidth=.3)+geom_vline(xintercept=c(-1,1),linetype=3,linewidth=.3)+ggrepel::geom_text_repel(data=labelled,aes(label=gene_symbol),size=2.7,seed=20261031,max.overlaps=Inf,min.segment.length=0,show.legend=FALSE)+scale_color_manual(values=c('HS higher'='#D94A38','LS higher'='#4B8DB9','Not selected'='#9C9C9C'),breaks=c('LS higher','HS higher','Not selected'),labels=c(sprintf('LS higher (n = %d)',down),sprintf('HS higher (n = %d)',up),'Not selected'))+labs(x='log2 fold change (HS / LS)',y='-log10(P-value)',color=NULL,title='Exploratory protein abundance differences',subtitle=sprintf('%d tested | Candidates: %d | HS higher: %d | LS higher: %d\nRaw P < 0.05 and |log2FC| >= 1; candidates selected without BH correction.',nrow(d),n,up,down))+theme_classic(base_size=8,base_family='Arial')+theme(plot.title=element_text(size=9,face='bold'),plot.subtitle=element_text(size=7),legend.position='bottom',legend.text=element_text(size=7))
grDevices::cairo_pdf(file.path(b,'figures/Fig2B_protein_volcano_nominalP05_FC2.pdf'),width=150/25.4,height=115/25.4,family='Arial');print(p);dev.off()
fwrite(d,file.path(b,'source_data/Fig2B_volcano_nominalP05_FC2.tsv'),sep='\t',na='NA')
fwrite(data.table(tested=nrow(d),candidates=n,HS_higher=up,LS_higher=down,BH_FDR05=sum(d$significant_FDR)),file.path(b,'audit/Fig2B_nominalP_FC2_counts.tsv'),sep='\t')
cat('PASS: 200 candidates, 149 HS higher, 51 LS higher; original FDR figure untouched.\n')

cat('ALL PROTEOMICS FIGURES COMPLETED\n')
