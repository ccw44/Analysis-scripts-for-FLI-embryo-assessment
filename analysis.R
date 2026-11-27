# Protein-level reanalysis from the vendor quantitative matrix (HS minus LS).
#   Rscript analysis.R results/proteomics
# Reads  <base>/inputs/Result_N0011/...  and <base>/inputs/sample_metadata.tsv
#!/usr/bin/env Rscript
# Protein-level reanalysis. Plot contract: evidence grid; all six independent embryos;
# FDR-aware volcano as main evidence, exploratory heatmaps explicitly labelled.
args <- commandArgs(TRUE); stopifnot(length(args)==1)
base <- normalizePath(args[1], mustWork=TRUE)
suppressPackageStartupMessages({library(data.table); library(limma); library(ggplot2); library(patchwork)})
set.seed(20261031)
options(stringsAsFactors=FALSE)
for (d in c('audit','statistics','figures','source_data','logs')) dir.create(file.path(base,d),recursive=TRUE,showWarnings=FALSE)
input <- file.path(base,'inputs','Result_N0011')
read <- function(p) fread(file.path(input,p))
put <- function(x,p) fwrite(as.data.table(x),file.path(base,p),sep='\t',na='NA')
samples <- c('H1','H2','H3','L1','L2','L3'); hi <- 1:3; lo <- 4:6
raw <- read('01_Raw_Data/RawFile.txt'); filt <- read('03_Quality_Control/Filter_Data.txt')
norm <- read('04_Clean_Data/Normalized_Data.txt'); imp <- read('04_Clean_Data/Normalized_Impute_Data.txt')
old <- read('06_SubGroup_Analysis/H_L/Compare_SubGroup_H-L_Result.txt')
vol <- read('06_SubGroup_Analysis/H_L/01_DEG/Volcano_DEG_2_0.05.txt')
stopifnot(ncol(raw)==9,identical(names(filt)[-1],samples), !anyDuplicated(raw[[1]]), !anyDuplicated(norm[[1]]))
stopifnot(identical(filt[[1]],norm[[1]]),identical(norm[[1]],imp[[1]]))
ids <- norm[[1]]; rx <- as.matrix(raw[,4:9]); rownames(rx)<-raw[[1]]; colnames(rx)<-samples
fx <- as.matrix(filt[,-1]); rownames(fx)<-ids
nx <- as.matrix(norm[,-1]); rownames(nx)<-ids
ix <- as.matrix(imp[,-1]); rownames(ix)<-ids
rr <- rx[ids,,drop=FALSE]
stopifnot(identical(is.na(rr),fx==0),max(abs(rr[!is.na(rr)]-fx[!is.na(rr)]))<1e-7)
stopifnot(max(abs(ix[nx>0]-nx[nx>0]))<1e-7,all(nx>=0),all(ix>0))
scales <- sapply(1:6,function(j)median(nx[fx[,j]>0,j]/fx[fx[,j]>0,j]))
for(j in 1:6)stopifnot(max(abs(nx[,j]-fx[,j]*scales[j]))<1e-7)
nx[nx==0]<-NA_real_; lx<-log2(nx); rownames(lx)<-ids; colnames(lx)<-samples
sx<-as.data.table(read('02_Stat_Info/Stat_Summary.txt'))
put(sx,'audit/instrument_QC_as_supplied.tsv')
meta<-fread(file.path(base,'inputs/sample_metadata.tsv'));stopifnot(identical(meta$sample_id,samples))
qc<-data.table(sample_id=samples,group=meta$group,raw_detected=colSums(!is.na(rx)),filtered_detected=colSums(!is.na(nx)),filtered_missing=colSums(is.na(nx)),filtered_missing_pct=100*colMeans(is.na(nx)),raw_median=apply(rx,2,median,na.rm=TRUE),normalized_median=apply(nx,2,median,na.rm=TRUE),normalization_factor=scales)
put(qc,'audit/sample_QC.tsv')
mask<-data.table(protein_group=ids,as.data.table(!is.na(nx)));put(mask,'source_data/observed_mask.tsv')
put(data.table(protein_group=ids,as.data.table(lx)),'source_data/observed_log2_normalized.tsv')
removed<-raw[!raw[[1]]%in%ids,1:3];put(removed,'audit/legacy_filtered_out_unverified_reason.tsv')
annotations<-data.table(protein_group=raw[[1]],gene_symbol=raw[[2]],n_stripped_sequences=raw[[3]])
put(annotations,'audit/protein_annotation.tsv')
nh<-rowSums(!is.na(lx[,hi]));nl<-rowSums(!is.na(lx[,lo]));keep<-nh>=2 & nl>=2; complete<-complete.cases(lx)
put(data.table(protein_group=ids,n_HS=nh,n_LS=nl,eligible_main=keep,complete_six=complete),'audit/eligibility.tsv')
det<-data.table(protein_group=ids,n_HS=nh,n_LS=nl,HS_mean_log2=rowMeans(lx[,hi],na.rm=TRUE),LS_mean_log2=rowMeans(lx[,lo],na.rm=TRUE))
det<-det[(n_HS>=2 & n_LS<2)|(n_LS>=2 & n_HS<2)]
det[,note:='Detection-pattern candidate only; no low-value imputation or conventional significance assigned']
put(det,'statistics/detection_pattern_candidates.tsv')
design<-model.matrix(~0+factor(meta$group,levels=c('LS','HS')));colnames(design)<-c('LS','HS')
contrast<-makeContrasts(HS-LS,levels=design)
pre<-lmFit(lx[keep,,drop=FALSE],design)
finite<-is.finite(pre$Amean)&is.finite(pre$sigma)&pre$sigma>0
diag<-cor.test(pre$Amean[finite],log2(pre$sigma[finite]^2),method='spearman',exact=FALSE)
trend<-abs(unname(diag$estimate))>=0.2 && diag$p.value<0.05
put(data.table(spearman_rho=unname(diag$estimate),p=diag$p.value,trend_selected=trend,rule='abs(rho)>=0.2 and P<0.05; selected from residual variance, not differential P values'),'audit/variance_trend_decision.tsv')
writeLines(c('Independent protein-group analysis, HS minus LS; no pairing.',
 'Use supplied median-normalized per-sample linear intensities after confirming fixed scaling and observation preservation.',
 'Confirmed zero means missing; log2; no imputation for primary inference.',
 'Primary filter: at least two observed values in each group; BH across all finite primary tests.',
 paste('limma eBayes robust=TRUE; trend=',trend,'; chosen by pre-test variance rule.'),
 'FDR<0.05 primary; abs(log2FC)>=1 subset separately. No data-driven sample exclusion.',
 'Unknown: original identification filters and exact normalization implementation, sample provenance/stage/batch, instrument and search settings.',
 'This is a post-review reanalysis specification, not a prospective preregistration.'),file.path(base,'audit/locked_analysis_specification.txt'))
fit_matrix<-function(mat,label,trend_use=trend) {
 f<-lmFit(mat,design);f<-contrasts.fit(f,contrast);f<-eBayes(f,trend=trend_use,robust=TRUE)
 z<-topTable(f,coef=1,number=Inf,sort.by='none',confint=TRUE,adjust.method='BH')
 out<-data.table(protein_group=rownames(z),log2FC=z$logFC,average_log2=z$AveExpr,CI95_low=z$CI.L,CI95_high=z$CI.R,moderated_t=z$t,Pvalue=z$P.Value,padj=z$adj.P.Val,B=z$B)
 out[,moderated_SE:=(CI95_high-CI95_low)/(2*qt(.975,df=f$df.total))]
 out[,status:=ifelse(is.finite(Pvalue),'tested','not_estimable')]
 out[,significant_FDR:=is.finite(padj)&padj<.05]
 out[,significant_FDR_FC2:=significant_FDR & abs(log2FC)>=1]
 out[,direction:=ifelse(log2FC>0,'Higher in HS',ifelse(log2FC<0,'Higher in LS','Equal'))]
 out<-merge(out,annotations,by='protein_group',all.x=TRUE,sort=FALSE)
 if(label!='gene_level'){
  k<-match(out$protein_group,ids);out[,':='(n_HS_observed=nh[k],n_LS_observed=nl[k])]
 }
 put(out,paste0('statistics/',label,'_all.tsv'))
 put(out[significant_FDR==TRUE],paste0('statistics/',label,'_FDR05.tsv'))
 put(out[significant_FDR_FC2==TRUE],paste0('statistics/',label,'_FDR05_FC2.tsv'))
 list(table=out,fit=f)
}
main<-fit_matrix(lx[keep,,drop=FALSE],'main')
cc<-fit_matrix(lx[complete,,drop=FALSE],'sensitivity_complete')
knn<-fit_matrix(log2(ix[keep,,drop=FALSE]),'sensitivity_supplied_KNN')
rawlog<-log2(rr); centers<-apply(rawlog,2,median,na.rm=TRUE);recenter<-sweep(rawlog,2,centers,'-')+median(centers)
med<-fit_matrix(recenter[keep,,drop=FALSE],'sensitivity_observed_median_centering')
cvrows<-lapply(seq_along(ids),function(i){data.table(protein_group=ids[i],CV_HS=sd(nx[i,hi],na.rm=TRUE)/mean(nx[i,hi],na.rm=TRUE),CV_LS=sd(nx[i,lo],na.rm=TRUE)/mean(nx[i,lo],na.rm=TRUE))});put(rbindlist(cvrows),'audit/group_CV.tsv')
old[,BH_recomputed:=p.adjust(Pvalue,'BH')]
put(old,'audit/legacy_statistics_BH_recomputed.tsv')
legacy85<-vol[threshold!='Not'];new<-main$table
legacy85[,legacy_label:=threshold]
cmp<-merge(legacy85[,.(protein_group=ProteinGroups,legacy_gene=Genes,legacy_label,legacy_log2FC=FC_AB,legacy_P=Pvalue,legacy_invalid_Padjust=Padjust)],new,by='protein_group',all.x=TRUE)
put(cmp,'statistics/legacy85_tracking.tsv')
sens<-rbindlist(lapply(list(Main=main,Complete=cc,KNN=knn,MedianRecentering=med),function(x){z<-x$table;data.table(n=nrow(z),tested=sum(z$status=='tested'),FDR05=sum(z$significant_FDR),FDR05_FC2=sum(z$significant_FDR_FC2),HS_higher_FDR=sum(z$significant_FDR&z$log2FC>0),LS_higher_FDR=sum(z$significant_FDR&z$log2FC<0),min_padj=min(z$padj,na.rm=TRUE))}),idcol='analysis')
put(sens,'statistics/analysis_summary.tsv')
compare<-Reduce(function(a,b)merge(a,b,by='protein_group',all.x=TRUE),list(new[,.(protein_group,main_log2FC=log2FC,main_padj=padj)],cc$table[,.(protein_group,complete_log2FC=log2FC,complete_padj=padj)],knn$table[,.(protein_group,KNN_log2FC=log2FC,KNN_padj=padj)],med$table[,.(protein_group,recenter_log2FC=log2FC,recenter_padj=padj)]))
put(compare,'statistics/sensitivity_protein_comparison.tsv')
mapping<-merge(data.table(protein_group=ids),annotations,by='protein_group',all.x=TRUE,sort=FALSE)
mapping[,mapping_status:=ifelse(is.na(gene_symbol)|gene_symbol=='','missing_gene',ifelse(grepl('[;,|]',gene_symbol),'ambiguous_multi_gene','single_symbol'))]
put(mapping,'audit/gene_mapping_initial.tsv')
saveRDS(list(log2=lx,keep=keep,design=design,contrast=contrast,main=new,trend=trend,annotation=annotations,mapping=mapping,legacy85=legacy85),file.path(base,'statistics/analysis_state.rds'))
stopifnot(all(new$padj>=0 & new$padj<=1,na.rm=TRUE))
stopifnot(max(abs(new$padj-p.adjust(new$Pvalue,'BH')),na.rm=TRUE)<1e-10)
expected<-rowMeans(lx[keep,hi,drop=FALSE],na.rm=TRUE)-rowMeans(lx[keep,lo,drop=FALSE],na.rm=TRUE)
stopifnot(max(abs(new$log2FC-expected[match(new$protein_group,names(expected))]),na.rm=TRUE)<1e-8)
put(data.table(check=c('raw_NA_matches_legacy_zero','normalization_is_samplewise_fixed_scaling','KNN_observed_values_unchanged','BH_range_and_recomputation','contrast_matches_HS_minus_LS'),pass=TRUE),'audit/statistical_assertions.tsv')
# R-only publication exports. Data points have no aggregate uncertainty; protein CIs are in source tables.
theme_set(theme_classic(base_size=8,base_family='Arial')+theme(plot.title=element_text(size=9,face='bold'),plot.subtitle=element_text(size=7),legend.position='bottom',legend.text=element_text(size=7),strip.text=element_text(size=8)))
save_plot<-function(p,name,w=183,h=125){
 grDevices::cairo_pdf(file.path(base,'figures',paste0(name,'.pdf')),width=w/25.4,height=h/25.4,family='Arial');print(p);dev.off()
}
colors<-c(HS='#B46B33',LS='#39799B')
missing<-ggplot(qc,aes(sample_id,filtered_missing_pct,fill=group))+geom_col(width=.65)+scale_fill_manual(values=colors)+labs(x=NULL,y='Missing protein groups (%)',title='a  Missing observations',fill=NULL)
long<-function(mat,stage){d<-as.data.table(mat,keep.rownames='protein_group');setnames(d,names(d)[1],'protein_group');d<-melt(d,id.vars='protein_group',variable.name='sample_id',value.name='log2_intensity',na.rm=TRUE);d[,stage:=stage];d[,group:=ifelse(substr(sample_id,1,1)=='H','HS','LS')];d}
dis<-rbind(long(log2(rr),'Before normalization'),long(lx,'After normalization'));put(dis,'source_data/QC_intensity_long.tsv')
distplot<-ggplot(dis,aes(sample_id,log2_intensity,fill=group))+geom_boxplot(outlier.shape=NA,width=.65)+facet_wrap(~stage,ncol=1)+scale_fill_manual(values=colors)+labs(x=NULL,y='Observed log2 intensity',title='b  Intensity distributions',fill=NULL)
pc<-prcomp(t(lx[complete,,drop=FALSE]),center=TRUE,scale.=FALSE);pv<-100*pc$sdev^2/sum(pc$sdev^2)
pcs<-data.table(sample_id=samples,group=meta$group,pc$x[,1:3]);put(pcs,'source_data/PCA_complete_observations.tsv')
put(data.table(PC=names(pv),variance_pct=pv),'source_data/PCA_variance.tsv')
pca<-ggplot(pcs,aes(PC1,PC2,color=group,label=sample_id))+geom_point(size=2)+ggrepel::geom_text_repel(size=2.7,seed=20261031,max.overlaps=Inf)+scale_color_manual(values=colors)+labs(x=sprintf('PC1 (%.1f%%)',pv[1]),y=sprintf('PC2 (%.1f%%)',pv[2]),title='c  PCA without imputation',color=NULL)
corr<-cor(lx[complete,,drop=FALSE]);put(data.table(sample_id=samples,as.data.table(corr)),'source_data/correlation_complete.tsv')
cd<-as.data.table(as.table(corr));setnames(cd,c('sample1','sample2','Pearson_r'))
heat<-ggplot(cd,aes(sample1,sample2,fill=Pearson_r))+geom_tile()+geom_text(aes(label=sprintf('%.2f',Pearson_r)),size=2.5)+scale_fill_gradient(low='#F4F3EE',high='#39799B')+coord_equal()+labs(x=NULL,y=NULL,title='d  Correlation without imputation',fill='Pearson r')
save_plot((missing|pca)/(distplot|heat),'Supplementary_protein_QC',183,180)
vdata<-copy(new);vdata[,display:='Not FDR significant'];vdata[significant_FDR==TRUE,display:=ifelse(log2FC>0,'FDR significant, HS higher','FDR significant, LS higher')]
put(vdata,'source_data/Fig2B_volcano.tsv')
vplot<-ggplot(vdata,aes(log2FC,-log10(pmax(padj,1e-300)),color=display))+geom_point(size=.65,alpha=.6)+geom_hline(yintercept=-log10(.05),linetype=2,linewidth=.3)+geom_vline(xintercept=c(-1,1),linetype=3,linewidth=.3)+scale_color_manual(values=c('Not FDR significant'='#9C9C9C','FDR significant, HS higher'='#B46B33','FDR significant, LS higher'='#39799B'))+labs(x='log2 fold change (HS / LS)',y='-log10 BH-adjusted P',color=NULL,title='Protein abundance differences',subtitle=sprintf('%d tested protein groups; FDR < 0.05: %d (%d HS higher, %d LS higher)',nrow(new),sum(new$significant_FDR),sum(new$significant_FDR&new$log2FC>0),sum(new$significant_FDR&new$log2FC<0)))
save_plot(vplot,'Fig2B_protein_volcano',140,105)
vp<-ggplot(vdata,aes(log2FC,-log10(pmax(Pvalue,1e-300)),color=display))+geom_point(size=.6,alpha=.5)+scale_color_manual(values=c('Not FDR significant'='#9C9C9C','FDR significant, HS higher'='#B46B33','FDR significant, LS higher'='#39799B'))+labs(x='log2 fold change (HS / LS)',y='-log10 unadjusted P',color=NULL,title='Protein abundance differences',subtitle='Colors denote BH FDR < 0.05; vertical axis shows unadjusted P')
save_plot(vp,'Supplementary_volcano_unadjusted_axis',140,105)
vd<-data.table(average_log2=pre$Amean,log2_residual_variance=log2(pre$sigma^2))
varianceplot<-ggplot(vd,aes(average_log2,log2_residual_variance))+geom_point(size=.5,alpha=.25)+geom_smooth(method='loess',se=FALSE,color='#39799B',linewidth=.5)+labs(x='Mean observed log2 intensity',y='Log2 residual variance',title='Variance trend diagnostic',subtitle=sprintf('Spearman rho = %.3f; trend = %s; robust = TRUE',diag$estimate,trend))
save_plot(varianceplot,'Supplementary_variance_trend',120,90)
# All 85 legacy candidates are displayed in pages of at most 30; missing stays blank.
oldids<-legacy85$ProteinGroups;stopifnot(all(oldids%in%ids));oldids<-oldids[order(legacy85$Genes,legacy85$ProteinGroups)]
for (page in seq_len(ceiling(length(oldids)/30))){
 chosen<-oldids[seq((page-1)*30+1,min(page*30,length(oldids)))];z<-lx[chosen,,drop=FALSE];z<-t(apply(z,1,function(v)(v-mean(v,na.rm=TRUE))/sd(v,na.rm=TRUE)));colnames(z)<-samples;rownames(z)<-chosen
 hd<-as.data.table(z,keep.rownames='protein_group');hd<-melt(hd,id.vars='protein_group',variable.name='sample_id',value.name='row_z')
 labels<-annotations$gene_symbol[match(chosen,annotations$protein_group)];labels[is.na(labels)]<-chosen[is.na(labels)];labsmap<-setNames(paste0(labels,' [',chosen,']'),chosen)
 hd[,protein_group:=factor(protein_group,levels=rev(chosen))];put(hd,paste0('source_data/legacy_heatmap_page',page,'.tsv'))
 hp<-ggplot(hd,aes(sample_id,protein_group,fill=row_z))+geom_tile()+scale_y_discrete(labels=labsmap)+scale_fill_gradient2(low='#39799B',mid='white',high='#B46B33',na.value='#DCDCDC',name='Row z')+labs(x=NULL,y=NULL,title=sprintf('Legacy exploratory candidates, page %d',page),subtitle='Unadjusted-P selected in original analysis; gray denotes missing')+theme(axis.text.y=element_text(size=6),axis.text.x=element_text(size=8))
 save_plot(hp,paste0('Supplementary_legacy_heatmap_',page),150,160)
}
writeLines(capture.output(sessionInfo()),file.path(base,'logs/sessionInfo_statistics.txt'))
cat('MAIN ANALYSIS COMPLETED\n');print(sens)
