# GO / KEGG over-representation and GSEA on the gene-collapsed matrix.
#   Rscript enrichment.R results/proteomics
# Needs the analysis_state.rds written by analysis.R
#!/usr/bin/env Rscript
# R-only figures; full eligible tested gene-set results, including P=1 terms.
args<-commandArgs(TRUE);stopifnot(length(args)==1);base<-normalizePath(args[1])
suppressPackageStartupMessages({library(data.table);library(limma);library(AnnotationDbi);library(org.Mm.eg.db);library(GO.db);library(fgsea);library(ggplot2)})
set.seed(20261031);options(stringsAsFactors=FALSE)
for(d in c('enrichment','reference'))dir.create(file.path(base,d),showWarnings=FALSE)
put<-function(x,p)fwrite(as.data.table(x),file.path(base,p),sep='\t',na='NA')
s<-readRDS(file.path(base,'statistics/analysis_state.rds'))
map<-s$mapping[mapping_status=='single_symbol']
ann<-as.data.table(AnnotationDbi::select(org.Mm.eg.db,keys=unique(map$gene_symbol),keytype='SYMBOL',columns='ENTREZID'))
ann<-unique(ann[!is.na(ENTREZID)]);ann[,n_entrez:=uniqueN(ENTREZID),by=SYMBOL]
map<-merge(map,ann,by.x='gene_symbol',by.y='SYMBOL',all.x=TRUE,allow.cartesian=TRUE)
map[,entrez_status:=ifelse(is.na(ENTREZID),'unmapped',ifelse(n_entrez!=1,'ambiguous_symbol_to_multiple_entrez','unique_entrez'))]
put(map,'audit/gene_mapping_entrez.tsv')
valid<-map[entrez_status=='unique_entrez' & protein_group%in%s$main$protein_group]
universe<-sort(unique(valid$ENTREZID));stopifnot(length(universe)>100)
gene_mat<-t(vapply(universe,function(g){ii<-valid[ENTREZID==g,protein_group];apply(s$log2[ii,,drop=FALSE],2,function(v)median(v,na.rm=TRUE))},numeric(6)))
rownames(gene_mat)<-universe;colnames(gene_mat)<-colnames(s$log2)
put(data.table(ENTREZID=universe,as.data.table(gene_mat)),'source_data/gene_collapsed_log2.tsv')
gf<-lmFit(gene_mat,s$design);ok<-is.finite(gf$sigma)&gf$sigma>0
diag<-cor.test(gf$Amean[ok],log2(gf$sigma[ok]^2),method='spearman',exact=FALSE)
gtrend<-abs(diag$estimate)>=.2 & diag$p.value<.05
put(data.table(gene_trend=gtrend,rho=unname(diag$estimate),P=diag$p.value),'audit/gene_variance_trend.tsv')
gf<-eBayes(contrasts.fit(gf,s$contrast),trend=gtrend,robust=TRUE)
gt<-as.data.table(topTable(gf,number=Inf,sort.by='none',confint=TRUE),keep.rownames='ENTREZID');setnames(gt,'adj.P.Val','padj')
sym<-valid[,.(symbols=paste(sort(unique(gene_symbol)),collapse=';'),protein_groups=paste(sort(unique(protein_group)),collapse=';')),by=ENTREZID]
gt<-merge(gt,sym,by='ENTREZID',sort=FALSE);put(gt,'statistics/gene_level_all.tsv')
ranks<-setNames(gt$t,gt$ENTREZID);ranks<-sort(ranks[is.finite(ranks)],decreasing=TRUE)
put(data.table(ENTREZID=names(ranks),moderated_t=as.numeric(ranks)),'source_data/GSEA_ranked_genes.tsv')
go<-as.data.table(AnnotationDbi::select(org.Mm.eg.db,keys=universe,keytype='ENTREZID',columns=c('GOALL','ONTOLOGYALL')))
go<-unique(go[!is.na(GOALL),.(ENTREZID,term=GOALL,family=paste0('GO_',ONTOLOGYALL))])
goids<-unique(go$term);got<-as.data.table(AnnotationDbi::select(GO.db,keys=goids,keytype='GOID',columns='TERM'));got<-unique(got[,.(term=GOID,Description=TERM)])
put(go,'reference/GO_memberships_tested_universe.tsv');put(got,'reference/GO_term_names.tsv')
writeLines(c(paste('org.Mm.eg.db',packageVersion('org.Mm.eg.db')),paste('GO.db',packageVersion('GO.db')),paste('retrieval_time_UTC',format(Sys.time(),tz='UTC',usetz=TRUE)),capture.output(AnnotationDbi::metadata(org.Mm.eg.db))),file.path(base,'reference/database_versions.txt'))
fetch<-function(url,path){
 if(!file.exists(path)){
  tmp<-paste0(path,'.download');status<-download.file(url,tmp,method='libcurl',quiet=TRUE)
  stopifnot(status==0,file.info(tmp)$size>100);stopifnot(file.rename(tmp,path))
 }
 stopifnot(file.info(path)$size>100)
}
kegg_ok<-TRUE
tryCatch({
 fetch('https://rest.kegg.jp/link/pathway/mmu',file.path(base,'reference/KEGG_mmu_gene_pathway.tsv'))
 fetch('https://rest.kegg.jp/list/pathway/mmu',file.path(base,'reference/KEGG_mmu_names.tsv'))
},error=function(e){kegg_ok<<-FALSE;writeLines(conditionMessage(e),file.path(base,'logs/KEGG_download_error.txt'))})
membership<-go;terms<-got
if(kegg_ok){
 k<-fread(file.path(base,'reference/KEGG_mmu_gene_pathway.tsv'),header=FALSE);setnames(k,c('ENTREZID','term'));k[,ENTREZID:=sub('^mmu:','',ENTREZID)];k[,term:=sub('^path:','',term)];k[,family:='KEGG'];k<-k[ENTREZID%in%universe]
 kt<-fread(file.path(base,'reference/KEGG_mmu_names.tsv'),header=FALSE,sep='\t');setnames(kt,c('term','Description'));kt[,term:=sub('^path:','',term)]
 terms<-rbind(got,kt);membership<-rbind(go,k)
}
membership<-unique(membership);put(membership,'reference/all_gene_set_memberships.tsv');put(terms,'reference/all_term_names.tsv')
writeLines(c('GO/KEGG mapping is gene level with a measured/tested universe; not all mouse genes.',
 'ORA universe is measured/tested genes with annotation in the respective family; foreground intersects that universe. BH separately for GO BP, CC, MF and KEGG and each list/direction; terms with 2-500 measured members.',
 'GSEA BH separately within each family; 15-500 measured genes, signed moderated t, fgseaMultilevel, sampleSize=101, nPermSimple=1000, eps=0, seed=20261031, nproc=1.',
 'Duplicate protein groups per uniquely mapped gene collapsed by within-sample median log2 intensity, followed by a new limma fit.',
 'Raw-P legacy candidate ORA is exploratory; primary ORA uses FDR-significant proteins mapped to genes.',
 paste('KEGG retrieval successful:',kegg_ok)),file.path(base,'audit/enrichment_specification.txt'))
old85<-s$legacy85$ProteinGroups
lists<-list(primary_up=s$main[significant_FDR & log2FC>0,protein_group],primary_down=s$main[significant_FDR & log2FC<0,protein_group],legacy85_all=old85,legacy85_up=s$legacy85[FC_AB>0,ProteinGroups],legacy85_down=s$legacy85[FC_AB<0,ProteinGroups])
empty_ora<-data.table(family=character(),term=character(),Description=character(),Count=integer(),foreground_n=integer(),set_size=integer(),background_n=integer(),GeneRatio=character(),BgRatio=character(),fold_enrichment=numeric(),Pvalue=numeric(),padj=numeric(),genes=character())
ora_status<-list();ora_results<-list()
for(nm in names(lists)){
 genes_all<-unique(valid[protein_group%in%lists[[nm]],ENTREZID]);mapped<-data.table(list=nm,ENTREZID=genes_all);put(mapped,paste0('source_data/ORA_',nm,'_foreground.tsv'))
 for(fam in sort(unique(membership$family))){
  m<-membership[family==fam];background<-sort(unique(m$ENTREZID));genes<-intersect(genes_all,background)
  put(data.table(ENTREZID=background),paste0('source_data/ORA_',fam,'_background.tsv'))
  put(data.table(ENTREZID=genes),paste0('source_data/ORA_',nm,'_',fam,'_foreground.tsv'))
  sets<-split(m$ENTREZID,m$term);sets<-lapply(sets,unique);sets<-sets[lengths(sets)>=2 & lengths(sets)<=500]
  if(!length(genes)||!length(sets)){
   z<-copy(empty_ora);status<-'NOT_APPLICABLE_empty_foreground_or_eligible_sets'
  }else{
   M<-length(background);n<-length(genes)
   z<-rbindlist(lapply(names(sets),function(id){a<-sets[[id]];hit<-intersect(genes,a);b<-length(a);cnt<-length(hit);data.table(family=fam,term=id,Count=cnt,foreground_n=n,set_size=b,background_n=M,GeneRatio=paste0(cnt,'/',n),BgRatio=paste0(b,'/',M),fold_enrichment=(cnt/n)/(b/M),Pvalue=phyper(cnt-1,b,M-b,n,lower.tail=FALSE),genes=paste(hit,collapse=';'))}))
   z[,padj:=p.adjust(Pvalue,'BH')];z<-merge(z,terms,by='term',all.x=TRUE,sort=FALSE);setorder(z,padj,Pvalue);status<-'COMPLETED'
  }
  put(z,paste0('enrichment/ORA_',nm,'_',fam,'_all.tsv'));ora_results[[paste(nm,fam)]]<-z
  ora_status[[length(ora_status)+1]]<-data.table(analysis='ORA',list=nm,family=fam,status=status,input_protein_groups=length(lists[[nm]]),mapped_foreground_genes=length(genes),tested_terms=nrow(z),FDR05=sum(z$padj<.05))
 }
}
gsresults<-list();gsstatus<-list()
for(fam in sort(unique(membership$family))){
 m<-membership[family==fam];sets<-lapply(split(m$ENTREZID,m$term),unique);sets<-sets[lengths(sets)>=15 & lengths(sets)<=500]
 set.seed(20261031);cat('GSEA',fam,'sets',length(sets),'\n');flush.console()
 g<-as.data.table(fgsea::fgseaMultilevel(pathways=sets,stats=ranks,minSize=15,maxSize=500,sampleSize=101,nPermSimple=1000,eps=0,nproc=1))
 g[,family:=fam];g[,leading_edge_entrez:=vapply(leadingEdge,paste,collapse=';',FUN.VALUE=character(1))]
 le<-rbindlist(lapply(seq_len(nrow(g)),function(i){if(!length(g$leadingEdge[[i]]))return(NULL);data.table(family=fam,term=g$pathway[i],ENTREZID=g$leadingEdge[[i]])}))
 le<-merge(le,gt[,.(ENTREZID,symbols,logFC,t,P.Value,padj)],by='ENTREZID',all.x=TRUE,sort=FALSE)
 put(le,paste0('enrichment/GSEA_',fam,'_leading_edge.tsv'))
 g[,leadingEdge:=NULL];setnames(g,'pathway','term');g<-merge(g,terms,by='term',all.x=TRUE,sort=FALSE);setorder(g,padj,pval)
 put(g,paste0('enrichment/GSEA_',fam,'_all.tsv'));put(g[padj<.05],paste0('enrichment/GSEA_',fam,'_FDR05.tsv'));gsresults[[fam]]<-g
 gsstatus[[length(gsstatus)+1]]<-data.table(analysis='GSEA',family=fam,status='COMPLETED',tested_terms=nrow(g),FDR05=sum(g$padj<.05,na.rm=TRUE),min_padj=min(g$padj,na.rm=TRUE))
}
put(rbindlist(ora_status),'enrichment/ORA_status.tsv');put(rbindlist(gsstatus),'enrichment/GSEA_status.tsv')
saveRDS(list(universe=universe,membership=membership,terms=terms,ORA=ora_results,GSEA=gsresults,gene_statistics=gt),file.path(base,'enrichment/enrichment_state.rds'))
# Existing historical tables are archived separately, not silently replaced.
for(f in c('GO_ONTOLOGY','GO_GSEA_ONTOLOGY','KEGG_Pathway','KEGG_GSEA_Pathway')){
 hist<-fread(file.path(base,'inputs/Result_N0011/06_SubGroup_Analysis/H_L/03_Annotation',paste0(f,'.txt')))
 put(hist,paste0('audit/legacy_',f,'.tsv'))
}
theme_set(theme_classic(base_size=8,base_family='Arial')+theme(plot.title=element_text(size=9,face='bold'),plot.subtitle=element_text(size=7)))
save_plot<-function(p,name,w=183,h=130){
 grDevices::cairo_pdf(file.path(base,'figures',paste0(name,'.pdf')),width=w/25.4,height=h/25.4,family='Arial');print(p);dev.off()
}
if(kegg_ok){
 z<-ora_results[['legacy85_all KEGG']];z<-z[Count>0];z<-head(z[order(padj,Pvalue)],15)
 if(nrow(z)){
  z[,label:=vapply(Description,function(s)paste(strwrap(sub(' - Mus musculus.*','',s),width=42),collapse='\n'),character(1))];z[,label:=factor(label,levels=rev(unique(label)))];z[,ratio:=Count/foreground_n]
  put(z,'source_data/Fig2C_exploratory_KEGG.tsv')
  p<-ggplot(z,aes(ratio,label,size=Count,color=padj))+geom_point()+scale_color_gradient(low='#B46B33',high='#39799B',name='BH adjusted P')+scale_size_continuous(range=c(1.8,4.5))+labs(x='Gene ratio in mapped legacy candidates',y=NULL,title='Exploratory KEGG analysis of legacy candidates',subtitle=sprintf('Unadjusted-P selected proteins; tested-universe background; %d displayed terms pass FDR < 0.05',sum(z$padj<.05)))
  save_plot(p,'Fig2C_exploratory_KEGG',183,145)
 }
 g<-gsresults[['KEGG']];show<-head(g[order(padj,pval)],15);show[,label:=vapply(Description,function(s)paste(strwrap(sub(' - Mus musculus.*','',s),width=42),collapse='\n'),character(1))];show[,label:=factor(label,levels=rev(unique(label)))];put(show,'source_data/KEGG_GSEA_display.tsv')
 p<-ggplot(show,aes(NES,label,color=padj,size=size))+geom_vline(xintercept=0,linetype=3,linewidth=.3)+geom_point()+scale_color_gradient(low='#B46B33',high='#39799B',name='BH adjusted P')+labs(x='NES (positive: higher in HS)',y=NULL,title='KEGG GSEA using all ranked measured genes',subtitle='Terms ranked by adjusted P; full table includes all tested terms')
 save_plot(p,'Supplementary_KEGG_GSEA',183,145)
}
writeLines(capture.output(sessionInfo()),file.path(base,'logs/sessionInfo_enrichment.txt'))
cat('ENRICHMENT COMPLETED\n');print(rbindlist(gsstatus))
