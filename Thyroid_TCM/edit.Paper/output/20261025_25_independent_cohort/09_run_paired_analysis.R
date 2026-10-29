# Internal reproducible analysis; run from repository root.
suppressPackageStartupMessages(library(limma))
root <- 'edit.Paper/output/20261025_25_independent_cohort'
matfile <- file.path(root, '02_GSE33630_series_matrix.txt.gz')
gplfile <- 'edit.Paper/output/20261024_25_transcriptome_validation/20261024_GPL570.annot.gz'
get_row <- function(lines, prefix) {
  i <- which(startsWith(lines, prefix))[1]
  if (is.na(i)) stop('Missing metadata row: ', prefix)
  read.delim(text=lines[i], header=FALSE, sep='\t', quote='"', check.names=FALSE,
             stringsAsFactors=FALSE, fill=TRUE, comment.char='')
}
ml <- readLines(gzfile(matfile), warn=FALSE)
ib <- which(ml == '!series_matrix_table_begin')
ie <- which(ml == '!series_matrix_table_end')
if (length(ib)!=1 || length(ie)!=1) stop('Matrix table boundary missing/multiple')
titles <- as.character(get_row(ml,'!Sample_title')[1,-1])
gsm <- as.character(get_row(ml,'!Sample_geo_accession')[1,-1])
sources <- as.character(get_row(ml,'!Sample_source_name_ch1')[1,-1])
if (length(titles)!=length(gsm)) stop('Metadata column count mismatch')
# Preserve all submitted sample characteristics for manual traceability.
charrows <- which(startsWith(ml,'!Sample_characteristics_ch1'))
charlist <- lapply(charrows,function(i) {z<-read.delim(text=ml[i],header=FALSE,sep='\t',quote='\"',check.names=FALSE,stringsAsFactors=FALSE,fill=TRUE,comment.char=''); as.character(z[1,-1])})
char_keys <- vapply(charlist,function(x) sub(':.*$','',x[1]),'')
characteristics <- as.data.frame(do.call(cbind,charlist),stringsAsFactors=FALSE,check.names=FALSE)
if(ncol(characteristics)) names(characteristics) <- make.unique(char_keys)
rownames(characteristics) <- gsm
# Source metadata identifies matched PTC and non-tumoral thyroid tissues.
patient <- ifelse(grepl('^[0-9]+[TN]$',titles),sub('[TN]$','',titles),NA_character_)
status <- ifelse(grepl('^[0-9]+T$',titles),'PTC',ifelse(grepl('^[0-9]+N$',titles),'Normal thyroid','Exclude'))
stopifnot(sum(status=='PTC')==49,sum(status=='Normal thyroid')==45,sum(status=='Exclude')==11)
manifest <- data.frame(GSM=gsm, title=titles, source_name=sources, patient_id=patient, tissue_group=status,
                       include=status!='Exclude' & !is.na(patient), stringsAsFactors=FALSE)
manifest <- cbind(manifest, characteristics[gsm,,drop=FALSE])
write.csv(manifest, file.path(root,'GSE33630_02_sample_manifest_all.csv'), row.names=FALSE)
# Retain only primary PTC and normal thyroid; collapse all tissue blocks within patient x group.
use <- which(manifest$include)
if (!length(use)) stop('No eligible samples parsed')
raw <- read.delim(text=paste(ml[(ib+1):(ie-1)], collapse='\n'), sep='\t', header=TRUE,
                  check.names=FALSE, quote='"', comment.char='', stringsAsFactors=FALSE)
if (ncol(raw)-1 != length(gsm)) stop('Expression and sample metadata columns disagree')
if (!identical(colnames(raw)[-1],gsm)) stop('Expression sample columns do not exactly match GSM order')
colnames(raw)[1] <- 'probe_id'
expr0 <- as.matrix(data.frame(lapply(raw[-1], as.numeric), check.names=FALSE))
rownames(expr0) <- raw$probe_id
expr0 <- expr0[,use,drop=FALSE]
meta <- manifest[use,,drop=FALSE]
# GPL annotation; only rows with one unambiguous gene symbol are retained.
gl <- readLines(gzfile(gplfile), warn=FALSE)
gb <- which(gl == '!platform_table_begin'); ge <- which(gl == '!platform_table_end')
if (length(gb)!=1 || length(ge)!=1) stop('GPL table boundary missing/multiple')
gpl <- read.delim(text=paste(gl[(gb+1):(ge-1)],collapse='\n'), sep='\t', header=TRUE,
                  check.names=FALSE, quote='"', comment.char='', stringsAsFactors=FALSE)
symcol <- grep('^Gene symbol$', names(gpl), value=TRUE)
if (length(symcol)!=1) stop('Cannot find unique GPL Gene symbol column')
ann <- data.frame(probe_id=as.character(gpl$ID), symbol=trimws(as.character(gpl[[symcol]])))
ann$single <- !is.na(ann$symbol) & ann$symbol!='' & !ann$symbol %in% c('---','NA','na') & !grepl('///|//|;|,',ann$symbol) & !grepl('|',ann$symbol,fixed=TRUE) & !grepl('[',ann$symbol,fixed=TRUE) & !grepl(']',ann$symbol,fixed=TRUE)
ann <- ann[ann$single & !duplicated(ann$probe_id),]
idx <- match(rownames(expr0),ann$probe_id)
keep <- !is.na(idx)
exprp <- expr0[keep,,drop=FALSE]
symbol <- ann$symbol[idx[keep]]
# mean log2 expression over probes per gene, the prespecified probe rule.
exprg <- rowsum(exprp, group=symbol, reorder=FALSE) / as.vector(table(factor(symbol, levels=unique(symbol))))
audit_gene <- unique(symbol)[1]; audit_manual <- colMeans(exprp[symbol==audit_gene,,drop=FALSE]); audit_computed <- as.numeric(exprg[audit_gene,]); stopifnot(isTRUE(all.equal(unname(audit_manual),unname(audit_computed),tolerance=1e-10))); write.csv(data.frame(gene=audit_gene,max_abs_difference=max(abs(audit_manual-audit_computed))),file.path(root,'GSE33630_15_probe_mean_manual_check.csv'),row.names=FALSE)
probe_counts <- as.data.frame(table(symbol),stringsAsFactors=FALSE); names(probe_counts)<-c('gene','retained_probe_n')
write.csv(data.frame(input_probe_n=nrow(expr0),single_symbol_probe_n=sum(keep),dropped_probe_n=sum(!keep),mapped_gene_n=nrow(exprg)),file.path(root,'GSE33630_11_probe_mapping_qc.csv'),row.names=FALSE)
write.csv(probe_counts,file.path(root,'GSE33630_12_gene_probe_counts.csv'),row.names=FALSE)
# Convert sample-level replicates to patient x tissue means; patient-level units thereafter.
keys <- paste(meta$patient_id,meta$tissue_group,sep='|')
uk <- unique(keys)
exprc <- sapply(uk,function(k) rowMeans(exprg[,keys==k,drop=FALSE],na.rm=TRUE))
colnames(exprc) <- uk
cm <- meta[match(uk,keys),c('patient_id','tissue_group')]
cm$block_n <- as.integer(table(factor(keys,levels=uk)))
cm$source_GSMs <- vapply(uk,function(k) paste(meta$GSM[keys==k],collapse=';'), '')
write.csv(cm,file.path(root,'GSE33630_13_patient_group_samples.csv'),row.names=FALSE)
rownames(cm) <- uk
# paired patients only for limma patient-block fixed effects.
tab <- table(cm$patient_id,cm$tissue_group)
paired <- rownames(tab)[tab[,'PTC']>0 & tab[,'Normal thyroid']>0]
sel <- cm$patient_id %in% paired
exprc <- exprc[,sel,drop=FALSE]; cm <- cm[sel,,drop=FALSE]
cm$tissue_group <- factor(cm$tissue_group, levels=c('Normal thyroid','PTC'))
if (any(!is.finite(exprc))) stop('Non-finite expression values present')
if (any(apply(exprc,1,var)==0)) exprc <- exprc[apply(exprc,1,var)>0,,drop=FALSE]
stopifnot(length(unique(cm$patient_id))==44,nrow(cm)==88)
if (length(unique(cm$patient_id)) < 3) stop('Fewer than 3 paired patients; do not run inferential model')
manifest$included_in_paired_analysis <- manifest$include & manifest$patient_id %in% paired & manifest$GSM %in% meta$GSM
write.csv(manifest,file.path(root,'GSE33630_02_sample_manifest_all.csv'),row.names=FALSE)
write.csv(data.frame(patient_id=cm$patient_id,tissue_group=cm$tissue_group,block_n=cm$block_n,source_GSMs=cm$source_GSMs),file.path(root,'GSE33630_13_patient_group_samples.csv'),row.names=FALSE)
write.csv(data.frame(gene=rownames(exprc),exprc,check.names=FALSE),file.path(root,'GSE33630_14_patient_group_expression.csv'),row.names=FALSE)
design <- model.matrix(~ patient_id + tissue_group, data=cm)
resid_df <- nrow(design)-qr(design)$rank
if (qr(design)$rank < ncol(design) || resid_df < 1) stop('Design matrix is rank deficient or has no residual degrees of freedom')
write.csv(data.frame(sample=rownames(design),design,check.names=FALSE),file.path(root,'GSE33630_17_design_matrix.csv'),row.names=FALSE)
write.csv(data.frame(n_patients=length(unique(cm$patient_id)),n_patient_group_observations=nrow(cm),design_rank=qr(design)$rank,residual_df=resid_df,log2_min=min(exprc),log2_max=max(exprc)),file.path(root,'GSE33630_18_model_qc.csv'),row.names=FALSE)
quantiles <- t(apply(exprc,2,quantile,probs=c(0,.25,.5,.75,1),na.rm=TRUE))
write.csv(data.frame(sample=rownames(quantiles),tissue_group=cm$tissue_group,quantiles,check.names=FALSE),file.path(root,'GSE33630_19_expression_distribution.csv'),row.names=FALSE)
if (min(exprc) < -10 || max(exprc)>50) stop('fRMA matrix value range is implausible for log2 normalized expression; inspect before analysis')
fit <- eBayes(lmFit(exprc, design), robust=TRUE)
coef <- grep('tissue_groupPTC',colnames(design))
if(length(coef)!=1) stop('PTC contrast not found')
deg <- topTable(fit,coef=coef,number=Inf,sort.by='none',adjust.method='BH')
deg$gene <- rownames(deg); deg$DE_threshold <- deg$adj.P.Val<0.05 & abs(deg$logFC)>=1
deg <- deg[,c('gene','logFC','AveExpr','t','P.Value','adj.P.Val','B','DE_threshold')]
write.csv(deg,file.path(root,'GSE33630_03_all_measured_gene_DE.csv'),row.names=FALSE)
# Fixed candidate sets and thyroid-disease source list from historical workbook artifacts.
sets <- list(XYT=read.csv('v20241015/output/02_common_gene_XYT.csv')$Com_Gene,
             XKC=read.csv('v20241015/output/02_common_gene_XKC.csv')$Com_Gene,
             ZBM=read.csv('v20241015/output/02_common_gene_ZBM.csv')$Com_Gene)
# Disease sensitivity set is the full historical GeneCards thyroid tumor list (Gene Symbol).
disease <- read.csv('v20241015/input/GeneCards-SearchResults_Thyroid_tumor.csv',check.names=FALSE,stringsAsFactors=FALSE)
dis <- unique(trimws(as.character(disease[['Gene Symbol']])))
measured <- deg$gene; degenes <- deg$gene[deg$DE_threshold]
res <- list(); coverage <- list()
for(nm in names(sets)) {
  s0 <- unique(trimws(sets[[nm]])); s <- intersect(s0,measured)
  for(bgname in c('all_measured','disease_measured')) {
    bg <- if(bgname=='all_measured') measured else intersect(measured,dis)
    sg <- intersect(s,bg); dx <- intersect(degenes,bg)
    a <- length(intersect(sg,dx)); b <- length(sg)-a; c <- length(dx)-a; d <- length(bg)-a-b-c
    ft <- if(length(sg)>0 && length(bg)>0) fisher.test(matrix(c(a,b,c,d),nrow=2),alternative='greater') else NULL
    res[[length(res)+1]] <- data.frame(set=nm,background=bgname,background_n=length(bg),set_n=length(sg),DE_in_set=a,
       odds_ratio=if(is.null(ft)) NA_real_ else unname(ft$estimate),p_value=if(is.null(ft)) NA_real_ else ft$p.value)
  }
  coverage[[nm]] <- data.frame(set=nm,input_n=length(s0),measured_n=length(s),unmeasured_n=length(setdiff(s0,measured)),
     DE_n=length(intersect(s,degenes)),measured_genes=paste(s,collapse=';'),DE_genes=paste(intersect(s,degenes),collapse=';'))
}
er <- do.call(rbind,res); er$FDR_6tests <- p.adjust(er$p_value,method='BH')
write.csv(er,file.path(root,'GSE33630_04_candidate_set_enrichment.csv'),row.names=FALSE)
write.csv(do.call(rbind,coverage),file.path(root,'GSE33630_05_candidate_coverage.csv'),row.names=FALSE)
write.csv(do.call(rbind,lapply(unique(c(unlist(sets),'CYP3A4','CACNA1H','POLB')),function(g) {z<-deg[deg$gene==g,,drop=FALSE]; if(!nrow(z)) data.frame(gene=g,logFC=NA,AveExpr=NA,t=NA,P.Value=NA,adj.P.Val=NA,B=NA,DE_threshold=NA) else z})),
          file.path(root,'GSE33630_06_candidate_gene_DE.csv'),row.names=FALSE)
# QC plots over patient x tissue samples.
pdf(file.path(root,'GSE33630_07_QC_PCA.pdf'),width=7,height=6)
pc <- prcomp(t(exprc),center=TRUE,scale.=FALSE)
plot(pc$x[,1],pc$x[,2],col=ifelse(cm$tissue_group=='PTC','#b33b35','#2878a5'),pch=19,
     xlab=paste0('PC1 (',round(100*summary(pc)$importance[2,1],1),'%)'),
     ylab=paste0('PC2 (',round(100*summary(pc)$importance[2,2],1),'%)'))
text(pc$x[,1],pc$x[,2],labels=cm$patient_id,pos=3,cex=.55)
legend('topright',legend=levels(cm$tissue_group),col=c('#2878a5','#b33b35'),pch=19)
dev.off()
pdf(file.path(root,'GSE33630_08_volcano.pdf'),width=7,height=6)
plot(deg$logFC,-log10(pmax(deg$adj.P.Val,1e-300)),pch=16,cex=.45,col=ifelse(deg$DE_threshold,'#b33b35','#999999'),
     xlab='log2 fold change (PTC vs normal thyroid)',ylab='-log10 BH FDR')
abline(v=c(-1,1),h=-log10(.05),lty=2)
dev.off()
writeLines(capture.output(sessionInfo()),file.path(root,'GSE33630_09_session_info.txt'))
writeLines(c('Source matrix: GSE33630, submitter RMA-normalized log2 expression (affy). No additional log transform or between-array normalization.',
'44 matched patients; 5 unpaired PTC, 1 unpaired normal and 11 ATC excluded. Patient ID is the numeric sample-title prefix; T and N denote matched tissue types.',
'RMA values aggregated by single-gene probe annotations; same GPL570 annotation and candidate sets as the first cohort.',
'Paired limma model with robust eBayes. BH<0.05 and absolute log2FC>=1; six one-sided Fisher tests with BH; exact two-sided confidence intervals computed separately.'),file.path(root,'GSE33630_methods.txt'))
# Exact two-sided confidence intervals; preserve one-sided primary P values.
er$CI_low<-NA_real_;er$CI_high<-NA_real_
for(i in seq_len(nrow(er))) {
 bg<-if(er$background[i]=='all_measured') measured else intersect(measured,dis)
 sg<-intersect(sets[[er$set[i]]],bg);dx<-intersect(degenes,bg)
 a<-length(intersect(sg,dx));b<-length(sg)-a;c<-length(dx)-a;d<-length(bg)-a-b-c
 z<-fisher.test(matrix(c(a,b,c,d),nrow=2),alternative='two.sided')
 er$CI_low[i]<-z$conf.int[1];er$CI_high[i]<-z$conf.int[2]
}
write.csv(er,file.path(root,'GSE33630_enrichment_with_CI.csv'),row.names=FALSE)
first<-read.csv('edit.Paper/output/20261024_25_transcriptome_validation/20261024_03_all_measured_gene_DE.csv')
both<-merge(first,deg,by='gene',suffixes=c('_first','_second'))
both$same_direction<-sign(both$logFC_first)==sign(both$logFC_second)
write.csv(both,file.path(root,'GSE33630_cross_cohort_all_genes.csv'),row.names=FALSE)
summary<-lapply(c('all_measured',names(sets)),function(nm){
 z<-if(nm=='all_measured')both else both[both$gene %in% sets[[nm]],]
 f<-z[z$DE_threshold_first,]
 data.frame(set=nm,n=nrow(z),spearman_rho=cor(z$logFC_first,z$logFC_second,method='spearman'),first_DE=nrow(f),first_DE_same_direction=sum(f$same_direction),first_DE_replicated_same_direction=sum(f$DE_threshold_second & f$same_direction))
})
write.csv(do.call(rbind,summary),file.path(root,'GSE33630_concordance_summary.csv'),row.names=FALSE)
print(er);print(do.call(rbind,summary));print(deg[deg$gene %in% c('POLB','CACNA1H','CYP3A4'),])
clinical<-read.csv(file.path(root,'11_verified_paired_manifest.csv'),stringsAsFactors=FALSE)
patient_status<-setNames(clinical$Status,clinical$patient)
strata_res<-list()
for(exposure in c('Ex','NE')) {
 idx<-patient_status[as.character(cm$patient_id)]==exposure
 submeta<-droplevels(cm[idx,]);submeta$patient_id<-factor(submeta$patient_id)
 design_sub<-model.matrix(~patient_id+tissue_group,data=submeta)
 fit_sub<-eBayes(lmFit(exprc[,idx],design_sub),robust=TRUE)
 tab_sub<-topTable(fit_sub,coef=grep('tissue_groupPTC',colnames(design_sub)),number=Inf,sort.by='none',adjust.method='BH')
 de_sub<-rownames(tab_sub)[tab_sub$adj.P.Val<.05 & abs(tab_sub$logFC)>=1]
 write.csv(data.frame(gene=rownames(tab_sub),tab_sub),file.path(root,paste0('GSE33630_',exposure,'_all_gene_DE.csv')),row.names=FALSE)
 for(nm in names(sets))for(bgname in c('all_measured','disease_measured')) {
  bg<-if(bgname=='all_measured')measured else intersect(measured,dis)
  sg<-intersect(sets[[nm]],bg);dx<-intersect(de_sub,bg)
  a<-length(intersect(sg,dx));b<-length(sg)-a;c<-length(dx)-a;d<-length(bg)-a-b-c
  test<-fisher.test(matrix(c(a,b,c,d),nrow=2),alternative='greater')
  strata_res[[length(strata_res)+1]]<-data.frame(stratum=exposure,n_patients=length(unique(submeta$patient_id)),set=nm,background=bgname,background_n=length(bg),background_DE=length(dx),set_n=length(sg),DE_in_set=a,odds_ratio=unname(test$estimate),p_value=test$p.value)
 }
}
strata<-do.call(rbind,strata_res);strata$FDR_12tests<-p.adjust(strata$p_value,'BH');write.csv(strata,file.path(root,'GSE33630_exposure_sensitivity.csv'),row.names=FALSE)
print(strata)
