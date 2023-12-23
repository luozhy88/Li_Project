# Li_Project
## NCBI_SRP388727
[-]文章来源  
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8005713/  
Gut Microbiome Alterations in Patients With Thyroid Nodules；  
Saliva microbiome changes in thyroid cancer and thyroid nodules patients；  
[-]数据来源  
https://trace.ncbi.nlm.nih.gov/Traces/sra?study=SRP388727  
[-]分组  
thyroid cancer patient;	healthy conrtrol;	thyroid nodules patient  
[-]bash
/data/zhiyu/software/nfcore/V2.7.0/nextflow run nf-core/ampliseq -r 2.7.0 -profile docker    --input_folder  /data/zhiyu/rawdata/16S/NCBI_SRP388727   --FW_primer ACTCCTACGGGAGGCAGCA   --RV_primer GGACTACHVGGGTWTCTAAT  --email 479321347@qq.com --picrust  --kraken2_ref_taxonomy silva=138  --extension "/*_{1,2}.fastq.gz" --outdir results


31：/data3/zhiyu/pipelines/TOFU_16S/NCBI/NCBI_SRP388727/results/Kraken/merge_data

### 进展
已经完成SRP388727该项目的phyloseq生成  

## NCBI_SRP151288
[-]文章来源    
Alterations in the gut microbiota and metabolite profiles of thyroid carcinoma patients  
https://onlinelibrary.wiley.com/doi/10.1002/ijc.32007  
Bacterial genomic DNA was amplified with primers 341F (CCTAYGGGRBGCASCAG) and 806R (GGACTACNNGGGTATCTAAT) specific to the V3-V4 hypervariable regions of the 16S rRNA gene.
[-]数据来源    
https://trace.ncbi.nlm.nih.gov/Traces/sra?study=SRP151288  
[-]bash
/data/zhiyu/software/nfcore/V2.7.0/nextflow run nf-core/ampliseq -r 2.7.0 -profile docker    --input_folder  /data_bk/zhiyu/Rawdata/16S/NCBI_SRP151288/rawdata   --FW_primer CCTAYGGGRBGCASCAG --RV_primer GGACTACNNGGGTATCTAAT --email 479321347@qq.com --picrust  --kraken2_ref_taxonomy silva=138  --extension "/*_{1,2}.fastq.gz" --outdir results

/data/zhiyu/software/nfcore/V2.7.0/nextflow run nf-core/ampliseq -r 2.7.0 -profile docker    --input_folder  /data_bk/zhiyu/Rawdata/16S/NCBI_SRP151288/rawdata   --FW_primer ACTGATGGACTACTTGG  --RV_primer GGCTACCCTATGGGDK  --email 479321347@qq.com --picrust  --kraken2_ref_taxonomy silva=138  --extension "/*_{1,2}.fastq.gz" --outdir results

[-]分组  
Control;TC  

31：/data3/zhiyu/pipelines/TOFU_16S/NCBI/NCBI_SRP151288/results/Kraken/merge_data



