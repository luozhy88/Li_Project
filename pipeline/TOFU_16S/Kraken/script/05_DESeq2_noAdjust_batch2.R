library(DESeq2)
library(phyloseq)
library(microbiome)
library(metagMisc)
library(openxlsx)
rm(list=ls())

# Create a directory named "Deseq2"
dir.create("../output/differential/Deseq2_noAdjust", recursive = TRUE)
path <- "../output/differential/Deseq2_noAdjust/"


# Read the phyloseq object from a file and filter out taxa with zero counts
phy <- readRDS("../../Kraken/input/output/Batch2_phy_meta.rds")
# Meta=meta(phy)

BATCH="SRP388727"
pseq = phyloseq::subset_samples(phy,grepl(BATCH,Batch )  &grepl("HC|TC",GROUP ))
Meta=meta(pseq)


###############################################################
############### empty values for ranks replaced with "undefined"##
###############################################################
# Extract the tax_table from the phyloseq object
taxa <- as.data.frame(tax_table(pseq))

# Replace empty or NA values with "undefined"
taxa[is.na(taxa)] <- "undefined"
taxa[taxa == ""] <- "undefined"

# Convert the modified data frame back to a tax_table format
tax_table_modified <- tax_table(as.matrix(taxa))

# Assign the modified tax_table back to the phyloseq object
tax_table(pseq) <- tax_table_modified
###############################################################
###############################################################

pseq <- phyloseq::prune_taxa(taxa_sums(pseq) > 0, pseq)
pseq <- phyloseq::prune_samples(sample_sums(pseq) > 0, pseq)


###### group1 ---> group

# names(sample_data(pseq))[7] <- "group"



# GROUP_VAR <- paste0("sex+age+", "group")
 GROUP_VAR <- paste0("GROUP")

## here we refer to the main interest even though one has to add several covariates like above
MAIN <- trimws(GROUP_VAR, whitespace = ".*\\+")

### here I rename all names (unreliable at all) with ASVs. But we must know they are ASVs not Species
## better to make sure there are no improper names (empty) at Genus level
phyloseq::tax_table(pseq)[, "Species"] <-sprintf("%s%04d", "ASV", seq_along(taxa_names(pseq)))

############# split the phyloseq into a list based on the variable of interest
#######
agegut_split <- metagMisc::phyloseq_sep_variable(pseq, MAIN)
############ make pairwise comparison
pair <- combn(names(agegut_split), 2)


for (i in 1:(length(pair) / 2)) {
  ## new phyloseq
  #i=1
  print(i)
  pseq.sub <- phyloseq::merge_phyloseq(agegut_split[[pair[1, i]]], agegut_split[[pair[2, i]]])


  pseq.temp <- core(pseq.sub, detection = 0.000001, prevalence = .1)

  ##
  for (my_rank in c("Species", "Genus", "Family", "Class","Order","Phylum")) {
    #my_rank="Species"
    print(my_rank)
    ## this function works different from the original phyloseq::tax_glom
    pseq.new <- microbiome::aggregate_taxa(pseq.temp, my_rank)
    fil.names <- paste0(unique(microbiome::meta(pseq.new)[, MAIN])[1], "_vs_", unique(microbiome::meta(pseq.new)[, MAIN])[2])

    ## run DESeq2 using the following step, do not run directly with DESeq () function
    gutdds <- phyloseq_to_deseq2(pseq.new, as.formula(paste0("~ ", GROUP_VAR)))
    gutdds <- estimateSizeFactors(gutdds, type = "poscounts")
    gutdds <- estimateDispersions(gutdds, fitType = "parametric")
    gutdds <- nbinomWaldTest(gutdds)

    # we can also use the following one, no big difference at all from the above classic DESeq2 way
    # gutdds <- DESeq2::DESeq(gutdds,
    #                         sfType = "poscounts",
    #                         fitType = "glmGamPoi",
    #                         test = "LRT",
    #                         reduced = ~ 1 )
    ##
    REF <- unique(microbiome::meta(pseq.new)[, MAIN])[1]
    notREF <- unique(microbiome::meta(pseq.new)[, MAIN])[2]
    CONTRAST <- c(MAIN, notREF, REF)
    ##
    sum_res_deseq <- results(gutdds, contrast = CONTRAST)
    sum_res <- as.data.frame(sum_res_deseq)
    ## we can take only the significantly differential taxa
    sig_deseq <- subset(sum_res, padj < 0.05)
    # Adding taxonomic labels
    taxa_info <- data.frame(tax_table(pseq.new))
    sum_res <- merge(sum_res, taxa_info, by = 0)
    sig_deseq <- merge(sig_deseq, taxa_info, by = 0)

    ##
    fil.names.new <- paste0(fil.names, "_***", REF, "_as_reference***",BATCH,"_")
    ## save the ../output for each comparison, reference shown on the file name
    openxlsx::write.xlsx(sum_res, file = paste0(path, fil.names.new, my_rank, "__all_DESeq2_cal.xlsx"))
    

    ## reformat the ../output for barplot
    sig_deseq$P.adj <-
      ifelse(sig_deseq$padj <= 0.05,
        "< 0.05",
        ifelse(sig_deseq$padj <= 0.1 & sig_deseq$padj > 0.05,
          "0.05 - 0.1",
          "0.1 - 0.2 "
        )
      )
    ##
    sig_deseq$group <- ifelse(sig_deseq$log2FoldChange < 0, REF, notREF)
    
    ### do not need to plot if there were no significantly differential taxa
    if (nrow(sig_deseq) > 0) {
      ##
      sig_deseq$taxa <- sig_deseq$Row.names

      ## one can also plot the results based on log fold change LFC but colored by padj range
      ggpubr::ggbarplot(sig_deseq,
        x = "taxa", y = "log2FoldChange",
        fill = "group",
        color = "white",
        palette = c("#00A4D3","#DD498D"),
        sort.val = "asc",
        sort.by.sex = FALSE,
        x.text.angle = 90,
        ylab = "log2FoldChange",
        xlab = "",
        rotate = TRUE,
        title = paste0("[", notREF, "]", "-", "[", REF, "]"),
        ggtheme = theme_minimal()
      ) 
      
      ## save the plot
      ##
      ggsave(filename = paste0(path, fil.names.new, my_rank, "__DESeq2_cal_barplot.pdf"), width = 6, height = nrow(sig_deseq) * 0.1 + 2)
      openxlsx::write.xlsx(sig_deseq, file = paste0(path, fil.names.new, my_rank, "__sig_DESeq2_cal.xlsx"))
      # ### we could also color the bar by family or other rank levels if needed
      # ## one can also plot the results based on lfc.
      # ggpubr::ggbarplot(sig_deseq,
      #   x = "taxa", y = "log2FoldChange",
      #   fill = "Family",
      #   color = "white",
      #  # palette = unname(see::metro_colors(1:25)),
      #   sort.val = "desc",
      #   sort.by.sex = FALSE,
      #   x.text.angle = 90,
      #   ylab = "log2FoldChange",
      #   xlab = "",
      #   rotate = TRUE,
      #   title = paste0("[", notREF, "]", "-", "[", REF, "]"),
      #   ggtheme = theme_minimal()
      # )
      # ##
      # ## save the plot
      # ggsave(filename = paste0(path, fil.names.new, my_rank, "__DESeq2_cal_barplot_FamilyColor.pdf"), width = 12, height = nrow(sig_deseq) * 0.1 + 2)
    }
  }
}

