library(LinDA)
library(phyloseq)
library(microbiome)
library(metagMisc)
library(openxlsx)
rm(list=ls())

# Create a directory named "LinDA"
dir.create("../output/differential/LinDA_noAdjust", recursive = TRUE)
path <- "../output/differential/LinDA_noAdjust/"



# Read the phyloseq object from a file and filter out taxa with zero counts
phy <- readRDS("../../Kraken/input/output/Batch2_phy_meta.rds")
# Meta=meta(phy)
phy


BATCH="SRP388727"
pseq = phyloseq::subset_samples(phy,grepl(BATCH,Batch )  &grepl("HC|TC",GROUP ))
pseq
Meta=meta(pseq)

#sample_data(pseq)$age <- scale(sample_data(pseq)$age)
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

#names(sample_data(pseq))[1] <- "group"




## GROUP_VAR <- paste0("sex+age+", "group")
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
  print(i)
  ## new phyloseq
  pseq.sub <- phyloseq::merge_phyloseq(agegut_split[[pair[1, i]]], agegut_split[[pair[2, i]]])
  
  
  pseq.temp <- core(pseq.sub, detection = 0.000001, prevalence = .1)
  
  ##
  for (my_rank in c("Species", "Genus", "Family", "Class","Order","Phylum")) {
    print(my_rank)
    ## this function works different from the original phyloseq::tax_glom
    pseq.new <- microbiome::aggregate_taxa(pseq.temp, my_rank)
    fil.names <- paste0(unique(microbiome::meta(pseq.new)[, MAIN])[1], "_vs_", unique(microbiome::meta(pseq.new)[, MAIN])[2])
    
    ## otu table and meta table
    pseq.new.otu <- microbiome::abundances(pseq.new)
    pseq.new.meta <- microbiome::meta(pseq.new)
    
    ### run Linda
    pseq.new.linda <- linda(
      otu.tab = pseq.new.otu,
      meta = pseq.new.meta,
      formula = paste0("~", GROUP_VAR),
      type = "count",
      adaptive = TRUE,
      imputation = FALSE,
      pseudo.cnt = 0.5,
      corr.cut = 0.1,
      p.adj.method = "BH",
      alpha = 0.1,
      prev.cut = 0,
      lib.cut = 1,
      winsor.quan = NULL,
      n.cores = 4
    )
    
    ###
    ## organize the results into a pretty dataframe for plotting.
    sum_res <- pseq.new.linda$output[[names(pseq.new.linda$output)[grepl(GROUP_VAR, names(pseq.new.linda$output))]]]
    
    if (length(sum_res) > 0) {
      
    
    ## we can take only the significantly differential taxa
    sig_linda <- subset(sum_res, reject == TRUE)
    
    # Adding taxonomic labels
    taxa_info <- data.frame(tax_table(pseq.new))
    sum_res <- merge(sum_res, taxa_info, by = 0)
    sig_linda <- merge(sig_linda, taxa_info, by = 0)
    
    ## take the reference level
    notREF <- gsub(pattern = MAIN, replacement = "", names(pseq.new.linda[["../output"]]))
    REF <- gsub(pattern = paste0("_vs_", notREF), replacement = "", fil.names)
    
    ##
    fil.names.new <- paste0(fil.names, "_***", REF, "_as_reference***",BATCH,"_")
    ## save the ../output for each comparison, reference shown on the file name
    openxlsx::write.xlsx(sum_res, file = paste0(path, fil.names.new, my_rank, "__all_linda2_cal.xlsx"))
    
    
    ## reformat the ../output for barplot
    sig_linda$P.adj <-
      ifelse(sig_linda$padj <= 0.05,
             "< 0.05",
             ifelse(sig_linda$padj <= 0.1 & sig_linda$padj > 0.05,
                    "0.05 - 0.1",
                    "0.1 - 0.2 "
             )
      )
    ##
    sig_linda$group <- ifelse(sig_linda$log2FoldChange < 0, REF, notREF)
    
    ### do not need to plot if there were no significantly differential taxa
    if (nrow(sig_linda) > 0) {
      ##
      sig_linda$taxa <- sig_linda$Row.names
      openxlsx::write.xlsx(sig_linda, file = paste0(path, fil.names.new, my_rank, "__sig_linda2_cal.xlsx"))
      ## one can also plot the results based on log fold change LFC but colored by padj range
      ggpubr::ggbarplot(sig_linda,
                        x = "taxa", y = "log2FoldChange",
                        fill = "group",
                        color = "white",
                        palette = c("#00A4D3","#DD498D"),
                        # palette = c("#DD498D","#00A4D3"),
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
      ggsave(filename = paste0(path, fil.names.new, my_rank, "__linda2_cal_barplot.pdf"), width = 8, height = nrow(sig_linda) * 0.1 + 2)
      
      # ### we could also color the bar by family or other rank levels if needed
      # ## one can also plot the results based on lfc.
      # ggpubr::ggbarplot(sig_linda,
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
      # ggsave(filename = paste0(path, fil.names.new, my_rank, "__linda2_cal_barplot_FamilyColor.pdf"), width = 12, height = nrow(sig_linda) * 0.1 + 2)
    }
  }
}

}