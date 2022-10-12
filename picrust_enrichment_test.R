
results <- data.frame(matrix(ncol = 5, nrow= 1))

colnames(results) <- c('Pathway_ID', 'Phase_Comparison', 'p_value', 'Odds_ratio', 'barcode')

results$p_value <- as.numeric(results$p_value)


#ko maps, pathway ids and pathway description all obtained from files within mapfiles PICRUST2 https://github.com/picrust/picrust2/tree/master/picrust2/default_files/pathway_mapfiles
komap <- read.delim("~/KEGG_maps/KEGG_map.tsv", colClasses = "character") #table of all picrust ko with class, subclass, pathway id, pathway name and ko name
pathway_id <- read.delim("~/KEGG_pathways_picrust2_v2.4.1/KEGG_pathways_to_KO_v2.tsv" ) #table of all ko_ids and pathway_ids which the ko is in
pathway_descrip <- read.delim("~/KEGG_pathways_picrust2_v2.4.1/KEGG_pathways_info.txt") #table of all pathwayids and pathway names


for(barcode in c("16sv1", "16sv4")){
  for (phase1 in c("SemiPristine", "Eutrophic", "Pesticide" )){
    for(phase2 in c("Eutrophic", "Pesticide", "Recovery" )){
      if(phase1 == phase2){
        next
      }else if((phase1 == "Pesticide") & (phase2 == "Eutrophic")){
        next
      }
      
      phasecombo <- paste0(phase1, phase2)
      stage <- paste0(barcode, phasecombo)
      print(stage)
      for(pathway in unique(pathway_id$pathway_id) ){
      
        
        sig_ko_file <- paste0(barcode,"/ancom_ko_data/ancom_", phase1, "_", phase2, ".tsv" )
        ancom_ko <- read.delim(sig_ko_file)
        names(ancom_ko)[names(ancom_ko) == 'X'] <- 'ko_id'
        sig_ko <- ancom_ko$ko_id[ancom_ko$Reject.null.hypothesis == "True"]
        
        kos_in_path <- pathway_id$ko_id[pathway_id$pathway_id == pathway]
        num_kos_in_path <- length(unique(kos_in_path))
        num_ko_not_in_path <- length(unique(pathway_id$ko_id)) - num_kos_in_path
        
        sig_kos_in_path <- sig_ko[which(sig_ko %in% kos_in_path)]
        num_sig_kos_in_path <- length(sig_kos_in_path)
        
        num_nonsig_kos_in_path <- num_kos_in_path - num_sig_kos_in_path
        num_sig_kos_not_in_path <- length(unique(sig_ko)) - num_sig_kos_in_path
        
        num_total_ko <- length(unique(pathway_id$ko_id))
        num_nonsig_ko_not_in_path <- num_total_ko - (num_sig_kos_in_path + num_sig_kos_not_in_path + num_nonsig_kos_in_path)
          
        fish_table <-  data.frame(diff_abun = c(num_sig_kos_in_path,num_sig_kos_not_in_path  ),
                                  no_sig_diff = c(num_nonsig_kos_in_path, num_nonsig_ko_not_in_path))
          
        fish_res <- fisher.test(fish_table, alternative = "greater")
        add <- c(pathway, phasecombo, fish_res$p.value, fish_res$estimate, barcode)
        results <- rbind(results, add)
        
      }
    }
  }
}

forloopresults <- results
results <- na.omit(results)
write.table(results, file = "enriched_pathways_obsgreater_usingpicrustKEGGpath.txt")

results$p_value <- as.numeric(results$p_value) 
res_pval <- results[results$p_value <= 0.05,]

res_pval$Odds_ratio <- as.numeric(res_pval$Odds_ratio)
write.table(res_pval, file = "significantly_enriched_pathways_obsgreater_usingpicrustKEGGpath.txt")

