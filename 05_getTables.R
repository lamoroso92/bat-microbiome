#### Script - Count tables ####

# Setup

library(tidyverse)

rm(list=ls())

getwd()
setwd("/data/dgpinheiro/MetaBat/")
getwd()

my_path = paste0(getwd(), "/")

# Functions

sample_groups = read.table("scripts/Sample_Details_LibrariesSequenced-01092021.txt")[2:3]

metadata = read.table("metadata_20220331.tsv", sep = "\t", header = T, quote = "")

names(sample_groups) = c("samples", "groups")

my_groups = unique(sample_groups$group)
my_groups

master_table = function(group){
  
  tmp = read.delim(paste0(my_path, "output/kraken_005//", group, "/", group, "_labels.txt"))
  names(tmp) = c("ID", "taxonomy")
  
  tmp = separate(tmp, taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "; ")
  
  tmp$Kingdom = ifelse(tmp$Kingdom == "k__", "k__Unclassified", tmp$Kingdom)
  tmp$Phylum = ifelse(tmp$Phylum == "p__", "p__Unclassified", tmp$Phylum)
  tmp$Class = ifelse(tmp$Class == "c__", "c__Unclassified", tmp$Class)
  tmp$Order = ifelse(tmp$Order == "o__", "o__Unclassified", tmp$Order)
  tmp$Family = ifelse(tmp$Family == "f__", "f__Unclassified", tmp$Family)
  tmp$Genus = ifelse(tmp$Genus == "g__", "g__Unclassified", tmp$Genus)
  tmp$Species = ifelse(tmp$Species == "s__", "s__Unclassified", tmp$Species)
  
  samples = sample_groups %>% filter(groups == group) %>% select(samples) %>% arrange(samples) %>% unlist()
  
  for (sample in samples) {
    
    tmp_sample = read.delim(paste0(my_path, "output/kallisto/", group, "/", sample, "/", "Kb_", group, "_", sample, "/", "abundance.tsv"))[,c(1,4)]
    names(tmp_sample) = c("ID", sample)
    
    tmp = merge(tmp, tmp_sample, by = "ID", sort = F)
    
  }
  
  tmp[,9:length(tmp)] = sapply(tmp[,9:length(tmp)], function(x) round(x))
  
  tmp_count_samples = as.data.frame(sapply(tmp[,9:length(tmp)], function(x) ifelse( x > 0, 1, 0)))
  tmp_count_samples = tmp_count_samples %>% mutate("count_samples" = rowSums(.[1:dim(tmp_count_samples)[2]]))
  
  tmp_count_sum = tmp %>% mutate("count_sum" = rowSums(.[9:length(tmp)]))
  tmp_final = cbind(tmp[,1:8], "count_sum" = tmp_count_sum$count_sum, "count_sample" = tmp_count_samples$count_samples, select(tmp, 9:length(tmp)))
  
  tmp_final = filter(tmp_final, count_sample > 0)
  
  return(tmp_final)
  
}

# "SpinturnicidaeMite" "MacronyssidaeMite" "BatFly" "BatSwab" "TickLarvae" "Undetermined"

master_bf = master_table("BatFly")
master_bs = master_table("BatSwab")
master_sm = master_table("SpinturnicidaeMite")
master_mm = master_table("MacronyssidaeMite")
master_tl = master_table("TickLarvae")

classifications = function(tb) {
  
  tmp = data.frame("Kingdom" = NA, "Phylum" = NA, "Class" = NA, "Order" = NA, "Family" = NA, "Genus" = NA, "Species" = NA)
  
  tmp[1,"Kingdom"] = length(which(tb$Kingdom != "k__Unclassified"))
  tmp[1,"Phylum"] = length(which(tb$Phylum != "p__Unclassified"))
  tmp[1,"Class"] = length(which(tb$Class != "c__Unclassified"))
  tmp[1,"Order"] = length(which(tb$Order != "o__Unclassified"))
  tmp[1,"Family"] = length(which(tb$Family != "f__Unclassified"))
  tmp[1,"Genus"] = length(which(tb$Genus != "g__Unclassified"))
  tmp[1,"Species"] = length(which(tb$Species != "s__Unclassified"))

  tmp[2,"Kingdom"] = length(which(tb$Kingdom == "k__Unclassified"))
  tmp[2,"Phylum"] = length(which(tb$Phylum == "p__Unclassified"))
  tmp[2,"Class"] = length(which(tb$Class == "c__Unclassified"))
  tmp[2,"Order"] = length(which(tb$Order == "o__Unclassified"))
  tmp[2,"Family"] = length(which(tb$Family == "f__Unclassified"))
  tmp[2,"Genus"] = length(which(tb$Genus == "g__Unclassified"))
  tmp[2,"Species"] = length(which(tb$Species == "s__Unclassified"))
  
  row.names(tmp) = c("classified", "unclassified")
 
  return(tmp)
}

count_classifications_bf = classifications(master_bf)
count_classifications_bs = classifications(master_bs)
count_classifications_sm = classifications(master_sm)
count_classifications_mm = classifications(master_mm)
count_classifications_tl = classifications(master_tl)

filt_master = function(tb) {
  tmp = filter(tb, Phylum != "p__Unclassified")
  return(tmp)
}

master_bf_filt = filt_master(master_bf)
master_bs_filt = filt_master(master_bs)
master_sm_filt = filt_master(master_sm)
master_mm_filt = filt_master(master_mm)
master_tl_filt = filt_master(master_tl)

collaps_master = function(tb) {
  tmp = tb %>% select(-c(ID, count_sample))
  tmp2 = aggregate(. ~ Kingdom + Phylum + Class + Order + Family + Genus + Species, data = tmp, FUN = sum)
  
  tmp_count_samples = as.data.frame(sapply(tmp2[,9:length(tmp2)], function(x) ifelse( x > 0, 1, 0)))
  tmp_count_samples = tmp_count_samples %>% mutate("count_samples" = rowSums(.[1:dim(tmp_count_samples)[2]]))
  
  tmp_final = cbind(tmp2[,1:7], "count_sum" = tmp2$count_sum, "count_sample" = tmp_count_samples$count_samples, select(tmp2, 9:length(tmp2)))
  
  tmp_final = filter(tmp_final, count_sample > 1)
  
  return(tmp_final)
}

master_bf_collaps = collaps_master(master_bf_filt)
master_bs_collaps = collaps_master(master_bs_filt)
master_sm_collaps = collaps_master(master_sm_filt)
master_mm_collaps = collaps_master(master_mm_filt)
master_tl_collaps = collaps_master(master_tl_filt)

unique(c(master_bf$Kingdom, master_bs$Kingdom, master_mm$Kingdom, master_sm$Kingdom, master_tl$Kingdom))

all_phylum = unique(c(paste0(master_bf$Kingdom, "; ", master_bf$Phylum), 
                      paste0(master_bs$Kingdom, "; ", master_bs$Phylum),
                      paste0(master_mm$Kingdom, "; ", master_mm$Phylum),
                      paste0(master_sm$Kingdom, "; ", master_sm$Phylum), 
                      paste0(master_tl$Kingdom, "; ", master_tl$Phylum)))

sort(grep(pattern = "k__Eukaryota", x = all_phylum, value = T))

get_tax_tables = function(master_table){
  
  tmp = master_table %>% arrange(., desc(count_sum), desc(count_sample))
  
  tb_vir = tmp %>% filter(Kingdom == "k__Viruses") %>% arrange(., desc(count_sum), desc(count_sample))

  tb_prok = tmp %>% filter(Kingdom %in% c("k__Bacteria", "k__Archaea")) %>% arrange(., desc(count_sum), desc(count_sample))
  
  tb_euk = tmp %>% filter(Kingdom == "k__Eukaryota") %>% arrange(., desc(count_sum), desc(count_sample))
  
  tb_euk_filt = tmp %>% filter(Kingdom == "k__Eukaryota") %>% filter(Phylum %in% c("p__Apicomplexa", "p__Euglenozoa", "p__Nematoda", "p__Platyhelminthes", "p__Ascomycota", "p__Basidiomycota", "p__Chytridiomycota", "p__Glomeromycota", "p__Mucoromycota", "p__Microsporidia", "p__Zygomycota", "p__Oomycota")) %>% arrange(., desc(count_sum), desc(count_sample))
  
  return(list("All" = tmp,
              "Viruses" = tb_vir, 
              "Prokaryote" = tb_prok,
              "Eukaryota" = tb_euk,
              "Euk_filt" = tb_euk_filt))
}

tb_final_bf = get_tax_tables(master_bf_collaps)
tb_final_bs = get_tax_tables(master_bs_collaps)
tb_final_sm = get_tax_tables(master_sm_collaps)
tb_final_mm = get_tax_tables(master_mm_collaps)
tb_final_tl = get_tax_tables(master_tl_collaps)

names(my_groups) = c("sm", "mm", "bf", "bs", "tl", "un")
my_groups

dir.create(paste0("output/count_tables/"))

for (i in names(my_groups)[names(my_groups) != "un"]) {

  tmp = get(paste0("tb_final_", i))
  
  master = get(paste0("master_", i))
  
  sample_ids = data.frame("sample_id" = names(master)[11:length(master)])
  
  meta = merge(sample_ids, metadata, by = "sample_id", all.x = T) %>% 
    rename("NAME" = "sample_id")
    
  dir.create(paste0("output/count_tables/", unname(my_groups[i])))
  
  write.table(x = master, file = paste0("output/count_tables/", unname(my_groups[i]), "/", "00_master_table_", i, ".txt"), quote = F, sep = "\t", row.names = F, col.names = T)
  write.table(x = meta, file = paste0("output/count_tables/", unname(my_groups[i]), "/", "02_metadata_", i, ".txt"), quote = F, sep = "\t", row.names = F, col.names = T)
  
  if (nrow(tmp[["All"]]) > 0) {
    
    all.count = cbind(data.frame("ID" = paste0("ID_", seq_along(tmp[["All"]][[1]]))), tmp[["All"]][, 10:ncol(tmp[["All"]])]) %>% rename("NAME" = "ID")
    all.taxon = cbind(data.frame("ID" = paste0("ID_", seq_along(tmp[["All"]][[1]]))), tmp[["All"]][, 1:7]) %>% rename("TAXONOMY" = "ID")
    
    write.table(x = all.count, file = paste0("output/count_tables/", unname(my_groups[i]), "/", "01_all_count_", i, ".txt"), quote = F, sep = "\t", row.names = F, col.names = T)
    write.table(x = all.taxon, file = paste0("output/count_tables/", unname(my_groups[i]), "/", "03_all_tax_", i, ".txt"), quote = F, sep = "\t", row.names = F, col.names = T)
    
  }
  
  if (nrow(tmp[["Viruses"]]) > 0) {

  vir.count = cbind(data.frame("ID" = paste0("ID_", seq_along(tmp[["Viruses"]][[1]]))), tmp[["Viruses"]][, 10:ncol(tmp[["Viruses"]])]) %>% rename("NAME" = "ID")
  vir.taxon = cbind(data.frame("ID" = paste0("ID_", seq_along(tmp[["Viruses"]][[1]]))), tmp[["Viruses"]][, 1:7]) %>% rename("TAXONOMY" = "ID") %>% mutate("Kingdom" = str_replace(Kingdom, "Viruses", "Bacteria"))
  
  write.table(x = vir.count, file = paste0("output/count_tables/", unname(my_groups[i]), "/", "01_vir_count_", i, ".txt"), quote = F, sep = "\t", row.names = F, col.names = T)
  write.table(x = vir.taxon, file = paste0("output/count_tables/", unname(my_groups[i]), "/", "03_vir_tax_", i, ".txt"), quote = F, sep = "\t", row.names = F, col.names = T)
  
  }
  
  if (nrow(tmp[["Prokaryote"]]) > 0) {
  
  prok.count = cbind(data.frame("ID" = paste0("ID_", seq_along(tmp[["Prokaryote"]][[1]]))), tmp[["Prokaryote"]][, 10:ncol(tmp[["Prokaryote"]])]) %>% rename("NAME" = "ID")
  prok.taxon = cbind(data.frame("ID" = paste0("ID_", seq_along(tmp[["Prokaryote"]][[1]]))), tmp[["Prokaryote"]][, 1:7]) %>% rename("TAXONOMY" = "ID") %>% mutate("Kingdom" = str_replace(Kingdom, "Archaea", "Bacteria"))
  
  write.table(x = prok.count, file = paste0("output/count_tables/", unname(my_groups[i]), "/", "01_prok_count_", i, ".txt"), quote = F, sep = "\t", row.names = F, col.names = T)
  write.table(x = prok.taxon, file = paste0("output/count_tables/", unname(my_groups[i]), "/", "03_prok_tax_", i, ".txt"), quote = F, sep = "\t", row.names = F, col.names = T)
  
  }
  
  if (nrow(tmp[["Eukaryota"]]) > 0) {
  
  euk.count = cbind(data.frame("ID" = paste0("ID_", seq_along(tmp[["Eukaryota"]][[1]]))), tmp[["Eukaryota"]][, 10:ncol(tmp[["Eukaryota"]])]) %>% rename("NAME" = "ID")
  euk.taxon = cbind(data.frame("ID" = paste0("ID_", seq_along(tmp[["Eukaryota"]][[1]]))), tmp[["Eukaryota"]][, 1:7]) %>% rename("TAXONOMY" = "ID") %>% mutate("Kingdom" = str_replace(Kingdom, "Eukaryota", "Bacteria"))
  
  write.table(x = euk.count, file = paste0("output/count_tables/", unname(my_groups[i]), "/", "01_euk_count_", i, ".txt"), quote = F, sep = "\t", row.names = F, col.names = T)
  write.table(x = euk.taxon, file = paste0("output/count_tables/", unname(my_groups[i]), "/", "03_euk_tax_", i, ".txt"), quote = F, sep = "\t", row.names = F, col.names = T)
  
  }
  
  if (nrow(tmp[["Euk_filt"]]) > 0) {
  
  eukf.count = cbind(data.frame("ID" = paste0("ID_", seq_along(tmp[["Euk_filt"]][[1]]))), tmp[["Euk_filt"]][, 10:ncol(tmp[["Euk_filt"]])]) %>% rename("NAME" = "ID")
  eukf.taxon = cbind(data.frame("ID" = paste0("ID_", seq_along(tmp[["Euk_filt"]][[1]]))), tmp[["Euk_filt"]][, 1:7]) %>% rename("TAXONOMY" = "ID") %>% mutate("Kingdom" = str_replace(Kingdom, "Eukaryota", "Bacteria"))
  
  write.table(x = eukf.count, file = paste0("output/count_tables/", unname(my_groups[i]), "/", "01_eukf_count_", i, ".txt"), quote = F, sep = "\t", row.names = F, col.names = T)
  write.table(x = eukf.taxon, file = paste0("output/count_tables/", unname(my_groups[i]), "/", "03_eukf_tax_", i, ".txt"), quote = F, sep = "\t", row.names = F, col.names = T)
  
  }
  
}

save.image("output/count_tables/count_tables.RData")

# 2022-04-04

master_bf_c = master_bf_collaps %>% mutate("taxonomy" = paste(Kingdom, Phylum, Class, Order, Family, Genus, Species, sep = "; ")) %>% select(-c(Kingdom:count_sample))
master_bs_c = master_bs_collaps %>% mutate("taxonomy" = paste(Kingdom, Phylum, Class, Order, Family, Genus, Species, sep = "; ")) %>% select(-c(Kingdom:count_sample))
master_sm_c = master_sm_collaps %>% mutate("taxonomy" = paste(Kingdom, Phylum, Class, Order, Family, Genus, Species, sep = "; ")) %>% select(-c(Kingdom:count_sample))
master_mm_c = master_mm_collaps %>% mutate("taxonomy" = paste(Kingdom, Phylum, Class, Order, Family, Genus, Species, sep = "; ")) %>% select(-c(Kingdom:count_sample))
master_tl_c = master_tl_collaps %>% mutate("taxonomy" = paste(Kingdom, Phylum, Class, Order, Family, Genus, Species, sep = "; ")) %>% select(-c(Kingdom:count_sample))

master_master = data.frame("taxonomy" = NA)

master_master = merge(master_master, master_bf_c, by = "taxonomy", all = T)
master_master = merge(master_master, master_bs_c, by = "taxonomy", all = T)
master_master = merge(master_master, master_sm_c, by = "taxonomy", all = T)
master_master = merge(master_master, master_mm_c, by = "taxonomy", all = T)
master_master = merge(master_master, master_tl_c, by = "taxonomy", all = T)
master_master = master_master[-nrow(master_master),]

master_master = 
  separate(master_master, col = "taxonomy", into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "; ", remove = T) %>% 
  replace(is.na(.), 0) %>% 
  mutate("count_sum" = rowSums(.[8:length(.)])) %>% 
  select(Kingdom:Species, count_sum, matches("\\d$"))

collaps_master2 = function(tb) {
  tmp = tb
  tmp2 = aggregate(. ~ Kingdom + Phylum + Class + Order + Family + Genus + Species, data = tmp, FUN = sum)
  
  tmp_count_samples = as.data.frame(sapply(tmp2[,9:length(tmp2)], function(x) ifelse( x > 0, 1, 0)))
  tmp_count_samples = tmp_count_samples %>% mutate("count_samples" = rowSums(.[1:dim(tmp_count_samples)[2]]))
  
  tmp_final = cbind(tmp2[,1:7], "count_sum" = tmp2$count_sum, "count_sample" = tmp_count_samples$count_samples, select(tmp2, 9:length(tmp2)))
  
  tmp_final = filter(tmp_final, count_sample > 1)
  
  return(tmp_final)
}

master_all = collaps_master2(master_master)

tb_final_all = get_tax_tables(master_collaps)

for (i in "all") {
  
  tmp = get(paste0("tb_final_", i))
  
  master = get(paste0("master_", i))
  
  sample_ids = data.frame("sample_id" = names(master)[11:length(master)])
  
  meta = merge(sample_ids, metadata, by = "sample_id", all.x = T) %>% 
    rename("NAME" = "sample_id")
  
  dir.create(paste0("output/count_tables/", i))
  
  write.table(x = master, file = paste0("output/count_tables/", i, "/", "00_master_table_", i, ".txt"), quote = F, sep = "\t", row.names = F, col.names = T)
  write.table(x = meta, file = paste0("output/count_tables/", i, "/", "02_metadata_", i, ".txt"), quote = F, sep = "\t", row.names = F, col.names = T)
  
  if (nrow(tmp[["All"]]) > 0) {
    
    all.count = cbind(data.frame("ID" = paste0("ID_", seq_along(tmp[["All"]][[1]]))), tmp[["All"]][, 10:ncol(tmp[["All"]])]) %>% rename("NAME" = "ID")
    all.taxon = cbind(data.frame("ID" = paste0("ID_", seq_along(tmp[["All"]][[1]]))), tmp[["All"]][, 1:7]) %>% rename("TAXONOMY" = "ID")
    
    write.table(x = all.count, file = paste0("output/count_tables/", i, "/", "01_all_count_", i, ".txt"), quote = F, sep = "\t", row.names = F, col.names = T)
    write.table(x = all.taxon, file = paste0("output/count_tables/", i, "/", "03_all_tax_", i, ".txt"), quote = F, sep = "\t", row.names = F, col.names = T)
    
  }
  
  if (nrow(tmp[["Viruses"]]) > 0) {
    
    vir.count = cbind(data.frame("ID" = paste0("ID_", seq_along(tmp[["Viruses"]][[1]]))), tmp[["Viruses"]][, 10:ncol(tmp[["Viruses"]])]) %>% rename("NAME" = "ID")
    vir.taxon = cbind(data.frame("ID" = paste0("ID_", seq_along(tmp[["Viruses"]][[1]]))), tmp[["Viruses"]][, 1:7]) %>% rename("TAXONOMY" = "ID") %>% mutate("Kingdom" = str_replace(Kingdom, "Viruses", "Bacteria"))
    
    write.table(x = vir.count, file = paste0("output/count_tables/", i, "/", "01_vir_count_", i, ".txt"), quote = F, sep = "\t", row.names = F, col.names = T)
    write.table(x = vir.taxon, file = paste0("output/count_tables/", i, "/", "03_vir_tax_", i, ".txt"), quote = F, sep = "\t", row.names = F, col.names = T)
    
  }
  
  if (nrow(tmp[["Prokaryote"]]) > 0) {
    
    prok.count = cbind(data.frame("ID" = paste0("ID_", seq_along(tmp[["Prokaryote"]][[1]]))), tmp[["Prokaryote"]][, 10:ncol(tmp[["Prokaryote"]])]) %>% rename("NAME" = "ID")
    prok.taxon = cbind(data.frame("ID" = paste0("ID_", seq_along(tmp[["Prokaryote"]][[1]]))), tmp[["Prokaryote"]][, 1:7]) %>% rename("TAXONOMY" = "ID") %>% mutate("Kingdom" = str_replace(Kingdom, "Archaea", "Bacteria"))
    
    write.table(x = prok.count, file = paste0("output/count_tables/", i, "/", "01_prok_count_", i, ".txt"), quote = F, sep = "\t", row.names = F, col.names = T)
    write.table(x = prok.taxon, file = paste0("output/count_tables/", i, "/", "03_prok_tax_", i, ".txt"), quote = F, sep = "\t", row.names = F, col.names = T)
    
  }
  
  if (nrow(tmp[["Eukaryota"]]) > 0) {
    
    euk.count = cbind(data.frame("ID" = paste0("ID_", seq_along(tmp[["Eukaryota"]][[1]]))), tmp[["Eukaryota"]][, 10:ncol(tmp[["Eukaryota"]])]) %>% rename("NAME" = "ID")
    euk.taxon = cbind(data.frame("ID" = paste0("ID_", seq_along(tmp[["Eukaryota"]][[1]]))), tmp[["Eukaryota"]][, 1:7]) %>% rename("TAXONOMY" = "ID") %>% mutate("Kingdom" = str_replace(Kingdom, "Eukaryota", "Bacteria"))
    
    write.table(x = euk.count, file = paste0("output/count_tables/", i, "/", "01_euk_count_", i, ".txt"), quote = F, sep = "\t", row.names = F, col.names = T)
    write.table(x = euk.taxon, file = paste0("output/count_tables/", i, "/", "03_euk_tax_", i, ".txt"), quote = F, sep = "\t", row.names = F, col.names = T)
    
  }
  
  if (nrow(tmp[["Euk_filt"]]) > 0) {
    
    eukf.count = cbind(data.frame("ID" = paste0("ID_", seq_along(tmp[["Euk_filt"]][[1]]))), tmp[["Euk_filt"]][, 10:ncol(tmp[["Euk_filt"]])]) %>% rename("NAME" = "ID")
    eukf.taxon = cbind(data.frame("ID" = paste0("ID_", seq_along(tmp[["Euk_filt"]][[1]]))), tmp[["Euk_filt"]][, 1:7]) %>% rename("TAXONOMY" = "ID") %>% mutate("Kingdom" = str_replace(Kingdom, "Eukaryota", "Bacteria"))
    
    write.table(x = eukf.count, file = paste0("output/count_tables/", i, "/", "01_eukf_count_", i, ".txt"), quote = F, sep = "\t", row.names = F, col.names = T)
    write.table(x = eukf.taxon, file = paste0("output/count_tables/", i, "/", "03_eukf_tax_", i, ".txt"), quote = F, sep = "\t", row.names = F, col.names = T)
    
  }
  
}

save.image("output/count_tables/count_tables.RData")
