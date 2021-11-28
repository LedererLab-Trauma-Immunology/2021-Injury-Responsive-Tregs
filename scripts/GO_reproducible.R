#my gene is not reproducible, so we try again
library(GO.db)
library(org.Mm.eg.db)

mapped_genes <- mappedkeys(org.Mm.egSYMBOL)
entrez_symbol_dic <- as.list(org.Mm.egSYMBOL[mapped_genes])

entrez_gene_list <- c()

gene_list <- sapply(gene_list, Cap)

for (i in 1:length(gene_list)){
  sym <- gene_list[i]
  index <- grep(sym,entrez_symbol_dic)
  if (length(index) == 1){
    entrez_gene_list <- c(entrez_gene_list,names(entrez_symbol_dic[index]))
  }
}


entrez_GOid_dic <- as.list(org.Mm.egGO[entrez_gene_list])

entrez_GOterm_dic <- list()

for (i in 1:length(entrez_GOid_dic)){
  entrez_GOterm_dic[[names(entrez_GOid_dic[i])]] <- Term(GOTERM)[names(entrez_GOid_dic[[i]])]
}

genes_func_map <- f_genes_func_map(entrez_GOterm_dic,f_upset_dic(f_upset_funcs()))
listInput <- f_listInput(genes_func_map)
upset(fromList(listInput), order.by = "freq",mainbar.y.label = 'Gene Count: Intersection',sets.x.label = 'Gene Count: Gene Ontology',nsets = 9,set_size.scale_max = 1000,set_size.show = TRUE,set_size.angles = 0,text.scale = 1.9)
