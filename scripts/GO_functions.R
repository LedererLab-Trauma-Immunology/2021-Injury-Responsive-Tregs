#Brandon Hancock

f_getBP <- function(res,gene){
  BP <- res[which(res$query == gene),]$go.BP[[1]]
  return(BP)
}

f_getMF <- function(res,gene){
  MF <- res[which(res$query == gene),]$go.MF[[1]]
  return(MF)
}

f_getCC <- function(res,gene){
  CC <- res[which(res$query == gene),]$go.CC[[1]]
  return(CC)
}


#input mygene querymany, output datframe with terms and count number for each term
f_term_counts <- function(res){
  
  funcs <- res$go.BP
  term_counts <- data.frame(term = c(),count = c())
  
  for (i in 1:length(funcs)){
    for (j in 1:length(funcs[[i]]$term)){
      if (is.null(funcs[i][[1]])){
        next
      }
      t <- funcs[[i]]$term[j]
      if (! t %in% term_counts$term){
        term_counts <- rbind(term_counts,data.frame(term = t,count = 1))
      } else{
        term_counts[which(term_counts$term == t),]$count <- term_counts[which(term_counts$term == t),]$count + 1
      }
    }
  }
  sorted_term_counts <- term_counts[order(term_counts$count, decreasing = TRUE),]
  return(sorted_term_counts)
}

#input term_counts and search terms, output filtered term counts 
f_grep_term_counts <- function(term_counts,sterm){
  grep_term_counts <- data.frame(term = c(),count = c())
  
  for (i in 1:length(term_counts$term)){
    t <- term_counts$term[[i]]
    if (any(str_contains(t,sterm,ignore.case = TRUE))){
      grep_term_counts <- rbind(grep_term_counts,data.frame(term = t,count = term_counts$count[[i]]))
    } 
  }  
  return(grep_term_counts)
}

f_genes_by_term <- function(res,sterm){
  genes_by_term <- c()
  for (i in 1:length(res$query)){
    if (length(res$go.BP[[i]]$term) == 0){
      next
    }
    for (j in 1:length(res$go.BP[[i]]$term)){
      t <- res$go.BP[[i]]$term[[j]]
      if (any(t %in% sterm)){
        genes_by_term <- c(genes_by_term,res$query[i])
      }
    }
  }
  genes_by_term <- genes_by_term[!duplicated(genes_by_term)]
  return(genes_by_term)
}

f_gene_FC <- function(gene_fc_csv,gene_list){
  gene_fc_csv <- gene_fc_csv[which(gene_fc_csv$gene_symbol %in% gene_list),]
  gene_FC <- gene_fc_csv$logFC
  names(gene_FC) <- gene_fc_csv$gene_symbol
  return(gene_FC)
}

topDiffGenes <- function(allScore) {
  return(allScore < 0.0001)
}

topDiffGenes_2pointo <- function(allScore) {
  shigh_slow_csv <- read_csv("C:/Users/bhanc/Dropbox (Partners HealthCare)/Harvard CyTof/for Brandon/Sally/Sham low vs Sham high.csv")
  gene_fc_csv <- shigh_slow_csv
  
  gene_fc_csv_pval <- gene_fc_csv[gene_fc_csv$P.Value < 0.01,]
  gene_fc_csv_logFC <- gene_fc_csv_pval[gene_fc_csv_pval$logFC > 1 |  gene_fc_csv_pval$logFC < -1,]
  gene_list <- gene_fc_csv_logFC$P.Value
  names(gene_list) <- gene_fc_csv_logFC$gene_symbol
  return(gene_list)
}

f_Goid_fisher <- function(res,diffexp_csv){
  geneList <- diffexp_csv$P.Value
  names(geneList) <- res$entrezgene
  
  #will it work using query instead of ezgene?
  #names(geneList) <- res$query

  ezgene_to_GO <- list()
  for (i in 1:length(res$entrezgene)){
    ezgene_to_GO[[res$entrezgene[[i]]]] <- res$go.BP[[i]]$id
  }
  #will it work using query instead of ezgene?
  gene_to_GO <- list()
  for (i in 1:length(res$query)){
    gene_to_GO[[res$query[[i]]]] <- res$go.BP[[i]]$id
  }

  GO_obj <- new("topGOdata", ontology = "BP",allGenes = geneList,annot = annFUN.gene2GO, gene2GO = ezgene_to_GO, geneSelectionFun = topDiffGenes)
  #GO_obj <- new("topGOdata", ontology = "BP",allGenes = geneList,annot = annFUN.gene2GO, gene2GO = gene_to_GO, geneSelectionFun = topDiffGenes_2pointo)
  
  
  resultFisher <- runTest(GO_obj, algorithm = "classic", statistic = "fisher")
  GOid_fisher <- score(resultFisher)[score(resultFisher) < 0.01]
  return(GOid_fisher)
}


f_genes_term_map_fisher <- function(res,gene_list,GOid_fisher){
  genes_term_map_fisher <- list()
  #all_terms_fish <- c()
  for (i in 1:length(gene_list)){
    gene_terms_fish <- c()
    res_row <- res[which(res$query == gene_list[[i]]),]$go.BP[[1]]
    if (length(res_row$id != 0)){
      for (j in 1:length(res_row$id)){
        if (res_row$id[[j]] %in% names(GOid_fisher)){
          gene_terms_fish <- c(gene_terms_fish, res_row$term[[j]])
          #all_terms_fish <- c(all_terms_fish, res_row$term[[j]] )
       }
      }
    }
  genes_term_map_fisher[[gene_list[[i]]]] <- gene_terms_fish
  }
  return(genes_term_map_fisher)
}

f_upset_funcs <- function(){
  #upset_funcs <- c('Phosphorylation','Immune Response','Cell Cycle','Transcription','Signaling','Cell Movement','Translation','No Gene Ontology Data')
  upset_funcs <- c('Phosphorylation','Immune Response','Cell Cycle','Transcription','Signaling','Cell Movement','Translation')
  return(upset_funcs)
}

f_upset_dic <- function(upset_funcs){
  upset_dic <- hash()
  #upset_dic[[upset_funcs[[1]]]] <- c('cytokine production','interleukin','BMP signaling pathway','cytokine','chemokine','interferon')
  upset_dic[[upset_funcs[[1]]]] <- c('Phosphorylation','MAPK')
  upset_dic[[upset_funcs[[2]]]] <- c('immune response','inflammatory response','defense response', 't cell','T cell activation','leukocyte','lymphocyte',' t cell','mast cell','b cell','monocyte','interleukin','response to bacterium','cellular response to lipopolysaccharide','immune system','toll-like receptor')
  upset_dic[[upset_funcs[[3]]]] <- c('cell cycle','cell proliferation','cell differentiation','cell division','MAPK cascade','ERK1 and ERK2','chromosome segregation','Ras protein signal transduction','DNA replication','spindle organization','chromosome organization')
  upset_dic[[upset_funcs[[4]]]] <- c('transcription','gene expression')
  upset_dic[[upset_funcs[[5]]]] <- c('Intracellular Signal','signal transduction','signaling',upset_dic[['Cytokine']],upset_dic[['Phosphorylation']],upset_dic[['Cytokine']],'GTPase')
  upset_dic[[upset_funcs[[6]]]] <- c('taxis', 'chemotaxis', 'migration', 'motility')
  upset_dic[[upset_funcs[[7]]]] <- c('translation')
  #upset_dic[[upset_funcs[[8]]]] <- c('No Gene Ontology Data')
  
  #upset_dic[[upset_funcs[[3]]]] <- c('taxis', 'chemotaxis', 'migration', 'motility')
  #upset_dic[[upset_funcs[[4]]]] <- c('cell cycle','cell proliferation','cell differentiation','cell division','MAPK cascade','ERK1 and ERK2','chromosome segregation','Ras protein signal transduction','DNA replication','regulation of growth','chromosome organization')
  #upset_dic[[upset_funcs[[5]]]] <- c('immune response','inflammatory response','defense response', upset_dic[['T Cell']],'leukocyte','lymphocyte',' t cell','mast cell','b cell','monocyte','interleukin-2','response to bacterium','cellular response to lipopolysaccharide','immune system')
  #upset_dic[[upset_funcs[[6]]]] <- c('Intracellular Signal','signal transduction','signaling')
  #upset_dic[[upset_funcs[[7]]]] <- c('apoptotic','induced cell death')
  #upset_dic[[upset_funcs[[8]]]] <- c('metabolic process')
  #upset_dic[[upset_funcs[[9]]]] <- c('Phosphorylation')
  
  
  #
  #upset_dic[[upset_funcs[[4]]]] <- c('immune response','inflammatory response','defense response', upset_dic[['T Cell']],'leukocyte','lymphocyte',' t cell','mast cell','b cell','monocyte','interleukin-2','response to bacterium','cellular response to lipopolysaccharide','immune system')
  #upset_dic[[upset_funcs[[5]]]] <- c('Intracellular Signal','signal transduction','signaling','GTPase')

  
  #upset_dic[[upset_funcs[[9]]]] <- c('transcription','gene expression')
  return(upset_dic)
}

check_func <- function(func,terms,upset_dic){
  for (j in 1:length(terms)){
    t <- terms[j]
    if (any(str_contains(t,func,ignore.case = TRUE))){
      return(TRUE)
    } 
  }
  return(FALSE)
}

f_genes_func_map <- function(genes_term_map,upset_dic){
  genes_func_map <- list()
  upset_funcs <- keys(upset_dic)
  count = 0
  for (i in 1:length(genes_term_map)){
    gene_upset_funcs <- c()
    for (j in 1:length(upset_funcs)){
      sterm <- upset_dic[[upset_funcs[[j]]]]
      if (check_func(sterm,genes_term_map[[i]])){
        gene_upset_funcs <- c(gene_upset_funcs, upset_funcs[[j]])
      }
    }
    if (length(gene_upset_funcs) == 0){
      print(genes_term_map[i])
      count = count + 1
      print(count)
      #print('huh')
    }
    #if(length(gene_upset_funcs) == 0){
    #  gene_upset_funcs <- c("Uncategorised")
    #}
    genes_func_map[[names(genes_term_map)[i]]] <- gene_upset_funcs
  }
  return(genes_func_map)
}

f_gene_func_FC <- function(gene_fc_csv,genes_func_map,func){
  gene_list <- c()
  gene_list_FC <- c()
  for (i in 1:length(genes_func_map)){
    g <- names(genes_func_map)[[i]]
    if (func %in% genes_func_map[[i]]){
      gene_list <- c(gene_list,g)
      gene_list_FC <- c(gene_list_FC,gene_fc_csv[which(gene_fc_csv$gene_symbol == g),]$logFC)
    }
  }
  names(gene_list_FC) <- gene_list
  return(gene_list_FC)
}

f_listInput <- function(genes_func_map){
  listInput <- list()
  upset_funcs <- c(f_upset_funcs(),'Uncategorised')
  for (i in 1:length(upset_funcs)){
    listInput[upset_funcs[[i]]] <- c()
  } 
  for (i in 1:length(upset_funcs)){
    for (j in 1:length(genes_func_map)){
      if (upset_funcs[[i]] %in% genes_func_map[[j]]){
        listInput[[upset_funcs[[i]]]] <- c(listInput[[upset_funcs[[i]]]], names(genes_func_map)[j])
      }
    }
  }
  return(listInput)
}

genes_for_heatmap <- function(){
  #innate immune response, TH1 type cytokine, just immune response, adaptive or innate 
  c('IL23R','TYROBP','TXK','ITK','LEF1')
  
  
  #TH2 type cytokine, anti inflammatory 
  c('IL10','CD83','TIGIT','HAVCR2','PARP14','IL18')
  c('GATA3') # transcription factor the promotes IL4,IL5,IL13 from Th2 helper cells, need for development of th2 cells 
}

#terms for innate immune response, th1 type cytokine
c('positive regulation of interleukin-6 production','positive regulation of T-helper 1 type immune response','innate immune response','negative regulation of interleukin-10 production')

#terms for th2 type cytokine 
c('positive regulation of interleukin-10 production','interleukin-10 production')

#anti inflammatory cytokines: IL4, IL10,IL11, IL13 - focus on IL4, 10, 13


#IL4 terms
c('cellular response to interleukin-4','response to interleukin-4','positive regulation of interleukin-4-mediated signaling pathway')
c('positive regulation of interleukin-4 production','interleukin-4 production')
c('negative regulation of interleukin-4 production','negative regulation of interleukin-4-mediated signaling pathway')

#NFIL3 binds DNA and is required for transcription of il3, probably ignore it
#HAVCR2 causes supression of TH1 and Th17 responses 
#PARP14 anti-inflamatory response of macrophages 
#IL18 facilitates type 1 and 2 responses 

#-----------------------
#IL6,IL1B -> Type 1 -> Th17

f_gene_MF_list <- function(gene_list,res){
  MF_row = data.frame(Molecular_Function = c())
  for (j in 1:length(gene_list)){
    MF_terms <- f_getMF(res,gene_list[[j]])$term
    if (str_contains(MF_terms,'cytokine')){
      MF_row <- rbind(MF_row,data.frame(Molecular_Function = 'Cytokine'))
    } else if (str_contains(MF_terms,'transcription')){
      MF_row <- rbind(MF_row,data.frame(Molecular_Function = 'Transcription'))
    } else {
      MF_row <- rbind(MF_row,data.frame(Molecular_Function = 'other'))
    }
  }
  row.names(MF_row) <- gene_list
  return(MF_row)
} 

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

f_gene_list_heat <- function(){
  gene_list_heat <- c('IL23R','TYROBP','TXK','ITK','LEF1','IL10','CD83','TIGIT','HAVCR2','PARP14','IL18','GATA3','CCR5')
  gene_list_heat <- c(gene_list_heat,'GZMB', 'FGL2','S100A11','CASP1','CASP3','CCNB2','KLRG1','LGALS3','LTB4R1','NCAM1')
  gene_list_heat <- c(gene_list_heat,'CCR1','CCR2','CCR3','CCR4','CCR6','CCR7','CCR8','CCR9','CCR10','CXCR3', 'CXCR6')
  gene_list_heat <- c(gene_list_heat,'IL1RL1','IL1R2','IL3RA','IL7R','IL9R','IL10RA','IL12RB1','IL13RA1')
  gene_list_heat <- c(gene_list_heat,'FOXP3','ITGAE','CD44','ICOS')
  return(gene_list_heat)
}

f_gene_Th_type <- function(gene_list){
  gene_Th_type <- c()
  Th_dic <- hash()
  Th_dic[['IL23R']] <- 'Type 1'
  Th_dic[['TYROBP']] <- 'Type 1'
  Th_dic[['TXK']] <- 'Type 1'
  Th_dic[['ITK']] <- 'Type 1'
  Th_dic[['LEF1']] <- 'Type 1'
  Th_dic[['IL10']] <- 'Type 2'
  Th_dic[['CD83']] <- 'Type 2'
  Th_dic[['TIGIT']] <- 'Type 2'
  Th_dic[['HAVCR2']] <- 'Type 2'
  Th_dic[['PARP14']] <- 'Unclear'
  Th_dic[['IL18']] <- 'Type 1'
  Th_dic[['GATA3']] <- 'Type 2'
  Th_dic[['CCR5']] <- 'Unclear' #Th1 cells preferentially express 
  Th_dic[['GZMB']] <- 'Unclear'
  Th_dic[['FGL2']] <- 'Type 2' #inhibits immune response
  Th_dic[['S100A11']] <- 'Unclear' #cell motility, tubulin polymerization 
  Th_dic[['CASP1']] <- 'Type 1' #proinflammatory 
  Th_dic[['CASP3']] <- 'Unclear' 
  Th_dic[['CCNB2']] <- 'Unclear' #cell cycle
  Th_dic[['KLRG1']] <- 'Type 2'
  Th_dic[['LGALS3']] <- 'Unclear'
  Th_dic[['LTB4R1']] <- 'Unclear'
  Th_dic[["NCAM1"]] <- 'Unclear' #motility 
  Th_dic[["CCR1"]] <- 'Unclear'
  Th_dic[["CCR2"]] <- 'Unclear'
  Th_dic[["CCR3"]] <- 'Unclear' #Th2 cells preferentially express
  Th_dic[["CCR4"]] <- 'Unclear' #Th2 cells preferentially express
  Th_dic[["CCR6"]] <- 'Unclear' #did not look up CCR# after 4
  Th_dic[["CCR7"]] <- 'Unclear' #
  Th_dic[["CCR8"]] <- 'Unclear'
  Th_dic[["CCR9"]] <- 'Unclear'
  Th_dic[["CCR10"]] <- 'Unclear'
  Th_dic[["CXCR3"]] <- 'Type 1' #recruitment of inflammatory cells 
  Th_dic[["CXCR6"]] <- 'Unclear'
  Th_dic[["IL1RL1"]] <- 'Unclear'
  Th_dic[["IL1R2"]] <- 'Unclear'
  Th_dic[["IL3RA"]] <- 'Unclear'
  Th_dic[["IL7R"]] <- 'Unclear' #didn't check
  Th_dic[["IL9R"]] <- 'Unclear' 
  Th_dic[["IL10RA"]] <- 'Type 2'
  Th_dic[["IL12RB1"]] <- 'Unclear'
  Th_dic[['IL13RA1']] <- 'Unclear'
  Th_dic[['ITGAE']] <- 'Type 2'
  Th_dic[['FOXP3']] <- 'Type 2'
  Th_dic[['CD44']] <- 'Unclear'
  Th_dic[['ICOS']] <- 'Unclear'
  
  for (i in 1:length(gene_list)){
    gene_Th_type <- c(gene_Th_type, Th_dic[[gene_list[[i]]]])
  }
  return(gene_Th_type) 
}

#cytokines <- rownames(treg)[grep("Ccl|Cxc|Ifn|Il|Tnf|Cd40lg|Fasl|Cd70", rownames(treg))]

f_diff_cytokines <- function(diff_gene_list){
  diff_cytokines <- gene_list[grep("^CCL|^CXC|^IFN|^IL|^TNF|CD40LG|FASL|CD70|TGFB|^CD[[:digit:]]",diff_gene_list)]
  return(diff_cytokines)
}

f_treg_genes <- function(){
  treg_genes <- c("TGFB","IL10","ENTPD1","NT5E","LAG3",'TIGIT','CTLA4','ITGAE','KLRG1','ICOS','IL10RA','FGL2','HAVCR2','CD83')
  return(treg_genes)
}

f_treg_annot <- function(genes){
  treg_annot <- c()
  for (i in 1:length(genes)){
    if(genes[i] %in% f_treg_genes()){
      treg_annot <- c(treg_annot,'Treg Activity') 
    } else {
      treg_annot <- c(treg_annot,"Other")
    }
  }
  names(treg_annot) <- genes
  return(treg_annot)
}

f_transcription_genes <- function(res,genes){
  tgenes <- c()
  for (i in 1:length(genes)){
    if (length(grep('transcription factor',f_getMF(res,genes[i])$term)) > 0){
      tgenes <- c(tgenes,i)
    }
  }
  return(genes[tgenes])
}

f_genes_term_map <- function(res,genes){
  genes_term_map <- list()
  #count <- 0
  for (i in 1:length(genes)){
    #gene_terms <- c(f_getBP(res,genes[i])$term, f_getMF(res,genes[i])$term)
    gene_terms <- c(f_getBP(res,genes[i])$term)
    genes_term_map[[genes[i]]] <- gene_terms
    #if(length(gene_terms) == 0){
    #  genes_term_map[[genes[i]]] <- c("No Gene Ontology Data")
    #}
  }
  return(genes_term_map)
}

f_get_var <- function(vsd){
  var_list <- c()
  gene_ids <- row.names(vsd)
  for (i in 1:length(gene_ids)){
    gene_row <- as.vector(vsd[gene_ids[i],]) #as vector is a disapointment
    gene_vec <- c()
    for (j in 1:length(gene_row)){
      gene_vec <- c(gene_vec,gene_row[[j]])
    }
    var_list[gene_ids[i]] <- var(gene_vec)
  }
  return(var_list)
}
