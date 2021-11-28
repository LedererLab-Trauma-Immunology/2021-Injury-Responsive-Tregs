
# https://www.biostars.org/p/383217/

add_clonotype <- function(tcr_folder, seurat_obj){
  tcr <- read.csv(paste(tcr_folder,"filtered_contig_annotations.csv", sep=""))
  
  # Remove the -1 at the end of each barcode.
  # Subsets so only the first line of each barcode is kept,
  # as each entry for given barcode will have same clonotype.
  
  #BH "this last step is for non aggregated data"
  #tcr$barcode <- gsub("-1", "", tcr$barcode)
  tcr <- tcr[!duplicated(tcr$barcode), ]
  
  # Only keep the barcode and clonotype columns. 
  # We'll get additional clonotype info from the clonotype table.
  tcr <- tcr[,c("barcode", "raw_clonotype_id")]
  names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"
  
  # Clonotype-centric info.
  clono <- read.csv(paste(tcr_folder,"clonotypes.csv", sep=""))
  
  # Slap the AA sequences onto our original table by clonotype_id.
  tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa")])
  
  # Reorder so barcodes are first column and set them as rownames.
  tcr <- tcr[, c(2,1,3)]
  rownames(tcr) <- tcr[,1]
  tcr[,1] <- NULL
  
  # Add to the Seurat object's metadata.
  clono_seurat <- AddMetaData(object=seurat_obj, metadata=tcr)
  
  #BH - group clones into expanded and non expanded 
  clone_counts <- f_clone_counts(clono_seurat)
  cell_state <- data.frame("Cell State" = c())
  clone_total <- data.frame("Clone Total" = c())
  pb <- progress_bar$new(total = length(clono_seurat$clonotype_id))
  for (i in 1:length(clono_seurat$clonotype_id)){
    clonotype_id <- clono_seurat$clonotype_id[[i]]
    cell <- names(clono_seurat$clonotype_id)[[i]]
    pb$tick()
    if (is.na(clonotype_id)){
      cell_state <- rbind(cell_state, data.frame("Cell State" = 'Non Expanded'))
      clone_total <- rbind(clone_total, data.frame("Clone Total" = 0))
      next
    }
    count <- clone_counts[which(clone_counts$clones == clonotype_id),]$count
    if (count > 1){
      cell_state <- rbind(cell_state, data.frame("Cell State" = 'Expanded'))
    } else {
      cell_state <- rbind(cell_state, data.frame("Cell State" = 'Non Expanded'))
    }
    clone_total <- rbind(clone_total, data.frame("Clone Total" = count))
  }
  rownames(cell_state) <- names(clono_seurat$clonotype_id)
  rownames(clone_total) <- names(clono_seurat$clonotype_id)
  clono_seurat <- AddMetaData(object = clono_seurat, metadata = cell_state)
  clono_seurat <- AddMetaData(object = clono_seurat, metadata = clone_total)
  
  return(clono_seurat)
}

#checking unique clonotypes 

f_clone_counts <- function(seurat_object){
  #clones <- seurat_object@meta.data$cdr3s_aa
  clones <- seurat_object@meta.data$clonotype_id
  clone_counts <- data.frame(clones = c(),count = c())
  for (i in 1:length(clones)){
    c <- clones[[i]]
    if (! clones[[i]] %in% clone_counts$clones){
      clone_counts <- rbind(clone_counts,data.frame(clones = clones[[i]],count = 1))
    } else {
      clone_counts[which(clone_counts$clone == c),]$count <- clone_counts[which(clone_counts$clone == c),]$count + 1
    }
  }
  clone_counts <- clone_counts[order(clone_counts$count, decreasing = TRUE),]
  return(clone_counts)
}

f_add_tsne <- function(seurat_object){
  tsne_fig <- data.frame(tsne_fig = c())
  pb <- progress_bar$new(total = length(seurat_object$seurat_clusters))
  for (i in 1:length(seurat_object$seurat_clusters)){
    pb$tick()
    ex <- seurat_object$Cell.State[[i]]
    if (ex == 'Expanded'){
      tsne_fig <- rbind(tsne_fig,data.frame(tsne_fig = ex))
    } else {
      cha <- paste()
      tsne_fig <- rbind(tsne_fig,data.frame(tsne_fig = str_pad(seurat_object$seurat_clusters[[i]],width = 2,side = 'left',pad = '0')))
    }
  }
  rownames(tsne_fig) <- names(seurat_object$seurat_clusters)
  seurat_object <- AddMetaData(object = seurat_object, metadata = tsne_fig)
  print("All Done :)")
  return(seurat_object)
}

f_add_act_pheno <- function(s_obj, act_clusters){
  act_pheno <- data.frame(act_pheno = c())
  pb <- progress_bar$new(total = length(s_obj$seurat_clusters))
  for (i in 1:length(s_obj$seurat_clusters)){
    pb$tick()
    clust <- s_obj$seurat_clusters[[i]]
    if (clust %in% act_clusters){
      act_pheno <- rbind(act_pheno, data.frame(act_pheno = 'Expanded Phenotype'))
    } else {
      act_pheno <- rbind(act_pheno, data.frame(act_pheno = 'Unexpanded Phenotype'))
    }
  }
  rownames(act_pheno) <- names(s_obj$seurat_clusters)
  s_obj <- AddMetaData(object = s_obj, metadata = act_pheno)
}

f_add_sample <- function(seurat_object,sample_names = c('Uninjured','7D after Injury')){
  sample <- data.frame(sample = c())
  pb <- progress_bar$new(total = length(seurat_object$orig.ident))
  for (i in 1:length(seurat_object$orig.ident)){
    pb$tick()
    cell <- names(seurat_object$orig.ident[i])
    if (grepl('-1',cell,fixed = TRUE)){
      sample <- rbind(sample, data.frame(sample = sample_names[1]))
    }
    else if (grepl('-2',cell,fixed = TRUE)){
      sample <- rbind(sample, data.frame(sample = sample_names[2]))
    }
  }
  row.names(sample) <- names(seurat_object$orig.ident)
  seurat_object <- AddMetaData(object = seurat_object, metadata = sample)
  print("All Done :)")
  return(seurat_object)
}

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}