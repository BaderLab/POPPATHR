# Retrieve edge overlap data to reconstruct gene set clusters (i.e. pathways)
edge <- commandsPOST("edge get attribute")
edge_name <- as.data.frame(unlist(lapply(edge, '[', "name")))
edge_name <- strsplit(as.character(edge_name[,1]), " \\(Geneset_Overlap_Dataset 1) ")

# Initialize values
len <- length(edge_name)
ks <- NULL
rem <- NULL

for (i in 1:len) {
  int_k <- i+1
  if (int_k > len) {
    cat("*Merged pathway list generated\n")
    break
  }
  for (k in int_k:len) {
    if (length(which(edge_name[[i]] %in% edge_name[[k]])) >= 1) {
      edge_name[[i]] <- union(edge_name[[i]], edge_name[[k]])
    }
    ks <- c(ks, k)
  }
  for (j in ks) {
    if (length(which(edge_name[[i]] %in% edge_name[[j]])) >= 1) {
      rem <- c(rem, j)
    }
  }
  if (length(rem)) {
    edge_name <- edge_name[-rem]
  }
  len <- length(edge_name)
  rem <- NULL
  ks <- NULL
}

pathways <- edge_name

# Retrieve node data to get un-annotated gene sets (single pathways)
node <- commandsPOST("node get attribute")
node_name <- as.data.frame(unlist(lapply(node, '[', "name")))
node_name <- as.character(node_name[,1])

# Create crude pathway names by selecting top 3 words in pathway group
# NOTE potential future update of Cytoscape API to retrieve nodes in AutoAnnotate clusters
single_path <- node_name[-which(node_name %in% unlist(pathways))]
pathway_groups <- c(pathways, single_path)

for (x in seq_along(pathway_groups)) {
  path <- gsub("\\%.*", "", pathway_groups[[x]])
  path_split <- unlist(strsplit(path, " "))
  # Remove small words (to omit words such as "of", "the", "as", etc)
  if (length(path_split) > 3) {
    small_words <- which(nchar(path_split) <= 3)
    if (length(small_words)) {
      path_split <- path_split[-small_words]
    }
    top <- sort(table(path_split), decreasing=TRUE)[1:3]
    names(pathway_groups)[x] <- paste(names(top), collapse="_")
  }
  else {
    names(pathway_groups)[x] <- paste(path_split, collapse="_")
  }
}

# Write out rda file with pathway groupings
save(pathway_groups, file=sprintf("%s/pathway_groups.rda", outDir))
