library(Giotto)
library(stringr)

args =  commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
    stop("Input data is is not provided.", call.=FALSE)
} else {
    inFile = args[1]
}

require(Seurat)
sobject <- readRDS(inFile)

fName = basename(inFile)
stopifnot(grepl(".result.rds",fName))

idBlock <- gsub(".result.rds", "", fName)
resId = str_split(idBlock,"_")[[1]][1]
sId = str_split(idBlock,"_")[[1]][2]

if (length(args) > 1) {
    annFile <- args[2]
    print("Loading annitation...")
    #annInfo <- readRDS(annFile)
    #annInfo <- annInfo[ (annInfo$orig.ident == sId) & (annInfo$sample == resId),]
    #stopifnot( nrow(annInfo) == ncol(sobject) )
    annInfo <- read.delim(annFile)
    stopifnot( rownames(annInfo) == rownames(sobject@meta.data))
    cId = 1
    print(paste("Target column:",colnames(annInfo)[cId] ))
} else {
    annInfo <- NULL
}
#

print(paste("Processing",resId,sId))

# matrix
exp = Seurat::GetAssayData(object = sobject, slot = "counts")

# PCA, UMAP (only fits for merged blocks)
pca_coord = as.matrix(Seurat::Embeddings(object = sobject, reduction = "pca"))
pca_load = as.matrix(Seurat::Loadings(object = sobject, reduction = "pca"))
pca_eig = sapply(Seurat::Stdev(sobject, reduction = "pca"), function(x) x ^ 2) 
colnames(pca_coord) = gsub(x = colnames(pca_coord), pattern = "PC_", replacement = "Dim.")  
colnames(pca_load) = gsub(x = colnames(pca_load), pattern = "PC_", replacement = "Dim.") 

pca = list(type = "cells",
           spat_unit = "cell",
           name = "pca", # issue here
           reduction_method = 'pca',
           coordinates = pca_coord,
           misc =list(eigenvalues = pca_eig, loadings = pca_load))

umap_coord = as.matrix(Seurat::Embeddings(object = sobject, reduction = "umap"))

colnames(umap_coord) = gsub(x = colnames(umap_coord), pattern = "UMAP_", replacement = "Dim.") 

umap = list(type = "cells",
            spat_unit = "cell",
            name = "umap",
            reduction_method = 'umap',
            coordinates = umap_coord,
            misc = NULL)

dim_reduc = list(pca,umap)

# spatial coords - only points required, cell borders are integrated in different manner

spat_coord = Seurat::GetTissueCoordinates(sobject)

#spat_coord = cbind(rownames(spat_coord), data.frame(spat_coord, row.names=NULL))
#colnames(spat_coord) = c("cell_ID", "sdimy", "sdimx")
#spat_loc = spat_coord
#spat_coord =  data.frame(spat_coord, row.names=NULL)
rownames(spat_coord) <- NULL
colnames(spat_coord) = c( "sdimy", "sdimx","cell_ID") 
summary(colnames(exp) == spat_coord$cell_ID)
# control order
spat_loc = spat_coord[,c("sdimx","sdimy","cell_ID")] # use cell centers


### integrate Seurat object into Giotto

# CNV clusters
temp_dir = paste0("cellposeRes/dfGiottoCnvRes/giotto_",resId,"_",sId,"/")
# CNV clusters with assignments
#temp_dir = paste0("cellposeRes/dfGiottoCnvABCRes/giotto_",resId,"_",sId,"/")



print(temp_dir)
# no clue what is the issue  : have to do this
dir.create(temp_dir, showWarnings = FALSE)


myinstructions = createGiottoInstructions(save_dir = temp_dir,
                                          save_plot = TRUE, 
                                          show_plot = TRUE)


gMb = createGiottoObject(exp,
                         spatial_locs = spat_loc,
                         dimension_reduction = dim_reduc,
                         instructions = myinstructions)

# normalization
gMb <- normalizeGiotto(gobject = gMb, verbose = T)
gMb <- addStatistics(gobject = gMb)
# Q: adjust normalization? is it helpful?
gMb <- adjustGiottoMatrix(gobject = gMb, 
                          expression_values = c('normalized'),
                          covariate_columns = c('nr_genes', 'total_expr'))

#screePlot(gMb, ncp = 20)
#plotPCA(gMb)
#plotUMAP(gobject =gMb)


# annotate cell types?


# Default : use initial clusters
if (is.null(annInfo)) {
    gMb <- addCellMetadata(gobject = gMb,new_metadata = sobject@meta.data$seurat_cluster, vector_name = "cluster")
}  else {
    print("Assign annotation...")
    gMb <- addCellMetadata(gobject = gMb,new_metadata = as.factor(annInfo[,cId]), vector_name = "cluster")
}



# visualization  per cluster

metadata <- pDataDT(gMb)
cTypes <- summary(metadata$cluster)

fName = paste0(temp_dir,sId,".spat_dim_clusters.pdf")
pdf(fName, width=12, height = 6)

p <- spatDimPlot(gMb,dim_point_size=2,
            cell_color = "cluster",  #select_cell_groups = '4',
            plot_alignment = "horizontal", spat_point_size = 2)


for (cType in names(cTypes)) {
    p <- spatDimPlot(gMb,dim_point_size=2,
            cell_color = "cluster",  select_cell_groups = cType,
            plot_alignment = "horizontal", spat_point_size = 2)
}

dev.off()



gMb <- createNearestNetwork(gobject = gMb, dimensions_to_use = 1:5, k = 5)



# initial UMAP
plotUMAP( gMb ,cell_color = "cluster", show_NN_network = T) 
# spatial view
spatPlot(gMb, return_plot = T, point_alpha = 0.8, cell_color = 'cluster', point_size=2)

gMb = createSpatialNetwork(gobject = gMb, minimum_k = 2, 
                           maximum_distance_delaunay = 400)
gMb = createSpatialNetwork(gobject = gMb, minimum_k = 2, 
                           method = 'kNN', k = 10)
showNetworks(gMb)

fName = paste0(temp_dir, sId,".delaunay_network.pdf")
pdf(fName,width = 6, height = 4)
spatPlot(gobject = gMb, show_network = T,
         network_color = 'blue', spatial_network_name = 'Delaunay_network',
         point_size = 1, cell_color = 'cluster')
dev.off()



## spatial co-expression

km_spatialgenes = binSpect(gMb) 

spatGenePlot(gMb, expression_values = 'scaled', 
             genes = km_spatialgenes[1:4]$genes,
             point_shape = 'border', point_border_stroke = 0.1,
             show_network = F, network_color = 'lightgrey', point_size = 1,
             cow_n_col = 2)


# 1. calculate spatial correlation scores 
#ext_spatial_genes = km_spatialgenes$genes
ext_spatial_genes = NULL # full
spat_cor_netw_DT = detectSpatialCorGenes(gMb,
                                         method = 'network', 
                                         spatial_network_name = 'Delaunay_network',
                                         subset_genes = ext_spatial_genes)

# 2. cluster correlation scores
# k= ?
spat_cor_netw_DT = clusterSpatialCorGenes(spat_cor_netw_DT, 
                                          name = 'spat_netw_clus', k = 6)
fName = paste0(temp_dir, sId,".spatial_heatmap.pdf")
pdf(fName,width=10,height = 15)
heatmSpatialCorGenes(gMb, spatCorObject = spat_cor_netw_DT, 
                     use_clus_name = 'spat_netw_clus',show_row_names = T)
dev.off()

netw_ranks = rankSpatialCorGroups(gMb, 
                                  spatCorObject = spat_cor_netw_DT, 
                                  use_clus_name = 'spat_netw_clus')
top_netw_spat_cluster = showSpatialCorGenes(spat_cor_netw_DT, 
                                            use_clus_name = 'spat_netw_clus',
                                            selected_clusters = 6, 
                                            show_top_genes = 1)

cluster_genes_DT = showSpatialCorGenes(spat_cor_netw_DT, 
                                       use_clus_name = 'spat_netw_clus',
                                       show_top_genes = 1)
cluster_genes = cluster_genes_DT$clus; names(cluster_genes) = cluster_genes_DT$gene_ID

gMb = createMetagenes(gMb, gene_clusters = cluster_genes, name = 'cluster_metagene')
fName = paste0(temp_dir,sId,".spatial_clusters.png")
png(fName,width=1000,height = 600)
spatCellPlot(gMb,show_plot = F,
             spat_enr_names = 'cluster_metagene',
             cell_annotation_values = netw_ranks$clusters,
             point_size = 1.5, cow_n_col = 3)
dev.off()



### cell neighborhood: 

set.seed(seed = 42)
cell_proximities = cellProximityEnrichment(gobject = gMb,
                                           cluster_column = 'cluster',
                                           spatial_network_name = 'Delaunay_network',
                                           adjust_method = 'fdr',
                                           number_of_simulations = 1000)


cellProximityBarplot(gobject =  gMb, 
                     CPscore = cell_proximities, 
                     min_orig_ints = 5, min_sim_ints = 5)

fName = paste0(temp_dir,sId,".cell_proximity_heatmap.png")
png(fName, width = 600, height = 400)

cellProximityHeatmap(gobject = gMb, CPscore = cell_proximities, 
                     order_cell_types = T, scale = T,
                     color_breaks = c(-1.5, 0, 1.5), 
                     color_names = c('blue', 'white', 'red'))
dev.off()


fName = paste0(temp_dir,idBlock,".cell_proximity_heatmap.pdf")
pdf(fName, width = 6, height = 4)

cellProximityHeatmap(gobject = gMb, CPscore = cell_proximities, 
                     order_cell_types = T, scale = T,
                     color_breaks = c(-1.5, 0, 1.5), 
                     color_names = c('blue', 'white', 'red'))
dev.off()

#cellProximityNetwork(gobject = gMb, CPscore = cell_proximities, 
#                     remove_self_edges = T, only_show_enrichment_edges = T)

cellProximityNetwork(gobject = gMb, CPscore = cell_proximities,
                     remove_self_edges = F, self_loop_strength = 0.3,
                     only_show_enrichment_edges = F,
                     rescale_edge_weights = T,
                     node_size = 8,
                     edge_weight_range_depletion = c(1, 2),
                     edge_weight_range_enrichment = c(2,5))


fName = paste0(temp_dir,idBlock,".cell_proximity_network.adjusted.pdf")
pdf(fName, width=8,height=6)
cellProximityNetwork(gobject = gMb, CPscore = cell_proximities,
                     remove_self_edges = F, self_loop_strength = 0.3,
                     only_show_enrichment_edges = F,
                     rescale_edge_weights = T,
                     node_size = 8,
                     edge_width_range = c(1, 1),
                     edge_weight_range_depletion = c(1, 2),
                     edge_weight_range_enrichment = c(2,5),
                     save_plot=F)
dev.off()


saveRDS(gMb,file=paste0(temp_dir, resId,"_",sId,".giotto_res.rds"))


