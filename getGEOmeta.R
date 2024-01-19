library(GEOquery)
library(data.table)
library(rhdf5)
library(ggplot2)
library(ggsci)
library(rhdf5)
library(DESeq2)
library(biomaRt)
library(AnnotationHub)
library(dplyr)
library(DGEobj.utils)
library(preprocessCore)
require(biomaRt)
library(factoextra)
library(FactoMineR)
library(ggsankey)
library(Seurat)
library(SCpubr)
options(scipen = 200)
# to identify pan cancer cell line RNAseq samples acquired from SRA database
setwd("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/pan_cell_line")
list=read.csv("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/pan_cell_line/info.csv")
list_clean=subset(list,list$TaxID=="9606")
list_clean=subset(list_clean,list_clean$LibraryStrategy=="RNA-Seq")
list_clean=subset(list_clean,grepl("GSM",list_clean$SampleName))
list_clean=subset(list_clean,list_clean$LibrarySource=="TRANSCRIPTOMIC")
index=duplicated(list_clean$SampleName)
list_clean=list_clean[!index,]
colnames(list_clean)[48]="Key"
write.csv(list_clean,"Cleaned_list.csv")
# define the GEOquert function to dowload meta data
download_metadata <- function(gsm_list) {
  # Create a character vector of GSM accession numbers
  gsm_ids <- as.character(gsm_list)
  
  # Download metadata for each GSM accession number
  metadata_list <- lapply(gsm_ids, function(gsm_id) {
    # Use tryCatch to handle any errors that may occur during downloading
    tryCatch(
      {
        # Download metadata for the specified GSM accession number
        raw=getGEO(gsm_id, destdir = ".", getGPL = FALSE,GSEMatrix = F)
        metadata=raw@header[c("geo_accession","library_selection","library_strategy","organism_ch1",
                              "platform_id","molecule_ch1","source_name_ch1","title")]
        metadata["series_id_1"]=raw@header["series_id"][[1]][1]
        metadata["series_id_2"]=raw@header["series_id"][[1]][2]
        metadata["Cell_type_1"]=raw@header["characteristics_ch1"][[1]][1]
        metadata["Cell_type_2"]=raw@header["characteristics_ch1"][[1]][2]
        metadata["Cell_type_3"]=raw@header["characteristics_ch1"][[1]][3]
        metadata["Cell_type_4"]=raw@header["characteristics_ch1"][[1]][4]
        metadata["Cell_type_5"]=raw@header["characteristics_ch1"][[1]][5]
        metadata["Cell_type_6"]=raw@header["characteristics_ch1"][[1]][6]
        metadata["growth_protocol_ch1"]=raw@header["growth_protocol_ch1"][[1]][1]
        metadata["description_1"]=raw@header["description"][[1]][1]
        metadata["description_2"]=raw@header["description"][[1]][2]
        metadata["description_3"]=raw@header["description"][[1]][3]
        # Return the metadata
        return(metadata)
      },
      error = function(e) {
        # Print an error message if there's an issue with downloading
        cat(paste("Error downloading metadata for GSM:", gsm_id, "\n"))
        return(NULL)
      }
    )
  })
  
  # Filter out NULL values (failed downloads)
  metadata_list <- metadata_list[!sapply(metadata_list, is.null)]
  
  # Return the list of metadata
  return(metadata_list)
}
# start to retrive data
metadata_list <- download_metadata(list_clean$SampleName)
combined_metadata <- rbindlist(lapply(metadata_list, as.data.table), fill = TRUE)
list_clean$ID=list_clean$SampleName
combined_metadata$ID=combined_metadata$geo_accession
meta_result=merge(list_clean,combined_metadata,by="ID")
# result turns out a lot of cell lines were not included because the names in CCLE collection are the not most used names
# start a second search from SRA to include more samples
# first add more PCa cell lines LNCaP, LNCaP is because the CCLE give a too long name
# remove the "NCI" letter from CCLE cell line names
ccle=read.csv("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/pan_cell_line/ccle_original.csv")
extra=ccle$StrippedCellLineName[grep("NCI",ccle$StrippedCellLineName)]
extra=gsub("NCI","",extra)
extra=c(extra,"A375","A673","NT2/D1","LNCAP","C42","C42B")
write.csv(extra,"extra_cellline_info.csv")
# all SRA records are downloaded using the bash file
# now retrive GEO meta data again for the extra samples
list_extra=read.csv("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/pan_cell_line/info_extra.csv")
list_extra=subset(list_extra,list_extra$TaxID=="9606")
list_extra=subset(list_extra,list_extra$LibraryStrategy=="RNA-Seq")
list_extra=subset(list_extra,grepl("GSM",list_extra$SampleName))
list_extra=subset(list_extra,list_extra$LibrarySource=="TRANSCRIPTOMIC")
index=duplicated(list_extra$SampleName)
list_extra=list_extra[!index,]
colnames(list_extra)[48]="Key"
write.csv(list_extra,"Cleaned_list_extra.csv")
#
metadata_list_extra <- download_metadata(list_extra$SampleName)
combined_metadata_extra <- rbindlist(lapply(metadata_list_extra, as.data.table), fill = TRUE)
list_extra$ID=list_extra$SampleName
combined_metadata_extra$ID=combined_metadata_extra$geo_accession
meta_result_extra=merge(list_extra,combined_metadata_extra,by="ID")
# now combine 2 reseach results
identical(colnames(meta_result),colnames(meta_result_extra))
final=rbind(meta_result,meta_result_extra)
# now screen out false positive samples
# Function to remove special characters "_", "-"
remove_special_chars <- function(x) gsub("[_-]", "", x)
final=subset(final,final$library_strategy=="RNA-Seq")
final$Cell_type_1_mod=remove_special_chars(final$Cell_type_1)
final$Cell_type_2_mod=remove_special_chars(final$Cell_type_2)
final$Cell_type_3_mod=remove_special_chars(final$Cell_type_3)
final$Cell_type_4_mod=remove_special_chars(final$Cell_type_4)
final$Cell_type_5_mod=remove_special_chars(final$Cell_type_5)
final$Cell_type_6_mod=remove_special_chars(final$Cell_type_6)
final$growth_protocol_ch1=remove_special_chars(final$Cell_type_1)
# this is to screen the key cell type that maches the identified samples' cell type
# Initialize a logical vector to store matches
matches <- logical(nrow(final))
# Loop through columns 60 to 75 and check for matches
cell_type_columns <- colnames(final)[c(56,57,60:65,67:75)]
# Create a new column 'celltye_check' by concatenating values from columns 60 to 75
final$celltye_check <- apply(final[, cell_type_columns], 1, function(row) paste0(row, collapse = " "))
matches <- logical(nrow(final))
for (row_index in 1:nrow(final)) {
  matches[row_index]=grepl(final$Key[row_index],final$celltye_check[row_index],ignore.case = T)
  
}
final$match=matches
# now start to clean
final_clean=subset(final,final$match)
rm(final)
final_clean=subset(final_clean,final_clean$Key!="clone")
write.csv(final_clean,"cleaned_allGEO_meta.csv")
# now start to retrive ARCHS4 data using the GEO accession number
destination_file = "/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/pan_cell_line/human_gene_v2.2.h5"
# Selected samples to be extracted
samp = final_clean$ID
# Retrieve information from compressed data
samples = h5read(destination_file, "meta/samples/geo_accession")
genes = h5read(destination_file, "meta/genes/symbol")
# Identify columns to be extracted
sample_locations = which(samples %in% samp)
# extract gene expression from compressed data
# too big to load together
# load partically
#load 1:10000

chunk_size <- 1000  # Adjust this value based on your available memory
expression_1 <- h5read(destination_file, name = "data/expression", index = list(sample_locations[1:chunk_size], 1:length(genes)))
expression_1 <- t(expression_1)  # Transpose the first chunk
cat("Reading chunk:", 1, "-", chunk_size, "\n")

# Loop to read and cbind more chunks
for (start_index in seq(chunk_size + 1, length(sample_locations), by = chunk_size)) {
  end_index <- min(start_index + chunk_size - 1, length(sample_locations))
  
  # Read data for the current chunk
  current_chunk <- h5read(destination_file, name = "data/expression", index = list(sample_locations[start_index:end_index], 1:length(genes)))
  
  # Transpose the current chunk before combining
  current_chunk <- t(current_chunk)
  
  # Combine the current chunk with the existing data
  expression_1 <- cbind(expression_1, current_chunk)
  
  # Print information about the current chunk
  cat("Reading chunk:", start_index, "-", end_index, "\n")
}

H5close()
gc()
count=as.data.frame(expression_1)
rm(expression_1)
colnames(count)=samples[sample_locations]
count$Symbol=genes
index=duplicated(count$Symbol)
count=count[!index,]
rownames(count)=count$Symbol
count=count[,-92504]
saveRDS(count,"count_matrix.rds")
# I indeitified some samples even have GAPDH gene count = 0, lets remove these low quality samples
samples_to_remove=names(count[, which(count["GAPDH",] < 10)])
count=count[,setdiff(colnames(count),samples_to_remove)]
samples_to_remove=names(count[, which(count["ACTB",] < 10)])
count=count[,setdiff(colnames(count),samples_to_remove)]
# also identify the metadata of collected samples
index=duplicated(final_clean$ID)
table(index)
final_clean=final_clean[!index,]
rownames(final_clean)=final_clean$ID
annotiation=final_clean[colnames(count),]
annotiation=annotiation[,-c(71:79)]
annotiation=annotiation[,-c(1,4:12,19,20,23:26,33:49)]
# filter out the cell lines with too less samples
annotiation$Key=gsub("NCI","",annotiation$Key)
freq=table(annotiation$Key)
annotiation=annotiation[annotiation$Key %in% names(freq[freq > 5]),]
# load CCLE annotiation file which also manually added the 293T and HEK293
ccle=read.csv("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/pan_cell_line/ccle_original.csv")
ccle$StrippedCellLineName=gsub("NCI","",ccle$StrippedCellLineName)
rownames(ccle)=ccle$StrippedCellLineName
cell_line_list=unique(annotiation$Key)
ccle=ccle[cell_line_list,]
ccle$Key=rownames(ccle)
annotiation=merge(annotiation,ccle,all.x = T)
write.csv(annotiation,"PCTA_annotiation.csv")
# generate statistics for the meta data
ggplot(annotiation, aes(x = "", fill =OncotreeLineage)) +
  geom_bar(width = 1, stat = "count", show.legend = FALSE, color = "black") +
  coord_polar("y") +
  theme_void() +
  scale_fill_ordinal()
ggsave("PCTA_cancer_type_pie_no_legend.tiff",width = 30,height = 30,units = "cm",dpi = 1200,compression="lzw")
ggplot(annotiation, aes(x = "", fill =OncotreeLineage)) +
  geom_bar(width = 1, stat = "count", show.legend = T, color = "black") +
  coord_polar("y") +
  theme_void() +
  scale_fill_ordinal()
ggsave("PCTA_cancer_type_pie_with_legend.tiff",width = 40,height = 30,units = "cm",dpi = 1200,compression="lzw")
# get statistics
length(table(annotiation$Key))
length(table(annotiation$OncotreeLineage))
# PCA analysis
# remove the low expression genes
count=count[,c(annotiation$ID)]
#
smallestGroupSize=5
keep=rowSums(count >= 10) >= smallestGroupSize
count=count[keep,]
saveRDS(count,"PCTA_count.rds")
# check if these is any sample with all genes with 0 count
#all.zero=apply(count, 2, function(x) all(x==0))
#table(all.zero)
# load to deseq2
coldata=annotiation[,c("Key","ID","OncotreeLineage","OncotreePrimaryDisease","OncotreeSubtype","OncotreeCode")]
dds=DESeqDataSetFromMatrix(countData = as.matrix(count),colData = coldata,design = ~Key)
# 
rm(count)
gc()
vsd=vst(dds, blind=T,nsub=500)


# now load the transcript level data
annotiation=read.csv("PCTA_annotiation.csv")
# now start to retrive ARCHS4 data using the GEO accession number
destination_file = "/Users/siyuan/Downloads/human_transcript_v2.2.h5"
# Selected samples to be extracted
samp = annotiation$ID
# Retrieve information from compressed data
samples = h5read(destination_file, "meta/samples/geo_accession")
tx=h5read(destination_file, "meta/transcripts")
genes = h5read(destination_file, "meta/transcripts/ensembl_id")
# Identify columns to be extracted
sample_locations = which(samples %in% samp)
# prepare for converting transcript counts to TPM
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org")
transcript_info <- getBM(attributes = c("ensembl_transcript_id", "transcript_length","hgnc_symbol"),filters = "ensembl_transcript_id", values = genes,mart = ensembl)
transcript_info_refine=transcript_info
index=duplicated(transcript_info_refine$ensembl_transcript_id)
transcript_info_refine=transcript_info_refine[!index,]
rownames(transcript_info_refine)=transcript_info_refine$ensembl_transcript_id
common_tx=intersect(transcript_info_refine$ensembl_transcript_id,genes)
transcript_info_refine=transcript_info_refine[common_tx,]
lengths=transcript_info_refine$transcript_length
geneids=transcript_info_refine$hgnc_symbol
#load 1:500
chunk_size <- 500  # Adjust this value based on your available memory
expression_1 <- h5read(destination_file, name = "data/expression", index = list(sample_locations[1:chunk_size], 1:length(genes)))
expression_1 <- t(expression_1)  # Transpose the first chunk
rownames(expression_1)=genes
expression_1=expression_1[common_tx,]
identical(rownames(expression_1),rownames(transcript_info_refine))
expression_1_tpm=convertCounts(countsMatrix = expression_1,geneLength = lengths,unit = "TPM")
expression_1_tpm=aggregate(expression_1_tpm, by = list(geneids), FUN = sum)
rm(expression_1)
cat("Reading chunk:", 1, "-", chunk_size, "\n")
# Loop to read and cbind more chunks
for (start_index in seq(chunk_size + 1, length(sample_locations), by = chunk_size)) {
  end_index <- min(start_index + chunk_size - 1, length(sample_locations))
  # Read data for the current chunk
  current_chunk <- h5read(destination_file, name = "data/expression", index = list(sample_locations[start_index:end_index], 1:length(genes)))
  # Transpose the current chunk before combining
  current_chunk <- t(current_chunk)
  rownames(current_chunk)=genes
  current_chunk=current_chunk[common_tx,]
  identical(rownames(current_chunk),rownames(transcript_info_refine))
  current_chunk_tpm=convertCounts(countsMatrix = current_chunk,geneLength = lengths,unit = "TPM")
  current_chunk_tpm=aggregate(current_chunk_tpm, by = list(geneids), FUN = sum)
  # Combine the current chunk with the existing data
  expression_1_tpm <- cbind(expression_1_tpm, current_chunk_tpm[,-1])
  # Print information about the current chunk
  cat("Reading chunk:", start_index, "-", end_index, "\n")
}
H5close()
gc()
#saveRDS(expression_1_tpm,"PCTA_transcript_tpm.rds")
expression_1_tpm=expression_1_tpm[-1,]
rownames(expression_1_tpm)=expression_1_tpm$Group.1
expression_1_tpm=expression_1_tpm[,-1]
colnames(expression_1_tpm)=samples[sample_locations]
# round the TPM matrix to save the space
# the space is enough, no need to rounded the data at this step
#rounded_expression_1_tpm <- as.data.table(lapply(expression_1_tpm, function(x) round(x, 0)))
#rounded_expression_1_tpm=as.data.table(lapply(rounded_expression_1_tpm,as.integer))
#rownames(rounded_expression_1_tpm)=rownames(expression_1_tpm)
annotiation=subset(annotiation,!(is.na(annotiation$OncotreeLineage)))


# normalize genes expression TPM values
expression_1_tpm=expression_1_tpm[,annotiation$ID]
gene=rownames(expression_1_tpm)
#matrix=as.matrix(expression_1_tpm)
#rownames(matrix)=rownames(expression_1_tpm)
saveRDS(expression_1_tpm,"PCTA_transcript_tpm.rds")
rm(expression_1_tpm)
gc()
# further clean the samples
# remove OncotreeLineage = NA samples
# normalize genes expression TPM values
data=readRDS("PCTA_transcript_tpm.rds")
mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "https://useast.ensembl.org")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup <- getBM(
  mart = mart,
  attributes = c(
    "hgnc_symbol",
    "entrezgene_id",
    "ensembl_gene_id",
    "gene_biotype"),
  filter = "hgnc_symbol",
  values = rownames(data),
  uniqueRows=TRUE)
index=duplicated(annotLookup$hgnc_symbol)
annotLookup=annotLookup[!index,]
gene_keep=c(annotLookup[which(annotLookup$gene_biotype=="protein_coding"),1],annotLookup[which(annotLookup$gene_biotype=="lncRNA"),1])
gene_keep=intersect(gene_keep,rownames(data))
data=data[gene_keep,]
saveRDS(data,"PCTA_tpm_protein_coding_only.rds")
# start to generate .csv files
gene_symbols=rownames(data)
output_dir <- "/Users/siyuan/Dropbox"
# Create the output directory if it doesn't exist
for (gene_symbol in gene_symbols) {
  # Extract the row corresponding to the current gene symbol
  gene_data <- data[gene_symbol, , drop = FALSE]  # Ensure it's a data.frame
  
  # Transpose the data
  gene_data_transposed <- t(gene_data)
  
  # Create a file path based on the gene symbol and the specified output directory
  file_path <- file.path(output_dir, paste0(gene_symbol, ".csv"))
  
  # Write the transposed data to a CSV file without row names
  write.csv(round(gene_data_transposed,digits = 2), file = file_path, row.names = FALSE)
  
  # Print the saved file path
  cat("Saved file:", file_path, "\n")
}




gene=rownames(data)
rname=rownames(data)
cname=colnames(data)
chunk_size <- 1000
for (i in seq(1, ncol(data), chunk_size)) {
  # Get the indices for the current chunk
  chunk_indices <- i:min(i + chunk_size - 1, ncol(data))
  
  # Round the selected columns in the current chunk
  data[, chunk_indices] <- lapply(data[, chunk_indices, drop = FALSE], function(x) round(x, 0))
  
  # Print a progress report
  cat("Processed columns:", min(chunk_indices), "-", max(chunk_indices), "\n")
}
matrix=as.matrix(data)
rm(data)
gc()
# round TPM values to save storage
rownames(matrix)=rname
colnames(matrix)=cname
#matrix=log2(matrix+1)
#matrix = normalize.quantiles(matrix)
#saveRDS(matrix,"PCTA_tpm_rounded.rds")
#gc()
# now start to generate H5 file
h5_file="PCTA_tpm.h5"
h5createFile(h5_file)
h5write(matrix,h5_file,"expression",level=6)
h5write(gene,h5_file,name = "gene")
h5ls('PCTA_tpm.h5')
# write sample meta files for PCTA app
annotiation=read.csv("PCTA_annotiation.csv")
#annotiation=annotiation[,c("Key","ID","OncotreeLineage","title)]
annotiation=subset(annotiation,!(is.na(annotiation$OncotreeLineage)))
annotiation$OncotreeLineage=gsub("Tool","HEK293/293T",annotiation$OncotreeLineage)
annotiation$OncotreeLineage=factor(annotiation$OncotreeLineage,levels = c("HEK293/293T",setdiff(unique(annotiation$OncotreeLineage),c("HEK293/293T","Other")),"Other"))
annotiation$Key=factor(annotiation$Key,levels = c("HEK293","293T",setdiff(unique(annotiation$Key),c("HEK293","293T"))))
pcta_meta=annotiation[,c("Key","ID","OncotreeLineage","title")]
write.csv(annotiation,"PCTA_meta.csv")
saveRDS(pcta_meta,"meta.rds")

# this is for testing
query = gene[1]
# Retrieve information from compressed data
genes = h5read("PCTA_tpm.h5", "gene")
# Identify columns to be extracted
gene_locations = which(genes %in% query)
gene_expression= h5read("PCTA_tpm.h5", "expression", index=list(gene_locations,1:ncol(data)))



# get DropBox tokens
# it is important to follow this https://github.com/karthik/rdrop2/issues/201
library(rdrop2)
#drop_auth()
#token <- drop_auth(new_user = T)
#saveRDS(token, file = "token.rds")


# try to use Seurat package to visulize the big count matrix 
data=readRDS("PCTA_tpm_protein_coding_only.rds")
rnames=rownames(data)
cnames=colnames(data)
rm(data)
count=readRDS("PCTA_count.rds")
count=count[rnames,cnames]
pcta_meta=readRDS("meta.rds")
rownames(pcta_meta)=pcta_meta$ID
identical(rownames(pcta_meta),colnames(count))
count=as.matrix(count)
count=na.omit(count)
# load meta data
annotiation=read.csv("PCTA_annotiation.csv")
annotiation=annotiation[,c("Key","ID","OncotreeLineage","series_id_1","Platform")]
annotiation=subset(annotiation,!(is.na(annotiation$OncotreeLineage)))
annotiation$OncotreeLineage=gsub("Tool","HEK293/293T",annotiation$OncotreeLineage)
annotiation$OncotreeLineage=factor(annotiation$OncotreeLineage,levels = c("HEK293/293T",setdiff(unique(annotiation$OncotreeLineage),c("HEK293/293T","Other")),"Other"))
annotiation$Key=factor(annotiation$Key,levels = c("HEK293","293T",setdiff(unique(annotiation$Key),c("HEK293","293T"))))
rownames(annotiation)=annotiation$ID
annotiation=na.omit(annotiation[colnames(count),])
sc=CreateSeuratObject(counts = count,assay = "RNA",meta.data = annotiation)
rm(count)
# start to visulize the metadata
do_AlluvialPlot(sample = sc,
                        first_group = "Key",
                        last_group = "OncotreeLineage",
                        fill.by = "OncotreeLineage",legend.position = "none",
                font.size = 10,repel = T,grid.color = NA)
ggsave("PCTA_sanky.tiff",width = 50,height = 100,dpi = 600,units = "cm",compression="lzw")
# now to generate PCA and UMAP analysis
gc()
sc <- NormalizeData(sc)
gc()
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)
gc()
all.genes <- rownames(sc)
sc <- ScaleData(sc, features = all.genes)
gc()
sc <- RunPCA(sc, features = VariableFeatures(object = sc))

# filter the marker genes
markers=read.csv("Cancer_type_markers.csv")
markers=subset(markers,(markers$pct.1>0.8)&(markers$avg_log2FC>4))
annotiation=read.csv("PCTA_annotiation.csv")
annotiation=annotiation[,c("Key","ID","OncotreeLineage","series_id_1","Platform","DepmapModelType")]
annotiation=subset(annotiation,!(is.na(annotiation$OncotreeLineage)))
annotiation$OncotreeLineage=gsub("Tool","HEK293/293T",annotiation$OncotreeLineage)
annotiation$OncotreeLineage=factor(annotiation$OncotreeLineage,levels = c("HEK293/293T",setdiff(unique(annotiation$OncotreeLineage),c("HEK293/293T","Other")),"Other"))
annotiation$Key=factor(annotiation$Key,levels = c("HEK293","293T",setdiff(unique(annotiation$Key),c("HEK293","293T"))))
rownames(annotiation)=annotiation$ID
tpm=readRDS("PCTA_tpm_protein_coding_only.rds")
genes=unique(markers$gene)
tpm=tpm[genes,]
tpm_subsets <- lapply(split(annotiation$ID, annotiation$DepmapModelType), function(ids) tpm[, ids, drop = FALSE])
tpm_subsets[[1]]=NULL # this is to remove the 293T samples
markers_refine=markers
markers_refine$ave=0
markers_refine$std=0
for (i in 1:nrow(markers)) {
  gene=markers[i,"gene"]
  type=markers[i,"cluster"]
  expression=log2(as.numeric(tpm_subsets[[type]][gene,])+1)
  ave=mean(expression)
  std=sd(expression)
  markers_refine[i,"ave"]=ave
  markers_refine[i,"std"]=std
}
write.csv(markers_refine,"cancer_type_markers_with_average.csv")
markers_refine_final=subset(markers_refine,markers_refine$ave>5)
markers_refine_final=subset(markers_refine_final,markers_refine_final$std<3)
write.csv(markers_refine_final,"cancer_type_markers_selected.csv")



# select marker genes for manuscript
markers=read.csv("cancer_type_markers_selected.csv")

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5
markers_top <- as.data.frame(top5)

# Group by cluster and concatenate gene values
result <- markers_top %>%
  group_by(cluster) %>%
  summarise(gene = paste(gene, collapse = ", "))
write.csv(result,"markers_for_table.csv")
fullname=read.csv("fullnames.csv")
colnames(result)=c("Cancer Type", "Marker genes")
fullname=fullname[,-1]
colnames(fullname)=c("Cancer Type","Full Name")
table=merge(result,fullname,by="Cancer Type")
table=table[,c(1,3,2)]
write.csv(table,"Table1_markers.csv",row.names = F)
# get PCTA meta info
annotiation=read.csv("PCTA_meta.csv")
cl=annotiation[!duplicated(annotiation$StrippedCellLineName),]
cl=cl[,c("StrippedCellLineName","DepmapModelType","OncotreeLineage")]

new_table <- cl %>%
  group_by(DepmapModelType, OncotreeLineage) %>%
  summarize(
    StrippedCellLineName = paste(StrippedCellLineName, collapse = ", ")
  ) %>%
  select(DepmapModelType, OncotreeLineage, StrippedCellLineName)
write.csv(new_table,"PCTA_cell_line_info.csv")




# Sankey plot
annotiation=read.csv("PCTA_annotiation.csv")
annotiation=annotiation[,c("Key","ID","OncotreeLineage","series_id_1","Platform")]
annotiation=subset(annotiation,!(is.na(annotiation$OncotreeLineage)))
annotiation$OncotreeLineage=gsub("Tool","HEK293/293T",annotiation$OncotreeLineage)
annotiation$OncotreeLineage=factor(annotiation$OncotreeLineage,levels = c("HEK293/293T",setdiff(unique(annotiation$OncotreeLineage),c("HEK293/293T","Other")),"Other"))
annotiation$Key=factor(annotiation$Key,levels = c("HEK293","293T",setdiff(unique(annotiation$Key),c("HEK293","293T"))))

sanky=annotiation[,c("Key","series_id_1","OncotreeLineage")] %>%
  make_long(Key, series_id_1, OncotreeLineage)
ggplot(sanky, aes(x = x, 
               next_x = next_x, 
               node = node, 
               next_node = next_node,
               fill = factor(node))) +
  geom_sankey() +
  theme_sankey(base_size = 16)+
  theme(legend.position = "none")
ggsave("PCTA_sanky.tiff",width = 10,height = 10,dpi = 1200,units = "cm",compression="lzw")


ggplot(sanky, aes(x = x, next_x = next_x, node = node, next_node = next_node,fill= factor(node),label = node)) +
  geom_alluvial(flow.alpha = .6) +
  scale_fill_viridis_d() +
  theme_alluvial(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) 










