#this script processes files obtained from Watson for Drug Discovery
#to perform drug repositioning analysis for Psychiatric and Neurological
#disorder

#####libraries to load####
library(igraph)
library(readr)
library(dplyr)
options(stringsAsFactors = F)

####set working directory####
dir <- "~/Área de Trabalho/GitHub_projects/PND_drug_repo/"
setwd(dir)

####Raw files - import dis_FNG_FNDrugs file####
workingdir <- "data/wdd_files/"
df <- list.files(path = paste0(workingdir),pattern = ".csv",full.names = T) %>%
  lapply(read_csv) %>% 
  bind_rows

df <- as.data.frame(df)
df <- df[!duplicated(df[,c(3,6)]),]
df <- df %>%
  filter((Confidence >=50 | Confidence == -1) & Documents > 1)

####Fig 1A: get gene-PND and drug-PND networks####
df_drugs <- df[df$`Source type`=="DRUG",]
wdd_all_drugs <- unique(df_drugs$`Source Display Name`)
save(wdd_all_drugs,file = "data/01_wdd_all_drugs.RData")

df_genes <- df[df$`Source type`=="GENE",]
wdd_all_genes <- unique(df_genes$`Source Display Name`)
save(wdd_all_genes,file = "data/02_wdd_all_genes.RData")

####Fig 1B: get exclusive coexpressed genes####
#make nodes and edges table to calculate degree
nodes <- data.frame(Id=c(unique(df_genes$`Source Display Name`),
                         unique(df_genes$`Target Display Name`)),
                    Class=c(rep("GENE",1588),rep("CONDITION",9)))
edges <- df_genes %>%
  select(`Source Display Name`,`Target Display Name`) %>%
  rename(Source=`Source Display Name`,Target=`Target Display Name`)

#make an igraph object and calculate degree of each node
graph <- graph_from_data_frame(d = edges,directed = F)
degree <- degree(graph)
nodes$degree <- degree(graph)

#save files to be viwed in gephi and analyzed with other scripts
#(Fig. 2 and Fig. S1)
write.csv(edges,file="data/03_gene_PNDs_edges.csv",quote = T,row.names = F)
write.csv(nodes,file="data/03_gene_PNDs_nodes.csv",quote = T,row.names = F)

#keep only genes with degree = 1
exclusive_genes <- nodes$Id[nodes$Class=="GENE" & nodes$degree==1]
exclusive_genes_edges <- edges %>%
  filter(Source %in% exclusive_genes)

save(exclusive_genes,file = "data/04_exclusive_genes.RData")

#Load Gandal coexpression information (genes coexpressed in PNDs)
gandal <- read.csv("data/gandal_2018a.csv",dec = ",")
gandal_genes <- gandal %>%
  filter(Module.name!="CD0") %>%
  pull(external_gene_id)

#get exclusive coexpressed genes
exclusive_genes_coex <- exclusive_genes[exclusive_genes %in% gandal_genes]
exclusive_genes_coex_edges <- edges %>%
  filter(Source %in% exclusive_genes_coex)
diseases <- unique(edges$Target)

exclusive_genes_coex_edges <- df_genes %>%
  filter(`Source Display Name` %in% exclusive_genes_coex) %>%
  select(`Source Display Name`,`Target Display Name`,Confidence,`Document IDs`) %>%
  rename(Source=`Source Display Name`,Target=`Target Display Name`)

#Save file with all exclusive coexpressed genes in PNDs (Table S1)
save(exclusive_genes_coex_edges,file="data/05_exclusive_genes_coex_edges.RData")
write.csv(exclusive_genes_coex_edges,file = "data/Table_S1_exclusive_coexpressed_genes.csv",row.names = F)

####Fig. 1C:Get WDD searches with exclusive coexpresed genes####
genes_files <- paste0("data/first_neighbor_genes_all/",exclusive_genes_coex,"_50perc_2doc_50y.csv")

#drug-gene relations from raw source
df_drug_genes <- list.files(path = "data/first_neighbor_genes_all",pattern = ".csv",full.names = T) %>%
  .[. %in% genes_files] %>%
  lapply(read_csv) %>% 
  bind_rows

all_drug_gene <- df_drug_genes %>%
  filter(Confidence >= 50 & Documents >1 | Documents==-1) %>%
  filter(`Source type`=="DRUG" & `Target type`=="GENE") %>%
  select(`Source Display Name`,`Target Display Name`,Confidence,`Document IDs`) %>%
  rename(Source=`Source Display Name`,Target=`Target Display Name`)

#drugs that affect exclusive genes coex
exclusive_genes_drugs <- all_drug_gene %>%
  filter(Target %in% exclusive_genes_coex)

length(unique(exclusive_genes_drugs$Source))

write.table(exclusive_genes_coex,file="data/05_exclusive_genes_coex.txt",quote = F,row.names = F,col.names = F)
save(exclusive_genes_drugs,file="data/06_exclusive_coex_genes_drugs.RData")  

####Fig. 1D: remove known PND drugs and drugs that affect > 1 target####
#remove known drugs
exclusive_genes_drugs <- exclusive_genes_drugs %>%
  filter(!Source %in% unique(df_drugs$`Source Display Name`))

length(unique(exclusive_genes_drugs$Source))

#remove drugs that affect more than one gene
ndrugs <- length(unique(exclusive_genes_drugs$Source))
ngenes <-length(unique(exclusive_genes_drugs$Target))

#make a nodes object to receive the degree of each nodes
drug_nodes <- data.frame(Id=c(unique(exclusive_genes_drugs$Source),
                              unique(exclusive_genes_drugs$Target)),
                         Class=c(rep("DRUG",ndrugs),rep("GENE",ngenes))
                         )

#make igraph object of the drug-gene network
drug_graph <- graph_from_data_frame(d = exclusive_genes_drugs[,-c(3,4)],directed = F)

#calculate drug-gene network degree
drug_nodes$degree <- degree(drug_graph)

#keep only drugs that affect one gene
final_drugs <- drug_nodes$Id[drug_nodes$Class=="DRUG" & drug_nodes$degree==1]
length(unique(final_drugs))

####Fig. 1E: create table for drug prioritization and selection####
#create table with all WDD gene-disease references
gene_disease_references <- df_genes %>%
  select(`Source Display Name`,`Target Display Name`,`Document IDs`) %>%
  rename(Gene=`Source Display Name`,PND=`Target Display Name`,gene_PND_ref=`Document IDs`)

#create table with all WDD drug-gene references
drug_gene_references <- exclusive_genes_drugs %>%
  filter(Source %in% final_drugs) %>%
  select(Source,Target,`Document IDs`) %>%
  rename(Drug=Source,Gene=Target,drug_gene_ref=`Document IDs`)

#create drug-gene-disease table with WDD references for prioritization
drug_gene_disease_reference <- merge(drug_gene_references,gene_disease_references,
                                  by.x="Gene",by.y="Gene") %>%
  select(Drug,Gene,PND,gene_PND_ref,drug_gene_ref) %>%
  arrange(PND)

final_genes <- unique(drug_gene_disease_curate$Gene)
length(unique(final_genes))

#save file to be used in the drug prioritization steps (Table S2)
write.csv(drug_gene_disease_reference,file="data/Table_S2_potential_novel_PND_drugs.csv",quote = T,row.names = F)
save(drug_gene_disease_reference,file="data/07_drug_gene_disease_final.RData")