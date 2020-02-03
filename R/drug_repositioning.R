#libraries to load
library(igraph)
library(readr)
library(dplyr)
options(stringsAsFactors = F)

dir <- "~/Área de Trabalho/GitHub_projects/PND_drug_repo/"
setwd(dir)

#import dis_FNG_FNDrugs file
workingdir <- "data/evolution_files/"
df <- list.files(path = paste0(workingdir),pattern = ".csv",full.names = T) %>%
  lapply(read_csv) %>% 
  bind_rows
df <- as.data.frame(df)
df <- df[!duplicated(df[,c(3,6)]),] #remove duplicated edges
df <- df[(df$Confidence > 50 | df$Confidence == -1) & df$Documents > 1,]

#gene-PND and drug-PND networks
df_drugs <- df[df$`Source type`=="DRUG",]
df_genes <- df[df$`Source type`=="GENE",]
gene_disease_edges <- df_genes %>%
  select(`Source Display Name`,`Target Display Name`) %>%
  rename(Source=`Source Display Name`,Target=`Target Display Name`)

#all drug-gene interactions
all_drug_gene <- read.csv("data/drugs_all_50percent_2documents.csv")

#gandal info
gandal <- read.csv("data/gandal_2018a.csv",dec = ",")
gandal_genes <- gandal %>%
  filter(Module.name!="CD0") %>%
  pull(external_gene_id)

gandal_genes_all <- unique(gandal$external_gene_id)

#get exclusive genes
nodes <- data.frame(Id=c(unique(df_genes$`Source Display Name`),
                         unique(df_genes$`Target Display Name`)),
                    Class=c(rep("GENE",1588),rep("CONDITION",9)))
edges <- df_genes %>%
  select(`Source Display Name`,`Target Display Name`) %>%
  rename(Source=`Source Display Name`,Target=`Target Display Name`)

graph <- graph_from_data_frame(d = edges,directed = F)
degree <- degree(graph)
nodes$degree <- degree(graph)

exclusive_genes <- nodes$Id[nodes$Class=="GENE" & nodes$degree==1]
exclusive_genes_coex <- exclusive_genes[exclusive_genes %in% gandal_genes]

#####GET ALL WDD SEARCHES WITH GENES AND MAKE DRUG_GENE DATAFRAME###############
genes_files <- paste0("/home/thomaz/first_neighbor_genes_all/",exclusive_genes_coex,"_50perc_2doc_50y.csv")

#drug-gene relations from raw source
df_drug_genes <- list.files(path = "/home/thomaz/first_neighbor_genes_all",pattern = ".csv",full.names = T) %>%
  .[. %in% genes_files] %>%
  lapply(read_csv) %>% 
  bind_rows
  
all_drug_gene <- df_drug_genes %>%
  filter(Confidence >= 50 & Documents >1 | Documents==-1) %>%
  filter(`Source type`=="DRUG" & `Target type`=="GENE") %>%
  select(`Source Display Name`,`Target Display Name`) %>%
  rename(Source=`Source Display Name`,Target=`Target Display Name`)


#drugs that affect exclusive genes coex
exclusive_genes_drugs <- all_drug_gene %>%
  filter(Target %in% exclusive_genes_coex)
  
length(unique(exclusive_genes_drugs$Source))

#remove known PND drugs
exclusive_genes_drugs <- exclusive_genes_drugs %>%
  filter(!Source %in% unique(df_drugs$`Source Display Name`))

length(unique(exclusive_genes_drugs$Source))

#remove drugs that affect more than one gene
# exclusive_genes_drugs <- all_drug_gene %>%
#   filter(Source %in% unique(exclusive_genes_drugs$Source)) %>%
#   filter(Target %in% unique(df_genes$`Source Display Name`))

#remove drugs that affect more than one gene
length(unique(exclusive_genes_drugs$Source))
length(unique(exclusive_genes_drugs$Target))

#make a nodes object to receive the degree of each nodes
drug_nodes <- data.frame(Id=c(unique(exclusive_genes_drugs$Source),
                              unique(exclusive_genes_drugs$Target)),
                         Class=c(rep("DRUG",907),rep("GENE",253))
                         )

#make igraph object of the drug-gene network
drug_graph <- graph_from_data_frame(d = exclusive_genes_drugs,directed = F)

#calculate drug-gene network degree
drug_nodes$degree <- degree(drug_graph)

#keep only drugs that affect one gene (among all genes in the original gene-PND network)
final_drugs <- drug_nodes$Id[drug_nodes$Class=="DRUG" & drug_nodes$degree==1]
length(unique(final_drugs))

final_genes <- unique(drug_nodes$Id[drug_nodes$Class=="GENE"])
length(unique(final_genes))

#create object to curate each drug-gene-disease relationship
drug_gene_disease_curate <- exclusive_genes_drugs %>%
  filter(Source %in% final_drugs) %>%
  select(Source,Target)

drug_gene_disease_curate <- merge(drug_gene_disease_curate,gene_disease_edges,
                                  by.x="Target",by.y="Source") %>%
  select(Source,Target,Target.y)

colnames(drug_gene_disease_curate) <- c("DRUG","GENE","DISEASE")

#save file outside of R to curate mannually
write.csv(drug_gene_disease_curate,file = "~/Área de Trabalho/Papers_Helder/PND_drug_repo/Rebutal/Files/drug_gene_disease_curate.csv",row.names = F)

exclusive_genes_coex_refereces <- df_genes %>%
  filter(`Source Display Name` %in% exclusive_genes_coex) %>%
  select(`Source Display Name`,`Target Display Name`,Confidence,`Document IDs`) %>%
  rename(Source=`Source Display Name`,Target=`Target Display Name`)

write.csv(exclusive_genes_coex_refereces,file="data/Table_S1_exclusive_genes_references.csv",quote = T,row.names = F)

final_gene_disease_references <- df_genes %>%
  filter(`Source Display Name` %in% final_genes) %>%
  select(`Source Display Name`,`Target Display Name`,Confidence,`Document IDs`) %>%
  rename(Source=`Source Display Name`,Target=`Target Display Name`)

#####GET ALL WDD SEARCHES WITH GENES AND MAKE DRUG_GENE DATAFRAME###############
genes_files <- paste0("/home/thomaz/first_neighbor_genes_all/",exclusive_genes_coex,"_50perc_2doc_50y.csv")

#drug-gene relations from raw source
df_drug_genes <- list.files(path = "~/first_neighbor_genes_all",pattern = ".csv",full.names = T) %>%
  .[. %in% genes_files] %>%
  lapply(read_csv) %>% 
  bind_rows

all_drug_gene <- df_drug_genes %>%
  filter(Confidence >= 50 & Documents >1 | Documents==-1) %>%
  filter(`Source type`=="DRUG" & `Target type`=="GENE") %>%
  select(`Source Display Name`,`Target Display Name`,Confidence,`Document IDs`) %>%
  rename(Source=`Source Display Name`,Target=`Target Display Name`)

################################################################################

final_drug_gene_references <- all_drug_gene %>%
  filter(Source %in% final_drugs & Target %in% final_genes)

final_drug_gene_disease_references <- merge(final_drug_gene_references,final_gene_disease_references,
                                            by.x="Target",by.y="Source")

final_drug_gene_disease_references <- final_drug_gene_disease_references %>%
  select(Source,Target,Target.y,`Document IDs.y`,`Document IDs.x`) %>%
  rename(drug=Source,gene=Target,PND=Target.y,
         gene_dis_ref=`Document IDs.y`,
         drug_gene_ref=`Document IDs.x`) %>%
  arrange(PND)

write.csv(final_drug_gene_disease_references,file="data/Table_S2_drug_gene_disease_references.csv",quote = T,row.names = F)
