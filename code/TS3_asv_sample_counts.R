#==============================================================================
# Sample counts per ASV

## Read in files
CoCos_false  <- readRDS('data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_FALSE.rds')
CoCos_true  <- readRDS('data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_TRUE.rds')
CoCos_pseudo <- readRDS('data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_pseudo.rds')

RS_pooled     <- readRDS('data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_TRUE.rds')
RS_independent <- readRDS('data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_FALSE.rds')
RS_pseudo     <- readRDS('data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_pseudo.rds')

NW_pooled     <- readRDS('data/phyloseq_objects/NWWA_16S_phyloseq_nt_TRUE.rds')
NW_independent <- readRDS('data/phyloseq_objects/NWWA_16S_phyloseq_nt_FALSE.rds')
NW_pseudo     <- readRDS('data/phyloseq_objects/NWWA_16S_phyloseq_nt_pseudo.rds')

## Cocos Island dataset
# Count number of samples each ASV occurs in
# pooled
total_counts <- as_tibble(rowSums(CoCos_true@otu_table))
ASV_id <- rownames(CoCos_true@otu_table)
total_counts <- data.frame(ASV_id, total_counts)
colnames(total_counts) <- c('ASV_id', 'Total_counts')

# Sample_count <- as.data.frame(apply(NW_independent@otu_table,1,function(x) sum(x > 0)))
sample_count <- as_tibble(apply(CoCos_true@otu_table,1,function(x) sum(x > 0)))
ASV_id <- rownames(CoCos_true@otu_table)
sample_count <- data.frame(ASV_id, sample_count)
colnames(sample_count) <- c('ASV_id', 'Sample_count')

asv_seqs <- as.data.frame(tax_table(CoCos_true))
asv_seqs <- as_tibble(asv_seqs$ASV_sequence)
ASV_id <- rownames(CoCos_true@tax_table)
asv_seqs <- data.frame(ASV_id, asv_seqs)
colnames(asv_seqs) <- c('ASV_id', 'ASV_seqs')

Cocos_true_asvs <- merge(total_counts, sample_count, by = "ASV_id")
Cocos_true_asvs <- merge(Cocos_true_asvs, asv_seqs, by = "ASV_id")
# Cocos_true_asvs_1 <- Cocos_true_asvs %>%
#   filter(Sample_count == 1)
# head(Cocos_true_asvs_1)

# independent
total_counts <- as_tibble(rowSums(CoCos_false@otu_table))
ASV_id <- rownames(CoCos_false@otu_table)
total_counts <- data.frame(ASV_id, total_counts)
colnames(total_counts) <- c('ASV_id', 'Total_counts')

# Sample_count <- as.data.frame(apply(NW_independent@otu_table,1,function(x) sum(x > 0)))
sample_count <- as_tibble(apply(CoCos_false@otu_table,1,function(x) sum(x > 0)))
ASV_id <- rownames(CoCos_false@otu_table)
sample_count <- data.frame(ASV_id, sample_count)
colnames(sample_count) <- c('ASV_id', 'Sample_count')

asv_seqs <- as.data.frame(tax_table(CoCos_false))
asv_seqs <- as_tibble(asv_seqs$ASV_sequence)
ASV_id <- rownames(CoCos_false@tax_table)
asv_seqs <- data.frame(ASV_id, asv_seqs)
colnames(asv_seqs) <- c('ASV_id', 'ASV_seqs')


Cocos_false_asvs <- merge(total_counts, sample_count, by = "ASV_id")
Cocos_false_asvs <- merge(Cocos_false_asvs, asv_seqs, by = "ASV_id")
Cocos_false_asvs_1 <- Cocos_false_asvs %>%
  filter(Sample_count == 1)
length(Cocos_false_asvs_1$Sample_count) #97
# plyr::count(NW_false_asvs$Sample_count == 2)

length(intersect(Cocos_false_asvs_1$ASV_seqs, Cocos_true_asvs$ASV_seqs)) #57
asvs <- Reduce(intersect, list(Cocos_false_asvs_1$ASV_seqs, Cocos_true_asvs$ASV_seqs))
Cocos_true_asvs_inc = Cocos_true_asvs[Cocos_true_asvs$ASV_seqs %in% asvs,]
Cocos_false_asvs_inc = Cocos_false_asvs_1[Cocos_false_asvs_1$ASV_seqs %in% asvs,]

df <- as_tibble(merge(Cocos_false_asvs_inc, Cocos_true_asvs_inc, by = "ASV_seqs"))
ASV_inc_Cocos <- df %>%
  rename(Total_count_false = Total_counts.x) %>%
  rename(Total_count_true = Total_counts.y) %>%
  rename(Sample_count_false = Sample_count.x) %>%
  rename(Sample_count_true = Sample_count.y) %>%
  select(Total_count_false, Sample_count_false, Total_count_true, Sample_count_true, ASV_seqs)
head(ASV_inc_Cocos)   
length(ASV_inc_Cocos$ASV_seqs) # 57
write_csv2(ASV_inc_Cocos, file = "figures/ASV_inc_Cocos.csv")




## Rowley Shoals dataset
# Count number of samples each ASV occurs in
# pooled
total_counts <- as_tibble(rowSums(RS_pooled@otu_table))
ASV_id <- rownames(RS_pooled@otu_table)
total_counts <- data.frame(ASV_id, total_counts)
colnames(total_counts) <- c('ASV_id', 'Total_counts')

sample_count <- as_tibble(apply(RS_pooled@otu_table,1,function(x) sum(x > 0)))
ASV_id <- rownames(RS_pooled@otu_table)
sample_count <- data.frame(ASV_id, sample_count)
colnames(sample_count) <- c('ASV_id', 'Sample_count')

asv_seqs <- as.data.frame(tax_table(RS_pooled))
asv_seqs <- as_tibble(asv_seqs$ASV_sequence)
ASV_id <- rownames(RS_pooled@tax_table)
asv_seqs <- data.frame(ASV_id, asv_seqs)
colnames(asv_seqs) <- c('ASV_id', 'ASV_seqs')

RS_true_asvs <- merge(total_counts, sample_count, by = "ASV_id")
RS_true_asvs <- merge(RS_true_asvs, asv_seqs, by = "ASV_id")
# RS_true_asvs_1 <- RS_true_asvs %>%
#   filter(Sample_count == 1)
head(RS_true_asvs)

# independent
total_counts <- as_tibble(rowSums(RS_independent@otu_table))
ASV_id <- rownames(RS_independent@otu_table)
total_counts <- data.frame(ASV_id, total_counts)
colnames(total_counts) <- c('ASV_id', 'Total_counts')

sample_count <- as_tibble(apply(RS_independent@otu_table,1,function(x) sum(x > 0)))
ASV_id <- rownames(RS_independent@otu_table)
sample_count <- data.frame(ASV_id, sample_count)
colnames(sample_count) <- c('ASV_id', 'Sample_count')

asv_seqs <- as.data.frame(tax_table(RS_independent))
asv_seqs <- as_tibble(asv_seqs$ASV_sequence)
ASV_id <- rownames(RS_independent@tax_table)
asv_seqs <- data.frame(ASV_id, asv_seqs)
colnames(asv_seqs) <- c('ASV_id', 'ASV_seqs')


RS_false_asvs <- merge(total_counts, sample_count, by = "ASV_id")
RS_false_asvs <- merge(RS_false_asvs, asv_seqs, by = "ASV_id")
RS_false_asvs_1 <- RS_false_asvs %>%
  filter(Sample_count == 1)
length(RS_false_asvs_1$Sample_count) #388
# plyr::count(NW_false_asvs$Sample_count == 2)

length(intersect(RS_false_asvs_1$ASV_seqs, RS_true_asvs$ASV_seqs)) #159
asvs <- Reduce(intersect, list(RS_false_asvs_1$ASV_seqs, RS_true_asvs$ASV_seqs))
RS_true_asvs_inc = RS_true_asvs[RS_true_asvs$ASV_seqs %in% asvs,]
RS_false_asvs_inc = RS_false_asvs_1[RS_false_asvs_1$ASV_seqs %in% asvs,]

df <- as_tibble(merge(RS_false_asvs_inc, RS_true_asvs_inc, by = "ASV_seqs"))
ASV_inc_RS <- df %>%
  rename(Total_count_false = Total_counts.x) %>%
  rename(Total_count_true = Total_counts.y) %>%
  rename(Sample_count_false = Sample_count.x) %>%
  rename(Sample_count_true = Sample_count.y) %>%
  select(Total_count_false, Sample_count_false, Total_count_true, Sample_count_true, ASV_seqs)
head(ASV_inc_RS)   
length(ASV_inc_RS$ASV_seqs) #159
write_csv2(ASV_inc_RS, file = "figures/ASV_inc_RS.csv")

## Northwest dataset
# Count number of samples each ASV occurs in
# pooled
# total_counts <- as.data.frame(rowSums(NW_independent@otu_table))
total_counts <- as_tibble(rowSums(NW_pooled@otu_table))
ASV_id <- rownames(NW_pooled@otu_table)
total_counts <- data.frame(ASV_id, total_counts)
colnames(total_counts) <- c('ASV_id', 'Total_counts')

# Sample_count <- as.data.frame(apply(NW_independent@otu_table,1,function(x) sum(x > 0)))
sample_count <- as_tibble(apply(NW_pooled@otu_table,1,function(x) sum(x > 0)))
ASV_id <- rownames(NW_pooled@otu_table)
sample_count <- data.frame(ASV_id, sample_count)
colnames(sample_count) <- c('ASV_id', 'Sample_count')

asv_seqs <- as.data.frame(tax_table(NW_pooled))
asv_seqs <- as_tibble(asv_seqs$ASV_sequence)
ASV_id <- rownames(NW_pooled@tax_table)
asv_seqs <- data.frame(ASV_id, asv_seqs)
colnames(asv_seqs) <- c('ASV_id', 'ASV_seqs')

NW_true_asvs <- merge(total_counts, sample_count, by = "ASV_id")
NW_true_asvs <- merge(NW_true_asvs, asv_seqs, by = "ASV_id")
NW_true_asvs_1 <- NW_true_asvs %>%
  filter(Sample_count == 1)
head(NW_true_asvs)

# independent
# total_counts <- as.data.frame(rowSums(NW_independent@otu_table))
total_counts <- as_tibble(rowSums(NW_independent@otu_table))
ASV_id <- rownames(NW_independent@otu_table)
total_counts <- data.frame(ASV_id, total_counts)
colnames(total_counts) <- c('ASV_id', 'Total_counts')

# Sample_count <- as.data.frame(apply(NW_independent@otu_table,1,function(x) sum(x > 0)))
sample_count <- as_tibble(apply(NW_independent@otu_table,1,function(x) sum(x > 0)))
ASV_id <- rownames(NW_independent@otu_table)
sample_count <- data.frame(ASV_id, sample_count)
colnames(sample_count) <- c('ASV_id', 'Sample_count')

asv_seqs <- as.data.frame(tax_table(NW_independent))
asv_seqs <- as_tibble(asv_seqs$ASV_sequence)
ASV_id <- rownames(NW_independent@tax_table)
asv_seqs <- data.frame(ASV_id, asv_seqs)
colnames(asv_seqs) <- c('ASV_id', 'ASV_seqs')


NW_false_asvs <- merge(total_counts, sample_count, by = "ASV_id")
NW_false_asvs <- merge(NW_false_asvs, asv_seqs, by = "ASV_id")
NW_false_asvs_1 <- NW_false_asvs %>%
  filter(Sample_count == 1)
length(NW_false_asvs_1$Sample_count) #225
# plyr::count(NW_false_asvs$Sample_count == 2)


length(intersect(NW_false_asvs_1$ASV_seqs, NW_true_asvs$ASV_seqs)) #50
asvs <- Reduce(intersect, list(NW_false_asvs_1$ASV_seqs, NW_true_asvs$ASV_seqs))
NW_true_asvs_inc = NW_true_asvs[NW_true_asvs$ASV_seqs %in% asvs,]
NW_false_asvs_inc = NW_false_asvs_1[NW_false_asvs_1$ASV_seqs %in% asvs,]

df <- as_tibble(merge(NW_false_asvs_inc, NW_true_asvs_inc, by = "ASV_seqs"))
ASV_inc_NW <- df %>%
  rename(Total_count_false = Total_counts.x) %>%
  rename(Total_count_true = Total_counts.y) %>%
  rename(Sample_count_false = Sample_count.x) %>%
  rename(Sample_count_true = Sample_count.y) %>%
  select(Total_count_false, Sample_count_false, Total_count_true, Sample_count_true, ASV_seqs)
head(ASV_inc_NW)   
length(ASV_inc_NW$ASV_seqs)
write_csv2(ASV_inc_NW, file = "figures/ASV_inc_NW.csv")

