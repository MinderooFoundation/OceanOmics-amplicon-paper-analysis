#=========================================================================
# Comparing the number of reads per sample across all voyages in the study
# 
# Seb Rauschert
#=========================================================================

library(tidyverse)

# Loading all data
CoCosV10I_16S     <- read_table('data/read_qc/Sample_statistics_CocosV1_16S_filtered.txt')
CoCosV10I_MiFish  <- read_table('data/read_qc/Sample_statistics_CocosV1_MiFish_filtered.txt')
# CoCosV10II_16S    <- read_table('data/read_qc/Sample_statistics_CoCosV10II_16S.txt')
# CoCosV10II_MiFish <- read_table('data/read_qc/Sample_statistics_CoCosV10II_MiFish.txt')
RSV5_16S          <- read_table('data/read_qc/Sample_statistics_RSV5_16S_filtered.txt')
RSV5_MiFish       <- read_table('data/read_qc/Sample_statistics_RSV5_MiFish_filtered.txt')
KWEST_16S         <- read_table('data/read_qc/Sample_statistics_KWEST_16S.txt')


# Filtering out the controls
CoCosV10I_16S %>%
  filter(!grepl('NTC|WC_|DI_|BC_|EB_', file)) %>%
  select(file, num_seqs, `Q30(%)`) %>% 
  mutate(data_source = 'Cocos',
         sequencer = 'NextSeq 2000') -> CoCosV10I_16S

CoCosV10I_MiFish %>%
  filter(!grepl('NTC|WC_|DI_|BC_|EB_', file)) %>%
  select(file, num_seqs, `Q30(%)`) %>% 
  mutate(data_source = 'CoCos VI',
         sequencer = 'NextSeq 2000') -> CoCosV10I_MiFish

RSV5_16S %>%
  filter(!grepl('NTC|WC_|DI_|BC_|EB_|MT1_', file)) %>%
  select(file, num_seqs, `Q30(%)`) %>% 
  mutate(data_source = 'Rowley Shoals',
         sequencer = 'NextSeq 2000') -> RSV5_16S

RSV5_MiFish %>%
  filter(!grepl('NTC|WC_|DI_|BC_|EB_', file))  %>%
  select(file, num_seqs, `Q30(%)`) %>% 
  mutate(data_source = 'Rowley Shoals MiFish',
         sequencer = 'NextSeq 2000') -> RSV5_MiFish

KWEST_16S %>%
  filter(!grepl('Extbl|FC_', file))  %>%
  select(file, num_seqs, `Q30(%)`) %>% 
  mutate(data_source = 'NW WA',
         sequencer = 'MiSeq') -> KWEST_16S


#========================================================================================
# Number of reads

rbind(CoCosV10I_16S,
      CoCosV10I_MiFish,
      RSV5_16S,
      RSV5_MiFish,
      KWEST_16S) %>%
  ggplot(aes(x = data_source, y = num_seqs, fill = sequencer)) +
  geom_boxplot(alpha = 0.8) +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ylab("Numbner of Sequences (R1 only)") +
  xlab('') +
  scale_y_continuous(labels = scales::comma) 


# QC (Q 30) of reads

rbind(CoCosV10I_16S,
      CoCosV10I_MiFish,
      RSV5_16S,
      RSV5_MiFish,
      KWEST_16S) %>%
  ggplot(aes(x = data_source, y = `Q30(%)`, fill = sequencer)) +
  geom_boxplot(alpha = 0.8) +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ylab("Q30 (% of read)") +
  xlab('') +
  scale_y_continuous(labels = scales::comma) 


#============================================================================================
# Make figure for paper


# par(mfrow=c(1,2))

library(cowplot)

## 16S only
# Number of reads
number <- rbind(CoCosV10I_16S,
      RSV5_16S,
      KWEST_16S) %>%
  ggplot(aes(x = factor(data_source, level=c("Cocos", "Rowley Shoals", "NW WA")), y = num_seqs, fill = data_source)) +
  geom_jitter() +
  geom_boxplot(alpha = 0.8) +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ylab("Numbner of Sequences (R1 only)") +
  xlab('') +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#D55E00", "#00FFFF", "#FF00FF")) +
  scale_y_continuous(labels = scales::comma) 

number

# QC (Q 30) of reads
QC30 <- rbind(CoCosV10I_16S,
      RSV5_16S,
      KWEST_16S) %>%
  ggplot(aes(x = factor(data_source, level=c("Cocos", "Rowley Shoals", "NW WA")), y = `Q30(%)`, fill = data_source)) +
  geom_jitter() +
  geom_boxplot(alpha = 0.8) +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ylab("Q30 (% of read)") +
  xlab('') +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#D55E00", "#00FFFF", "#FF00FF")) +
  scale_y_continuous(labels = scales::comma) 

plot_grid(number, QC30, labels = NULL, label_size = 12)

# labels = c("B.1", "B.2")
#ggsave("figures/QC_plot_16S.pdf", QC_plot)
#save_plot("figures/QC_plot_16S.pdf", QC_plot, ncol = 2)



## MiFish only for supp materials

# Number of reads
number <- rbind(CoCosV10I_MiFish,
      RSV5_MiFish) %>%
  ggplot(aes(x = data_source, y = num_seqs, fill = sequencer)) +
  geom_boxplot(alpha = 0.8) +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ylab("Numbner of Sequences (R1 only)") +
  xlab('') +
  theme(legend.position = c(0.2, 0.8)) +
  scale_y_continuous(labels = scales::comma) 

# QC (Q 30) of reads
QC30 <- rbind(CoCosV10I_MiFish,
              RSV5_MiFish) %>%
  ggplot(aes(x = data_source, y = `Q30(%)`, fill = sequencer)) +
  geom_boxplot(alpha = 0.8) +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ylab("Q30 (% of read)") +
  xlab('') +
  scale_y_continuous(labels = scales::comma) 

QC_plot_MiFish <- plot_grid(number, QC30, labels = "AUTO")
