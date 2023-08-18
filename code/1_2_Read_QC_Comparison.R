#=========================================================================
# Comparing the number of reads per sample across all voyages in the study
# 
# Seb Rauschert
#=========================================================================

library(tidyverse)
library(cowplot)

# Loading all data
CoCosV10I_16S     <- read_table('data/read_qc/Sample_statistics_CoCosV10I_16S.txt')
RSV5_16S          <- read_table('data/read_qc/Sample_statistics_RS21AUG_16S.txt')
KWEST_16S         <- read_table('data/read_qc/Sample_statistics_KWEST_16S.txt')


# Filtering out the controls
CoCosV10I_16S %>%
  filter(!grepl('NTC|WC_|DI_|BC_|EB_', file)) %>%
  select(file, num_seqs, `Q30(%)`) %>% 
  mutate(data_source = 'Cocos Islands',
         Sequencer = 'NextSeq 2000') -> CoCosV10I_16S

RSV5_16S %>%
  filter(!grepl('NTC|WC_|DI_|BC_|EB_|MT1_', file)) %>%
  select(file, num_seqs, `Q30(%)`) %>% 
  mutate(data_source = 'Rowley Shoals',
         Sequencer = 'NextSeq 2000') -> RSV5_16S

KWEST_16S %>%
  filter(!grepl('Extbl|FC_', file))  %>%
  select(file, num_seqs, `Q30(%)`) %>% 
  mutate(data_source = 'NW WA',
         Sequencer = 'MiSeq') -> KWEST_16S

#============================================================================================
# Make figure for paper

# Number of reads
number <- rbind(CoCosV10I_16S,
      RSV5_16S,
      KWEST_16S) %>%
  ggplot(aes(x = factor(data_source, level=c("Cocos Islands", "Rowley Shoals", "NW WA")), y = num_seqs, fill = data_source)) +
  geom_jitter() +
  geom_boxplot(alpha = 0.8) +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ylab("Numbner of Sequences (R1 only)") +
  xlab('') +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#D55E00", "#00FFFF", "#FF00FF")) +
  scale_y_continuous(labels = scales::comma) 

plot_grid(map_plot, number, labels = NULL, label_size = 16)

median(CoCosV10I_16S$num_seqs)
median(RSV5_16S$num_seqs)
median(KWEST_16S$num_seqs)

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

