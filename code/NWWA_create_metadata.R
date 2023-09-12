# Create metadata for KWEST
# The number in the sample name is the site ID
library(tidyverse)


samples <- unlist(lapply(list.files('01-demultiplexed/16S'), function(x) str_remove(x, ".fastq")))
samples <- c(samples[samples != "Controls"], unlist(lapply(list.files('01-demultiplexed/16S/Controls'), function(x) str_remove(x, ".fastq"))))
samples <- unlist(lapply(samples, function(x) (str_split(x, '~')[[1]][1])))
sites <- unlist(lapply(samples, function(x) parse_number(str_split(x, '~')[[1]][1])))

tibble(`Sample ID` = samples,
           assay = '16S',
           site = sites) %>%
  write_csv('06-report/KWest_metadata.csv')

parse_number('10a')
