
# Load required packages
require(XLConnect)
require(magrittr)
require(tidyr)
require(dplyr)
require(ggplot2)
library(stringr)

# Move to wanted directory

setwd("~/Mis proyectos con R/data_example")

# Import files

data.mock <- readWorksheetFromFile("Mock_sRNA.xlsx", sheet = "mock")
data.virus <- readWorksheetFromFile("virus_sRNA.xlsx", sheet = "virus", startRow = 2, endRow = 302, startCol = 2, endCol = 14)

head(data.mock); head(data.virus)

reads <- read.table("Reads_by_time_by_sample.txt")

reads %>%
  select(RNA, pos, strand, ends_with("Mock")) %>%
  gather(key = Time, value = Reads, -RNA, -pos, -strand) %>%
  mutate(Time = str_split_fixed(Time, "_", 3) [,2]) %>%
  mutate(Time = as.numeric(unlist(strsplit(Time, "dpi")))) -> reads.mock

reads %>%
  select(RNA, pos, strand, ends_with("TYLCV")) %>%
  gather(key = Time, value = Reads, -RNA, -pos, -strand) %>%
  mutate(Time = str_split_fixed(Time, "_", 3) [,2]) %>%
  mutate(Time = as.numeric(unlist(strsplit(Time, "dpi")))) -> reads.virus

data.mock <- mutate(data.mock, Sample = "Mock")
data.virus <- mutate(data.virus, Sample = "TYLCV")

data.mock.r <- inner_join(data.mock, reads.mock, by = c("miRNA_aligned_fragment" = "RNA"))
data.virus.r <- inner_join(data.virus, reads.virus, by = c("miRNA_aligned_fragment" = "RNA"))



data.full <- full_join(data.mock.r, data.virus.r)

data.full %>%
  group_by(Sample, Time) %>%
  summarize(Number_of_pairs = n()) %>%
  arrange(Sample, Time)




data.full %>%  # Welcome to piping!
  select(miRNA_Acc., Target_Acc., miRNA_end, Sample, pos, strand, Time, Reads) %>%
  filter(Reads >= 2000) %>%
  arrange(desc(Reads)) %>%
  group_by(Sample, Time) %>%
  slice(1:5) 

# Let's plot our data
data.full %>% 
  select(miRNA_Acc., Target_Acc., miRNA_end, Sample, pos, strand, Time, Reads) %>%
  filter(Reads >= 2000) %>%
  mutate(Reads = ifelse(test = strand == "-", yes = -Reads, no = Reads)) %>%
  mutate(miRNA_end = as.character(miRNA_end)) %>%
  ggplot(mapping = aes(x = pos, y = Reads, color = miRNA_end)) -> plot
plot <- plot + facet_grid(Time~Sample) + geom_point() + guides(fill = guide_legend())
ggsave(plot = plot, filename = "FilteredsRNA_byTime_bySample.png")


write.table(data.full, file = "Combined_data.txt")

excel.full <- loadWorkbook("Combined_data.xlsx", create = TRUE)
createSheet(excel.full, "Combined_data")
writeWorksheet(excel.full, 
               data.full, 
               sheet = "Combined_data", 
               startRow = 1, 
               startCol = 1)
saveWorkbook(excel.full)