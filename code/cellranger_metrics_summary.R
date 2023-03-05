library(tidyverse)
metrics_summary_SRR17655997 <- read_csv("/data/PRJNA798401/cellranger/SRR17655997/outs/metrics_summary.csv") %>% mutate(Run = "SRR17655997")
metrics_summary_SRR17655998 <- read_csv("/data/PRJNA798401/cellranger/SRR17655998/outs/metrics_summary.csv") %>% mutate(Run = "SRR17655998")
metrics_summary_SRR17655999 <- read_csv("/data/PRJNA798401/cellranger/SRR17655999/outs/metrics_summary.csv") %>% mutate(Run = "SRR17655999")
metrics_summary_SRR17656000 <- read_csv("/data/PRJNA798401/cellranger/SRR17656000/outs/metrics_summary.csv") %>% mutate(Run = "SRR17656000")
metrics_summary_SRR17656001 <- read_csv("/data/PRJNA798401/cellranger/SRR17656001/outs/metrics_summary.csv") %>% mutate(Run = "SRR17656001")
metrics_summary_SRR17656002 <- read_csv("/data/PRJNA798401/cellranger/SRR17656002/outs/metrics_summary.csv") %>% mutate(Run = "SRR17656002")
metrics_summary_SRR17656003 <- read_csv("/data/PRJNA798401/cellranger/SRR17656003/outs/metrics_summary.csv") %>% mutate(Run = "SRR17656003")
metrics_summary_SRR17656004 <- read_csv("/data/PRJNA798401/cellranger/SRR17656004/outs/metrics_summary.csv") %>% mutate(Run = "SRR17656004")
metrics_summary_SRR17656005 <- read_csv("/data/PRJNA798401/cellranger/SRR17656005/outs/metrics_summary.csv") %>% mutate(Run = "SRR17656005")
metrics_summary_SRR17656006 <- read_csv("/data/PRJNA798401/cellranger/SRR17656006/outs/metrics_summary.csv") %>% mutate(Run = "SRR17656006")
metrics_summary_SRR17656007 <- read_csv("/data/PRJNA798401/cellranger/SRR17656007/outs/metrics_summary.csv") %>% mutate(Run = "SRR17656007")
metrics_summary_SRR17656008 <- read_csv("/data/PRJNA798401/cellranger/SRR17656008/outs/metrics_summary.csv") %>% mutate(Run = "SRR17656008")
metrics_summary_SRR17656009 <- read_csv("/data/PRJNA798401/cellranger/SRR17656009/outs/metrics_summary.csv") %>% mutate(Run = "SRR17656009")
metrics_summary_SRR17656010 <- read_csv("/data/PRJNA798401/cellranger/SRR17656010/outs/metrics_summary.csv") %>% mutate(Run = "SRR17656010")
metrics_summary_SRR17656011 <- read_csv("/data/PRJNA798401/cellranger/SRR17656011/outs/metrics_summary.csv") %>% mutate(Run = "SRR17656011")
metrics_summary_SRR17656012 <- read_csv("/data/PRJNA798401/cellranger/SRR17656012/outs/metrics_summary.csv") %>% mutate(Run = "SRR17656012")
metrics_summary <-
  bind_rows(
    metrics_summary_SRR17655997,
    metrics_summary_SRR17655998,
    metrics_summary_SRR17655999,
    metrics_summary_SRR17656000,
    metrics_summary_SRR17656001,
    metrics_summary_SRR17656002,
    metrics_summary_SRR17656003,
    metrics_summary_SRR17656004,
    metrics_summary_SRR17656005,
    metrics_summary_SRR17656006,
    metrics_summary_SRR17656007,
    metrics_summary_SRR17656008,
    metrics_summary_SRR17656009,
    metrics_summary_SRR17656010,
    metrics_summary_SRR17656011,
    metrics_summary_SRR17656012)

metrics_summary |>
  select("Estimated Number of Cells", "Run")

write_tsv(metrics_summary, here("metrics_summary.tsv"))

