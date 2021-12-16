if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R") # current version is 0.99.20

install.packages('svglite')

require(qiime2R)
library(tidyverse)
library(data.table)
library(reshape2)

ExtractEvenness <- function(sampDepth=9500) {
  D <- read_qza(paste("C:\\Users\\jamsh\\evenness_vector_",sampDepth,".qza", sep=""))
  
  DD <- D$data %>% rownames_to_column("SampleID")

  vars = c(DD[DD$SampleID %like% "HH", ]$pielou_evenness, 
           DD[DD$SampleID %like% "SH", ]$pielou_evenness,
           DD[DD$SampleID %like% "SS", ]$pielou_evenness)
  
  groups = c(rep("HH", sum(DD$SampleID %like% "HH")),
             rep("SH", sum(DD$SampleID %like% "SH")),
             rep("SS", sum(DD$SampleID %like% "SS")))
  
  
  data=data.frame(vars, groups)
  return(data)
}

dat9000 <- ExtractEvenness(6000)
dat9000$depth <- as.character(6000)
dat9500 <- ExtractEvenness(9500)
dat9500$depth <- as.character(9500)
dat12000 <- ExtractEvenness(14000)
dat12000$depth <- as.character(14000)

dat <- do.call("rbind", list(dat9000,dat9500,dat12000))

ggplot(dat, aes(x=groups, y=vars, fill=depth)) + 
  geom_boxplot() + 
  scale_y_continuous(name="Pielou’s Evenness") + 
  scale_fill_discrete(name = "Sampling Depth") +
  theme(legend.position="top")

ggsave("../Evenness.svg", height=3, width=5)
