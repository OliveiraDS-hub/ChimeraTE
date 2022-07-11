library(ggplot2)
library(scales)
library(tidyr)
library(dplyr)

## TE-exonized
data = read.table("TE-exon-freq.tsv", sep = "\t", header = F)
data<-rename(data,c(TE_ID="V1", freq="V2"))

pdf("TE-exonized.pdf")
data %>% 
  count(TE_ID) %>%   
  mutate(perc = data$freq / sum(data$freq)) -> data2 

ggplot(data2, aes(x = reorder(TE_ID, -perc), y = perc))+
  geom_bar(stat = "identity", fill = "#FF3000", colour="black")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=12, color = "black"),
        axis.title.y = element_text(size=12, color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5))+
  ggtitle("TE-exonized transcripts")+
  xlab("") + ylab("Frequency in %")
dev.off()

## TE-initiated 5' UTR
data = read.table("TE-init-freq.tsv", sep = "\t", header = F)
data<-rename(data,c(TE_ID="V1", freq="V2"))

pdf("TE-initiated.pdf")
data %>%
  count(TE_ID) %>%
  mutate(perc = data$freq / sum(data$freq)) -> data2

ggplot(data2, aes(x = reorder(TE_ID, -perc), y = perc))+
  geom_bar(stat = "identity", fill = "#FF3000", colour="black")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=12, color = "black"),
        axis.title.y = element_text(size=12, color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5))+
  ggtitle("TE-initiated transcripts - 5' UTR")+
  xlab("") + ylab("Frequency in %")
dev.off()

## TE-initiated UPSTREAM

data = read.table("TE-init-freq-UPSTREAM.tsv", sep = "\t", header = F)
data<-rename(data,c(TE_ID="V1", freq="V2"))

pdf("TE-initiated-UPSTREAM.pdf")
data %>%
  count(TE_ID) %>%
  mutate(perc = data$freq / sum(data$freq)) -> data2

ggplot(data2, aes(x = reorder(TE_ID, -perc), y = perc))+
  geom_bar(stat = "identity", fill = "#FF3000", colour="black")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=12, color = "black"),
        axis.title.y = element_text(size=12, color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5))+
  ggtitle("TE-initiated transcripts - UPSTREAM")+
  xlab("") + ylab("Frequency in %")
dev.off()

## TE-terminated 3' UTR

data = read.table("TE-term-freq.tsv", sep = "\t", header = F)
data<-rename(data,c(TE_ID="V1", freq="V2"))

pdf("TE-terminated.pdf")
data %>%
  count(TE_ID) %>%
  mutate(perc = data$freq / sum(data$freq)) -> data2

ggplot(data2, aes(x = reorder(TE_ID, -perc), y = perc))+
  geom_bar(stat = "identity", fill = "#FF3000", colour="black")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=12, color = "black"),
        axis.title.y = element_text(size=12, color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5))+
  ggtitle("TE-terminated transcripts")+
  xlab("") + ylab("Frequency in %")
dev.off()

## TE-terminated 3' UTR

data = read.table("TE-term-freq-DOWN.tsv", sep = "\t", header = F)
data<-rename(data,c(TE_ID="V1", freq="V2"))

pdf("TE-terminated-DOWN.pdf")
data %>%
  count(TE_ID) %>%
  mutate(perc = data$freq / sum(data$freq)) -> data2

ggplot(data2, aes(x = reorder(TE_ID, -perc), y = perc))+
  geom_bar(stat = "identity", fill = "#FF3000", colour="black")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=12, color = "black"),
        axis.title.y = element_text(size=12, color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5))+
  ggtitle("TE-terminated transcripts - DOWNSTREAM")+
  xlab("") + ylab("Frequency in %")
dev.off()






