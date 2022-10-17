#require(data.table)
library("data.table")
library(bit64)
library(bit)
library(readr)
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(tidyverse)
require(stats); require(graphics)
library(tidyr)
library(plyr)



list_vcf_files <- list.files(path = "/Volumes/homes/am21t144/IAV_bioinformatics/P+P_all_outputs",
                            recursive = TRUE,
                            pattern = "*.csv",
                            full.names = TRUE)

#exaple of importing only 1 and adding the filename to its column
vcf1 <- read.delim("/Volumes/homes/am21t144/IAV_bioinformatics/P+P_all_outputs/AM_124_minQ20.recode.csv")
View(vcf1)
vcf1$ID <- 124

#Alternative with lapply, but I can't fix the df because I don't know how to enclose the loop
# Lapply <- lapply(myFiles, read.delim, header=TRUE)
# names(Lapply) <- myFiles
# for(i in myFiles)
#   Lapply[[i]]$Source = i
# do.call(rbind, Lapply)

setwd("/Volumes/homes/am21t144/IAV_bioinformatics/P+P_all_outputs")
myFiles <- list.files(path="/Volumes/homes/am21t144/IAV_bioinformatics/P+P_all_outputs/", pattern = "*.csv$")
names(myFiles) <- basename(myFiles)
all <- ldply(myFiles, read.csv)
view(all)

pp_table <-separate(data = all, col = variant_id.POS.REF.ALTRO.AO, into = c("variant_id", "pos", "ref", "alt","countRef","countAlt"), sep = "\t")
dim(pp_table)#3489
view(pp_table)
tail(pp_table)
#pp_table <- write.table(pp_table, file = "/Volumes/homes/am21t144/IAV_bioinformatics/pp_table.csv", dec = ",",  sep = " ")

#count the number of rows by .id in pp_table

library(plyr)

new_counts <- count(pp_table, ".id")

#create new column with sum values
pp_table$countRef<- as.numeric(pp_table$countRef)
pp_table$countAlt <- as.numeric(pp_table$countAlt)
pp_table$total_counts <- pp_table$countRef + pp_table$countAlt
dim(pp_table)#3489
min(pp_table$total_counts)#2
max(pp_table$total_counts)#151646
length(which(pp_table$total_counts==2))#only 1
length(which(pp_table$total_counts<10))

#add column with allele frequencies
pp_table$snp_freq <- pp_table$countAlt*100/pp_table$total_counts
View(pp_table)

min(pp_table$snp_freq)#0.03
length(which(pp_table$snp_freq<0.1)) #only 4

#drop the id long name and use only numbers

drop_string <-gsub("AM_(+[0-9]+)_.*.*", "\\1", pp_table$.id)
tail(drop_string)

pp_table$new_id <- drop_string
pp_table$new_id <- as.numeric(as.character(pp_table$new_id))

#try to merge dta frames based on new_id
library(readxl)
list_samples_seq <- read_excel("/Volumes/homes/am21t144/IAV_bioinformatics/vcf_IAV_analysis/list_samples_seq.xlsx")
library(dplyr)

# Adding column based on other columns. Here I add replicate number and drug value into separate columns:
new_list_samples <- list_samples_seq %>%
  mutate(replicate = case_when(
    endsWith(Sample_name, "1") ~ "1",
    endsWith(Sample_name, "2") ~ "2",
    endsWith(Sample_name, "3") ~ "3"
  ))%>%  
  mutate(drug = case_when(grepl("0uM", Sample_name) ~ "0",
                           grepl("4uM", Sample_name, ignore.case = TRUE) ~"4"))

list <-as.data.frame(new_list_samples)
names(list)[names(list) == 'Sample#'] <- 'new_id'

new_pp_table <- merge(pp_table,list, by  = "new_id") 
View(new_pp_table)
#toss out those with less than 0.1%

pptable_filter1 <-new_pp_table[which(new_pp_table$snp_freq >= 0.1),]

library(dplyr)
fil_counts <- pptable_filter1%>%count("new_id")

fil_counts<-pptable_filter1%>% 
  group_by(new_id,Passage,replicate,drug) %>%
  dplyr::summarize(count = n())%>%na.omit()

  View(fil_counts)
  
#View(pptable_filter1)

ggplot(fil_counts, aes(x=Passage, y=count))+
  geom_line()+geom_point()+ ylim(0,150)+
  scale_x_continuous(breaks = seq(from = 1, to = 10, by = 1))+
  labs(title = "PB1+PA Control line and mutants", 
       subtitle= "SNPs >= 0.1% ", 
       x= "Passage number", y= "SNPs' counts")+
  geom_rect(data = subset(fil_counts, replicate %in% "3" & drug %in% "4"),
                          fill = NA, colour = "red", xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf) +
  facet_wrap(~drug+replicate)
 
ggsave(file="/Volumes/homes/am21t144/IAV_bioinformatics/P+P_analysis/filt_P+P.pdf", width=12, height=8, dpi=300)


#group replicate 1
fil_rep1_4um <- fil_counts[c('2','17','20', '24', '28','6','9'),]
fil_rep1_4um$pass <- c(1,3,5,6,7,8,9)

ggplot(un_rep1_4um, aes(x=pass, y=freq))+
  geom_line()+geom_point()+ ylim(0,150)+
  scale_x_continuous(breaks = seq(from = 1, to = 10, by = 1))+
  labs(title = "PB1+PA Replicate 1", 
       subtitle= "SNPs >= 1", 
       x= "Passage number", y= "Mutations counts")

ggsave(file="/Volumes/homes/am21t144/IAV_bioinformatics/P+P_analysis/fil_rep1_4uM_P+P.pdf", width=4, height=4, dpi=300)



#PLOTS

unfiltered_counts <- count(pp_table, ".id")

unfil_control_passage <- unfiltered_counts[c('1','16','23', '27', '5','12','13'),]
unfil_control_passage$pass <- c(1,3,6,7,8,9,10)

unfilt_control <- ggplot(unfil_control_passage, aes(x=pass, y=freq))+
  geom_line()+geom_point()+ ylim(0,150)+
  scale_x_continuous(breaks = seq(from = 1, to = 10, by = 1))+
  labs(title = "PB1+PA_Control line", 
       subtitle= "Total depth ", 
       x= "Passage number", y= "Mutations counts")

ggsave(file="/Volumes/homes/am21t144/IAV_bioinformatics/P+P_analysis/unfilt_controlP+P.pdf", width=4, height=4, dpi=300)

#group replicate 1
un_rep1_4um <- unfiltered_counts[c('2','17','20', '24', '28','6','9'),]
un_rep1_4um$pass <- c(1,3,5,6,7,8,9)

ggplot(un_rep1_4um, aes(x=pass, y=freq))+
  geom_line()+geom_point()+ ylim(0,150)+
  scale_x_continuous(breaks = seq(from = 1, to = 10, by = 1))+
  labs(title = "PB1+PA Replicate 1", 
       subtitle= "total depth ", 
       x= "Passage number", y= "Mutations counts")

ggsave(file="/Volumes/homes/am21t144/IAV_bioinformatics/P+P_analysis/un_rep1_4uM_P+P.pdf", width=4, height=4, dpi=300)


#group replicate 2

un_rep2_4um <- unfiltered_counts[c('3','18','21', '25', '29','7','10'),]
un_rep2_4um$pass <- c(1,3,5,6,7,8,9)

ggplot(un_rep2_4um, aes(x=pass, y=freq))+
  geom_line()+geom_point()+ ylim(0,150)+
  scale_x_continuous(breaks = seq(from = 1, to = 10, by = 1))+
  labs(title = "PB1+PA Replicate 2", 
       subtitle= "total depth ", 
       x= "Passage number", y= "Mutations counts")

ggsave(file="/Volumes/homes/am21t144/IAV_bioinformatics/P+P_analysis/un_rep2_4uM_P+P.pdf", width=4, height=4, dpi=300)

#group replicate 3 unfiltered

un_rep3_4um <- unfiltered_counts[c('4','19','22', '26', '30','8','11','14'),]
un_rep3_4um$pass <- c(1,3,5,6,7,8,9,10)

ggplot(un_rep3_4um, aes(x=pass, y=freq))+
  geom_line()+geom_point()+ ylim(0,150)+
  scale_x_continuous(breaks = seq(from = 1, to = 10, by = 1))+
  labs(title = "PB1+PA Replicate 3", 
       subtitle= "total depth ", 
       x= "Passage number", y= "Mutations counts")+
  theme(plot.background = element_rect(fill = "green"))

ggsave(file="/Volumes/homes/am21t144/IAV_bioinformatics/P+P_analysis/un_rep3_4uM_P+P.pdf", width=4, height=4, dpi=300)


#drop rows with a less than 100 counts

filtered_pp_table <- subset(pp_table, total_counts>=100) 
dim(filtered_pp_table) #1973

#count again the number of snps per samples
library(plyr)
filtered_new_counts <- count(filtered_pp_table, ".id")
View(filtered_new_counts)

#group by control lines
control_passage <- filtered_new_counts[c('1','16','23', '27', '5','12','13'),]
control_passage$pass <- c(1,3,6,7,8,9,10)

control <- ggplot(control_passage, aes(x=pass, y=freq))+
  geom_line()+geom_point()+ ylim(0,105)+
  scale_x_continuous(breaks = seq(from = 1, to = 10, by = 1))+
  labs(title = "PB1+PA_Control line", 
       subtitle= "filtered out <100 reads total depth ", 
       x= "Passage number", y= "Mutations counts")

ggsave(file="/Volumes/homes/am21t144/IAV_bioinformatics/P+P_analysis/controlP+P.pdf", width=4, height=4, dpi=300)

#group replicate 1
rep1_4um <- filtered_new_counts[c('2','17','20', '24', '28','6','9'),]
rep1_4um$pass <- c(1,3,5,6,7,8,9)

ggplot(rep1_4um, aes(x=pass, y=freq))+
  geom_line()+geom_point()+ ylim(0,105)+
  scale_x_continuous(breaks = seq(from = 1, to = 10, by = 1))+
  labs(title = "PB1+PA Replicate 1", 
       subtitle= "filtered out <100 reads total depth ", 
       x= "Passage number", y= "Mutations counts")

ggsave(file="/Volumes/homes/am21t144/IAV_bioinformatics/P+P_analysis/rep1_4uM_P+P.pdf", width=4, height=4, dpi=300)


#group replicate 2

rep2_4um <- filtered_new_counts[c('3','18','21', '25', '29','7','10'),]
rep2_4um$pass <- c(1,3,5,6,7,8,9)

ggplot(rep2_4um, aes(x=pass, y=freq))+
  geom_line()+geom_point()+ ylim(0,105)+
  scale_x_continuous(breaks = seq(from = 1, to = 10, by = 1))+
  labs(title = "PB1+PA Replicate 2", 
       subtitle= "filtered out <100 reads total depth ", 
       x= "Passage number", y= "Mutations counts")

ggsave(file="/Volumes/homes/am21t144/IAV_bioinformatics/P+P_analysis/rep2_4uM_P+P.pdf", width=4, height=4, dpi=300)

#group replicate 3

rep3_4um <- filtered_new_counts[c('4','19','22', '26', '30','8','11','14'),]
rep3_4um$pass <- c(1,3,5,6,7,8,9,10)

ggplot(rep3_4um, aes(x=pass, y=freq))+
  geom_line()+geom_point()+ ylim(0,105)+
  scale_x_continuous(breaks = seq(from = 1, to = 10, by = 1))+
  labs(title = "PB1+PA Replicate 3", 
       subtitle= "filtered out <100 reads total depth ", 
       x= "Passage number", y= "Mutations counts")+
 theme(plot.background = element_rect(fill = "green"))

ggsave(file="/Volumes/homes/am21t144/IAV_bioinformatics/P+P_analysis/rep3_4uM_P+P.pdf", width=4, height=4, dpi=300)


