###########################################################################################
# Author: Rodrigo Álvarez Pardo                                                           #
# Date: 09/10/21                                                                          #
# Description: Script to create a gtf file with the tRNAS fasta file from GtRNAdb, with   #
# the same format that the headers of the primary gtf file.                               #
###########################################################################################

########## LIBRARIES ##############

library(Biostrings)
library(tidyverse)
library(seqRFLP)
library(plyr)
library(naniar)
library(GenomicRanges)
library(rtracklayer)
library(remotes)
library(CAGEfightR)

#### VARIABLES AND DATA FRAME #####

setwd("G:/Mi unidad/Doctorado/Coding/Data")
getwd()
fastaFile <- readDNAStringSet("mm10-tRNAs.fa") # Insert fasta file in a variable
seq_name = names(fastaFile)                    # Insert sequences names in a variable
sequence = paste(fastaFile)                    # Insert sequences in a variable
df <- data.frame(seq_name, sequence)           # Insert names and sequences in a dataframe

###### Regular expressions ######

geneid<-str_match(seq_name,"t[A-Z]*\\-[A-Z]+[a-z]+\\-[A-Z]+\\-[0-9]+\\-[0-9]+|t[A-Z]*\\-[a-z][A-Z][a-z]+\\-[A-Z]+\\-[0-9]\\-[0-9]|t[A-Z]*\\-[A-Z][a-z][A-z]+\\-[A-Z]+\\-[0-9]\\-[0-9]") #Regex to extract the tRNA-id
geneid<-sapply(strsplit(geneid, '[, ]+'), function(x) toString(dQuote(x)))
dfid<-data.frame(geneid) #Create a dataframe with the gene id
#dfid[,1] <- dQuote(dfid[,1])#Adding double quotes to the geneid names
chr_number<-str_match(seq_name, "chr\\d+|chr\\w+") #Regex to extract the number of the chromosomes 
seq_length<-str_match(seq_name, "\\s\\d+\\s")
seq_length<-trimws(seq_length)#Regex to extract the lenght of the sequences
seq_length<-as.numeric(seq_length)
dflen<-data.frame(seq_length) #Create a dataframe with the sequence length
foward_strand<-str_match(seq_name, "\\+")#Regex to extract the symbol of foward strands
dffow<-data.frame(foward_strand)#Create a dataframe with the foward strands and a empty field in the reverse ones
dfstrand<-dffow %>% replace_na(list(foward_strand="-")) #Replacing NA values for minus sign
both_strands <- foward_strand %>% replace_na(list(foward_strand="-"))
both_strands<-as.character(both_strands)
#### Chromosomes numbers only #####

dfchr<-data.frame(chr_number) #Create dataframe from the list
chr_number<-gsub("chr","", chr_number) #Remove the "chr" word

####################################

line6<-':; gene_source' 
dfline6 <- data.frame(line6) 
n<-408
dfline6<-do.call("rbind", replicate(n, dfline6, simplify = FALSE))

line3<-'gene_id' 
dfline3 <- data.frame(line3) 
n<-408
dfline3<-do.call("rbind", replicate(n, dfline3, simplify = FALSE))

line5<-'; gene_version' 
dfline5 <- data.frame(line5) 
n<-408
dfline5<-do.call("rbind", replicate(n, dfline5, simplify = FALSE))


line9<-'; gene_name ' 
dfline9 <- data.frame(line9) 
n<-408
dfline9<-do.call("rbind", replicate(n, dfline9, simplify = FALSE))

line7<-'; gene_biotype' 
dfline7 <- data.frame(line7) 
n<-408
dfline7<-do.call("rbind", replicate(n, dfline7, simplify = FALSE))

if (file.exists("tRNA.txt")){ # Check the gtf file exits if it exits remove it
  file.remove("tRNA.txt")
}


gtf_primary<-rtracklayer::import("Mus_musculus.GRCm38.96_500.gtf")
gtf_primary_df=as.data.frame(gtf_primary)

dfstructure<-gtf_primary_df[1,]
dfstructure<-rbind(dfstructure, dfstructure[rep(1,407), ])

dfstructure$seqnames <- chr_number

dfstructure$start[dfstructure$start<3073254] <- 1

dfstructure$end <- seq_length

dfstructure$gene_id <- geneid
dfstructure$gene_id <- paste(dfline3$line3, dfstructure$gene_id)

dfstructure$gene_version<-dQuote(dfstructure$gene_version)
dfstructure$gene_version <- paste(dfline5$line5, dfstructure$gene_version)

dfstructure$gene_name <- geneid
dfstructure$gene_name<-paste(dfline9$line9, dfstructure$gene_name)

dfstructure$strand <- both_strands

dfstructure$gene_source<-dQuote(dfstructure$gene_source)
dfstructure$gene_source<-paste(dfline6$line6, dfstructure$gene_source)

dfstructure$gene_biotype<-dQuote(dfstructure$gene_biotype)
dfstructure$gene_biotype<-paste(dfline7$line7, dfstructure$gene_biotype)


file.create("tRNA.txt") # Create a txt file to transfer the database
write.table(dfstructure,file="tRNA.txt",row.names=FALSE,col.names = FALSE, na="", quote=F, sep ="\t")# the quote format allow to decide which column has doublequote
file.rename("tRNA.txt", "tRNA.gtf") # Changing the extension to gtf

