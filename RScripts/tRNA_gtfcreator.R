###########################################################################################
# Author: Rodrigo �lvarez Pardo                                                           #
# Date: 09/10/21                                                                          #
# Description: Script to create a gtf file with the tRNAS fasta file from GtRNAdb, with   #
# the same format that the headers of the primary gtf file.                               #
###########################################################################################

########## LIBRARIES ##############

library(Biostrings)
library(tidyr, help, pos = 2, lib.loc = NULL)
library(seqRFLP)
library(plyr)
library(stringr, help, pos = 2, lib.loc = NULL)
library(vscDebugger, help, pos = 2, lib.loc = NULL)
library(filesstrings, help, pos = 2, lib.loc = NULL)

#### VARIABLES AND DATA FRAME #####

setwd("~/Coding/Coding/Data/refgenome")
getwd()
fastaFile <- readDNAStringSet("mm10-tRNAs.fa") # Insert fasta file in a variable
seq_name = names(fastaFile)                    # Insert sequences names in a variable
sequence = paste(fastaFile)                    # Insert sequences in a variable
df <- data.frame(seq_name, sequence)           # Insert names and sequences in a dataframe

###### Regular expressions ######

geneid<-str_match(seq_name,"t[A-Z]*\\-[A-Z]+[a-z]+\\-[A-Z]+\\-[0-9]+\\-[0-9]+|t[A-Z]*\\-[a-z][A-Z][a-z]+\\-[A-Z]+\\-[0-9]\\-[0-9]|t[A-Z]*\\-[A-Z][a-z][A-z]+\\-[A-Z]+\\-[0-9]\\-[0-9]") #Regex to extract the tRNA-id
#geneid<-sapply(strsplit(geneid, '[, ]+'), function(x) toString(dQuote(x)))
dfid<-data.frame(geneid) #Create a dataframe with the gene id
geneid[,1] <- dQuote(geneid[,1])#Adding double quotes to the geneid names
chr_number<-str_match(seq_name, "chr\\d+|chr\\w+") #Regex to extract the number of the chromosomes 
seq_length<-str_match(seq_name, "\\s\\d+\\s")
seq_length<-trimws(seq_length)#Regex to extract the lenght of the sequences
dflen<-data.frame(seq_length) #Create a dataframe with the sequence length
foward_strand<-str_match(seq_name, "\\+")#Regex to extract the symbol of foward strands
dffow<-data.frame(foward_strand)#Create a dataframe with the foward strands and a empty field in the reverse ones
dfstrand<-dffow %>% replace_na(list(foward_strand="-", stringsAsFactors=FALSE)) #Replacing NA values for minus sign

#### Chromosomes numbers only #####

dfchr<-data.frame(chr_number) #Create dataframe from the list
dfchr$chr_number<-gsub("chr","",as.character(dfchr$chr_number)) #Remove the "chr" word

### Create tRNA GTF File ###

if (file.exists("tRNA.txt")){ # Check the gtf file exits if it exits remove it
  file.remove("tRNA.txt")
}

### Create and replicate exact phrase which is in a gtf by columns ###

n<-408

line1<-'havana' 
dfline1 <- data.frame(line1) 
dfline1<-do.call("rbind", replicate(n, dfline1, simplify = FALSE))

line1.1<-'gene' 
dfline1.1 <- data.frame(line1.1) 
dfline1.1<-do.call("rbind", replicate(n, dfline1.1, simplify = FALSE))

line2<-"." 
dfline2 <- data.frame(line2) 
dfline2<-do.call("rbind", replicate(n, dfline2, simplify = FALSE))

line3<-'gene_id' 
dfline3 <- data.frame(line3) 
dfline3<-do.call("rbind", replicate(n, dfline3, simplify = FALSE))

line4<-'";'
dfline4 <- data.frame(line4) 
dfline4<-do.call("rbind", replicate(n, dfline4, simplify = FALSE))

linex<-';'
dflinex <- data.frame(linex) 
dflinex<-do.call("rbind", replicate(n, dflinex, simplify = FALSE))

liney<-'"'
dfliney <- data.frame(liney) 
dfliney<-do.call("rbind", replicate(n, dfliney, simplify = FALSE))

line5<-'gene_version' 
dfline5 <- data.frame(line5) 
dfline5<-do.call("rbind", replicate(n, dfline5, simplify = FALSE))

line5.1<-'"1"' 
dfline5.1 <- data.frame(line5.1) 
dfline5.1<-do.call("rbind", replicate(n, dfline5.1, simplify = FALSE))

line6<-'gene_source' 
dfline6 <- data.frame(line6) 
dfline6<-do.call("rbind", replicate(n, dfline6, simplify = FALSE))

line6.1<-'"havana"' 
dfline6.1 <- data.frame(line6.1) 
dfline6.1<-do.call("rbind", replicate(n, dfline6.1, simplify = FALSE))

line7<-'gene_biotype' 
dfline7 <- data.frame(line7) 
dfline7<-do.call("rbind", replicate(n, dfline7, simplify = FALSE))

line7.1<-'"TEC"' 
dfline7.1 <- data.frame(line7.1) 
dfline7.1<-do.call("rbind", replicate(n, dfline7.1, simplify = FALSE))

line8<-'1' 
dfline8 <- data.frame(line8) 
dfline8<-do.call("rbind", replicate(n, dfline8, simplify = FALSE))

line9<-'gene_name' 
dfline9 <- data.frame(line9) 
dfline9<-do.call("rbind", replicate(n, dfline9, simplify = FALSE))


# Creation of a data frame with the different parts of the phrase of a gtf in columns #

df_gtf<<-data.frame(column1=dfchr,
                   column2=dfline1,
                   column3=dfline1.1,
                   column4=dfline8,
                   column5=dflen,
                   column6=dfline2,
                   column7=dfstrand,
                   column8=dfline2,
                   column9=dfline3,
                   column10=dfliney,
                   column11=dfid,
                   column12=dfline4,
                   column13=dfline5,
                   column14=dfline5.1,
                   column15=dflinex,
                   columnn16=dfline9,
                   column17=dfliney,
                   column18=dfid,
                   column19=dfline4,
                   column20=dfline6,
                   column21=dfline6.1,
                   column22=dflinex,
                   column23=dfline7,
                   column24=dfline7.1,
                   column25=dflinex
                   )

df_gtf$liney<- gsub(" ", "", paste(df_gtf$liney,df_gtf$geneid))
df_gtf$geneid<-NULL

df_gtf$liney<- gsub(" ", "", paste(df_gtf$liney,df_gtf$line4))
df_gtf$line4<-NULL

df_gtf$line3<- paste(df_gtf$line3,df_gtf$liney)
df_gtf$liney<-NULL

df_gtf$liney.1<- gsub(" ", "", paste(df_gtf$liney,df_gtf$geneid.1))
df_gtf$geneid.1<-NULL

df_gtf$liney.1<- gsub(" ", "", paste(df_gtf$liney.1,df_gtf$line4.1))
df_gtf$line4.1<-NULL

df_gtf$line9<- paste(df_gtf$line9,df_gtf$liney.1)
df_gtf$liney.1<-NULL

df_gtf$line7.1<- gsub(" ", "", paste(df_gtf$line7.1,df_gtf$line4.4))
df_gtf$line4.4<-NULL

df_gtf$line5.1<- gsub(" ", "", paste(df_gtf$line5.1,df_gtf$linex))
df_gtf$linex<-NULL

df_gtf$line6.1<- gsub(" ", "", paste(df_gtf$line6.1,df_gtf$linex.1))
df_gtf$linex.1<-NULL

df_gtf$line7.1<- gsub(" ", "", paste(df_gtf$line7.1,df_gtf$linex.2))
df_gtf$linex.2<-NULL

df_gtf$column9 <- apply( df_gtf[,c(9:16)] , 1 , paste , collapse = " ") 
df_gtf2 <- df_gtf[,c(1:8,17)]



### Creation of the gtf file ###

file.create("tRNA.txt") # Create a txt file to transfer the database
write.table(df_gtf2,file="tRNA.txt",row.names=FALSE,col.names = FALSE, quote=F, sep ="\t")# the quote format allow to decide which column has doublequote
file.rename("tRNA.txt", "tRNA.gtf") # Changing the extension to gtf
file.move("tRNA.gtf", "~/Coding/Coding/Data/annotation/")

Line_gtf<- readLines(paste("~/Coding/Coding/Data/annotation/","Mus_musculus.GRCm38.96.gtf", sep = ""),
                     n=6)
print(Line_gtf)

Line_trna<- readLines(paste("~/Coding/Coding/Data/annotation/","tRNA.gtf", sep = ""),
                      n=1)
print(Line_trna)

