###########################################################################################
# Author: Rodrigo Álvarez Pardo                                                           #
# Date: 07/11/21                                                                          #
# Description: Script to create a gtf file with the tRNAS fasta file from GtRNAdb, with   #
# the same format that the headers of the primary gtf file.                               #
###########################################################################################

############# LIBRARIES ##############

library(Biostrings)
library(stringr)
library(filesstrings)

######################################

setwd("~/Coding/Coding/Data/refgenome") # Set directory 
fasta_file <- readDNAStringSet("mm10-tRNAs_17.fa") # Import fastafile and keep in a variable
seq_names = names(fasta_file) # Keep the names of the sequences in a variable
df_fastanames <-data.frame(seq_names) # Keep the sequences names in a dataframe

########### Check GTF file ############

if (file.exists("tRNA.gtf")){ # Check the gtf file exits if it exits remove it
  file.remove("tRNA.gtf")
}

############### REGEX #################

df_length<-data.frame() # Empty dataframe to store sequence length
df_length<-str_extract(df_fastanames$seq_names, "(\\s\\d+\\s)") # Extract with Regex sequence length
print(df_length)

df_chr<-data.frame() # Empty dataframe to store chromosome numbers
df_chr<-gsub("chr", "", str_extract(df_fastanames$seq_names, 
                            "chr\\d+|chr\\w+"))  #  Regex to extract chromosome numbers and remove chr word
print(df_chr)                                    

df_coordinates1<-data.frame() # Empty dataframe to store first coordinate
df_coordinates1<-str_extract(df_fastanames$seq_names, "[0-9]{7,}(?=\\-)") # Regex to extract first coordinate
print(df_coordinates1)

df_coordinates2<-data.frame() # Empty dataframe to store second coordinate
df_coordinates2<-str_extract(df_fastanames$seq_names, "(?<=\\-)[0-9]{7,}") # Regex to extract second coordinate
print(df_coordinates2)

df_geneid<-data.frame() # Empty dataframe to store geneid
df_geneid<-gsub("Mus_musculus_", "",gsub(" .*$", "", df_fastanames$seq_names)) # Regex to extract geneid
print(df_geneid)

df_orientation<-data.frame() # Empty dataframe to store strand orientation
df_orientation<-gsub("(\\d+)|\\)|\\(|\\s","",
                                str_sub(df_fastanames$seq_names,106,120)) # Regex to extract orientation
print(df_orientation)

########### GTF Structure ###############

n <- nrow(df_fastanames)

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

df_gtf<<-data.frame(column1=df_chr,
                   column2=dfline1,
                   column3=dfline1.1,
                   column4=df_coordinates1,
                   column5=df_coordinates2,
                   column6=dfline2,
                   column7=df_orientation,
                   column8=dfline2,
                   column9=dfline3,
                   column10=dfliney,
                   column11=df_geneid,
                   column12=dfline4,
                   column13=dfline5,
                   column14=dfline5.1,
                   column15=dflinex,
                   columnn16=dfline9,
                   column17=dfliney,
                   column18=df_geneid,
                   column19=dfline4,
                   column20=dfline6,
                   column21=dfline6.1,
                   column22=dflinex,
                   column23=dfline7,
                   column24=dfline7.1,
                   column25=dflinex
                   )

print(df_gtf)


df_gtf$liney<- gsub(" ", "", paste(df_gtf$liney,df_gtf$column11))
df_gtf$column11<-NULL

df_gtf$liney<- gsub(" ", "", paste(df_gtf$liney,df_gtf$line4))
df_gtf$line4<-NULL

df_gtf$line3<- paste(df_gtf$line3,df_gtf$liney)
df_gtf$liney<-NULL

df_gtf$liney.1<- gsub(" ", "", paste(df_gtf$liney,df_gtf$column18))
df_gtf$column18<-NULL

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

print(df_gtf2)

### Creation of the gtf file ###

file.create("tRNA.txt") # Create a txt file to transfer the database
write.table(df_gtf2,file="tRNA.txt",row.names=FALSE,col.names = FALSE, quote=F, sep ="\t")# the quote format allow to decide which column has doublequote
file.rename("tRNA.txt", "tRNA_17.gtf") # Changing the extension to gtf

Line_gtf<- readLines(paste("~/Coding/Coding/Data/annotation/","Mus_musculus.GRCm38.96.gtf", sep = ""),
                     n=6)
print(Line_gtf)

Line_trna<- readLines(paste("~/Coding/Coding/Data/refgenome/","tRNA_17.gtf", sep = ""),
                      n=1)
print(Line_trna)

file.move("tRNA_17.gtf", "~/Coding/Coding/Data/annotation/", overwrite=TRUE)

