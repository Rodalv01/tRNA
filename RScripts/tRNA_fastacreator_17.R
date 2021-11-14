###########################################################################################
# Author: Rodrigo �lvarez Pardo                                                           #
# Date: 09/10/21                                                                          #
# Description: Script to create a fasta file with the tRNAS sequences, the tRNA sequences #
# with the same format that the headers of the primary fasta file.                        #
###########################################################################################

########## LIBRARIES ##############

library(Biostrings)
library(stringr)
library(seqRFLP)
library(filesstrings)

#### VARIABLES AND DATA FRAME #####

setwd("~/Coding/Coding/Data/refgenome") # Set directory 
getwd()
fastaFile <- readDNAStringSet("mm10-tRNAs_17.fa") # Insert fasta file in a variable
seq_name = names(fastaFile)                    # Insert sequences names in a variable
df_fastanames <-data.frame(seq_name)                # Dataframe with sequence information
sequence = paste(fastaFile)                    # Insert sequences in a variable
df_sequence <-data.frame(sequence)             # Dataframe with sequence information
df <- data.frame(seq_name, sequence)           # Insert names and sequences in a dataframe

########### Check fasta file #############

if (file.exists("tRNA_fasta.17.fa")){ # Check the gtf file exits if it exits remove it
  file.remove("tRNA_fasta.17.fa")
}

################ REGEX ##################

df_length<-data.frame() # Empty dataframe to store sequence length
df_length<-str_extract(df_fastanames$seq_name, "(\\s\\d+\\s)") # Extract with Regex sequence length
print(df_length)

df_chr<-data.frame() # Empty dataframe to store chromosome numbers
df_chr<-gsub("chr", "", str_extract(df_fastanames$seq_name, 
                            "chr\\d+|chr\\w+"))  #  Regex to extract chromosome numbers and remove chr word
print(df_chr)                                    

df_coordinates1<-data.frame() # Empty dataframe to store first coordinate
df_coordinates1<-str_extract(df_fastanames$seq_name, "[0-9]{7,}(?=\\-)") # Regex to extract first coordinate
print(df_coordinates1)

df_coordinates2<-data.frame() # Empty dataframe to store second coordinate
df_coordinates2<-str_extract(df_fastanames$seq_name, "(?<=\\-)[0-9]{7,}") # Regex to extract second coordinate
print(df_coordinates2)

df_geneid<-data.frame() # Empty dataframe to store geneid
df_geneid<-gsub("Mus_musculus_", "",gsub(" .*$", "", df_fastanames$seq_name)) # Regex to extract geneid
print(df_geneid)

df_orientation<-data.frame() # Empty dataframe to store strand orientation
df_orientation<-gsub("(\\d+)|\\)|\\(|\\s","",
                                str_sub(df_fastanames$seq_name,106,120)) # Regex to extract orientation
print(df_orientation)

##### Create tRNA fasta file ######

n<-nrow(df_fastanames)

line1<-"dna:chromosome" 
dfline1 <- data.frame(line1) 
dfline1<-do.call("rbind", replicate(n, dfline1, simplify = FALSE)) #Create a dataframe column with the same line of text 

line2<-":1:"
dfline2 <- data.frame(line2) 
dfline2<-do.call("rbind", replicate(n, dfline2, simplify = FALSE)) #Create a dataframe column with the same line of text
print(dfline2)

line3<-":1"
dfline3 <- data.frame(line3) 
dfline3<-do.call("rbind", replicate(n, dfline3, simplify = FALSE)) #Create a dataframe column with the same line of text
print(dfline3)

line4<-"REF"
dfline4 <- data.frame(line4) 
dfline4<-do.call("rbind", replicate(n, dfline4, simplify = FALSE)) #Create a dataframe column with the same line of text
print(dfline4)

line5<-"chromosome:GRCm38:" 
dfline5 <- data.frame(line5) 
dfline5<-do.call("rbind", replicate(n, dfline5, simplify = FALSE))



# Creation of a data frame with the different parts of the header in columns #

dffasta<<-data.frame(second_column =df_geneid, 
                     third_column = dfline1,
                     four_column= dfline5,
                     five_column = df_chr,
                     six_column=dfline2,
                     seven_column=df_length,
                     eigth_column=dfline3,
                     nine_column=dfline4
                     )


# Combine the columns of the data frame to create the header #

dffasta$line5<- paste(dffasta$line5,dffasta$five_column, sep = "")
dffasta$five_column<-NULL 
dffasta$line5<- paste(dffasta$line5,dffasta$line2, sep = "")
dffasta$line2<-NULL     
dffasta$line5<- gsub(" ", "", paste(dffasta$line5,dffasta$seven_column))
dffasta$seven_column<-NULL 
dffasta$line5<- paste(dffasta$line5,dffasta$line3, sep = "")
dffasta$line3<-NULL
dffasta$line5<- paste(dffasta$line5,dffasta$line4)
dffasta$line4<-NULL
dffasta$line1<-paste(dffasta$line1,dffasta$line5) 
dffasta$line5<-NULL
dffasta$second_column<-paste(dffasta$second_column,dffasta$line1) 
dffasta$line1<-NULL

dffinal<-data.frame(dffasta,sequence) # Create a dataframe with headers and sequences
print(dffinal)

# Creation of the fasta file of tRNAS with headers and nucleotide sequences

tRNA_fasta.17.fa = dataframe2fas(dffinal, file="tRNA_fasta.17.fa")

Line_fa<-readLines(paste("Mus_musculus.GRCm38.dna.primary_assembly.fa", sep = ""),
                     n=1)
print(Line_fa)

Line_trna<- readLines(paste("tRNA_fasta.17.fa", sep = ""),
                      n=1)
print(Line_trna)

