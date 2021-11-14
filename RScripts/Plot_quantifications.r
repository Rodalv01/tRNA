###########################################################################################
# Author: Rodrigo Álvarez Pardo                                                           #
# Date: 13/11/21                                                                          #
# Description: Script to compare with a scatter plot the differences between the quantifi #
# cations in the differents verisons of fasta tRNAs and gtf with coordinates              #
###########################################################################################

library(data.table, help, pos = 2, lib.loc = NULL)

setwd("/home/rodrigo/Coding/Coding/Data/Featurecounts")
list.files()

# Empty dataframes to store the quantifications #

df_t18.1<-data.frame()
df_t18<-data.frame()
df_t19<-data.frame()
df_tc19<-data.frame()

table_18.1 <- read.table('tRNA_Featurecounts_18.1.count')
df_t18.1<-data.frame(table_18.1)

setnames(df_t18.1, old=c('V1','V2','V3','V4','V5','V6','V7','V8','V9','V10',
                        'V11','V12','V13','V14','V15','V16','V17','V18','V19'),
                        new = c('geneid','tRNAid','start','end','strand','seq_length',
                        'S12','S13','S14','S15','S16','S17','S18','S19','S20','S21',
                        'S22','S23','S24'))