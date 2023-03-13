"""
Created April 2021
@author: Harikrishnan Ramadasan <harikrishnan@students.iisertirupati.ac.in>
- HR Essentialome analysis
"""

library(Biostrings)
library(splitstackshape)
library(stringr)
library(dplyr)
data("BLOSUM62")   

filenames <- Sys.glob("../data/*.fasta") #get a list of all sequence files

dfs<-list()
j<-0
# 
for (i in filenames){
  j<-j+1
  
  y=tools::file_path_sans_ext(i)
  print(y)
  t=as.list(strsplit(y,"/")[[1]])

  fastaFile = readAAStringSet(i)
  seq_name = as.vector(names(fastaFile))
  sequence = paste(fastaFile)
  df = data.frame(seq_name, sequence, stringsAsFactors = FALSE)

  df = as.data.frame(cSplit(df,"seq_name",sep = " "))

  df$seq_name_1 =  as.character(df$seq_name_1)
  reqd = as.vector(c("sequence","seq_name_1"))
  df2 = df[,reqd]
  
  # Calculate pid and store it in a list
  seqs = list()
  pids = list()
  for (i in 1:nrow(df2)){
    align = pairwiseAlignment(df2[i,1],df2[1,1],
                              type = "local",
                              substitutionMatrix = BLOSUM62,
                              gapOpening = 0,
                              gapExtension = 8)
    pcnt_id = pid(align, type = "PID2") #pid2 calculation
    seqs = c(seqs,df2$seq_name_1[i])
    pids = c(pids, pcnt_id)
    
  }
  
  df4 = as.data.frame(cbind(seq = unlist(seqs),pid = unlist(as.vector(pids))))
  dfs[[j]]<-df4

}
dff<-bind_rows(dfs, .id = "column_label")   # convert the output to a dataframe and write to a file
write.csv(dff,"../outputs/Pid2_output.csv",quote = FALSE,row.names = FALSE) #Write pid output to a file
print(dff)




