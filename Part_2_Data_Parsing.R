#Zach Walters

#Data Processing
require(reticulate)
require(stringr)

#load python script into R
source_python("Path/to/ms_automator_v2.py")

#run these python functions
get_formulas()
msms_structures()
get_spectra_matches()

#load annotation files
ms1_file <- read.csv("compound_formulas.txt", sep = "\t", header = T, stringsAsFactors = F)
ms1_file <- ms1_file[,1:10]
ms1_file <- as.matrix(ms1_file, ncol = 10)
ms2_file <- read.csv("ms2_output.csv", sep = "\t", header = T, stringsAsFactors = F)
ms2_file <- as.matrix(ms2_file[,c(1:14)], ncol = 14)
score = ms2_file[,8]
score[is.na(score)] <- 0
ms2_file[,8] <- score
ms2_file[is.na(ms2_file)] <- ""
ms3_file <- read.csv("spectra_matches.txt", sep = "\t", header = F, stringsAsFactors = F)
ms3_file = as.matrix(ms3_file)
ms3_file[is.na(ms3_file)] <- 0


#combine ms1 and msms annotation files
output = matrix(c(1:15),nrow = 1)
for(i in (1:nrow(ms1_file))){
  ms1_line = ms1_file[i,]
  txt = rep('',14)
  txt[1:9] = ms1_line[1:9]
  txt[15] = ms1_line[10]
  original_line = subset(data, as.numeric(data[,1]) == as.numeric(ms1_line[1]))
  txt[2] = original_line[1,3]

  ms3_lines = subset(ms3_file, as.numeric(ms3_file[,1]) == as.numeric(ms1_line[1]))
  if(nrow(ms3_lines) != 0){
    txt[c(10,11,12,13,14)] = ms3_lines[c(4,8,5,7,6)]  #original
  }

  ms2_lines = subset(ms2_file, as.numeric(ms2_file[,1]) == as.numeric(ms1_line[1]))
  if(nrow(ms2_lines)!=0){
    if(txt[13] == '' | as.numeric(txt[13]) < 1){
      txt[10:14] = ms2_lines[c(4,6,7,8,9)]
      if(ms2_lines[9] != ''){
        txt[4:6] = ms2_lines[10:12]
        txt[8] = ms2_lines[13]
        txt[9] = ms2_lines[5]
        txt[15] = ms2_lines[14]
      }
    }
  }
  output = rbind(output, txt)
}


output = output[-1,c(1:15)]
colnames(output) = c('id','mz','rt','formula','dmz','ms1_score','isotope_source','isotope_match','adduct','ms2_source','ms2_match','compound_name','compound_score','inchi','possible_compounds')
write.table(output,'data_table_for_nx.txt', row.names = F, quote = T,sep = "\t")
output2 = read.table('data_table_for_nx.txt',header = T,sep = '\t',stringsAsFactors = F)

#make sure the annotated output matches the aligned output
paired_data = subset(data, as.numeric(data[,1]) %in% as.numeric(output2[,1]))
refined_table = cbind(output2, paired_data[,seq(22, ncol(paired_data))])


#generate tables for each sample using intensity values
for(i in 16:ncol(refined_table)){
  name = colnames(refined_table)[i]
  output = cbind(refined_table[,1:15], refined_table[i])
  colnames(refined_table)[i] = name
  out_file = paste(name, 'for_shiny.txt')
  write.table(output,out_file,sep = '\t')
}


