#Zachary Walters

#automates compassxport, msdial and msfinder to perform compound identification on all features in a dataset
#working directory should be set to a project frledr, with mzML files in a subfolder labeled 'MZML'


#select an ionization mode for MS data acquisition. Should be 'Positive' or 'Negative' (Case Sensitive)
ion_mode <- 'Positive'

#Enter paths to MS-DIAL and MS-FINDER executables and parameter files here
msdial <- "Path/to/MSDIALConsoleApp.exe"
msdial_param <- "Path/to/MS-DIAL_parameter_file"
msfinder_loc <- "Path/to/MsFinderConsoleApp.exe"
msfinder_param <- "Path/to/MS-FINDER_parameter_file"

#(Optional) The path to CompassXport can be supplied if you are using Bruker-format datafiles. 
#Otherwise, data files need to be in MZML format
compassXport <- ''  #'Path/to/CompassXport.exe'



require(stringr)


#export all .d files in directory to mzML format
if(compassxport != ''){
  da_files = list.files(pattern = "\\.d")
  if(length(da_files > 0)){
  
    for(i in(da_files)){
      cmd <- paste(c(compassXport, "-a",i,"-mode 2"), collapse = " ")
      system(cmd)
    }
  }
}

#relocate mzML files to new folder
mzml_files = list.files(pattern = "mzML", recursive = T)
if(dir.exists("MZML")){
  mzml_files2= list.files(path = "MZML",pattern = "mzML")
}

if(!dir.exists("MZML")){
  dir.create("MZML")
  mzml_files2 = list.files(path = "MZML",pattern = "mzML")
}

if(length(mzml_files2) == 0){
  for(i in (1:length(mzml_files))){
    file_path <- paste(getwd(), mzml_files[i], sep = "/")
    file_path <- shQuote(file_path, type = "cmd")
    file_split <- str_split(mzml_files[i],"/")
    renumbered <- gsub("\\.",paste(c("___________",i,"."),collapse=""),file_split[[1]][2])
    outf <- paste(getwd(), "MZML",renumbered,sep = "/")
    outf <- shQuote(outf, type = "cmd")
    cmd <- paste(c("copy" , file_path, outf), collapse = " ")
    cmd <- str_replace_all(cmd, "/","\\\\")
    shell(cmd,mustWork = T,translate = F)
  }
}


#run msdial using mzML files as input
if(!dir.exists("MSDIAL")){
  dir.create("MSDIAL")
}
inpf <- shQuote(paste(getwd(),"MZML",sep = "/"))
outf <- shQuote(paste(getwd(),"MSDIAL",sep = "/"))
cmd <- paste(c(msdial, "lcmsdda -i",inpf, "-o", outf, "-m", msdail_param), collapse = " ")
cmd <- str_replace_all(cmd,"/","\\\\")
system(cmd)


#modify alignmnet file to be readable as a data frame:
msdial_f <- list.files(pattern = "\\.msdial", recursive = T)
af <- grepl("Align",msdial_f)
alignment <- subset(msdial_f, af == T)
categories <- subset(msdial_f, af == F)
file <- readLines(alignment)
writeLines(file[4:length(file)],'alignment.csv')

#load alignemnt file as data frame, separate out features with MS/MS data
data <- read.csv('alignment.csv', header = T, sep = "\t")
data_ms2 <- subset(data, data$MS.MS.included == "True")


#Create folder sfor MS-FINDER output
dir.create("MAT1")
dir.create("MAT2")
dir.create('MAT3')

# Set possible adducts based on ionization mode
if(ion_mode == 'Positive'){
  adduct.list = c('H',"NH4","K","Na")
}
if(ion_mode == "Negative"){
  adduct_type = c('M-H', "M-Cl",'M+Na-2H','M+Hac-H')
}

#function to generate MS-FINDER input files
mat_generator <- function(peak.list,directory,adduct.list){
  peak.list = as.matrix(peak.list)
  for(i in 1:nrow(peak.list)){
    peak <-  peak.list[i,]
    ms1 <- peak[20]
    np1 <- str_count(ms1, ":")
    ms1 <- gsub(" ","\n", ms1)
    ms1 <- gsub(":"," ", ms1)
    ms2 <- peak[21]
    np2 <- str_count(ms2,":")
    ms2 <- gsub(" ","\n", ms2)
    ms2 <- gsub(":"," ", ms2)
    name <- paste("NAME:", peak[1], sep = " ")
    rt <- paste("RETENTIONTIME:", peak[2], sep = " ")
    ion.type <- paste("IONTYPE: ",ion_mode, sep = "")
    ms.type1 <- "MSTYPE: MS1"
    num.peaks1 <- paste("NUMPEAKS:",np1,sep = " ")
    ms.type2 <- "MSTYPE: MS2"
    num.peaks2 <- paste("NUMPEAKS:",np2,sep = " ")
    precursor.mz <- paste("PRECURSORMZ:", peak[3], sep = " ")
    for(a in 1:length(adduct.list)){
      adduct <- adduct.list[[a]][1]
      mz = as.numeric(peak[3])# - as.numeric(adduct.list[[a]][2])
      precursor.mz <- paste("PRECURSORMZ:", mz, sep = " ")
      if(ion_mode== 'Positive'){
        precursor.type <- paste("PRECURSORTYPE:",paste(c("[M+",adduct,"]+"),collapse = ""),sep= " ")
      }
      if(ion_mode == 'Negative'){
        precursor.type <- paste("PRECURSORTYPE:",paste(c("[",adduct,"]-"),collapse = ""),sep= " ")
      }
      output <- paste(c(name,precursor.mz,precursor.type,rt,
                         ion.type,ms.type1,num.peaks1,ms1,
                         ms.type2,num.peaks2,ms2),collapse = "\n")
      file.name<- paste(c(getwd(),"/",directory,"/",paste(peak[1],adduct,sep = '_'),".MAT"),collapse = "")
      write.table(output, file.name,quote = F,
                  row.names = F,col.names = F)
    }
  }
  return(1)
}
mat_generator(data,"MAT1",adduct_type)
mat_generator(data_ms2,"MAT2",adduct_type)
mat_generator(data_ms2,"MAT3",adduct.list = c('H'))



#Generate batch file to run MS-FINDER on the previously-generated input files
dir.create("MSFINDER")
next_script <- shQuote("D:/data/Zach/Methods/Shiny/TableGenerationForShinyApp.R")
next_script <- paste("CMD BATCH", next_script, sep = " ")
input <- shQuote(paste(getwd(),"MAT1",sep = "/"))
output <- shQuote(paste(getwd(),"MSFINDER",sep = "/"))
cmd <- paste(c(msfinder_loc, "predict -i", 
               input,"-o",output,"-m",msfinder_param),collapse = " ")
cmd <- str_replace_all(cmd,"/","\\\\")

input <- shQuote(paste(getwd(),"MAT2",sep = "/"))
output <- shQuote(paste(getwd(),"MSFINDER",sep = "/"))
cmd2 <- paste(c(msfinder_loc, "predict -i", 
               input,"-o",output,"-m",msfinder_param),collapse = " ")
cmd2 <- str_replace_all(cmd2,"/","\\\\")

input <- shQuote(paste(getwd(),"MAT3",sep = "/"))
output <- shQuote(paste(getwd(),"MSFINDER",sep = "/"))
cmd3 <- paste(c(msfinder_loc, "mssearch -i", 
                input,"-o",output,"-m",msfinder_param),collapse = " ")
cmd3 <- str_replace_all(cmd3,"/","\\\\")
writeLines(paste(c("@echo on",cmd2,cmd3,cmd,
                   #"\"C:\\Program Files\\R\\R-3.5.0\\bin\\Rscript.exe\" \"D:\\data\\Zach\\Methods\\Shiny\\TableGenerationForShinyApp.R\"",
                  "pause"),
                 collapse = "\n"),"msfinder.bat")
shell.exec("msfinder.bat")
