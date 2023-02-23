#!/usr/bin/env Rscript

#Cancer Data Services - CatchERR.R

#This script will take a CDS metadata manifest file and try to blindly fix the most common errors before the validation step.

##################
#
# USAGE
#
##################

#Run the following command in a terminal where R is installed for help.

#Rscript --vanilla CDS-CatchERR.R --help


##################
#
# Env. Setup
#
##################

#List of needed packages
list_of_packages=c("readr","openxlsx","dplyr","stringi","readxl","janitor","optparse","tools")

#Based on the packages that are present, install ones that are required.
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
suppressMessages(if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org"))

#Load libraries.
suppressMessages(library(readr,verbose = F))
suppressMessages(library(readxl,verbose = F))
suppressMessages(library(openxlsx, verbose = F))
suppressMessages(library(dplyr, verbose = F))
suppressMessages(library(stringi,verbose = F))
suppressMessages(library(janitor,verbose = F))
suppressMessages(library(optparse,verbose = F))
suppressMessages(library(tools,verbose = F))


#remove objects that are no longer used.
rm(list_of_packages)
rm(new.packages)


##################
#
# Arg parse
#
##################

#Option list for arg parse
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file (.xlsx, .tsv, .csv)", metavar="character"),
  make_option(c("-t", "--template"), type="character", default=NULL, 
              help="dataset template file, CDS_submission_metadata_template.xlsx", metavar="character")
)

#create list of options and values for file input
opt_parser = OptionParser(option_list=option_list, description = "\nCDS-CatchERR v2.0.4")
opt = parse_args(opt_parser)

#If no options are presented, return --help, stop and print the following message.
if (is.null(opt$file)&is.null(opt$template)){
  print_help(opt_parser)
  cat("Please supply both the input file (-f) and template file (-t), CDS_submission_metadata_template.xlsx.\n\n")
  suppressMessages(stop(call.=FALSE))
}


#Data file pathway
file_path=file_path_as_absolute(opt$file)

#Template file pathway
template_path=file_path_as_absolute(opt$template)

###########
#
# File name rework
#
###########

#Rework the file path to obtain a file extension.
file_name=stri_reverse(stri_split_fixed(stri_reverse(basename(file_path)),pattern = ".", n=2)[[1]][2])
ext=tolower(stri_reverse(stri_split_fixed(stri_reverse(basename(file_path)),pattern = ".", n=2)[[1]][1]))
path=paste(dirname(file_path),"/",sep = "")

#Output file name based on input file name and date/time stamped.
output_file=paste(file_name,
                  "_CatchERR",
                  stri_replace_all_fixed(
                    str = Sys.Date(),
                    pattern = "-",
                    replacement = ""),
                  sep="")

output_file_index=paste(file_name,
                  "_index",
                  stri_replace_all_fixed(
                    str = Sys.Date(),
                    pattern = "-",
                    replacement = ""),
                  sep="")


NA_bank=c("NA","na","N/A","n/a","")

#Read in file with trim_ws=TRUE
if (ext == "tsv"){
  df=suppressMessages(read_tsv(file = file_path, trim_ws = TRUE, na=NA_bank, guess_max = 1000000, col_types = cols(.default = col_character())))
}else if (ext == "csv"){
  df=suppressMessages(read_csv(file = file_path, trim_ws = TRUE, na=NA_bank, guess_max = 1000000, col_types = cols(.default = col_character())))
}else if (ext == "xlsx"){
  df=suppressMessages(read_xlsx(path = file_path, trim_ws = TRUE, na=NA_bank, sheet = "Metadata", guess_max = 1000000, col_types = "text"))
}else{
  stop("\n\nERROR: Please submit a data file that is in either xlsx, tsv or csv format.\n\n")
}


#############
#
# Data frame manipulation
#
#############

#Start write out for log file
sink(paste(path,output_file,".txt",sep = ""))


df_tavs=suppressMessages(read_xlsx(path = template_path, sheet = "Terms and Value Sets"))
df_tavs=remove_empty(df_tavs,c('rows','cols'))

#Pull out the positions where the value set names are located
VSN=grep(pattern = FALSE, x = is.na(df_tavs$`Value Set Name`))

df_all_terms=list()

#for each instance of a value_set_name, note the position on the Terms and Value Sets page, create a list for each with all accepted values.
for (x in 1:length(VSN)){
  if (!is.na(VSN[x+1])){
    df_all_terms[as.character(df_tavs[VSN[x],1])] = as.vector(df_tavs[VSN[x]:(VSN[x+1]-1),3])
    
  }else{
    df_all_terms[as.character(df_tavs[VSN[x],1])] = as.vector(df_tavs[VSN[x]:dim(df_tavs)[1],3])
  }
}

#Enumerated Array properties
enum_arrays=c('therapeutic_agents',"treatment_type","study_data_types","morphology","primary_site","race")

#Use the list of all accepted values for each value_set_name, and compare that against the Metadata page and determine if the values, if present, match the accepted terms.
for (value_set_name in names(df_all_terms)){
  if (value_set_name %in% enum_arrays){
    unique_values=unique(df[value_set_name][[1]])
    unique_values=unique(trimws(unlist(stri_split_fixed(str = unique_values,pattern = ";"))))
    unique_values=unique_values[!is.na(unique_values)]
    if (length(unique_values)>0){
      if (!all(unique_values%in%df_all_terms[value_set_name][[1]])){
        for (x in 1:length(unique_values)){
          check_value=unique_values[x]
          if (!is.na(check_value)){
            if (!as.character(check_value)%in%df_all_terms[value_set_name][[1]]){
              cat(paste("ERROR: ",value_set_name," property contains a value that is not recognized: ", check_value,"\n",sep = ""))
              #NEW ADDITION to push the most correct value into the proper case
              if (tolower(as.character(check_value))%in%tolower(df_all_terms[value_set_name][[1]])){
                pv_pos=grep(pattern = TRUE, x = tolower(df_all_terms[value_set_name][[1]])%in%tolower(as.character(check_value)))
                value_pos=grep(pattern = TRUE, x = df[value_set_name][[1]]%in%as.character(check_value))
                df[value_pos,value_set_name]<-df_all_terms[value_set_name][[1]][pv_pos]
                cat(paste("\tThe value in ",value_set_name," was changed: ", check_value," ---> ",df_all_terms[value_set_name][[1]][pv_pos],"\n",sep = ""))
              }
            }
          }
          
        }
        
      }else{
        cat(paste("PASS:",value_set_name,"property contains all valid values.\n"))
      }
    }
  }else if (value_set_name%in%colnames(df)){
    unique_values=unique(df[value_set_name][[1]])
    unique_values=unique_values[!is.na(unique_values)]
    if (length(unique_values)>0){
      if (!all(unique_values%in%df_all_terms[value_set_name][[1]])){
        for (x in 1:length(unique_values)){
          check_value=unique_values[x]
          if (!is.na(check_value)){
            if (!as.character(check_value)%in%df_all_terms[value_set_name][[1]]){
              cat(paste("ERROR: ",value_set_name," property contains a value that is not recognized: ", check_value,"\n",sep = ""))
              #NEW ADDITION to push the most correct value into the proper case
              if (tolower(as.character(check_value))%in%tolower(df_all_terms[value_set_name][[1]])){
                pv_pos=grep(pattern = TRUE, x = tolower(df_all_terms[value_set_name][[1]])%in%tolower(as.character(check_value)))
                value_pos=grep(pattern = TRUE, x = df[value_set_name][[1]]%in%as.character(check_value))
                df[value_pos,value_set_name]<-df_all_terms[value_set_name][[1]][pv_pos]
                cat(paste("\tThe value in ",value_set_name," was changed:\n\t", check_value," ---> ",df_all_terms[value_set_name][[1]][pv_pos],"\n",sep = ""))
              }
            }
          }
        }
      }else{
        cat(paste("PASS:",value_set_name,"property contains all valid values.\n"))
      }
    }
  }
}

##############
#
# Check and replace for non-UTF-8 characters
#
##############

cat("\n\nCertain characters do not handle being transformed into certain file types, due to this, the following characters were changed.\n----------\n")

#initialize table and then populate the key value pairs for the value that is present and what it needs to be converted to.
df_translations=tibble(value="",translation="")[0,]

#Expression to expand further pairs, copy the following line and add the new value pair.
df_translations=rbind(df_translations,tibble(value="™",translation="(TM)"))




#Grep through column and rows looking for values to change.
for (value in 1:dim(df_translations)[1]){
  df_value=df_translations[value,"value"]
  df_translate=df_translations[value,"translation"]
  colnums=grep(pattern = df_value, df)
  for (colnum in colnums){
    rownums=grep(pattern = df_value, df[,colnum][[1]])
    for (rownum in rownums){
      df[rownum,colnum]=stri_replace_all_fixed(str = df[rownum,colnum], pattern = df_value, replacement = df_translate)
    }
    cat(paste("\nWARNING: The following property, ",colnames(df[,colnum]),", contain values, ",df_value,", that can create issues when transforming, these were changed to : ",df_translate,"\n", sep = ""))
  }
}


##############
#
# ACL pattern check
#
##############

cat("\n\nThe value for ACL will be check to determine it follows the required structure, ['.*'].\n----------\n")

acl_check=unique(df$acl)

if (length(acl_check)>1){
  cat("ERROR: There is more than one ACL associated with this study and submission templte. Please only submit one ACL and corresponding data to a submission template.\n")
}else if(length(acl_check)==1){
  if (is.na(acl_check)){
    cat("ERROR: Please submit an ACL value to the 'acl' property.\n")
  }else if (!is.na(acl_check)){
    acl_test=grepl(pattern = "\\[\\'.*\\'\\]" , x= acl_check)
    if (!acl_test){
      acl_fix=paste("['",acl_check,"']", sep="")
      cat("The following ACL does not match the required structure, it will be changed:\n\t\t", acl_check, " ---> ", acl_fix,"\n",sep = "")
      df$acl=acl_fix
    }else if (acl_test){
      cat("The following ACL matches the required structure:\n\t\t", acl_check,"\n",sep = "")
    }
  }
}else{
  cat("ERROR: Something is wrong with the ACL value submitted in the acl property.\n")
}


##############
#
# Fix url paths
#
##############

#Fix urls if the url does not contain the file name but only the base url
#If full file path is in the url
for (bucket_loc in 1:dim(df)[1]){
  bucket_url=df$file_url_in_cds[bucket_loc]
  bucket_file=df$file_name[bucket_loc]
  #skip if bucket_url is NA (no associated url for file)
  if (!is.na(bucket_url)){
    #see if the file name is found in the bucket_url
    if (grepl(pattern = bucket_file,x = bucket_url)){
      #download file, run md5sum, copy to data frame and delete file
      file_name=basename(bucket_url)
      if (bucket_file!=file_name){
        cat(paste("\nERROR: There is an unresolvable issue with the file url for file: ",bucket_file,sep = ""))
      }
    }
    #if the file url has to be reworked to include the file with the base directory.
    else{
      if (substr(bucket_url,start = nchar(bucket_url),stop = nchar(bucket_url))=="/"){
        #fix the 'file_url_in_cds' section to have the full file location
        bucket_url_change=paste(bucket_url,bucket_file,sep = "")
        #double check changes made
        file_name=basename(bucket_url_change)
        if (bucket_file!=file_name){
          cat(paste("\nERROR: There is an unresolvable issue with the file url for file: ",bucket_file,sep = ""))
        }else{
          #if file name is still found, then change the url
          df$file_url_in_cds[bucket_loc]=bucket_url_change
          cat(paste("\nWARNING: The file location for the file, ", bucket_file,", has been changed:\n\t", bucket_url, " ---> ", bucket_url_change,sep = ""))
        }
      }else{
        #fix the 'file_url_in_cds' section to have the full file location
        bucket_url_change=paste(bucket_url,"/",bucket_file,sep = "")
        #double check changes made
        file_name=basename(bucket_url_change)
        if (bucket_file!=file_name){
          cat(paste("\nERROR: There is an unresolvable issue with the file url for file: ",bucket_file,sep = ""))
        }else{
          #if file name is still found, then change the url
          df$file_url_in_cds[bucket_loc]=bucket_url_change
          cat(paste("\nWARNING: The file location for the file, ", bucket_file,", has been changed:\n\t", bucket_url, " ---> ", bucket_url_change,sep = ""))
        }
      }
    }
  }
}


#close log file write out
sink()


###############
#
# Assign GUIDS to files
#
###############

cat("The file based nodes will now have a GUID assigned to each unique file.")

#Function to determine if operating system is OS is mac or linux, to run the UUID generation.
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

if ("GUID" %in% colnames(df)){
  df_index=df%>%
    select(GUID,file_url_in_cds, file_name, file_size, md5sum)
}else{
  df_index=df%>%
    select(file_url_in_cds, file_name, file_size, md5sum)%>%
    mutate(GUID=NA)
}

df_index=unique(df_index)
#For each unique file, apply a uuid to the GUID column. There is logic to handle this in both OSx and Linux, as the UUID call is different from R to the console.
pb=txtProgressBar(min=0,max=dim(df_index)[1],style = 3)
cat("\nGUID creation: \n", sep = "")

for (x in 1:dim(df_index)[1]){
  setTxtProgressBar(pb,x)
  if (is.na(df_index$GUID[x])){
    if (get_os()=="osx"){
      uuid=tolower(system(command = "uuidgen", intern = T))
    }else{
      uuid=system(command = "uuid", intern = T)
    }
    #Take the uuids in the GUID column and paste on the 'dg.4DFC/' prefix to create GUIDs for all the files.
    df_index$GUID[x]=paste("dg.4DFC/",uuid,sep = "")
  }
}


df_for_index=suppressMessages(left_join(df,df_index))


df_for_index=df_for_index%>%
  mutate(url = file_url_in_cds)%>%
  select(GUID, file_size, md5sum, url, acl, everything())


###############
#
# Write out
#
###############

#Write out file
if (ext == "tsv"){
  suppressMessages(write_tsv(df, file = paste(path,output_file,".tsv",sep = ""), na=""))
}else if (ext == "csv"){
  suppressMessages(write_csv(df, file = paste(path,output_file,".csv",sep = ""), na=""))
}else if (ext == "xlsx"){
  wb=openxlsx::loadWorkbook(file = file_path)
  openxlsx::deleteData(wb, sheet = "Metadata",rows = 1:(dim(df)[1]+1),cols=1:(dim(df)[2]+1),gridExpand = TRUE)
  openxlsx::writeData(wb=wb, sheet="Metadata", df)
  openxlsx::saveWorkbook(wb = wb,file = paste(path,output_file,".xlsx",sep = ""), overwrite = T)
}

write_tsv(x = df_for_index, file = paste(path,output_file_index,".tsv",sep = ""),na="")


cat(paste("\n\nProcess Complete.\n\nThe output file can be found here: ",path,"\n\n",sep = "")) 
