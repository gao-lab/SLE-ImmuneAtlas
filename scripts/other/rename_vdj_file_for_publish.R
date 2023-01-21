library(stringr)

files <- list.files('./data/10x_bcr/', full.names = T)
files <- files[grepl('csv$', files)]

rename_dict <- read_csv('./tmp/h5ad_rename.csv',col_names = F)

for (file in files) {
  print(file)
  name <- stringr::str_match(file, pattern ='\\/\\/[A-Z0-9]*_')
  name <- stringr::str_match(name, pattern ='[A-Z0-9]+')[[1]]
  new_name <- rename_dict[which(rename_dict$X1 == name),]$X2
  print(name)
  new_file <-  gsub(pattern = name, replacement = new_name, file)
  new_file <-  gsub(pattern = '10x_bcr', replacement = '10x_bcr_pub', new_file)
  print(new_file)
  
  file.copy(from=file, to=new_file, 
            overwrite = TRUE, recursive = FALSE, 
            copy.mode = TRUE)
}



files <- list.files('./data/10x_tcr/', full.names = T)
files <- files[grepl('csv$', files)]

rename_dict <- read_csv('./tmp/h5ad_rename.csv',col_names = F)

for (file in files) {
  print(file)
  name <- stringr::str_match(file, pattern ='\\/\\/[A-Z0-9]*_')
  name <- stringr::str_match(name, pattern ='[A-Z0-9]+')[[1]]
  new_name <- rename_dict[which(rename_dict$X1 == name),]$X2
  pattern <- paste0('\\/',file)
  new_file <-  gsub(pattern = name, replacement = new_name, file)
  new_file <-  gsub(pattern = '10x_tcr', replacement = '10x_tcr_pub', new_file)
  print(new_file)
  
  file.copy(from=file, to=new_file, 
            overwrite = TRUE, recursive = FALSE, 
            copy.mode = TRUE)
}

