# FUNCTIONS

getRawPath <- function(filename){
  
  path <- c(getwd(), "/data/raw/", filename)
  return(paste(path, collapse = ""))

}

getPlotPath <- function(filename, folder = ""){
  
  curDir = getwd()            # file path of the current working directory
  subDir = "/output/plots"    # file path of the plots directory
  
  if(exists('folder')){
    
    subDir = paste(subDir, "/", folder, sep = "")  
    
    if(!dir.exists(file.path(curDir, subDir))){
      dir.create(file.path(curDir, subDir), showWarnings = FALSE)
      path <- paste(curDir, subDir, "/", filename, sep = "")
      
    } else {
      
      path <- paste(curDir, subDir, "/", filename, sep = "")
      
    }
    
  } else {
    
    path <- paste(curDir, subDir, "/", filename, sep = "")
    
  }
  
  return(path)
  
}