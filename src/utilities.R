getRawPath <- function(filename){
  
  path <- paste(getwd(), "/data/raw/", filename, sep = "")
  
  # cat(path, "\n")
  
  if(!file.exists(path)){
    
    stop(paste("The file", filename, "you're looking for in data/raw doesn't exist.", sep = " "))
  
  } else {
    
    return(path)

  }

}

getTempPath <- function(filename){
  
  path <- paste(getwd(), "/data/temp/", filename, sep = "")
  
  if(!file.exists(path)){
    
    stop(paste("The file", filename, "you're looking for in data/temp doesn't exist.", sep = " "))
    
  } else {
    
    return(path)
    
  }
  
}

getProcessedPath <- function(filename){
  
  path <- paste(getwd(), "/data/processed/", filename, sep = "")
  
  if(!file.exists(path)){
    
    stop(paste("The file", filename, "you're looking for in data/processed doesn't exist.", sep = " "))
    
  } else {
    
    return(path)
    
  }
  
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

ciao <- function(){}