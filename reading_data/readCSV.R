readCSV = function(biodata){
  library(data.table)
  
  # out = list()
  
  # for(i in 1:length(biodata)){
    link = as.character(biodata)
    temp = fread(link, sep = ";")
    # temp = as.data.table(temp)
    
    out = list(temp)
  # }

  
  return(out)
}