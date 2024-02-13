# Saving a list of MOGP emulators
SaveMulti <- function(ems, filename){
  ell <- length(ems)
  dir.create(filename)
  for (i in 1:ell){
    tmp_em <- ems[[i]]
    save_ExUQmogp(tmp_em, filename = paste0(filename, '/', filename, '_', i))
  }
}

# Loading a list of MOGP emulators
LoadMulti <- function(filename){
  ell <- length(list.files(paste0(filename, '/'))) / 2
  tmp_list <- NULL
  for (i in 1:ell){
    tmp_list[[i]] <- load_ExUQmogp(filename = paste0(filename, '/', filename, '_', i))
  }
  return(tmp_list)
}