pls_vip <- function(models = models,
                    ncomp = ncomp) {
  
  #Collector
  vip <- data.table()
  
  #Number of models
  iterations <- length(models)
  
  # VIP
  for(i in 1:iterations) {
    
    #Select model
    object <- models[[i]]
    
    #Get VIP
    vip_model <- VIP(object, ncomp)
    vip <- rbind(vip, matrix(vip_model, nrow = 1))
    
  }
  
  colnames(vip) <- names(VIP(models[[1]], ncomp))
  vip <- cbind(iteration = 1:nrow(vip),
               vip)
  
  return(vip)
  
}
