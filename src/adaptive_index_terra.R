
adaptive_index <- function(RDA, K, env_pres, range = NULL, 
                           method = "loadings", scale_env, center_env){
  
  library(terra)
  
  # Environmental data to dataframe with coordinates
  var_env_proj_pres <- as.data.frame(
    env_pres[[row.names(RDA$CCA$biplot)]], 
    xy = TRUE, na.rm = TRUE
  )
  
  # Standardization of environmental variables
  var_env_proj_RDA <- as.data.frame(
    scale(
      var_env_proj_pres[ , -c(1,2)],
      center = center_env[row.names(RDA$CCA$biplot)],
      scale  = scale_env[row.names(RDA$CCA$biplot)]
    )
  )
  
  # Storage
  Proj_pres <- list()
  
  # --- Method 1: loadings ---
  if(method == "loadings"){
    for(i in 1:K){
      
      values <- apply(
        var_env_proj_RDA[, names(RDA$CCA$biplot[,i]), drop = FALSE],
        1,
        function(x) sum(x * RDA$CCA$biplot[,i])
      )
      
      df <- data.frame(
        x = var_env_proj_pres[,1],
        y = var_env_proj_pres[,2],
        z = values
      )
      
      ras_pres <- rast(df, type = "xyz", crs = crs(env_pres))
      names(ras_pres) <- paste0("RDA_pres_", i)
      
      Proj_pres[[i]] <- ras_pres
      names(Proj_pres)[i] <- paste0("RDA", i)
    }
  }
  
  # --- Method 2: predict ---
  if(method == "predict"){
    
    pred <- predict(
      RDA,
      var_env_proj_RDA[, names(RDA$CCA$biplot), drop = FALSE],
      type = "lc"
    )
    
    for(i in 1:K){
      
      df <- data.frame(
        x = var_env_proj_pres[,1],
        y = var_env_proj_pres[,2],
        z = pred[,i]
      )
      
      ras_pres <- rast(df, type = "xyz", crs = crs(env_pres))
      names(ras_pres) <- paste0("RDA_pres_", i)
      
      Proj_pres[[i]] <- ras_pres
      names(Proj_pres)[i] <- paste0("RDA", i)
    }
  }
  
  # --- Mask if range is provided ---
  if(!is.null(range)){
    Proj_pres <- lapply(Proj_pres, function(x) mask(x, range))
  }
  
  return(Proj_pres)
}