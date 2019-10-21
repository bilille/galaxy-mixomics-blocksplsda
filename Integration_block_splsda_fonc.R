# La fonction meanSpotRepl remplace les valeurs des spots répliqués par la
# moyenne de leurs intensités.

meanSpotRepl <-function(mat)
{
  ProbeName = colnames(mat)
  isDup = duplicated(ProbeName)
  dupNames = ProbeName[isDup]
  
  for (dups in unique(dupNames))
  {
    mat[dups,] = apply(mat[which(colnames(mat) == dups), ], 2, mean)
    
  }
  
  res = mat[, -which(isDup)] # On retire de la matrice mat toutes les spots qui sont répliqués.
  
  return(res)
  
}

# La fonction supprimerVaConst supprime de la matrice mat les variables constantes.

supprimerVaConst <-function(mat)
{
  name_mat = deparse(substitute(mat))
  
  cat(paste0("Pour ", name_mat, ", avant suppression des variables constantes, il y a ", dim(mat)[2], " variables."), "\n")
  
  indiceVaConst = sapply(1:dim(mat)[2], FUN = function(i){
    col_mat_i = mat[, i]
    res = all(col_mat_i == col_mat_i[1])
    
    return(res)  
  })
  
  if(length(which(indiceVaConst == FALSE)) != 0)
  {
    res = mat[, which(indiceVaConst == FALSE)]
    
  }else{
    res = mat
    
  }
  
  cat(paste0("Pour ", name_mat, ", après suppression des variables constantes, il reste ", dim(res)[2], " variables."), "\n")
  
  return(res)
}

# La fonction supprimerVaNaInd supprime les variables contenant des NA 
# (si on garde les NA dans ces variables, les composantes du bloc ne
# sont plus orthogonales).

supprimerVaNaInd <-function(mat)
{
  indiceVaNaInd = sapply(1:dim(mat)[2], FUN = function(i){
    col_mat_i = mat[, i]
    
    indNA = length(which(is.na(col_mat_i) == TRUE))
    
    if(indNA >= 1)
    {
      res = FALSE
      
    }else{
      res = TRUE
      
    }
    
    return(res)  
  })
  
  mat2 = mat[, which(indiceVaNaInd == TRUE)]
  name_mat = deparse(substitute(mat))
  
  cat(paste0("Pour ", name_mat, ", après suppression des variables contenant des NA, il reste ", dim(mat2)[2], " variables."), "\n")
  
  return(mat2)
  
}


# La fonction varAnnotation permet de fournir des informations sur les variables
# sélectionnées pour un design.

varAnnotation <-function(variablesSelect,
                         data_transcripto_col,
                         annot_metageno_caecum)
{
  res_variablesSelect = variablesSelect
  
  noms_blocks = sapply(1:length(variablesSelect), FUN = function(i){
    ch = strsplit(names(variablesSelect)[i], split = "_")[[1]]
    
    if(ch[1] == "resBio")
    {
      res = "resBio"
      
    }else{
      res = paste(ch[1:2], collapse = "_")
      
    }
    
    return(res)
  })
  
  ind_transcripto_col = which(noms_blocks == "transcripto_col")
  
  if(length(variablesSelect[[ind_transcripto_col]]) != 0)
  {
    
    
    varSelect_transcripto_colTemp = variablesSelect[[ind_transcripto_col]] 
    varSelect_transcripto_col = sapply(1:length(varSelect_transcripto_colTemp), FUN = function(i){
      ch = strsplit(varSelect_transcripto_colTemp[i], split = "_")[[1]]
      
      res = paste(ch[2:length(ch)], collapse = "_")
      
      return(res)
    })
    
    dataframe_varSelect_transcripto_col  = data.frame(ProbeName = varSelect_transcripto_col)
    
    dataframe_annot_transcripto_col = data.frame(ProbeName = data_transcripto_col$genes$ProbeName,
                                                 GeneName = data_transcripto_col$genes$GeneName,
                                                 Description = data_transcripto_col$genes$Description,
                                                 SystematicName = data_transcripto_col$genes$SystematicName)
    
    tab_transcripto_col = join(x = dataframe_varSelect_transcripto_col, y = dataframe_annot_transcripto_col,
                               type = "inner",
                               by = "ProbeName")
    
    res_variablesSelect[[ind_transcripto_col]] = tab_transcripto_col
    
    
  }
  
  ind_metageno_caecum = which(noms_blocks == "metageno_caecum")
  
  if(length(variablesSelect[[ind_metageno_caecum]]) != 0)
  {
    varSelect_metageno_caecum_Temp1 = variablesSelect[[ind_metageno_caecum]] 
    varSelect_metageno_caecum_Temp2 = sapply(1:length(varSelect_metageno_caecum_Temp1), FUN = function(i){
      ch = strsplit(varSelect_metageno_caecum_Temp1[i], split = "")[[1]]
      
      if(ch[1] == "X")
      {
        res = paste(ch[2:length(ch)], collapse = "")
        
      }else{
        res = varSelect_metageno_caecum_Temp1[i]
        
      }
      
      return(res)
    }) 
    
    
    dataframe_varSelect_metageno_caecum = data.frame(taxon = varSelect_metageno_caecum_Temp2)
    
    dataframe_annot_metageno_caecum = data.frame(annot_metageno_caecum)
    colnames(dataframe_annot_metageno_caecum)[1] = "taxon"
    
    tab_metageno_caecum  = join(x = dataframe_varSelect_metageno_caecum, y = dataframe_annot_metageno_caecum,
                                type = "inner",
                                by = "taxon")
    
    res_variablesSelect[[ind_metageno_caecum]] = tab_metageno_caecum
    
    
  }
  
  ind_metabo_S1 = which(noms_blocks == "metabo_S1")
  
  if(length(variablesSelect[[ind_metabo_S1]]) != 0)
  {
    res_variablesSelect[[ind_metabo_S1]] = data.frame(variable = variablesSelect[[ind_metabo_S1]])
    
  }
  
  ind_resBio = which(noms_blocks == "resBio")
  
  if(length(variablesSelect[[ind_resBio]]) != 0)
  {
    res_variablesSelect[[ind_resBio]] = data.frame(variable = variablesSelect[[ind_resBio]])
    
  }
  
  
  return(res_variablesSelect)
  
  
}


# La fonction varAnnotation_gene_6blocks permet de fournir des informations sur les variables
# sélectionnées pour un design.

varAnnotation_gene_6blocks <-function(variablesSelect,
                                      data_transcripto_col,
                                      data_transcripto_tae,
                                      annot_metageno_caecum,
                                      metavar_metaboLC_S1,
                                      metavar_resBio,
                                      metavar_cyto)
{
  res_variablesSelect = variablesSelect
  noms_blocks = names(variablesSelect)
  
  ind_transcripto_col = which(noms_blocks == "transcripto_col")
  
  if(length(ind_transcripto_col) != 0)
  {
    
    varSelect_transcripto_colTemp = variablesSelect[[ind_transcripto_col]]
    
    if(length(varSelect_transcripto_colTemp) != 0)
    {
      varSelect_transcripto_col = sapply(1:length(varSelect_transcripto_colTemp), FUN = function(i){
        variable_i = varSelect_transcripto_colTemp[i]
        res = gsub("Colon_", "", variable_i, fixed = TRUE)
        
        return(res)
      })
      
      dataframe_varSelect_transcripto_col  = data.frame(GeneName = varSelect_transcripto_col)
      
      dataframe_annot_transcripto_col = data.frame(ProbeName = data_transcripto_col$genes$ProbeName,
                                                   GeneName = data_transcripto_col$genes$GeneName,
                                                   Description = data_transcripto_col$genes$Description,
                                                   SystematicName = data_transcripto_col$genes$SystematicName)
      
      tab_transcripto_col = join(x = dataframe_varSelect_transcripto_col, y = dataframe_annot_transcripto_col,
                                 type = "inner",
                                 by = "GeneName")
      
      res_variablesSelect[[ind_transcripto_col]] = tab_transcripto_col
      
    }else{
      res_variablesSelect[[ind_transcripto_col]] = ""
      
    }
    
    
    
    
  }
  
  ind_transcripto_tae = which(noms_blocks == "transcripto_tae")
  
  if(length(ind_transcripto_tae) != 0)
  {
    
    varSelect_transcripto_taeTemp = variablesSelect[[ind_transcripto_tae]] 
    
    if(length(varSelect_transcripto_taeTemp) != 0)
    {
      varSelect_transcripto_tae = sapply(1:length(varSelect_transcripto_taeTemp), FUN = function(i){
        variable_i = varSelect_transcripto_taeTemp[i]
        res = gsub("TAE_", "", variable_i, fixed = TRUE)
        
        return(res)
      })
      
      dataframe_varSelect_transcripto_tae  = data.frame(GeneName = varSelect_transcripto_tae)
      
      dataframe_annot_transcripto_tae = data.frame(ProbeName = data_transcripto_tae$genes$ProbeName,
                                                   GeneName = data_transcripto_tae$genes$GeneName,
                                                   Description = data_transcripto_tae$genes$Description,
                                                   SystematicName = data_transcripto_tae$genes$SystematicName)
      
      tab_transcripto_tae = join(x = dataframe_varSelect_transcripto_tae, y = dataframe_annot_transcripto_tae,
                                 type = "inner",
                                 by = "GeneName")
      
      res_variablesSelect[[ind_transcripto_tae]] = tab_transcripto_tae
      
    }else{
      res_variablesSelect[[ind_transcripto_tae]] = ""
      
      
    }
    
    
    
    
  }
  
  
  ind_metageno_caecum = which(noms_blocks == "metageno_caecum")
  
  if(length(ind_metageno_caecum) != 0)
  {
    varSelect_metageno_caecum_Temp1 = variablesSelect[[ind_metageno_caecum]] 
    
    if(length(varSelect_metageno_caecum_Temp1) != 0)
    {
      varSelect_metageno_caecum_Temp2 = sapply(1:length(varSelect_metageno_caecum_Temp1), FUN = function(i){
        ch = strsplit(varSelect_metageno_caecum_Temp1[i], split = "")[[1]]
        
        if(ch[1] == "X")
        {
          res = paste(ch[2:length(ch)], collapse = "")
          
        }else{
          res = varSelect_metageno_caecum_Temp1[i]
          
        }
        
        return(res)
      }) 
      
      
      dataframe_varSelect_metageno_caecum = data.frame(taxon = varSelect_metageno_caecum_Temp2)
      
      dataframe_annot_metageno_caecum = data.frame(annot_metageno_caecum)
      colnames(dataframe_annot_metageno_caecum)[1] = "taxon"
      
      tab_metageno_caecum  = join(x = dataframe_varSelect_metageno_caecum, y = dataframe_annot_metageno_caecum,
                                  type = "inner",
                                  by = "taxon")
      
      res_variablesSelect[[ind_metageno_caecum]] = tab_metageno_caecum
      
    }else{
      res_variablesSelect[[ind_metageno_caecum]] = ""
      
      
    }
    
    
    
    
  }
  
  ind_metaboLC_S1 = which(noms_blocks == "metaboLC_S1")
  
  if(length(ind_metaboLC_S1) != 0)
  {
    if(length(variablesSelect[[ind_metaboLC_S1]]) != 0)
    {
      res_variablesSelect[[ind_metaboLC_S1]] = data.frame(variable = variablesSelect[[ind_metaboLC_S1]])
      
    }else{
      res_variablesSelect[[ind_metaboLC_S1]] = ""
      
    }
    
    
  }
  
  ind_resBio = which(noms_blocks == "resBio")
  
  if(length(ind_resBio) != 0)
  {
    if(length(variablesSelect[[ind_resBio]]) != 0)
    {
      res_variablesSelect[[ind_resBio]] = data.frame(variable = variablesSelect[[ind_resBio]])
      
    }else{
      res_variablesSelect[[ind_resBio]] = ""
      
    }
    
    
  }
  
  ind_cyto = which(noms_blocks == "cyto")
  
  if(length(ind_cyto) != 0)
  {
    if(length(variablesSelect[[ind_cyto]]) != 0)
    {
      res_variablesSelect[[ind_cyto]] = data.frame(variable = variablesSelect[[ind_cyto]])
      
    }else{
      res_variablesSelect[[ind_cyto]] = ""
      
    }
    
    
  }
  
  
  return(res_variablesSelect)
  
  
}

# La fonction varAnnotation_gene_6blocks permet de fournir des informations sur les variables
# sélectionnées pour un design.

varAnnotation_gene_6blocks <-function(variablesSelect,
                                      data_transcripto_col,
                                      data_transcripto_tae,
                                      annot_metageno_caecum,
                                      metavar_metaboLC_S1,
                                      metavar_resBio,
                                      metavar_cyto)
{
  res_variablesSelect = variablesSelect
  noms_blocks = names(variablesSelect)
  
  ind_transcripto_col = which(noms_blocks == "transcripto_col")
  
  if(length(ind_transcripto_col) != 0)
  {
    
    varSelect_transcripto_colTemp = variablesSelect[[ind_transcripto_col]]
    
    if(length(varSelect_transcripto_colTemp) != 0)
    {
      varSelect_transcripto_col = sapply(1:length(varSelect_transcripto_colTemp), FUN = function(i){
        variable_i = varSelect_transcripto_colTemp[i]
        res = gsub("Colon_", "", variable_i, fixed = TRUE)
        
        return(res)
      })
      
      dataframe_varSelect_transcripto_col  = data.frame(GeneName = varSelect_transcripto_col)
      
      dataframe_annot_transcripto_col = data.frame(ProbeName = data_transcripto_col$genes$ProbeName,
                                                   GeneName = data_transcripto_col$genes$GeneName,
                                                   Description = data_transcripto_col$genes$Description,
                                                   SystematicName = data_transcripto_col$genes$SystematicName)
      
      tab_transcripto_col = join(x = dataframe_varSelect_transcripto_col, y = dataframe_annot_transcripto_col,
                                 type = "inner",
                                 by = "GeneName")
      
      res_variablesSelect[[ind_transcripto_col]] = tab_transcripto_col
      
    }else{
      res_variablesSelect[[ind_transcripto_col]] = ""
      
    }
    
    
    
    
  }
  
  ind_transcripto_tae = which(noms_blocks == "transcripto_tae")
  
  if(length(ind_transcripto_tae) != 0)
  {
    
    varSelect_transcripto_taeTemp = variablesSelect[[ind_transcripto_tae]] 
    
    if(length(varSelect_transcripto_taeTemp) != 0)
    {
      varSelect_transcripto_tae = sapply(1:length(varSelect_transcripto_taeTemp), FUN = function(i){
        variable_i = varSelect_transcripto_taeTemp[i]
        res = gsub("TAE_", "", variable_i, fixed = TRUE)
        
        return(res)
      })
      
      dataframe_varSelect_transcripto_tae  = data.frame(GeneName = varSelect_transcripto_tae)
      
      dataframe_annot_transcripto_tae = data.frame(ProbeName = data_transcripto_tae$genes$ProbeName,
                                                   GeneName = data_transcripto_tae$genes$GeneName,
                                                   Description = data_transcripto_tae$genes$Description,
                                                   SystematicName = data_transcripto_tae$genes$SystematicName)
      
      tab_transcripto_tae = join(x = dataframe_varSelect_transcripto_tae, y = dataframe_annot_transcripto_tae,
                                 type = "inner",
                                 by = "GeneName")
      
      res_variablesSelect[[ind_transcripto_tae]] = tab_transcripto_tae
      
    }else{
      res_variablesSelect[[ind_transcripto_tae]] = ""
      
      
    }
    
    
    
    
  }
  
  
  ind_metageno_caecum = which(noms_blocks == "metageno_caecum")
  
  if(length(ind_metageno_caecum) != 0)
  {
    varSelect_metageno_caecum_Temp1 = variablesSelect[[ind_metageno_caecum]] 
    
    if(length(varSelect_metageno_caecum_Temp1) != 0)
    {
      varSelect_metageno_caecum_Temp2 = sapply(1:length(varSelect_metageno_caecum_Temp1), FUN = function(i){
        ch = strsplit(varSelect_metageno_caecum_Temp1[i], split = "")[[1]]
        
        if(ch[1] == "X")
        {
          res = paste(ch[2:length(ch)], collapse = "")
          
        }else{
          res = varSelect_metageno_caecum_Temp1[i]
          
        }
        
        return(res)
      }) 
      
      
      dataframe_varSelect_metageno_caecum = data.frame(taxon = varSelect_metageno_caecum_Temp2)
      
      dataframe_annot_metageno_caecum = data.frame(annot_metageno_caecum)
      colnames(dataframe_annot_metageno_caecum)[1] = "taxon"
      
      tab_metageno_caecum  = join(x = dataframe_varSelect_metageno_caecum, y = dataframe_annot_metageno_caecum,
                                  type = "inner",
                                  by = "taxon")
      
      res_variablesSelect[[ind_metageno_caecum]] = tab_metageno_caecum
      
    }else{
      res_variablesSelect[[ind_metageno_caecum]] = ""
      
      
    }
    
    
    
    
  }
  
  ind_metaboLC_S1 = which(noms_blocks == "metaboLC_S1")
  
  if(length(ind_metaboLC_S1) != 0)
  {
    if(length(variablesSelect[[ind_metaboLC_S1]]) != 0)
    {
      
      dataframe_varSelect_metaboLC_S1 = data.frame(Variable = variablesSelect[[ind_metaboLC_S1]])
      dataframe_metavar_metaboLC_S1 = data.frame(metavar_metaboLC_S1)
      
      
      dataframe_metavar_varSelect_metaboLC_S1 = join(x = dataframe_varSelect_metaboLC_S1, y = dataframe_metavar_metaboLC_S1,
                                                     type = "inner",
                                                     by = "Variable")
      
      res_variablesSelect[[ind_metaboLC_S1]] = dataframe_metavar_varSelect_metaboLC_S1
      
    }else{
      res_variablesSelect[[ind_metaboLC_S1]] = ""
      
    }
    
  }
  
  
  ind_resBio = which(noms_blocks == "resBio")
  
  if(length(ind_resBio) != 0)
  {
    if(length(variablesSelect[[ind_resBio]]) != 0)
    {
      
      dataframe_varSelect_resBio = data.frame(Variable = variablesSelect[[ind_resBio]])
      dataframe_metavar_resBio = data.frame(metavar_resBio)
      
      
      dataframe_metavar_varSelect_resBio = join(x = dataframe_varSelect_resBio, y = dataframe_metavar_resBio,
                                                type = "inner",
                                                by = "Variable")
      
      res_variablesSelect[[ind_resBio]] = dataframe_metavar_varSelect_resBio
      
    }else{
      res_variablesSelect[[ind_resBio]] = ""
      
    }
    
  }
  
  ind_cyto = which(noms_blocks == "cyto")
  
  if(length(ind_cyto) != 0)
  {
    if(length(variablesSelect[[ind_cyto]]) != 0)
    {
      
      dataframe_varSelect_cyto = data.frame(Variable = variablesSelect[[ind_cyto]])
      dataframe_metavar_cyto = data.frame(metavar_cyto)
      
      
      dataframe_metavar_varSelect_cyto = join(x = dataframe_varSelect_cyto, y = dataframe_metavar_cyto,
                                              type = "inner",
                                              by = "Variable")
      
      res_variablesSelect[[ind_cyto]] = dataframe_metavar_varSelect_cyto
      
    }else{
      res_variablesSelect[[ind_cyto]] = ""
      
    }
    
  }
  
  
  return(res_variablesSelect)
  
  
}

# La fonction varAnnotation_gene_7blocks permet de fournir des informations sur les variables
# sélectionnées pour un design.

varAnnotation_gene_7blocks <-function(variablesSelect,
                                      data_transcripto_col,
                                      data_transcripto_tae,
                                      annot_metageno_caecum,
                                      metavar_metaboLC_S1,
                                      metavar_resBio,
                                      metavar_cyto)
{
  res_variablesSelect = variablesSelect
  noms_blocks = names(variablesSelect)
  
  ind_transcripto_col = which(noms_blocks == "transcripto_col")
  
  if(length(ind_transcripto_col) != 0)
  {
    
    varSelect_transcripto_colTemp = variablesSelect[[ind_transcripto_col]]
    
    if(length(varSelect_transcripto_colTemp) != 0)
    {
      varSelect_transcripto_col = sapply(1:length(varSelect_transcripto_colTemp), FUN = function(i){
        variable_i = varSelect_transcripto_colTemp[i]
        res = gsub("Colon_", "", variable_i, fixed = TRUE)
        
        return(res)
      })
      
      dataframe_varSelect_transcripto_col  = data.frame(GeneName = varSelect_transcripto_col)
      
      dataframe_annot_transcripto_col = data.frame(ProbeName = data_transcripto_col$genes$ProbeName,
                                                   GeneName = data_transcripto_col$genes$GeneName,
                                                   Description = data_transcripto_col$genes$Description,
                                                   SystematicName = data_transcripto_col$genes$SystematicName)
      
      tab_transcripto_col = join(x = dataframe_varSelect_transcripto_col, y = dataframe_annot_transcripto_col,
                                 type = "inner",
                                 by = "GeneName")
      
      res_variablesSelect[[ind_transcripto_col]] = tab_transcripto_col
      
    }else{
      res_variablesSelect[[ind_transcripto_col]] = ""
      
    }
    
    
    
    
  }
  
  ind_transcripto_tae = which(noms_blocks == "transcripto_tae")
  
  if(length(ind_transcripto_tae) != 0)
  {
    
    varSelect_transcripto_taeTemp = variablesSelect[[ind_transcripto_tae]] 
    
    if(length(varSelect_transcripto_taeTemp) != 0)
    {
      varSelect_transcripto_tae = sapply(1:length(varSelect_transcripto_taeTemp), FUN = function(i){
        variable_i = varSelect_transcripto_taeTemp[i]
        res = gsub("TAE_", "", variable_i, fixed = TRUE)
        
        return(res)
      })
      
      dataframe_varSelect_transcripto_tae  = data.frame(GeneName = varSelect_transcripto_tae)
      
      dataframe_annot_transcripto_tae = data.frame(ProbeName = data_transcripto_tae$genes$ProbeName,
                                                   GeneName = data_transcripto_tae$genes$GeneName,
                                                   Description = data_transcripto_tae$genes$Description,
                                                   SystematicName = data_transcripto_tae$genes$SystematicName)
      
      tab_transcripto_tae = join(x = dataframe_varSelect_transcripto_tae, y = dataframe_annot_transcripto_tae,
                                 type = "inner",
                                 by = "GeneName")
      
      res_variablesSelect[[ind_transcripto_tae]] = tab_transcripto_tae
      
    }else{
      res_variablesSelect[[ind_transcripto_tae]] = ""
      
      
    }
    
    
    
    
  }
  
  
  ind_metageno_caecum = which(noms_blocks == "metageno_caecum")
  
  if(length(ind_metageno_caecum) != 0)
  {
    varSelect_metageno_caecum_Temp1 = variablesSelect[[ind_metageno_caecum]] 
    
    if(length(varSelect_metageno_caecum_Temp1) != 0)
    {
      varSelect_metageno_caecum_Temp2 = sapply(1:length(varSelect_metageno_caecum_Temp1), FUN = function(i){
        ch = strsplit(varSelect_metageno_caecum_Temp1[i], split = "")[[1]]
        
        if(ch[1] == "X")
        {
          res = paste(ch[2:length(ch)], collapse = "")
          
        }else{
          res = varSelect_metageno_caecum_Temp1[i]
          
        }
        
        return(res)
      }) 
      
      
      dataframe_varSelect_metageno_caecum = data.frame(taxon = varSelect_metageno_caecum_Temp2)
      
      dataframe_annot_metageno_caecum = data.frame(annot_metageno_caecum)
      colnames(dataframe_annot_metageno_caecum)[1] = "taxon"
      
      tab_metageno_caecum  = join(x = dataframe_varSelect_metageno_caecum, y = dataframe_annot_metageno_caecum,
                                  type = "inner",
                                  by = "taxon")
      
      res_variablesSelect[[ind_metageno_caecum]] = tab_metageno_caecum
      
    }else{
      res_variablesSelect[[ind_metageno_caecum]] = ""
      
      
    }
    
    
    
    
  }
  
  ind_metaboLC_S1 = which(noms_blocks == "metaboLC_S1")
  
  if(length(ind_metaboLC_S1) != 0)
  {
    if(length(variablesSelect[[ind_metaboLC_S1]]) != 0)
    {
      
      dataframe_varSelect_metaboLC_S1 = data.frame(Variable = variablesSelect[[ind_metaboLC_S1]])
      dataframe_metavar_metaboLC_S1 = data.frame(metavar_metaboLC_S1)
      
      
      dataframe_metavar_varSelect_metaboLC_S1 = join(x = dataframe_varSelect_metaboLC_S1, y = dataframe_metavar_metaboLC_S1,
                                                     type = "inner",
                                                     by = "Variable")
      
      res_variablesSelect[[ind_metaboLC_S1]] = dataframe_metavar_varSelect_metaboLC_S1
      
    }else{
      res_variablesSelect[[ind_metaboLC_S1]] = ""
      
    }
    
  }
  
  ind_metaboGC_S1 = which(noms_blocks == "metaboGC_S1")
  
  if(length(ind_metaboGC_S1) != 0)
  {
    if(length(variablesSelect[[ind_metaboGC_S1]]) != 0)
    {
      res_variablesSelect[[ind_metaboGC_S1]] = data.frame(Variable = variablesSelect[[ind_metaboGC_S1]])
      
    }else{
      res_variablesSelect[[ind_metaboGC_S1]] = ""
      
    }
    
    
  }
  
  ind_resBio = which(noms_blocks == "resBio")
  
  if(length(ind_resBio) != 0)
  {
    if(length(variablesSelect[[ind_resBio]]) != 0)
    {
      
      dataframe_varSelect_resBio = data.frame(Variable = variablesSelect[[ind_resBio]])
      dataframe_metavar_resBio = data.frame(metavar_resBio)
      
      
      dataframe_metavar_varSelect_resBio = join(x = dataframe_varSelect_resBio, y = dataframe_metavar_resBio,
                                                type = "inner",
                                                by = "Variable")
      
      res_variablesSelect[[ind_resBio]] = dataframe_metavar_varSelect_resBio
      
    }else{
      res_variablesSelect[[ind_resBio]] = ""
      
    }
    
  }
  
  ind_cyto = which(noms_blocks == "cyto")
  
  if(length(ind_cyto) != 0)
  {
    if(length(variablesSelect[[ind_cyto]]) != 0)
    {
      
      dataframe_varSelect_cyto = data.frame(Variable = variablesSelect[[ind_cyto]])
      dataframe_metavar_cyto = data.frame(metavar_cyto)
      
      
      dataframe_metavar_varSelect_cyto = join(x = dataframe_varSelect_cyto, y = dataframe_metavar_cyto,
                                              type = "inner",
                                              by = "Variable")
      
      res_variablesSelect[[ind_cyto]] = dataframe_metavar_varSelect_cyto
      
    }else{
      res_variablesSelect[[ind_cyto]] = ""
      
    }
    
  }
  
  
  return(res_variablesSelect)
  
  
}



# Integration 6 blocs norm sonde ------------------------------------------

# La fonction plotVarZoom permet de zoomer sur le cercle de corrélation et de récupérer les variables 
# contenues dans ce rectangle.

plotVarZoom <-function(res_block_splsda,
                       comp = 1:2,
                       blocks,
                       block_Y = NULL,
                       vec_col,
                       cutoff,
                       rad.in = 0.5,
                       min.X = -1,
                       max.X = 1,
                       min.Y = -1,
                       max.Y = 1,
                       cex = 0.7,
                       cex_legend = 0.8,
                       pos = c(1.2, 0),
                       pch = 20,
                       inset = c(-0.25, 0))
{
  if(class(res_block_splsda)[1] == "block.splsda")
  {
    circle = list()
    circle[[1]] = ellipse(0, levels = 1, t = 1)
    circle[[2]] = ellipse(0, levels = 1, t = rad.in)
    circle = data.frame(do.call("rbind", circle), "Circle" = c(rep("Main circle", 100), rep("Inner circle", 100)))
    
    MainCircle = circle[grep("Main circle", circle[, 3]), ]
    InnerCircle = circle[grep("Inner circle", circle[, 3]), ]
    
    
    if(length(blocks) > 1)
    {
      noms_bloc = names(res_block_splsda$variates)
      mat_comp1 = sapply(blocks, FUN = function(i){
        res = res_block_splsda$variates[[i]][, 1]
        
        return(res)
      })
      
      colnames(mat_comp1) = noms_bloc[blocks]
      mat_cor_comp1 = cor(mat_comp1)
      
      
      mat_comp2 = sapply(blocks, FUN = function(i){
        res = res_block_splsda$variates[[i]][, 2]
        
        return(res)
      })
      
      colnames(mat_comp2) = noms_bloc[blocks]
      mat_cor_comp2 = cor(mat_comp2)
      
    } # Fin if(length(blocks) > 1).
    
    # Pour chaque bloc, calcul des corrélations entre la première
    # composante et les variables sélectionnées et les corrélations entre
    # la deuxième composante et les variables sélectionnées. Pour la réponse,
    # calcul de la corrélation entre les variables de la réponse et la première composante
    # du premier bloc sélectionné et de la corrélation entre les variables de
    # la réponse et la deuxième composante du premier bloc sélectionné.
    
    liste_matCor_comp_var = list()
    
    varSelect_comp1 = selectVar(res_block_splsda,
                                comp = comp[1])
    
    varSelect_comp2 = selectVar(res_block_splsda,
                                comp = comp[2])
    
    vec_nom_blockEtReponse = c()
    blocksEtReponse = c(blocks, which(res_block_splsda$names$blocks == "Y"))
    
    for(i in 1:length(blocksEtReponse))
    {
      indice_block_i = blocksEtReponse[i]
      nom_blockEtReponse_i = res_block_splsda$names$blocks[indice_block_i]
      
      if(nom_blockEtReponse_i == "Y")
      {
        if(!is.null(block_Y))
        {
          block_i = block_Y
          
          comp1 = res_block_splsda$variates[[blocks[1]]][, comp[1]]
          comp2 = res_block_splsda$variates[[blocks[1]]][, comp[2]]
          
          vecCor_comp1_var = sapply(1:dim(block_i)[2], FUN = function(j){
            cor(comp1, block_i[, j], use = "complete.obs")
          })
          
          vecCor_comp2_var = sapply(1:dim(block_i)[2], FUN = function(j){
            cor(comp2, block_i[, j], use = "complete.obs")
            
          })
          
          
        }else{
          cat("La réponse n'est pas saisie comme paramètre d'entrée", "\n")
          
        }
        
        
        
      }else{
        comp1 = res_block_splsda$variates[[indice_block_i]][, comp[1]]
        comp2 = res_block_splsda$variates[[indice_block_i]][, comp[2]]
        
        varSelect_comp1_i = varSelect_comp1[[indice_block_i]][[1]]
        varSelect_comp2_i = varSelect_comp2[[indice_block_i]][[1]]
        varSelect_i = unique(c(varSelect_comp1_i, varSelect_comp2_i))
        
        block_i = res_block_splsda$X[[indice_block_i]][, varSelect_i]
        
        if(i == 1)
        {
          vecCor_comp1_var = sapply(1:dim(block_i)[2], FUN = function(j){
            cor(comp1, block_i[, j], use = "complete.obs")
          })
          
          vecCor_comp2_var = sapply(1:dim(block_i)[2], FUN = function(j){
            cor(comp2, block_i[, j], use = "complete.obs")
          })
          
        }else{
          signeCor_comp1 = sign(mat_cor_comp1[1, indice_block_i])
          
          vecCor_comp1_var = sapply(1:dim(block_i)[2], FUN = function(j){
            res = signeCor_comp1*cor(comp1, block_i[, j], use = "complete.obs")
            
            return(res)
          })
          
          signeCor_comp2 = sign(mat_cor_comp2[1, indice_block_i])
          
          vecCor_comp2_var = sapply(1:dim(block_i)[2], FUN = function(j){
            res = signeCor_comp2*cor(comp2, block_i[, j], use = "complete.obs")
            
            return(res)
          })
          
        }
        
        
      }
      
      matCor_comp_var = rbind(vecCor_comp1_var,
                              vecCor_comp2_var)
      
      colnames(matCor_comp_var) = colnames(block_i)        
      
      
      liste_matCor_comp_var[[i]] = matCor_comp_var
      
      vec_nom_blockEtReponse_i = rep(nom_blockEtReponse_i, dim(block_i)[2])
      vec_nom_blockEtReponse = c(vec_nom_blockEtReponse, vec_nom_blockEtReponse_i)
      
    } # Fin for(i in 1:length(blocks)).
    
    
    matCor_Allcomp_Allvar = t(Reduce(cbind, liste_matCor_comp_var))
    
    # indice permet de récupérer les variables de chaque bloc fortement corrélées avec soit
    # la première composante ou la deuxième composante dans une partie du cercle de corrélation
    # et de récupérer la variable réponse dans une partie du cercle de corrélation. 
    
    indice = sapply(1:dim(matCor_Allcomp_Allvar)[1], FUN = function(k){
      cor1 = matCor_Allcomp_Allvar[k, 1]
      cor2 = matCor_Allcomp_Allvar[k, 2]
      blockEtReponse_k = vec_nom_blockEtReponse[k]
      
      if(blockEtReponse_k == "Y")
      {
        cond2 = cor1 > min.X & cor1 < max.X & cor2 > min.Y & cor2 < max.Y
        
      }else{
        cond1 = abs(cor1) > cutoff | abs(cor2) > cutoff
        
        cond2 = cor1 > min.X & cor1 < max.X & cor2 > min.Y & cor2 < max.Y & cond1
        
      }
      
      
      return(cond2)
    })
    
    matCor_Allcomp_AllvarSelect = matCor_Allcomp_Allvar[indice, , drop = FALSE]
    varSelect = rownames(matCor_Allcomp_AllvarSelect)
    
    dataframe_Cor_Allcomp_Allvar = data.frame(cbind(rownames(matCor_Allcomp_Allvar),
                                                    vec_nom_blockEtReponse,
                                                    matCor_Allcomp_Allvar))
    
    colnames(dataframe_Cor_Allcomp_Allvar) = c("variable",
                                               "bloc",
                                               "cor_comp1_var",
                                               "cor_comp2_var")
    
    dataframe_Cor_Allcomp_Allvar[, 1:2] = apply(dataframe_Cor_Allcomp_Allvar[, 1:2], 2, as.character)
    
    dataframe_Cor_Allcomp_AllvarSelect = dataframe_Cor_Allcomp_Allvar[indice, ]
    
    
    # Tracé de la superposition des cerles de corrélation.
    
    plot(MainCircle[, 1], MainCircle[, 2],
         type = "l",
         xlab = paste0("composante ", comp[1]),
         ylab = paste0("composante ", comp[2]))
    
    points(InnerCircle[, 1], InnerCircle[, 2],
           type = "l")
    
    if(dim(matCor_Allcomp_AllvarSelect)[1] != 0)
    {
      
      nom_blockEtReponseSelect = unique(dataframe_Cor_Allcomp_AllvarSelect[, 2])
      indice_blockEtReponseSelect = sapply(1:length(nom_blockEtReponseSelect), FUN = function(i){
        ind = which(res_block_splsda$names$blocks == nom_blockEtReponseSelect[i])
        
        if(length(ind) != 0)
        {
          res = ind
          
        }else{
          res = which(res_block_splsda$names$blocks == "Y")
          
        }
        
        return(res)
      })
      vec_colSelect = vec_col[indice_blockEtReponseSelect]
      
      if(length(blocksEtReponse) == 1)
      {
        points(matCor_Allcomp_AllvarSelect[, 1], matCor_Allcomp_AllvarSelect[, 2],
               col  = NULL)
        
        text(matCor_Allcomp_AllvarSelect[, 1], matCor_Allcomp_AllvarSelect[, 2],
             labels = rownames(matCor_Allcomp_AllvarSelect),
             cex = cex,
             col = vec_colSelect[1])
        
      }else{
        
        nbVarSelect_bloc = cumsum(sapply(1:length(nom_blockEtReponseSelect), FUN = function(j){
          res = length(which(dataframe_Cor_Allcomp_AllvarSelect[, 2] == nom_blockEtReponseSelect[j]))
          
          return(res)
        }))
        
        for(i in 1:length(nbVarSelect_bloc))
        {
          if(i == 1)
          {
            indice1 = 1:nbVarSelect_bloc[1]
            
            if(length(indice1) != 0)
            {
              matCor_Allcomp_AllvarSelect2 = matCor_Allcomp_AllvarSelect[indice1, , drop = FALSE]
              
              points(matCor_Allcomp_AllvarSelect2[, 1], matCor_Allcomp_AllvarSelect2[, 2],
                     col  = NULL)
              
              text(matCor_Allcomp_AllvarSelect2[, 1], matCor_Allcomp_AllvarSelect2[, 2],
                   labels = rownames(matCor_Allcomp_AllvarSelect2),
                   cex = cex,
                   col = rep(vec_colSelect[i], dim(matCor_Allcomp_AllvarSelect2)[1]))
              
            }else{
              cat(paste0("Il n'y a de variables dans cette zone du cercle de corrélation pour le bloc ", nom_blockEtReponseSelect[i]), "\n")
              
              
            }
            
            
          }else{
            indice2 = (nbVarSelect_bloc[i - 1] + 1):nbVarSelect_bloc[i]
            
            if(length(indice2) != 0)
            {
              matCor_Allcomp_AllvarSelect2 = matCor_Allcomp_AllvarSelect[indice2, , drop = FALSE]
              
              points(matCor_Allcomp_AllvarSelect2[, 1], matCor_Allcomp_AllvarSelect2[, 2],
                     col  = NULL)
              
              text(matCor_Allcomp_AllvarSelect2[, 1], matCor_Allcomp_AllvarSelect2[, 2],
                   labels = rownames(matCor_Allcomp_AllvarSelect2),
                   cex = cex,
                   col = rep(vec_colSelect[i], dim(matCor_Allcomp_AllvarSelect2)[1]))
              
              
            }else{
              cat(paste0("Il n'y a de variables dans cette zone du cercle de corrélation pour le bloc ", nom_blockEtReponseSelect[i]), "\n")
              
            }
            
            
          }
          
        }
        
      }
      
      par(xpd = TRUE)
      legend(x = pos[1], y = pos[2],
             legend = nom_blockEtReponseSelect,
             pch = pch,
             col = vec_colSelect,
             cex = cex_legend,
             inset = inset)
      
      
    }else{
      cat("Il n'y a de variables dans cette zone du cercle de corrélation", "\n")
      
    }
    
    
    # Suppression de dataframe_Cor_Allcomp_Allvar, dataframe_Cor_Allcomp_AllvarSelect et varSelect de la variable
    # réponse.
    
    if(!is.null(block_Y))
    {
      for(i in 1:dim(block_Y)[2])
      {
        variableReponse = colnames(block_Y)[i]
        
        # dataframe_Cor_Allcomp_Allvar
        if(dim(dataframe_Cor_Allcomp_Allvar)[1] != 0)
        {
          ind_Allvar = which(dataframe_Cor_Allcomp_Allvar$variable == variableReponse)
          
          if(length(ind_Allvar) != 0)
          {
            dataframe_Cor_Allcomp_Allvar = dataframe_Cor_Allcomp_Allvar[- ind_Allvar, ]
            
          }
          
        }
        
        # dataframe_Cor_Allcomp_AllvarSelect
        if(dim(dataframe_Cor_Allcomp_AllvarSelect)[1] != 0)
        {
          ind_AllvarSelect = which(dataframe_Cor_Allcomp_AllvarSelect$variable == variableReponse)
          
          if(length(ind_AllvarSelect) != 0)
          {
            dataframe_Cor_Allcomp_AllvarSelect = dataframe_Cor_Allcomp_AllvarSelect[- ind_AllvarSelect, ]
            
          }
          
        }
        
        # varSelect
        if(length(varSelect) != 0)
        {
          ind_varSelect = which(varSelect == variableReponse)
          
          if(length(ind_varSelect) != 0)
          {
            varSelect = varSelect[- ind_varSelect]
            
          }
          
        }
        
        
        
      } # Fin for(i in 1:dim(block_Y)[2]).
      
    }
    
    
    res = list(dataframe_Cor_Allcomp_Allvar = dataframe_Cor_Allcomp_Allvar,
               dataframe_Cor_Allcomp_AllvarSelect = dataframe_Cor_Allcomp_AllvarSelect,
               varSelect = varSelect)
    
    return(res)
    
  }else{
    cat("Erreur : il ne s'agit pas de la sortie de la fonction block.splsda", "\n")
    
  }
  
}


# La fonction networkVarSelect permet de tracer un réseau pour les variables de
# certains blocs.

networkVarSelect <-function(object,
                            mat_Y,
                            comp = 1:2,
                            listeVar,
                            blocks,
                            cutoff = 0
)
{
  
  if(class(object)[1] == "block.splsda")
  {
    nomBlocs = names(listeVar)
    
    liste_XSelect = lapply(1:length(listeVar), FUN = function(i){
      nomBloc_i = names(listeVar)[i]
      ind1 = which(object$names$blocks == nomBloc_i)
      ind2 = colnames(object$X[[ind1]])%in%listeVar[[i]]
      res = object$X[[ind1]][, ind2, drop = FALSE]
      
      return(res)
    })
    names(liste_XSelect) = nomBlocs
    
    liste_matSelect = liste_XSelect
    
    if(!is.null(mat_Y))
    {
      liste_matSelect = c(liste_matSelect, list(mat_Y))
      names(liste_matSelect)[length(liste_matSelect)] = "Y"
      
    }
    
    # compute the similarity between var1 of block1 and var2 of block2.
    coord = list()
    
    for(k in 1:length(blocks))
    {
      nomBloc_k = names(liste_matSelect)[k]
      mat_k = liste_matSelect[[k]][, ,drop = FALSE]
      
      if(nomBloc_k == "Y")
      {
        coord[[k]] = cor(mat_k, object$variates[[blocks[1]]][, comp])
        
        
        
      }else{
        coord[[k]] = cor(mat_k, object$variates[[blocks[k]]][, comp])
        
        
      }
      
      if(dim(mat_k)[2] == 1)
      {
        coord[[k]] = as.matrix(coord[[k]],
                               nrow = 1)
        rownames(coord[[k]]) = colnames(mat_k)
        
      }
      
    } # Fin for(k in 1:length(blocks)).
    
    l = 1
    M_block = list()
    node.X1 = node.X2 = w = NULL
    
    
    for(j in 1:(length(blocks) - 1))
    {
      for(k in (j + 1):length(blocks))
      {
        M_block[[l]] = coord[[j]][, comp, drop = FALSE] %*% t(coord[[k]][, comp, drop = FALSE])
        
        X1 = rownames(coord[[j]])
        X2 = rownames(coord[[k]])
        
        rep.X1 = rep(X1, each = length(X2))
        rep.X2 = rep(X2, length(X1))
        
        node.X1= c(node.X1, rep.X1)
        node.X2 = c(node.X2, rep.X2)
        
        w = c(w, as.vector(t(M_block[[l]])))
        
        l = l + 1
        
      } # Fin for(k in (j + 1):length(blocks)).
      
    } # Fin for(j in 1:(length(blocks) - 1)).
    
    # nodes
    group = NULL
    temp = lapply(liste_matSelect, function(x) colnames(x))
    
    for (i in 1:length(temp))
    {
      group = c(group, rep(names(temp)[i], length(temp[[i]])))
      
    } # Fin for (i in 1:length(temp)).
    
    # nodes
    nodes = data.frame(name = unlist(temp),
                       group = group)
    
    # gR
    relations = data.frame(from = node.X1,
                           to = node.X2,
                           weight = w)
    
    idx = (abs(w) >= cutoff)
    relations = relations[idx, , drop = FALSE]
    
    gR = graph.data.frame(relations,
                          directed = FALSE,
                          vertices = nodes)
    
    
    block.var.names = sapply(1:length(liste_matSelect), FUN = function(i){
      res = colnames(liste_matSelect[[i]])
      
      return(res)
    })
    V(gR)$label = unlist(block.var.names)
    
    gR = delete.vertices(gR, which(degree(gR) == 0))
    
    
    res = list(gR = gR)
    
    l = 1
    for (i in 1:(length(blocks)-1))
    {
      for (j in (i + 1):length(blocks))
      {
        res[paste("M", names(liste_matSelect)[i], names(liste_matSelect)[j], sep="_")] = list(M_block[[l]])
        l = l + 1
        
      } # Fin for (j in (i + 1):length(blocks)).
      
    } # Fin for (i in 1:(length(blocks)-1).
    
    res$cutoff = cutoff
    
    return(res)
    
  } # Fin if(class(object)[1] == "block.splsda").
  
  
  
}

# La fonction varAnnotation_6blocks permet de fournir des informations sur les variables
# sélectionnées pour un design.

varAnnotation_6blocks <-function(variablesSelect,
                                 data_transcripto_col,
                                 data_transcripto_tae,
                                 annot_metageno_caecum)
{
  res_variablesSelect = variablesSelect
  noms_blocks = names(variablesSelect)
  
  ind_transcripto_col = which(noms_blocks == "transcripto_col")
  
  if(length(variablesSelect[[ind_transcripto_col]]) != 0)
  {
    
    varSelect_transcripto_colTemp = variablesSelect[[ind_transcripto_col]] 
    varSelect_transcripto_col = sapply(1:length(varSelect_transcripto_colTemp), FUN = function(i){
      ch = strsplit(varSelect_transcripto_colTemp[i], split = "_")[[1]]
      
      res = paste(ch[2:length(ch)], collapse = "_")
      
      return(res)
    })
    
    dataframe_varSelect_transcripto_col  = data.frame(ProbeName = varSelect_transcripto_col)
    
    dataframe_annot_transcripto_col = data.frame(ProbeName = data_transcripto_col$genes$ProbeName,
                                                 GeneName = data_transcripto_col$genes$GeneName,
                                                 Description = data_transcripto_col$genes$Description,
                                                 SystematicName = data_transcripto_col$genes$SystematicName)
    
    tab_transcripto_col = join(x = dataframe_varSelect_transcripto_col, y = dataframe_annot_transcripto_col,
                               type = "inner",
                               by = "ProbeName")
    
    res_variablesSelect[[ind_transcripto_col]] = tab_transcripto_col
    
    
  }
  
  ind_transcripto_tae = which(noms_blocks == "transcripto_tae")
  
  if(length(variablesSelect[[ind_transcripto_tae]]) != 0)
  {
    varSelect_transcripto_taeTemp = variablesSelect[[ind_transcripto_tae]] 
    varSelect_transcripto_tae = sapply(1:length(varSelect_transcripto_taeTemp), FUN = function(i){
      ch = strsplit(varSelect_transcripto_taeTemp[i], split = "_")[[1]]
      
      res = paste(ch[2:length(ch)], collapse = "_")
      
      return(res)
    })
    
    dataframe_varSelect_transcripto_tae  = data.frame(ProbeName = varSelect_transcripto_tae)
    
    dataframe_annot_transcripto_tae = data.frame(ProbeName = data_transcripto_tae$genes$ProbeName,
                                                 GeneName = data_transcripto_tae$genes$GeneName,
                                                 Description = data_transcripto_tae$genes$Description,
                                                 SystematicName = data_transcripto_tae$genes$SystematicName)
    
    tab_transcripto_tae = join(x = dataframe_varSelect_transcripto_tae, y = dataframe_annot_transcripto_tae,
                               type = "inner",
                               by = "ProbeName")
    
    res_variablesSelect[[ind_transcripto_tae]] = tab_transcripto_tae
    
  }
  
  
  ind_metageno_caecum = which(noms_blocks == "metageno_caecum")
  
  if(length(variablesSelect[[ind_metageno_caecum]]) != 0)
  {
    varSelect_metageno_caecum_Temp1 = variablesSelect[[ind_metageno_caecum]] 
    varSelect_metageno_caecum_Temp2 = sapply(1:length(varSelect_metageno_caecum_Temp1), FUN = function(i){
      ch = strsplit(varSelect_metageno_caecum_Temp1[i], split = "")[[1]]
      
      if(ch[1] == "X")
      {
        res = paste(ch[2:length(ch)], collapse = "")
        
      }else{
        res = varSelect_metageno_caecum_Temp1[i]
        
      }
      
      return(res)
    }) 
    
    
    dataframe_varSelect_metageno_caecum = data.frame(taxon = varSelect_metageno_caecum_Temp2)
    
    dataframe_annot_metageno_caecum = data.frame(annot_metageno_caecum)
    colnames(dataframe_annot_metageno_caecum)[1] = "taxon"
    
    tab_metageno_caecum  = join(x = dataframe_varSelect_metageno_caecum, y = dataframe_annot_metageno_caecum,
                                type = "inner",
                                by = "taxon")
    
    res_variablesSelect[[ind_metageno_caecum]] = tab_metageno_caecum
    
    
  }
  
  ind_metaboLC_S1 = which(noms_blocks == "metaboLC_S1")
  
  if(length(variablesSelect[[ind_metaboLC_S1]]) != 0)
  {
    res_variablesSelect[[ind_metaboLC_S1]] = data.frame(variable = variablesSelect[[ind_metaboLC_S1]])
    
  }
  
  ind_resBio = which(noms_blocks == "resBio")
  
  if(length(variablesSelect[[ind_resBio]]) != 0)
  {
    res_variablesSelect[[ind_resBio]] = data.frame(variable = variablesSelect[[ind_resBio]])
    
  }
  
  ind_cyto = which(noms_blocks == "cyto")
  
  if(length(variablesSelect[[ind_cyto]]) != 0)
  {
    res_variablesSelect[[ind_cyto]] = data.frame(variable = variablesSelect[[ind_cyto]])
    
  }
  
  
  return(res_variablesSelect)
  
  
}



# Integration 6 blocs norm gene -------------------------------------------


# La fonction matriceCorrelation_comp calcule la matrice de corrélation entre les comp[1]
# composantes de chaque bloc et la matrice de corrélation entre les comp[2] composantes de 
# chaque bloc.

matriceCorrelation_comp <-function(res_block_splsda,
                                   comp = 1:2)
{
  vec_blocksTemp = res_block_splsda$names$blocks
  ind_Y = which(vec_blocksTemp == "Y")
  vec_blocks = vec_blocksTemp[- ind_Y]
  vec_indice_blocks = sapply(1:length(vec_blocks), FUN = function(i){
    res = which(res_block_splsda$names$blocks == vec_blocks[i])
    
    return(res)
  })
  
  # Calcul de la matrice de corrélations entre les comp[1] composantes de chaque bloc
  # et de la matrice de corrélation entre les comp[2] composantes de chaque bloc.
  
  mat_comp1 = sapply(1:length(vec_indice_blocks), FUN = function(i){
    res = res_block_splsda$variates[[vec_indice_blocks[i]]][, comp[1]]
    
    return(res)
  })
  
  colnames(mat_comp1) = vec_blocks
  mat_cor_comp1 = cor(mat_comp1)
  
  mat_comp2 = sapply(1:length(vec_indice_blocks), FUN = function(i){
    res = res_block_splsda$variates[[vec_indice_blocks[i]]][, comp[2]]
    
    return(res)
  })
  
  colnames(mat_comp2) = vec_blocks
  mat_cor_comp2 = cor(mat_comp2)
  
  return(list(mat_cor_comp1 = mat_cor_comp1,
              mat_cor_comp2 = mat_cor_comp2))
}

# La fonction matCor détermine toutes les combinaisons possibles des blocs
# dont nous pouvons superposer les cercles de corrélation. Cette fonction calcule,
# pour chaque variable (variable d'un bloc ou une variable réponse), la corrélation de
# cette variable avec la première composante et la corrélation de cette variable avec la
# deuxième composante.


matCor <-function(res_block_splsda,
                  mat_cor_comp1,
                  mat_cor_comp2,
                  block_Y,
                  comp = 1:2,
                  cutoff_comp)
{
  
  noms_blocTemp = res_block_splsda$names$blocks
  ind_Y = which(noms_blocTemp == "Y")
  noms_bloc = noms_blocTemp[- ind_Y]
  
  # Détermination de toutes les combinaisons possibles de blocs dont nous pouvons superposer
  # les cercles de corrélation.
  
  blockSelect = unique(lapply(1:dim(mat_cor_comp1)[1], FUN = function(i){
    col_mat_cor_comp1_i = mat_cor_comp1[, i]
    col_mat_cor_comp2_i = mat_cor_comp2[, i]
    
    resultat = c()
    
    for(j in 1:length(col_mat_cor_comp1_i))
    {
      cond =   abs(col_mat_cor_comp1_i[j]) > cutoff_comp & abs(col_mat_cor_comp2_i[j]) > cutoff_comp  
      
      if(cond)
      {
        resultat = c(resultat, j)
        
      }
      
    }
    
    return(resultat)
    
  }))
  
  
  liste_noms_blocks = list()
  
  for(i in 1:length(blockSelect))
  {
    blockSelect_i = blockSelect[[i]]
    
    for(k in 1:length(blockSelect_i))
    {
      matComb = combn(blockSelect_i, m = k) 
      
      res = lapply(1:dim(matComb)[2], FUN = function(i){
        comb_i = matComb[, i]
        resultat = noms_bloc[comb_i]
        
        return(resultat)
      })
      
      for(j in 1:length(res))
      {
        liste_noms_blocks = c(liste_noms_blocks, list(res[[j]]))
        
      } # Fin for(j in 1:length(liste_noms_blocks)).
      
      
    } # Fin for(k in 1:length(blockSelect_i)).
    
  } # Fin for(i in 1:length(blockSelect)).
  
  liste_noms_blocks = unique(liste_noms_blocks)
  
  # Pour chaque bloc, calcul des corrélations entre la première
  # composante et les variables sélectionnées et les corrélations entre
  # la deuxième composante et les variables sélectionnées. Pour la réponse,
  # calcul de la corrélation entre les variables de la réponse et la première composante
  # du premier bloc sélectionné et de la corrélation entre les variables de
  # la réponse et la deuxième composante du premier bloc sélectionné.
  
  liste_matCor_comp_var_all = list()
  
  varSelect_comp1 = selectVar(res_block_splsda,
                              comp = comp[1])
  
  varSelect_comp2 = selectVar(res_block_splsda,
                              comp = comp[2])
  
  blockSelectEtReponse = list()
  
  for(i in 1:length(blockSelect))
  {
    blockSelectEtReponse[[i]] = c(blockSelect[[i]], which(res_block_splsda$names$blocks == "Y"))
    
    
  } # Fin for(i in 1:length(blockSelect)).
  
  
  for(i in 1:length(blockSelectEtReponse))
  {
    blockSelectEtReponse_i = blockSelectEtReponse[[i]]
    vec_nom_blockEtReponse = c()
    liste_matCor_comp_var = list()
    
    for(j in 1:length(blockSelectEtReponse_i))
    {
      indice_block_i_j = blockSelectEtReponse_i[j]
      nom_blockEtReponse_i_j = res_block_splsda$names$blocks[indice_block_i_j]
      
      if(nom_blockEtReponse_i_j == "Y")
      {
        if(!is.null(block_Y))
        {
          block_i_j = block_Y
          
          comp1 = res_block_splsda$variates[[blockSelectEtReponse_i[1]]][, comp[1]]
          comp2 = res_block_splsda$variates[[blockSelectEtReponse_i[1]]][, comp[2]]
          
          vecCor_comp1_var = sapply(1:dim(block_Y)[2], FUN = function(j){
            cor(comp1, block_Y[, j], use = "complete.obs")
          })
          
          vecCor_comp2_var = sapply(1:dim(block_Y)[2], FUN = function(j){
            cor(comp2, block_Y[, j], use = "complete.obs")
            
          })
          
          
        }else{
          cat("La réponse n'est pas saisie comme paramètre d'entrée", "\n")
          
        }
        
        
        
      }else{
        comp1 = res_block_splsda$variates[[indice_block_i_j]][, comp[1]]
        comp2 = res_block_splsda$variates[[indice_block_i_j]][, comp[2]]
        
        varSelect_comp1_i_j = varSelect_comp1[[indice_block_i_j]][[1]]
        varSelect_comp2_i_j = varSelect_comp2[[indice_block_i_j]][[1]]
        varSelect_i_j = unique(c(varSelect_comp1_i_j, varSelect_comp2_i_j))
        
        block_i_j = res_block_splsda$X[[indice_block_i_j]][, varSelect_i_j]
        
        if(i == 1)
        {
          vecCor_comp1_var = sapply(1:dim(block_i_j)[2], FUN = function(j){
            cor(comp1, block_i_j[, j], use = "complete.obs")
          })
          
          vecCor_comp2_var = sapply(1:dim(block_i_j)[2], FUN = function(j){
            cor(comp2, block_i_j[, j], use = "complete.obs")
          })
          
        }else{
          signeCor_comp1 = sign(mat_cor_comp1[1, indice_block_i_j])
          
          vecCor_comp1_var = sapply(1:dim(block_i_j)[2], FUN = function(j){
            res = signeCor_comp1*cor(comp1, block_i_j[, j], use = "complete.obs")
            
            return(res)
          })
          
          signeCor_comp2 = sign(mat_cor_comp2[1, indice_block_i_j])
          
          vecCor_comp2_var = sapply(1:dim(block_i_j)[2], FUN = function(j){
            res = signeCor_comp2*cor(comp2, block_i_j[, j], use = "complete.obs")
            
            return(res)
          })
          
        }
        
        
      }
      
      matCor_comp_var = rbind(vecCor_comp1_var,
                              vecCor_comp2_var)
      
      colnames(matCor_comp_var) = colnames(block_i_j)        
      
      
      liste_matCor_comp_var[[j]] = matCor_comp_var
      
      vec_nom_blockEtReponse_i_j = rep(nom_blockEtReponse_i_j, dim(block_i_j)[2])
      vec_nom_blockEtReponse = c(vec_nom_blockEtReponse, vec_nom_blockEtReponse_i_j)
      
    } # Fin for(j in 1:length(blockSelectEtReponse_i)).
    
    matCor_Allcomp_Allvar = t(Reduce(cbind, liste_matCor_comp_var))
    dataframe_Cor_Allcomp_Allvar = data.frame(cbind(rownames(matCor_Allcomp_Allvar),
                                                    vec_nom_blockEtReponse,
                                                    matCor_Allcomp_Allvar))
    
    colnames(dataframe_Cor_Allcomp_Allvar) = c("variable",
                                               "bloc",
                                               "cor_comp1_var",
                                               "cor_comp2_var")
    
    liste_matCor_comp_var_all[[i]] = dataframe_Cor_Allcomp_Allvar
    
    
  } # Fin for(i in 1:length(blocks)).
  
  
  names(liste_matCor_comp_var_all) = sapply(1:length(blockSelectEtReponse), FUN = function(i){
    blockSelectEtReponse_i = blockSelectEtReponse[[i]]
    blocks_i = blockSelectEtReponse_i[blockSelectEtReponse_i != ind_Y]
    nomsBlocks = paste(res_block_splsda$names$blocks[blocks_i], collapse = "_")
    
    return(nomsBlocks)
  })
  
  
  return(list(liste_matCor_comp_var_all = liste_matCor_comp_var_all,
              liste_noms_blocks = liste_noms_blocks))
  
}


# La fonction circleCorZoom permet superposer des cercles de corrélations et de zoomer un rectangle du cercle de corrélations.
# Cette fonction permet de récupérer les variables contenues dans ce rectangle.

circleCorZoom <-function(dataframe_Cor_Allcomp_Allvar,
                         mat_cor_comp1,
                         mat_cor_comp2,
                         nomsBlock,
                         comp = 1:2,
                         cutoff = 0.85,
                         cutoff_comp = 0.8,
                         min.X = -1,
                         max.X = 1,
                         min.Y = -1,
                         max.Y = 1,
                         vec_col = colorRampPalette(brewer.pal(9, "Spectral"))(length(unique(dataframe_Cor_Allcomp_Allvar$bloc))),
                         rad.in = 0.5,
                         cex = 0.7,
                         cex_legend = 0.8,
                         pos = c(1.2, 0),
                         pch = 20,
                         inset = c(-0.25, 0))
{
  
  # Nous vérifions que nous pouvons superposer les cercles de corrélation.
  
  blockSelect = unique(lapply(1:dim(mat_cor_comp1)[1], FUN = function(i){
    col_mat_cor_comp1_i = mat_cor_comp1[, i]
    col_mat_cor_comp2_i = mat_cor_comp2[, i]
    
    resultat = c()
    
    for(j in 1:length(col_mat_cor_comp1_i))
    {
      cond =   abs(col_mat_cor_comp1_i[j]) > cutoff_comp & abs(col_mat_cor_comp2_i[j]) > cutoff_comp  
      
      if(cond)
      {
        resultat = c(resultat, j)
        
      }
      
    }
    
    return(resultat)
    
  }))
  
  indice_block_nomsBlock = sapply(1:length(nomsBlock), FUN = function(i){
    res =  which(colnames(mat_cor_comp1) == nomsBlock[i])
    
    return(res)
  })
  
  boolean = FALSE
  i = 1
  
  while(i <=length(blockSelect) & !boolean)
  {
    blockSelect_i = blockSelect[[i]]
    cond = length(which(blockSelect_i%in%indice_block_nomsBlock == TRUE)) == length(nomsBlock)
    
    if(cond)
    {
      boolean = TRUE
      
    }
    
    i = i + 1
    
  } # Fin while(i <=length(blockSelect) & !boolean).
  
  varSelect = NULL
  
  if(!boolean)
  {
    cat(paste0("Erreur : les blocs : ", paste(nomsBlock, collapse = ", "), " ne peuvent pas être superposés.", "\n"))
    
  }else{
    nomsBlockEtReponse = c(nomsBlock, "Y")
    indice1 = dataframe_Cor_Allcomp_Allvar$bloc%in%nomsBlockEtReponse
    dataframe_Cor_Allcomp_Allvar2 = dataframe_Cor_Allcomp_Allvar[indice1, ]
    matCor_Allcomp_Allvar = apply(dataframe_Cor_Allcomp_Allvar2[, 3:4], 2, as.numeric)
    rownames(matCor_Allcomp_Allvar) = dataframe_Cor_Allcomp_Allvar2$variable
    
    # indice permet de récupérer les variables de chaque bloc fortement corrélées avec soit
    # la première composante ou la deuxième composante dans une partie du cercle de corrélation
    # et de récupérer la variable réponse.
    
    indice = sapply(1:dim(dataframe_Cor_Allcomp_Allvar2)[1], FUN = function(k){
      cor1 = matCor_Allcomp_Allvar[k, 1]
      cor2 = matCor_Allcomp_Allvar[k, 2]
      blockEtReponse_k = dataframe_Cor_Allcomp_Allvar2[k, 2]
      
      if(blockEtReponse_k == "Y")
      {
        cond2 = TRUE
        
      }else{
        cond1 = abs(cor1) > cutoff | abs(cor2) > cutoff
        
        cond2 = cor1 > min.X & cor1 < max.X & cor2 > min.Y & cor2 < max.Y & cond1
        
      }
      
      
      return(cond2)
    })
    
    dataframe_Cor_Allcomp_Allvar2[, 1:2] = apply(dataframe_Cor_Allcomp_Allvar2[, 1:2], 2, as.character)
    dataframe_Cor_Allcomp_Allvar2Select = dataframe_Cor_Allcomp_Allvar2[indice, ]
    varSelectTemp = dataframe_Cor_Allcomp_Allvar2Select$variable
    ind_Y = which(dataframe_Cor_Allcomp_Allvar2Select$bloc == "Y")
    if(length(ind_Y) != 0)
    {
      varSelect = varSelectTemp[- ind_Y]
      
    }else{
      varSelect = varSelectTemp
      
    }
    
    
    matCor_Allcomp_AllvarSelect = matCor_Allcomp_Allvar[indice, , drop = FALSE]
    
    
    
    # Tracé de la superposition des cerles de corrélation.
    circle = list()
    circle[[1]] = ellipse(0, levels = 1, t = 1)
    circle[[2]] = ellipse(0, levels = 1, t = rad.in)
    circle = data.frame(do.call("rbind", circle), "Circle" = c(rep("Main circle", 100), rep("Inner circle", 100)))
    
    MainCircle = circle[grep("Main circle", circle[, 3]), ]
    InnerCircle = circle[grep("Inner circle", circle[, 3]), ]
    
    plot(MainCircle[, 1], MainCircle[, 2],
         type = "l",
         xlab = paste0("composante ", comp[1]),
         ylab = paste0("composante ", comp[2]))
    
    points(InnerCircle[, 1], InnerCircle[, 2],
           type = "l")
    
    if(dim(matCor_Allcomp_AllvarSelect)[1] != 0)
    {
      
      nom_blockEtReponseSelect = unique(dataframe_Cor_Allcomp_Allvar2Select$bloc)
      indice_blockEtReponseSelect = sapply(1:length(nom_blockEtReponseSelect), FUN = function(i){
        ind = which(colnames(mat_cor_comp1) == nom_blockEtReponseSelect[i])
        
        if(length(ind) != 0)
        {
          res = ind
          
        }else{
          res = dim(mat_cor_comp1)[1] + 1
          
        }
        
        return(res)
      })
      vec_colSelect = vec_col[indice_blockEtReponseSelect]
      
      if(length(nom_blockEtReponseSelect) == 1)
      {
        points(matCor_Allcomp_AllvarSelect[, 1], matCor_Allcomp_AllvarSelect[, 2],
               col  = NULL)
        
        text(matCor_Allcomp_AllvarSelect[, 1], matCor_Allcomp_AllvarSelect[, 2],
             labels = rownames(matCor_Allcomp_AllvarSelect),
             cex = cex,
             col = vec_colSelect[1])
        
      }else{
        
        nbVarSelect_bloc = cumsum(sapply(1:length(nom_blockEtReponseSelect), FUN = function(j){
          res = length(which(dataframe_Cor_Allcomp_Allvar2Select$bloc == nom_blockEtReponseSelect[j]))
          
          return(res)
        }))
        
        for(i in 1:length(nbVarSelect_bloc))
        {
          if(i == 1)
          {
            indice1 = 1:nbVarSelect_bloc[1]
            
            if(length(indice1) != 0)
            {
              matCor_Allcomp_AllvarSelect2 = matCor_Allcomp_AllvarSelect[indice1, , drop = FALSE]
              
              points(matCor_Allcomp_AllvarSelect2[, 1], matCor_Allcomp_AllvarSelect2[, 2],
                     col  = NULL)
              
              text(matCor_Allcomp_AllvarSelect2[, 1], matCor_Allcomp_AllvarSelect2[, 2],
                   labels = rownames(matCor_Allcomp_AllvarSelect2),
                   cex = cex,
                   col = rep(vec_colSelect[i], dim(matCor_Allcomp_AllvarSelect2)[1]))
              
            }else{
              cat(paste0("Il n'y a de variables dans cette zone du cercle de corrélation pour le bloc ", nom_blockEtReponseSelect[i]), "\n")
              
              
            }
            
            
          }else{
            indice2 = (nbVarSelect_bloc[i - 1] + 1):nbVarSelect_bloc[i]
            
            if(length(indice2) != 0)
            {
              matCor_Allcomp_AllvarSelect2 = matCor_Allcomp_AllvarSelect[indice2, , drop = FALSE]
              
              points(matCor_Allcomp_AllvarSelect2[, 1], matCor_Allcomp_AllvarSelect2[, 2],
                     col  = NULL)
              
              text(matCor_Allcomp_AllvarSelect2[, 1], matCor_Allcomp_AllvarSelect2[, 2],
                   labels = rownames(matCor_Allcomp_AllvarSelect2),
                   cex = cex,
                   col = rep(vec_colSelect[i], dim(matCor_Allcomp_AllvarSelect2)[1]))
              
              
            }else{
              cat(paste0("Il n'y a de variables dans cette zone du cercle de corrélation pour le bloc ", nom_blockEtReponseSelect[i]), "\n")
              
            }
            
            
          }
          
        }
        
      }
      
      par(xpd = TRUE)
      legend(x = pos[1], y = pos[2],
             legend = nom_blockEtReponseSelect,
             pch = pch,
             col = vec_colSelect,
             cex = cex_legend,
             inset = inset)
      
      
    }else{
      cat("Il n'y a de variables dans cette zone du cercle de corrélation", "\n")
      
    }
    
    
  }
  
  return(varSelect)
  
  
}

# La fonction networkVariable permet de tracer un réseau pour les variables de
# certains blocs et une variable réponse.

networkVariable <-function(dataframe_Cor_Allcomp_Allvar,
                           vec_Var,
                           nomVar_block_Y,
                           comp = 1:2,
                           cutoff = 0
)
{
  # Pour chaque variable de vec_Var et la variable réponse, nous récupérons les noms des blocs,
  # les corrélations entre la variable et chaque composante.
  
  ind = which(dataframe_Cor_Allcomp_Allvar$variable%in%vec_Var == TRUE)
  ind_Var_block_Y = which(dataframe_Cor_Allcomp_Allvar$variable == nomVar_block_Y)
  
  dataframe_Cor_Allcomp_AllvarSelect = dataframe_Cor_Allcomp_Allvar[c(ind, ind_Var_block_Y), ]
  nomBlocs = unique(dataframe_Cor_Allcomp_AllvarSelect$bloc)
  
  if((length(ind) == length(vec_Var)) & length(ind_Var_block_Y) == 1)
  {
    if(length(nomBlocs) == 1)
    {
      cat("Erreur : il n'y a des variables d'un seul bloc dans listeVar", "\n")
      
    }else{
      
      liste_Cor_Allcomp_Allvar_bloc = lapply(1:length(nomBlocs), FUN = function(i){
        nomBloc_i = nomBlocs[i]
        indice = which(dataframe_Cor_Allcomp_AllvarSelect$bloc == nomBloc_i)
        resTemp = dataframe_Cor_Allcomp_AllvarSelect[indice, ]
        
        if(length(indice) == 1)
        {
          res = resTemp
          
        }else{
          res = resTemp
          
        }
        
        rownames(res) = resTemp$variable
        
        return(res)
      })
      names(liste_Cor_Allcomp_Allvar_bloc) = nomBlocs
      
      coord  = lapply(1:length(liste_Cor_Allcomp_Allvar_bloc), FUN = function(i){
        resTemp = liste_Cor_Allcomp_Allvar_bloc[[i]]
        res = apply(resTemp[, comp + 2, drop = FALSE], 2, as.numeric)
        
        if(dim(resTemp)[1] == 1)
        {
          res = matrix(res, nrow = 1)
          
        }
        
        rownames(res) = resTemp$variable
        
        return(res)
      })
      
      
      l = 1
      M_block = list()
      node.X1 = node.X2 = w = NULL
      
      # Calcul de la similarité.
      
      for(j in 1:(length(nomBlocs) - 1))
      {
        for(k in (j + 1):length(nomBlocs))
        {
          M_block[[l]] = coord[[j]][, drop = FALSE] %*% t(coord[[k]][, drop = FALSE])
          
          X1 = rownames(coord[[j]])
          X2 = rownames(coord[[k]])
          
          rep.X1 = rep(X1, each = length(X2))
          rep.X2 = rep(X2, length(X1))
          
          node.X1= c(node.X1, rep.X1)
          node.X2 = c(node.X2, rep.X2)
          
          w = c(w, as.vector(t(M_block[[l]])))
          
          l = l + 1
          
        } # Fin for(k in (j + 1):length(blocks)).
        
      } # Fin for(j in 1:(length(blocks) - 1)).
      
      # nodes
      group = NULL
      temp = lapply(1:length(liste_Cor_Allcomp_Allvar_bloc), function(i){
        res = liste_Cor_Allcomp_Allvar_bloc[[i]]$variable
        
        return(res)
      })
      names(temp) = names(liste_Cor_Allcomp_Allvar_bloc)
      
      for (i in 1:length(temp))
      {
        group = c(group, rep(names(temp)[i], length(temp[[i]])))
        
      } # Fin for (i in 1:length(temp)).
      
      nodes = data.frame(name = unlist(temp),
                         group = group)
      
      # gR
      relations = data.frame(from = node.X1,
                             to = node.X2,
                             weight = w)
      
      idx = (abs(w) >= cutoff)
      relations = relations[idx, , drop = FALSE]
      
      gR = graph.data.frame(relations,
                            directed = FALSE,
                            vertices = nodes)
      
      gR = delete.vertices(gR, which(degree(gR) == 0))
      
      res = list(gR = gR)
      
      l = 1
      for (i in 1:(length(nomBlocs)-1))
      {
        for (j in (i + 1):length(nomBlocs))
        {
          res[paste("M", names(liste_Cor_Allcomp_Allvar_bloc)[i], names(liste_Cor_Allcomp_Allvar_bloc)[j], sep="_")] = list(M_block[[l]])
          l = l + 1
          
        } # Fin for (j in (i + 1):length(blocks)).
        
      } # Fin for (i in 1:(length(blocks)-1).
      
      res$cutoff = cutoff
      
      return(res)
      
      
    }
    
    
    # Fin if((length(which(ind))= length(listeVar)) & length(ind_Var_block_Y) == 1).
  }else{
    if(length(ind) != length(vec_Var))
    {
      cat("Erreur  : les variables de vec_Var ne sont pas contenues dans dataframe_Cor_Allcomp_Allvar$variable.", "\n")
      
    }
    
    if(length(ind_Var_block_Y) == 1)
    {
      cat(paste0("Erreur  : la variable réponse ", nomVar_block_Y," n'est pas contenue dans dataframe_Cor_Allcomp_Allvar$variable."), "\n")
      
      
    }
    
  }
  
  
}



# Integration 7 blocs norm gene -------------------------------------------



# La fonction compute_cor_comp_var calcule, pour chaque variable d'un bloc, les corrélations
# entre la variable et les composantes sélectionnées par comp.

compute_cor_comp_var <-function(res_block_splsda,
                                comp = c(1:2))
{
  vec_blocksTemp = res_block_splsda$names$blocks
  ind_Y = which(vec_blocksTemp == "Y")
  vec_blocks1 = vec_blocksTemp[- ind_Y]
  vec_indice_blocks = sapply(1:length(vec_blocks1), FUN = function(i){
    res = which(res_block_splsda$names$blocks == vec_blocks1[i])
    
    return(res)
  })
  
  # Calcul, pour chaque variable d'un bloc, des corrélations de cette variable avec
  # les composantes dont les indices sont indiqués dans le vecteur comp.
  
  vec_blocks2 = c()
  liste_cor_comp_var_global = list()
  
  for(i in 1:length(vec_indice_blocks))
  {
    indice_blocks_i = vec_indice_blocks[i]
    block_i = res_block_splsda$names$blocks[indice_blocks_i]
    
    vec_varSelect_i = c()
    liste_comp_i = list()
    
    for(j in 1:length(comp))
    {
      indice_comp_j = comp[j]
      liste_comp_i[[j]] = res_block_splsda$variates[[indice_blocks_i]][, indice_comp_j]
      vec_varSelect_comp_j = selectVar(res_block_splsda,
                                       comp = indice_comp_j)[[indice_blocks_i]][[1]]
      
      vec_varSelect_i = c(vec_varSelect_i, vec_varSelect_comp_j)
      
    } # Fin for(j in 1:length(comp)).
    
    vec_varSelect_i = unique(vec_varSelect_i)
    
    mat_block_i = res_block_splsda$X[[indice_blocks_i]][, vec_varSelect_i]
    
    liste_cor_comp_var = list()
    
    for(j in 1:length(liste_comp_i))
    {
      liste_cor_comp_var[[j]] = sapply(1:dim(mat_block_i)[2], FUN = function(k){
        cor(liste_comp_i[[j]], mat_block_i[, k])
      })
      
    } # Fin for(j in 1:length(liste_comp_i)).
    
    mat_cor_Allcomp_Allvar = Reduce(cbind, liste_cor_comp_var)
    rownames(mat_cor_Allcomp_Allvar) = colnames(mat_block_i)
    
    vec_blocks_i = rep(block_i, dim(mat_block_i)[2])
    vec_blocks2 = c(vec_blocks2, vec_blocks_i)
    
    liste_cor_comp_var_global[[i]] = mat_cor_Allcomp_Allvar
    
    
  } # Fin for(i in 1:length(blocks)).
  
  
  mat_cor_comp_var_global = Reduce(rbind, liste_cor_comp_var_global)
  dataframe_cor_comp_var_global = data.frame(cbind(rownames(mat_cor_comp_var_global),
                                                   vec_blocks2,
                                                   mat_cor_comp_var_global))
  
  colnames(dataframe_cor_comp_var_global) = c("variable",
                                              "bloc",
                                              paste0("cor_var_comp", comp))
  
  dataframe_cor_comp_var_global[, 1:2] = apply(dataframe_cor_comp_var_global[, 1:2], 2, as.character)
  dataframe_cor_comp_var_global[, 3:dim(dataframe_cor_comp_var_global)[2]] = apply(dataframe_cor_comp_var_global[, 3:dim(dataframe_cor_comp_var_global)[2]], 2, as.numeric)
  
  
  return(dataframe_cor_comp_var_global)
  
}

# La fonction composanteColin renvoie une liste. Le ième élément de cette liste 
# contient les indices des blocs tels que, pour chaque paire de ces blocs, la 
# première composante du bloc 1 est fortement corrélée à la première composante
# du bloc2 en valeur absolue et la deuxième composante du bloc 1 est fortement 
# corrélée à la deuxième composante du bloc2 en valeur absolue.

composanteColin <-function(mat_cor_comp1,
                           mat_cor_comp2,
                           cutoff_comp)
{
  res = list()
  
  index = 1:dim(mat_cor_comp1)[1]
  i = 1
  compt = 1
  
  while(length(index) != 0 & compt <= dim(mat_cor_comp1)[1])
  {
    index_i = index[i]
    res[[compt]] = c(index_i)
    
    index2 = index[-i]
    
    if(length(index2) != 1)
    {
      
      for(j in index2)
      {
        if(length(res[[compt]]) == 1)
        {
          if(abs(mat_cor_comp1[j, index_i]) > cutoff_comp & abs(mat_cor_comp2[j, index_i]) > cutoff_comp)
          {
            res[[compt]] = c(res[[compt]], j)
            index = index[- c(which(index == j))]
            
          }
          
        }else{
          
          indice = sapply(1:length(res[[compt]]), FUN = function(k){
            index_k = res[[compt]][k]
            cond = abs(mat_cor_comp1[index_k, j]) > cutoff_comp & abs(mat_cor_comp2[index_k, j]) > cutoff_comp
            
            return(cond)
          })
          
          if(all(indice))
          {
            res[[compt]] = c(res[[compt]], j)
            index = index[- c(which(index == j))]
            
          }
          
          
        }
        
      }
      
    }else{
      res[[compt]] = index_i
      
    }
    
    index = index[- c(which(index == index_i))]
    compt = compt + 1
    
  } # Fin for(i in 1:dim(mat_cor_comp1)[2]).
  
  if(length(index) != 0)
  {
    for(i in 1:length(index))
    {
      res = c(res, list(index[i]))
      
    }
    
  }
  
  
  return(res)
}


# La fonction compute_blockSelect permet de déterminer toutes les combinaisons
# possibles des blocs qui peuvent être superposés dans le cercle de corrélations
# et dont les variables peuvent être présentes dans le réseau.

compute_blockSelect <-function(mat_cor_comp1,
                               mat_cor_comp2,
                               cutoff_comp)
{
  
  liste_vec_indice_blockSelect = composanteColin(mat_cor_comp1,
                                                 mat_cor_comp2,
                                                 cutoff_comp)
  
  vec_blocks = colnames(mat_cor_comp1)
  liste_vec_blocks = list()
  
  # liste_blocks est une liste contenant toutes les combinaisons possibles de blocs dont
  # on peut superposer les cercles de corrélations et dont les variables peuvent être
  # présentes dans le réseau.
  
  for(i in 1:length(liste_vec_indice_blockSelect))
  {
    vec_indice_blockSelect_i = liste_vec_indice_blockSelect[[i]]
    
    for(k in 1:length(vec_indice_blockSelect_i))
    {
      matComb = combn(vec_indice_blockSelect_i, m = k) 
      
      liste_vec_blocks_i = lapply(1:dim(matComb)[2], FUN = function(i){
        comb_i = matComb[, i]
        resultat = vec_blocks[comb_i]
        
        return(resultat)
      })
      
      for(j in 1:length(liste_vec_blocks_i))
      {
        liste_vec_blocks = c(liste_vec_blocks, list(liste_vec_blocks_i[[j]]))
        
      } # Fin for(j in 1:length(liste_noms_blocks)).
      
      
    } # Fin for(k in 1:length(blockSelect_i)).
    
  } # Fin for(i in 1:length(blockSelect)).
  
  liste_vec_blocks = unique(liste_vec_blocks)
  
  
  
  return(list(liste_vec_indice_blockSelect = liste_vec_indice_blockSelect,
              liste_vec_blocks = liste_vec_blocks))
  
}

# La fonction matCorEtBlockSelect permet de calculer la matrice de
# corrélation entre les comp[1] composantes de chaque bloc et la 
# matrice de corrélation entre les comp[2] composantes de chaque 
# bloc. Elle permet aussi de calculer, pour chaque variable d'un bloc,
# les corrélations entre les composantes comp[1] et comp[2] et cette variable 
# et de calculer toutes les combinaisons possibles des blocs pour lesquels nous
# pouvons superposer les cercles de corrélations.

matCorEtBlockSelect <-function(res_block_splsda,
                               cutoff_comp,
                               comp)
{
  liste_mat_cor_comp = matriceCorrelation_comp(res_block_splsda = res_block_splsda,
                                               comp = comp)  
  
  mat_cor_comp1 = liste_mat_cor_comp$mat_cor_comp1
  mat_cor_comp2 = liste_mat_cor_comp$mat_cor_comp2
  
  
  dataframe_cor_comp_var_global = compute_cor_comp_var(res_block_splsda = res_block_splsda,
                                                       comp = comp)
  
  
  liste_blockSelect = compute_blockSelect(mat_cor_comp1 = mat_cor_comp1,
                                          mat_cor_comp2 = mat_cor_comp2,
                                          cutoff_comp = cutoff_comp)
  
  liste_vec_indice_blockSelect = liste_blockSelect$liste_vec_indice_blockSelect
  liste_vec_blocks = liste_blockSelect$liste_vec_blocks
  
  return(list(mat_cor_comp1 = mat_cor_comp1,
              mat_cor_comp2 = mat_cor_comp2,
              dataframe_cor_comp_var_global = dataframe_cor_comp_var_global,
              liste_vec_indice_blockSelect = liste_vec_indice_blockSelect,
              liste_vec_blocks = liste_vec_blocks))
  
}

# La fonction addVariablesReponses permet de calculer, pour chaque
# variable réponse, la corrélation entre cette variable et la
# comp[1] composante et la corrélation entre cette variable et la comp[2]
# composante pour chaque groupe de blocs.

addVariablesReponses <-function(res_block_splsda,
                                dataframe_cor_comp_var_global,
                                liste_vec_indice_blockSelect,
                                mat_block_Y)
{
  # On récupère les indices des composantes utilisées pour calculer les
  # corrélations entre les variables des blocs et les composantes.
  
  comp = as.numeric(sapply(3:4, FUN = function(i){
    col_i = colnames(dataframe_cor_comp_var_global)[i]
    ch = strsplit(col_i, split = "_")[[1]]
    resTemp = ch[length(ch)]
    res = substring(resTemp, nchar(resTemp), nchar(resTemp))
    
    return(res)
  }))
  
  # On calcule, pour chaque variable réponse, pour chaque composante du premier
  # bloc du groupe de blocs, la corrélation entre cette variable réponse et la composante.
  
  liste_dataframe_cor_allcomp_varReponses = list()
  vec_groupe_blocks = c()
  
  for(i in 1:length(liste_vec_indice_blockSelect))
  {
    vec_indice_blockSelect_i = liste_vec_indice_blockSelect[[i]]
    indice_first_block_i = vec_indice_blockSelect_i[1]
    
    liste_comp_i = list()
    liste_cor_comp_var = list()
    
    for(j in 1:length(comp))
    {
      indice_comp_j = comp[j]
      comp_j = res_block_splsda$variates[[indice_first_block_i]][, indice_comp_j]
      
      liste_cor_comp_var[[j]] = sapply(1:dim(mat_block_Y)[2], FUN = function(k){
        cor(comp_j, mat_block_Y[, k])
      })
      
    } # Fin for(j in 1:length(liste_comp_i)).
    
    mat_cor_allcomp_varReponses = Reduce(cbind, liste_cor_comp_var)
    rownames(mat_cor_allcomp_varReponses) = colnames(mat_block_Y)
    
    dataframe_allcomp_varReponses = data.frame(colnames(mat_block_Y),
                                               rep("Y", dim(mat_block_Y)[2]),
                                               mat_cor_allcomp_varReponses)
    
    colnames(dataframe_allcomp_varReponses) = c("variable",
                                                "bloc",
                                                paste0("cor_var_comp", comp))
    
    
    liste_dataframe_cor_allcomp_varReponses[[i]] = dataframe_allcomp_varReponses
    
    groupe_blocks_i = res_block_splsda$names$blocks[vec_indice_blockSelect_i]
    vec_groupe_blocks = c(vec_groupe_blocks, paste(groupe_blocks_i, collapse = "-"))
    
    
  } # Fin for(in in 1:length(blockSelect)).
  
  names(liste_dataframe_cor_allcomp_varReponses) = vec_groupe_blocks
  
  liste_dataframe_cor_comp_var_global = list()
  
  for(i in 1:length(liste_vec_indice_blockSelect))
  {
    vec_indice_blockSelect_i = liste_vec_indice_blockSelect[[i]]
    groupe_blocks_i = res_block_splsda$names$blocks[vec_indice_blockSelect_i]
    indice_i = which(dataframe_cor_comp_var_global$bloc%in%groupe_blocks_i == TRUE)
    dataframe_cor_comp_var_global_indice_i = dataframe_cor_comp_var_global[indice_i, ]
    
    dataframe_cor_comp_varBlockEtVarRep_global = rbind(dataframe_cor_comp_var_global_indice_i,
                                                       liste_dataframe_cor_allcomp_varReponses[[i]])
    
    dataframe_cor_comp_varBlockEtVarRep_global[, 1:2] = apply(dataframe_cor_comp_varBlockEtVarRep_global[, 1:2], 2, as.character)
    dataframe_cor_comp_varBlockEtVarRep_global[, 3:dim(dataframe_cor_comp_varBlockEtVarRep_global)[2]] = apply(dataframe_cor_comp_varBlockEtVarRep_global[, 3:dim(dataframe_cor_comp_varBlockEtVarRep_global)[2]], 2, as.numeric)
    
    
    liste_dataframe_cor_comp_var_global[[i]] = dataframe_cor_comp_varBlockEtVarRep_global
    
    
  } # Fin for(i in 1:length(indice_blockSelect)).
  
  names(liste_dataframe_cor_comp_var_global) = vec_groupe_blocks
  
  
  return(liste_dataframe_cor_comp_var_global)
  
}

# La fonction addVariablesReponsesModified permet de calculer, pour chaque
# variable réponse, la corrélation entre cette variable et la
# comp[1] composante et la corrélation entre cette variable et la comp[2]
# composante pour chaque groupe de blocs.

addVariablesReponsesModified <-function(res_block_splsda,
                                        dataframe_cor_comp_var_global,
                                        liste_vec_indice_blockSelect,
                                        mat_block_Y)
{
  # On récupère les indices des composantes utilisées pour calculer les
  # corrélations entre les variables des blocs et les composantes.
  
  comp = as.numeric(sapply(3:4, FUN = function(i){
    col_i = colnames(dataframe_cor_comp_var_global)[i]
    ch = strsplit(col_i, split = "_")[[1]]
    resTemp = ch[length(ch)]
    res = substring(resTemp, nchar(resTemp), nchar(resTemp))
    
    return(res)
  }))
  
  # On calcule, pour chaque variable réponse, pour chaque composante du premier
  # bloc du groupe de blocs, la corrélation entre cette variable réponse et la composante.
  
  liste_dataframe_cor_allcomp_varReponses = list()
  vec_groupe_blocks = c()
  
  for(i in 1:length(liste_vec_indice_blockSelect))
  {
    vec_indice_blockSelect_i = liste_vec_indice_blockSelect[[i]]
    indice_first_block_i = vec_indice_blockSelect_i[1]
    
    liste_comp_i = list()
    liste_cor_comp_var = list()
    
    for(j in 1:length(comp))
    {
      indice_comp_j = comp[j]
      comp_j = res_block_splsda$variates[[indice_first_block_i]][, indice_comp_j]
      
      liste_cor_comp_var[[j]] = sapply(1:dim(mat_block_Y)[2], FUN = function(k){
        resTemp = comp_j%*%mat_block_Y[, k]
        res = resTemp/(norm(comp_j - mean(comp_j), "2")*norm(mat_block_Y[, k], "2"))
      })
      
    } # Fin for(j in 1:length(liste_comp_i)).
    
    mat_cor_allcomp_varReponses = Reduce(cbind, liste_cor_comp_var)
    rownames(mat_cor_allcomp_varReponses) = colnames(mat_block_Y)
    
    dataframe_allcomp_varReponses = data.frame(colnames(mat_block_Y),
                                               rep("Y", dim(mat_block_Y)[2]),
                                               mat_cor_allcomp_varReponses)
    
    colnames(dataframe_allcomp_varReponses) = c("variable",
                                                "bloc",
                                                paste0("cor_var_comp", comp))
    
    
    liste_dataframe_cor_allcomp_varReponses[[i]] = dataframe_allcomp_varReponses
    
    groupe_blocks_i = res_block_splsda$names$blocks[vec_indice_blockSelect_i]
    vec_groupe_blocks = c(vec_groupe_blocks, paste(groupe_blocks_i, collapse = "-"))
    
    
  } # Fin for(in in 1:length(blockSelect)).
  
  names(liste_dataframe_cor_allcomp_varReponses) = vec_groupe_blocks
  
  liste_dataframe_cor_comp_var_global = list()
  
  for(i in 1:length(liste_vec_indice_blockSelect))
  {
    vec_indice_blockSelect_i = liste_vec_indice_blockSelect[[i]]
    groupe_blocks_i = res_block_splsda$names$blocks[vec_indice_blockSelect_i]
    indice_i = which(dataframe_cor_comp_var_global$bloc%in%groupe_blocks_i == TRUE)
    dataframe_cor_comp_var_global_indice_i = dataframe_cor_comp_var_global[indice_i, ]
    
    dataframe_cor_comp_varBlockEtVarRep_global = rbind(dataframe_cor_comp_var_global_indice_i,
                                                       liste_dataframe_cor_allcomp_varReponses[[i]])
    
    dataframe_cor_comp_varBlockEtVarRep_global[, 1:2] = apply(dataframe_cor_comp_varBlockEtVarRep_global[, 1:2], 2, as.character)
    dataframe_cor_comp_varBlockEtVarRep_global[, 3:dim(dataframe_cor_comp_varBlockEtVarRep_global)[2]] = apply(dataframe_cor_comp_varBlockEtVarRep_global[, 3:dim(dataframe_cor_comp_varBlockEtVarRep_global)[2]], 2, as.numeric)
    
    
    liste_dataframe_cor_comp_var_global[[i]] = dataframe_cor_comp_varBlockEtVarRep_global
    
    
  } # Fin for(i in 1:length(indice_blockSelect)).
  
  names(liste_dataframe_cor_comp_var_global) = vec_groupe_blocks
  
  
  return(liste_dataframe_cor_comp_var_global)
  
}


# La fonction circleCor permet de zoomer sur un rectangle du cercle de corrélations
# et de récupérer les variables des blocs dans cette partie zoomée du cercle de 
# corrélations.

circleCor <-function(liste_dataframe_cor_comp_var_global,
                     liste_vec_indice_blockSelect,
                     mat_cor_comp1,
                     mat_cor_comp2,
                     vec_blocks,
                     nomsVarReponses,
                     cutoff = 0.85,
                     min.X = -1,
                     max.X = 1,
                     min.Y = -1,
                     max.Y = 1,
                     vec_col = colorRampPalette(brewer.pal(9, "Spectral"))(dim(mat_cor_comp1)[1] + 1),
                     rad.in = 0.5,
                     cex = 0.7,
                     cex_legend = 0.8,
                     pos = c(1.2, 0),
                     pch = 20,
                     inset = c(-0.25, 0))
{
  
  # On vérifie que nous pouvons superposer les cercles de corrélation.
  
  vec_indice_blocks = sapply(1:length(vec_blocks), FUN = function(i){
    res =  which(colnames(mat_cor_comp1) == vec_blocks[i])
    
    return(res)
  })
  
  boolean = FALSE
  i = 1
  
  while(i <= length(liste_vec_indice_blockSelect) & !boolean)
  {
    vec_indice_blockSelect_i = liste_vec_indice_blockSelect[[i]]
    cond = length(which(vec_indice_blockSelect_i%in%vec_indice_blocks == TRUE)) == length(vec_blocks)
    
    if(cond)
    {
      boolean = TRUE
      
    }
    
    i = i + 1
    
  } # Fin  while(i <= length(liste_vec_indice_blockSelect) & !boolean).
  
  varSelect = NULL
  
  if(!boolean)
  {
    stop(paste0("The blocks : ", paste(vec_blocks, collapse = ", "), " can not be superimposed."))
    
  }else{
    # On récupère le groupe de blocs auxquels appartient vec_blocks.
    
    indice_nomsBlock = sapply(1:length(liste_dataframe_cor_comp_var_global), FUN = function(i){
      name_iTemp = names(liste_dataframe_cor_comp_var_global)[i]
      name_i = strsplit(name_iTemp, split = "-")[[1]]
      res = all(vec_blocks%in%name_i)
      
      return(res)
    })
    
    dataframe_cor_comp_var_global = liste_dataframe_cor_comp_var_global[[which(indice_nomsBlock == TRUE)]]
    dataframe_cor_comp_var_globalTemp1 = dataframe_cor_comp_var_global
    
    indice_nomsVarReponses = all(nomsVarReponses%in%dataframe_cor_comp_var_globalTemp1$variable)
    
    if(!indice_nomsVarReponses)
    {
      stop("All the correlations between the response variables and the first component and the correlations 
            between the responses variables and the second component have not been computed.")
      
    }else{
      comp = as.numeric(sapply(3:4, FUN = function(i){
        col_i = colnames(dataframe_cor_comp_var_globalTemp1)[i]
        ch = strsplit(col_i, split = "_")[[1]]
        resTemp = ch[length(ch)]
        res = substring(resTemp, nchar(resTemp), nchar(resTemp))
        
        return(res)
      }))
      
      
      mat_cor_comp_var_globalTemp1 = t(sapply(1:dim(dataframe_cor_comp_var_globalTemp1)[1], FUN = function(i){
        dataframe_cor_comp_var_globalTemp1_i = dataframe_cor_comp_var_globalTemp1[i, ] 
        block_i = dataframe_cor_comp_var_globalTemp1_i$bloc
        cor1  = dataframe_cor_comp_var_globalTemp1_i[paste0("cor_var_comp", comp[1])]
        cor2  = dataframe_cor_comp_var_globalTemp1_i[paste0("cor_var_comp", comp[2])]
        
        if(block_i == "Y")
        {
          cor1_sign = cor1
          cor2_sign = cor2
          
        }else{
          indice_block_comp1_i =  which(colnames(mat_cor_comp1) == block_i)   
          indice_block_comp2_i =  which(colnames(mat_cor_comp2) == block_i)
          
          cor1_sign = sign(mat_cor_comp1[vec_indice_blocks[1], indice_block_comp1_i])*cor1
          cor2_sign = sign(mat_cor_comp2[vec_indice_blocks[1], indice_block_comp2_i])*cor2
          
        }
        
        res = c(dataframe_cor_comp_var_globalTemp1_i[1:2], cor1_sign, cor2_sign)
        
        return(res)
      }))
      
      dataframe_cor_comp_var_globalTemp2 = as.data.frame(mat_cor_comp_var_globalTemp1)
      dataframe_cor_comp_var_globalTemp2[, 1:2] = apply(dataframe_cor_comp_var_globalTemp2[, 1:2], 2, as.character)
      colnames(dataframe_cor_comp_var_globalTemp2) = colnames(dataframe_cor_comp_var_global)
      
      # Pour les variables de vec_blocks et les variables réponses nomsVarReponses, on récupère les corrélations
      # entre ces variables et les composantes. 
      
      indice1 = sapply(1:dim(dataframe_cor_comp_var_globalTemp2)[1], FUN = function(i){
        block_i = dataframe_cor_comp_var_globalTemp2$bloc[i]
        
        if(block_i == "Y")
        {
          res = dataframe_cor_comp_var_globalTemp2$variable[i]%in%nomsVarReponses
          
        }else{
          res = block_i%in%vec_blocks
          
        }
        
        return(res)
      })
      dataframe_cor_comp_var_global2 = dataframe_cor_comp_var_globalTemp2[indice1, ]
      mat_cor_comp_var_global2 = apply(dataframe_cor_comp_var_global2[ , 3:4], 2, as.numeric)
      rownames(mat_cor_comp_var_global2) = dataframe_cor_comp_var_global2$variable
      
      # indice permet de récupérer les variables de chaque bloc fortement corrélées avec soit
      # la première composante ou la deuxième composante dans un rectangle du cercle de corrélations
      # et de récupérer les variables réponses.
      
      indice2 = sapply(1:dim(dataframe_cor_comp_var_global2)[1], FUN = function(k){
        cor1 = mat_cor_comp_var_global2[k, 1]
        cor2 = mat_cor_comp_var_global2[k, 2]
        blockEtReponse_k = dataframe_cor_comp_var_global2[k, 2]
        
        if(blockEtReponse_k == "Y")
        {
          cond2 = TRUE
          
        }else{
          cond1 = abs(cor1) > cutoff | abs(cor2) > cutoff
          
          cond2 = cor1 > min.X & cor1 < max.X & cor2 > min.Y & cor2 < max.Y & cond1
          
        }
        
        return(cond2)
      })
      
      dataframe_cor_comp_var_global2Select = dataframe_cor_comp_var_global2[indice2, ]
      varSelectTemp = dataframe_cor_comp_var_global2Select$variable
      ind_Y = which(dataframe_cor_comp_var_global2Select$bloc == "Y")
      
      if(length(ind_Y) != 0)
      {
        varSelect = varSelectTemp[- ind_Y]
        
      }else{
        varSelect = varSelectTemp
        
      }
      
      mat_cor_comp_var_global2Select = mat_cor_comp_var_global2[indice2, , drop = FALSE]
      
      
      # Tracé de la superposition des cerles de corrélation.
      circle = list()
      circle[[1]] = ellipse(0, levels = 1, t = 1)
      circle[[2]] = ellipse(0, levels = 1, t = rad.in)
      circle = data.frame(do.call("rbind", circle), "Circle" = c(rep("Main circle", 100), rep("Inner circle", 100)))
      
      MainCircle = circle[grep("Main circle", circle[, 3]), ]
      InnerCircle = circle[grep("Inner circle", circle[, 3]), ]
      
      plot(MainCircle[, 1], MainCircle[, 2],
           type = "l",
           xlab = paste0("composante ", comp[1]),
           ylab = paste0("composante ", comp[2]))
      
      points(InnerCircle[, 1], InnerCircle[, 2],
             type = "l")
      
      if(dim(mat_cor_comp_var_global2Select)[1] != 0)
      {
        vec_blockEtReponseSelect = unique(dataframe_cor_comp_var_global2Select$bloc)
        indice_blockEtReponseSelect = sapply(1:length(vec_blockEtReponseSelect), FUN = function(i){
          ind = which(colnames(mat_cor_comp1) == vec_blockEtReponseSelect[i])
          
          if(length(ind) != 0)
          {
            res = ind
            
          }else{
            res = dim(mat_cor_comp1)[1] + 1
            
          }
          
          return(res)
        })
        vec_colSelect = vec_col[indice_blockEtReponseSelect]
        
        if(length(vec_blockEtReponseSelect) == 1)
        {
          points(mat_cor_comp_var_global2Select[, 1], mat_cor_comp_var_global2Select[, 2],
                 col  = NULL)
          
          text(mat_cor_comp_var_global2Select[, 1], mat_cor_comp_var_global2Select[, 2],
               labels = rownames(mat_cor_comp_var_global2Select),
               cex = cex,
               col = vec_colSelect[1])
          
        }else{
          
          nbVarSelect_bloc = cumsum(sapply(1:length(vec_blockEtReponseSelect), FUN = function(j){
            res = length(which(dataframe_cor_comp_var_global2Select$bloc == vec_blockEtReponseSelect[j]))
            
            return(res)
          }))
          
          for(i in 1:length(nbVarSelect_bloc))
          {
            if(i == 1)
            {
              indice_nbVar1 = 1:nbVarSelect_bloc[1]
              
              mat_cor_comp_var_global2Select2 = mat_cor_comp_var_global2Select[indice_nbVar1, , drop = FALSE]
              
              points(mat_cor_comp_var_global2Select2[, 1], mat_cor_comp_var_global2Select2[, 2],
                     col  = NULL)
              
              text(mat_cor_comp_var_global2Select2[, 1], mat_cor_comp_var_global2Select2[, 2],
                   labels = rownames(mat_cor_comp_var_global2Select2),
                   cex = cex,
                   col = rep(vec_colSelect[i], dim(mat_cor_comp_var_global2Select2)[1]))
              
              
            }else{
              indice_nbVar2 = (nbVarSelect_bloc[i - 1] + 1):nbVarSelect_bloc[i]
              
              
              mat_cor_comp_var_global2Select2 = mat_cor_comp_var_global2Select[indice_nbVar2, , drop = FALSE]
              
              points(mat_cor_comp_var_global2Select2[, 1], mat_cor_comp_var_global2Select2[, 2],
                     col  = NULL)
              
              text(mat_cor_comp_var_global2Select2[, 1], mat_cor_comp_var_global2Select2[, 2],
                   labels = rownames(mat_cor_comp_var_global2Select2),
                   cex = cex,
                   col = rep(vec_colSelect[i], dim(mat_cor_comp_var_global2Select2)[1]))
              
              
            }
            
          }
          
        }
        
        par(xpd = TRUE)
        legend(x = pos[1], y = pos[2],
               legend = vec_blockEtReponseSelect,
               pch = pch,
               col = vec_colSelect,
               cex = cex_legend,
               inset = inset)
        
        
      }else{
        warning("There is no variables in this rectangle of the correlation circle.")
        
      }
      
      
    }
    
    
    
  }
  
  return(varSelect)
  
  
}

tabVarSelect <-function(varSelect)
{
  
  
}

# La fonction compute_matSimilarity calcule, pour chaque groupe de blocs, les
# similarités entre le bloc1 et le bloc2 et entre les blocs et la réponse.

compute_matSimilarity <-function(liste_dataframe_cor_comp_var_global)
{
  comp = as.numeric(sapply(3:4, FUN = function(i){
    col_i = colnames(liste_dataframe_cor_comp_var_global[[1]])[i]
    ch = strsplit(col_i, split = "_")[[1]]
    resTemp = ch[length(ch)]
    res = substring(resTemp, nchar(resTemp), nchar(resTemp))
    
    return(res)
  }))
  
  liste_matSimilarity_group = list()
  
  # On calcule, pour chaque groupe de blocs, les matrices de similarités entre
  # chaque paire de blocs.
  
  for(i in 1:length(liste_dataframe_cor_comp_var_global))
  {
    dataframe_cor_comp_var_global_i = liste_dataframe_cor_comp_var_global[[i]]
    blocks_i = unique(dataframe_cor_comp_var_global_i$bloc)
    
    coord = lapply(1:length(blocks_i), FUN = function(j){
      ind_j = which(dataframe_cor_comp_var_global_i$bloc == blocks_i[j])
      dataframe_cor_comp_var_global_i_j = dataframe_cor_comp_var_global_i[ind_j, paste0("cor_var_comp", comp)]
      res = as.matrix(dataframe_cor_comp_var_global_i_j, drop = FALSE)
      rownames(res) = rownames(dataframe_cor_comp_var_global_i_j)
      
      return(res)
    })
    
    M_block = list()
    l = 1
    
    vec_blocks_i = c()
    
    for(j in 1:(length(blocks_i) - 1))
    {
      blocks_i_j = blocks_i[j]
      
      for(k in (j + 1):length(blocks_i))
      {
        blocks_i_k = blocks_i[k]
        
        M_block[[l]] = coord[[j]][, drop = FALSE] %*% t(coord[[k]][, drop = FALSE])
        rownames(M_block[[l]]) = rownames(coord[[j]])
        colnames(M_block[[l]]) = rownames(coord[[k]])
        
        blocks_j_k = paste(c(blocks_i_j, blocks_i_k), collapse = "-")
        vec_blocks_i = c(vec_blocks_i, blocks_j_k)
        
        l = l + 1
        
      } # Fin for(k in (j + 1):length(blocks)).
      
    } # Fin for(j in 1:(length(blocks) - 1)).
    names(M_block) = vec_blocks_i
    
    liste_matSimilarity_group[[i]] = M_block
    
  } # Fin for(i in 1:length(liste_dataframe_Cor_comp_var_global)).
  
  names(liste_matSimilarity_group) = names(liste_dataframe_cor_comp_var_global) 
  
  return(list(liste_matSimilarity_group = liste_matSimilarity_group,
              comp = comp))
  
}

# La fonction networkVariableSelect permet de tracer un réseau pour les variables de
# certains blocs et des variables réponses.

networkVariableSelect <-function(liste_matSimilarity_group,
                                 comp,
                                 res_block_splsda,
                                 cutoff_comp = 0.8,
                                 vec_varBlock,
                                 vec_varRep,
                                 cutoff = 0
)
{
  
  vec_varBlockEtReponse = c(vec_varBlock, vec_varRep)
  
  # Nous vérifions que nous pouvons créer un réseau pour les variables des 
  # blocs vec_Var_blockEtReponse.
  
  # Nous recherchons le groupe de blocs associé à vec_varBlock.
  
  indice_group_vecVar = sapply(1:length(liste_matSimilarity_group), FUN = function(i){
    liste_matSimilarity_group_i = liste_matSimilarity_group[[i]]
    boolean = FALSE
    j = 1
    
    while((j <= length(liste_matSimilarity_group_i))&!boolean)
    {
      matSimilarity_group_i_j = liste_matSimilarity_group_i[[j]]
      vec_var_block1 = rownames(matSimilarity_group_i_j)
      vec_var_block2 = colnames(matSimilarity_group_i_j)
      
      vec_var_block = c(vec_var_block1, vec_var_block2)
      
      if(any(vec_var_block%in%vec_varBlock))
      {
        boolean = TRUE
        
      }
      
      j = j + 1
      
    } # Fin while((j <= length(liste_matSimilarity_group_i))&!boolean).
    
    res = boolean
    
    return(res)
  }) 
  
  res = NULL
  
  if(length(which(indice_group_vecVar == TRUE)) >= 2)
  {
    cat("Erreur : les variables de vec_varBlock. doivent appartenir à un seul élément de la liste
        liste_res_matSimilarity_group$liste_matSimilarity_group.", "\n")
    
  }else{
    liste_matSimilarity = liste_matSimilarity_group[[which(indice_group_vecVar == TRUE)]]
    
    blocks_liste_matSimilarityTemp1 = sapply(1:length(liste_matSimilarity), FUN = function(i){
      noms_block1_block2_i = names(liste_matSimilarity)[i]
      ch = strsplit(noms_block1_block2_i, split = "-")[[1]]
      block1 = ch[1]
      block2 = ch[2]
      res = c(block1, block2)
      
      return(res)
    })
    
    blocks_liste_matSimilarity = unique(as.vector(blocks_liste_matSimilarityTemp1))
    
    indice_blocks_liste_matSimilarityTemp = sapply(1:length(blocks_liste_matSimilarity), FUN = function(i){
      res = which(res_block_splsda$names$blocks == blocks_liste_matSimilarity[i])
      
      return(res)
    })
    
    ind_Y = which(res_block_splsda$names$blocks == "Y")
    indice_blocks_liste_matSimilarity = indice_blocks_liste_matSimilarityTemp[indice_blocks_liste_matSimilarityTemp != ind_Y]
    
    boolean_pos_cor = TRUE
    
    # Nous vérifions que les ièmes composantes de chaque bloc sont fortement corrélées positivement.
    
    if(length(indice_blocks_liste_matSimilarity) == 1)
    {
      
      
    }else{
      for(i in 1:length(comp))
      {
        comp_i = comp[i]
        
        for(j in 1:(length(indice_blocks_liste_matSimilarity) - 1))
        {
          indice_blocks_liste_matSimilarity_j = indice_blocks_liste_matSimilarity[j]
          comp_indice_blocks_liste_matSimilarity_j = res_block_splsda$variates[[indice_blocks_liste_matSimilarity_j]][, comp_i]
          
          
          for(k in (j + 1):length(indice_blocks_liste_matSimilarity))
          {
            indice_blocks_liste_matSimilarity_k = indice_blocks_liste_matSimilarity[k]
            comp_indice_blocks_liste_matSimilarity_k = res_block_splsda$variates[[indice_blocks_liste_matSimilarity_k]][, comp_i]
            
            cor = cor(comp_indice_blocks_liste_matSimilarity_j, comp_indice_blocks_liste_matSimilarity_k)
            
            boolean_pos_cor = boolean_pos_cor & all(cor > cutoff_comp)
            
            
          } # Fin for(k in (j + 1):length(indice_blocks_liste_matSimilarity)).
          
        } # Fin for(j in 1:(length(indice_blocks_liste_matSimilarity) - 1)).
        
      } # Fin for(i in 1:length(comp)).
      
      
    }
    
    
    
    if(!boolean_pos_cor)
    {
      cat("Erreur : pour chaque paire de bloc, la ième composante de chaque bloc doivent être corrélées positivement afin de pouvoir créer 
          un réseau.", "\n")
      
    }else{
      # Nous récupérons les matrices de similarités associés aux variables de vec_varBlockEtReponse .
      
      liste_matSimilaritySelectTemp = lapply(1:length(liste_matSimilarity), FUN = function(j){
        matSimilarity_j = liste_matSimilarity[[j]]
        
        indice_row_matSimilarity_j  = which(rownames(matSimilarity_j)%in%vec_varBlockEtReponse == TRUE)
        indice_col_matSimilarity_j  = which(colnames(matSimilarity_j)%in%vec_varBlockEtReponse == TRUE)
        
        if((length(indice_row_matSimilarity_j) != 0) &(length(indice_col_matSimilarity_j) != 0))
        {
          res = matSimilarity_j[indice_row_matSimilarity_j, indice_col_matSimilarity_j, drop = FALSE]   
          
        }else{
          res = NA
          
        }
        
        return(res)
      })
      names(liste_matSimilaritySelectTemp) = names(liste_matSimilarity)
      
      indice_NA_liste_matSimilaritySelectTemp = sapply(1:length(liste_matSimilaritySelectTemp), FUN = function(i){
        liste_matSimilaritySelectTemp_i = liste_matSimilaritySelectTemp[[i]]
        
        if(is.matrix(liste_matSimilaritySelectTemp_i))
        {
          res = FALSE
          
        }else{
          if(is.na(liste_matSimilaritySelectTemp_i))
          {
            res = TRUE
            
          }else{
            res = FALSE
            
          }
          
        }
        
        return(res)
      })
      
      
      
      liste_matSimilaritySelect = liste_matSimilaritySelectTemp[!indice_NA_liste_matSimilaritySelectTemp]
      
      # Nous créons le réseau.
      
      w = c()
      node.X1 = c()
      node.X2 = c()
      vec_group = c()
      vec_nomsVar = c()
      
      for(i in 1:length(liste_matSimilaritySelect))
      {
        
        noms_bloc1_bloc2 = names(liste_matSimilaritySelect)[i]
        matSimilaritySelect_i = liste_matSimilaritySelect[[i]]
        
        X1 = rownames(matSimilaritySelect_i)
        X2 = colnames(matSimilaritySelect_i)
        
        rep.X1 = rep(X1, each = length(X2))
        rep.X2 = rep(X2, length(X1))
        
        node.X1= c(node.X1, rep.X1)
        node.X2 = c(node.X2, rep.X2)
        
        ch = strsplit(noms_bloc1_bloc2, split = "-")[[1]]
        nom_bloc1 = ch[1]
        nom_bloc2 = ch[2]
        vec_group = c(vec_group, c(rep(nom_bloc1, length(X1)), rep(nom_bloc2, length(X2))))
        vec_nomsVar = c(vec_nomsVar, c(X1, X2))
        
        w = c(w, as.vector(t(matSimilaritySelect_i)))
        
      } # Fin for(i in 1:length(liste_matSimilaritySelect)).
      
      dup = duplicated(vec_nomsVar)
      vec_nomsVar = vec_nomsVar[!dup]
      vec_group = vec_group[!dup]
      
      nodes = data.frame(name = vec_nomsVar,
                         group = vec_group)
      
      # gR
      relations = data.frame(from = node.X1,
                             to = node.X2,
                             weight = w)
      
      idx = (abs(w) >= cutoff)
      relations = relations[idx, , drop = FALSE]
      
      gR = graph.data.frame(relations,
                            directed = FALSE,
                            vertices = nodes)
      
      # On supprime les noeuds qui n'ont pas d'arêtes.
      gR = delete.vertices(gR, which(degree(gR) == 0))
      
      res = list(gR = gR)
      res$cutoff = cutoff
      
      
      
    }
    
    
  }
  
  return(res)
  
}

# La fonction networkVar permet de tracer un réseau pour certaines variables des
# blocs et des variables réponses.

networkVar <-function(liste_matSimilarity_group = liste_matSimilarity_group,
                      comp = comp,
                      res_block_splsda,
                      cutoff_comp = 0.8,
                      vec_varBlock,
                      vec_varRep,
                      vec_varBlockInteret = NULL,
                      cutoff = 0
)
{
  if(!is.null(vec_varBlockInteret))
  {
    Var_AllBlock_list = sapply(1:(length(res_block_splsda$names$blocks) - 1), FUN = function(i){
      res = colnames(res_block_splsda$X[[i]])
      
      return(res)
    })
    Var_AllBlock_vec = unlist(Var_AllBlock_list)
    
    isInterestVariableInVar_AllBlock_vec = vec_varBlockInteret%in%Var_AllBlock_vec
    
    if(length(which(isInterestVariableInVar_AllBlock_vec == FALSE)) != 0)
    {
      stop(paste0("The variables of interest ", paste(vec_varBlockInteret[which(isInterestVariableInVar_AllBlock_vec == FALSE)], collapse = ","), " are not
                   variables of a block."))
      
    }
    
  }
  

  
  
  
  
  if(!is.null(vec_varBlock) & !is.null(vec_varBlockInteret))
  {
    indice_varCom = which(vec_varBlock%in%vec_varBlockInteret == TRUE)
    
    if(length(indice_varCom) != 0)
    {
      vec_varBlock2 = vec_varBlock[- indice_varCom]
      
    }else{
      vec_varBlock2 = vec_varBlock
      
    }
    
    vec_varBlock3 = c(vec_varBlock2, vec_varBlockInteret)
    
  }else if(!is.null(vec_varBlock) & is.null(vec_varBlockInteret))
  {
    vec_varBlock3 = vec_varBlock
    
  }else if(is.null(vec_varBlock) & !is.null(vec_varBlockInteret))
  {
    vec_varBlock3 = vec_varBlockInteret
    
  }
  
  vec_var = c(vec_varBlock3, vec_varRep)
  
  
  
  # Nous vérifions que nous pouvons créer un réseau pour les variables des 
  # blocs vec_VarBlock.
  
  # Nous recherchons le groupe de blocs associé à vec_Var.
  
  indice_group_vecVar = sapply(1:length(liste_matSimilarity_group), FUN = function(i){
    liste_matSimilarity_group_i = liste_matSimilarity_group[[i]]
    boolean = FALSE
    j = 1
    
    while((j <= length(liste_matSimilarity_group_i))&!boolean)
    {
      matSimilarity_group_i_j = liste_matSimilarity_group_i[[j]]
      vec_var_block1 = rownames(matSimilarity_group_i_j)
      vec_var_block2 = colnames(matSimilarity_group_i_j)
      
      vec_var_block1_block2 = c(vec_var_block1, vec_var_block2)
      
      if(any(vec_var_block1_block2%in%vec_varBlock3))
      {
        boolean = TRUE
        
      }
      
      j = j + 1
      
    } # Fin while((j <= length(liste_matSimilarity_group_i))&!boolean).
    
    res = boolean
    
    return(res)
  }) 
  
  res = NULL
  
  if(length(which(indice_group_vecVar == TRUE)) >= 2)
  {
    stop("The variables of vec_var have to belong to only one element of
          liste_res_matSimilarity_group$liste_matSimilarity_group.")
    
  }else{
    liste_matSimilarity = liste_matSimilarity_group[[which(indice_group_vecVar == TRUE)]]
    
    blocks_liste_matSimilarityTemp1 = sapply(1:length(liste_matSimilarity), FUN = function(i){
      noms_block1_block2_i = names(liste_matSimilarity)[i]
      ch = strsplit(noms_block1_block2_i, split = "-")[[1]]
      block1 = ch[1]
      block2 = ch[2]
      res = c(block1, block2)
      
      return(res)
    })
    
    blocks_liste_matSimilarity = unique(as.vector(blocks_liste_matSimilarityTemp1))
    
    indice_blocks_liste_matSimilarityTemp = sapply(1:length(blocks_liste_matSimilarity), FUN = function(i){
      res = which(res_block_splsda$names$blocks == blocks_liste_matSimilarity[i])
      
      return(res)
    })
    
    ind_Y = which(res_block_splsda$names$blocks == "Y")
    indice_blocks_liste_matSimilarity = indice_blocks_liste_matSimilarityTemp[indice_blocks_liste_matSimilarityTemp != ind_Y]
    
    boolean_pos_cor = TRUE
    
    if(length(indice_blocks_liste_matSimilarity) == 1)
    {
      
      
    }else{
      for(i in 1:length(comp))
      {
        comp_i = comp[i]
        
        for(j in 1:(length(indice_blocks_liste_matSimilarity) - 1))
        {
          indice_blocks_liste_matSimilarity_j = indice_blocks_liste_matSimilarity[j]
          comp_indice_blocks_liste_matSimilarity_j = res_block_splsda$variates[[indice_blocks_liste_matSimilarity_j]][, comp_i]
          
          
          for(k in (j + 1):length(indice_blocks_liste_matSimilarity))
          {
            indice_blocks_liste_matSimilarity_k = indice_blocks_liste_matSimilarity[k]
            comp_indice_blocks_liste_matSimilarity_k = res_block_splsda$variates[[indice_blocks_liste_matSimilarity_k]][, comp_i]
            
            cor = cor(comp_indice_blocks_liste_matSimilarity_j, comp_indice_blocks_liste_matSimilarity_k)
            
            boolean_pos_cor = boolean_pos_cor & cor > cutoff_comp
            
            
          } # Fin for(k in (j + 1):length(indice_blocks_vec_VarBlock)).
          
        } # Fin for(j in 1:(length(indice_blocks_vec_VarBlock)) - 1).
        
      } # Fin for(i in 1:length(comp)).
      
      
    }
    
    
    
    if(!boolean_pos_cor)
    {
      stop("For each pair of blocks, the ith component of the first block
            and the ith component of the second block have to be positively correlated in order
            to create a network.")
      
    }else{
      liste_matSimilaritySelectTemp = lapply(1:length(liste_matSimilarity), FUN = function(j){
        matSimilarity_j = liste_matSimilarity[[j]]
        
        indice_row_matSimilarity_j  = which(rownames(matSimilarity_j)%in%vec_var == TRUE)
        indice_col_matSimilarity_j  = which(colnames(matSimilarity_j)%in%vec_var == TRUE)
        
        if((length(indice_row_matSimilarity_j) != 0) &(length(indice_col_matSimilarity_j) != 0))
        {
          res = matSimilarity_j[indice_row_matSimilarity_j, indice_col_matSimilarity_j, drop = FALSE]   
          
        }else{
          res = NA
          
        }
        
        return(res)
      })
      names(liste_matSimilaritySelectTemp) = names(liste_matSimilarity)
      
      indice_NA_liste_matSimilaritySelectTemp = sapply(1:length(liste_matSimilaritySelectTemp), FUN = function(i){
        liste_matSimilaritySelectTemp_i = liste_matSimilaritySelectTemp[[i]]
        
        if(is.matrix(liste_matSimilaritySelectTemp_i))
        {
          res = FALSE
          
        }else{
          if(is.na(liste_matSimilaritySelectTemp_i))
          {
            res = TRUE
            
          }else{
            res = FALSE
            
          }
          
        }
        
        return(res)
      })
      
      
      liste_matSimilaritySelect = liste_matSimilaritySelectTemp[!indice_NA_liste_matSimilaritySelectTemp]
      
      w = c()
      node.X1 = c()
      node.X2 = c()
      vec_group = c()
      vec_nomsVar = c()
      
      for(i in 1:length(liste_matSimilaritySelect))
      {
        
        noms_block1_block2 = names(liste_matSimilaritySelect)[i]
        matSimilaritySelect_i = liste_matSimilaritySelect[[i]]
        
        X1 = rownames(matSimilaritySelect_i)
        X2 = colnames(matSimilaritySelect_i)
        
        rep.X1 = rep(X1, each = length(X2))
        rep.X2 = rep(X2, length(X1))
        
        node.X1 = c(node.X1, rep.X1)
        node.X2 = c(node.X2, rep.X2)
        
        ch = strsplit(noms_block1_block2, split = "-")[[1]]
        nom_block1 = ch[1]
        nom_block2 = ch[2]
        vec_group = c(vec_group, c(rep(nom_block1, length(X1)), rep(nom_block2, length(X2))))
        vec_nomsVar = c(vec_nomsVar, c(X1, X2))
        
        w = c(w, as.vector(t(matSimilaritySelect_i)))
        
      } # Fin for(i in 1:length(liste_matSimilaritySelect)).
      
      dup = duplicated(vec_nomsVar)
      vec_nomsVar = vec_nomsVar[!dup]
      vec_group = vec_group[!dup]
      
      nodes = data.frame(name = vec_nomsVar,
                         group = vec_group)
      
      # gR
      relations = data.frame(from = node.X1,
                             to = node.X2,
                             weight = w)
      
      # idx
      if(!is.null(vec_varBlock) & !is.null(vec_varBlockInteret) & !is.null(vec_varRep))
      {
        
        idx = sapply(1:dim(relations)[1], FUN = function(i){
          node.X1_i = relations$from[i]
          node.X2_i = relations$to[i]
          
          if(node.X1_i%in%vec_varBlockInteret | node.X2_i%in%vec_varBlockInteret)
          {
            res = TRUE
            
          }else if(node.X1_i%in%vec_varRep | node.X2_i%in%vec_varRep){
            res = TRUE
            
          }else{
            res = abs(w)[i] >= cutoff
            
          }
          
          return(res)
        }) 
        
      }else if(!is.null(vec_varBlock) & is.null(vec_varBlockInteret) & !is.null(vec_varRep))
      {
        idx = sapply(1:dim(relations)[1], FUN = function(i){
          node.X1_i = relations$from[i]
          node.X2_i = relations$to[i]
          
          if(node.X1_i%in%vec_varRep | node.X2_i%in%vec_varRep){
            res = TRUE
            
          }else{
            res = abs(w)[i] >= cutoff
            
          }
          
          return(res)
        }) 
        
      }else if(is.null(vec_varBlockEtReponse) & !is.null(vec_varBlockInteret) & !is.null(vec_varRep))
      {
        idx = rep(TRUE, dim(relations)[1])
        
      }
      
      relations = relations[idx, , drop = FALSE]
      
      gR = graph.data.frame(relations,
                            directed = FALSE,
                            vertices = nodes)
      
      # On supprime les noeuds qui n'ont pas d'arêtes.
      gR = delete.vertices(gR, which(degree(gR) == 0))
      
      res = list(gR = gR)
      res$cutoff = cutoff
      
      return(res)
      
      
    }
    
    
  }
  
  
  
  
}









