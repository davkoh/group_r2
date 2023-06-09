
#----- Full simulation auxiliary function for comparison

full_simulation <- function(sim_conds,
                            sim_params,
                            smqoi,
                            ncores_simulation = 1,
                            seed = NULL,
                            path = NULL, 
                            special_name= "") {
  # Set seed for reproducibility.
  if (!is.null(seed)) {
    set.seed(seed)
  }
  # Global seeds
  seed_list <- sample(1000000000:.Machine$integer.max,
                      size = nrow(sim_conds)) 
  
  final_result <- vector(mode = "list", length = nrow(sim_conds))
  
  # Iterate over dataset configurations and combine the results
  for (i in seq_len(nrow(sim_conds))) {
    final_result[[i]] <- dataset_cond_sim(
      sim_cond = sim_conds[i, ],
      smqoi= smqoi,
      sim_params = sim_params,
      seed = seed_list[[i]],
      path = path,
      ncores = ncores_simulation,
      special_name=special_name
    )
  }
  
  final_result <- do.call(rbind, final_result)
  #final_result$global_seed <- seed
  
  if (!is.null(path)) {
    saveRDS(final_result, paste(path, paste0("full_sim_comp_result_", special_name ,".RDS"), sep = "/"))
  }
  return(final_result)
}




