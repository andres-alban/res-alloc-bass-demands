library(rpart)
library(tidyverse)
library(parallel)
library(assertable)

run_and_predict_P = function(){
  load("Data_Ex_1.RData")
  train_visits = arrange(visits[visits$training == 1,],cluster_id,cluster_visit_id)
  train_clusters = arrange(clusters[clusters$cluster_id %in% train_visits$cluster_id,],cluster_id)
  test_visits = arrange(visits[visits$training == 0,],cluster_id,cluster_visit_id)
  test_clusters = arrange(clusters[clusters$cluster_id %in% test_visits$cluster_id,],cluster_id)
  return_clusters = test_clusters
  seed = 6913
  
  estim_res = estimate_p(train_clusters = train_clusters,train_visits = train_visits,
                         return_clusters=return_clusters,
                         seed=seed)
  
  prediction = predict_p(test_clusters = test_clusters,test_visits=test_visits,
                         test_cluster_predictions = estim_res$test_cluster_predictions,
                         train_visits = train_visits)
  hist(prediction$visit_prediction - test_visits$demand,breaks=100)
  temp = true_parameters[true_parameters$cluster_id %in% test_clusters$cluster_id,]
  hist(prediction$cluster_prediction$r*120 - temp$r*120,breaks=100)
  hist(prediction$cluster_prediction$m - temp$m,breaks=100)
  hist(prediction$cluster_prediction$p - temp$p,breaks=100)
  hist(prediction$cluster_prediction$q - temp$q,breaks=100)
  hist(prediction$cluster_prediction$c - temp$c,breaks=100)
  save(estim_res, prediction, file = "Results_PGB_Ex_1.RData")
}



# Generate hyperparameter grid and test all combinations
test_parameters = function(train_clusters,train_visits,return_clusters = NULL,seed,keep_cores_free = 1){

  # For reproducibility
  set.seed(seed)
  
  # Generate grid
  hyper_grid = expand.grid(
    r0 = c(4/120),
    p0 = c(0.001),
    q0 = c(0.1),
    c0 = c(400),
    rho = c(15000),
    eps = c(0.2), 
    k = c(0.2),
    shrinkage = c(0.1),
    validation = c(0.2),
    n_trees = 8000,
    sample = c(0.6,0.8), 
    min_split = c(10,20,30), 
    max_depth = c(4,8,12), 
    cp = c(0.01),
    optimal_trees = 0,               # to dump results
    min_error = 0,                   # to dump results
    final_error = 0,                 # to dump results
    runtime = 0                      # to dump results
  )
  
  # run estimation
  start_time = Sys.time()
  estim_res_list = estimate_p(train_clusters=train_clusters,train_visits=train_visits,hyper_grid=hyper_grid,return_models=FALSE,return_clusters=return_clusters,seed = seed,keep_cores_free = keep_cores_free)
  end_time = Sys.time()
  
  end_time - start_time
  
  # add min training error and trees to grid
  for(i in 1:length(estim_res_list)){
    if(!is.null(estim_res_list[[i]])){
      hyper_grid$optimal_trees[i] = estim_res_list[[i]]$optimal_trees
      hyper_grid$min_error[i] = estim_res_list[[i]]$min_error
      hyper_grid$final_error[i] = estim_res_list[[i]]$final_error
      hyper_grid$runtime[i] = estim_res_list[[i]]$runtime
    }
  }
  
  # see the results for the tests
  hyper_grid %>% filter(runtime > 0) %>% arrange(min_error) %>% head(10)
}


estimate_p = function(train_clusters,train_visits,seed,hyper_grid=NULL,return_models=FALSE,return_clusters=NULL,keep_cores_free = 1,
                      r0=4/120,p0=0.001,q0=0.1,c0=400,rho=15000,eps=0.2,k=0.2,shrinkage=0.1,validation = 0,
                      n_trees=12000,sample=0.8,min_split=20,max_depth=4,cp=0.01){
  
  # Run grid in parallel, or run one iteration with standard settings otherwise
  if(!is.null(hyper_grid)){
    cluster = makeCluster(detectCores() - keep_cores_free,type = 'FORK')
    clusterSetRNGStream(cl = cluster, seed)
    estim_res = parApply(cluster,hyper_grid,1,function(x){out = NULL
                                                          try({ 
                                                            out = estimate_bass_ml(train_clusters = train_clusters,
                                                                                   train_visits = train_visits,
                                                                                   starting_constant = list(r0 = x["r0"],p0 = x["p0"],q0 = x["q0"],c0 = x["c0"]),
                                                                                   penalty_params = list(rho = x["rho"], eps = x["eps"], k = x["k"]),
                                                                                   tree_params = list(shrinkage=x["shrinkage"],validation=x["validation"],n_trees=x["n_trees"],sample=x["sample"],min_split=x["min_split"],max_depth=x["max_depth"],cp=x["cp"]),
                                                                                   optim_params = list(tol=1e-16,max_iter=500),
                                                                                   return_models = return_models,
                                                                                   return_clusters=return_clusters)
                                                          },silent=TRUE)
                                                          out})
    stopCluster(cluster)
  } else{
    set.seed(seed)
    estim_res = estimate_bass_ml(train_clusters = train_clusters,
                                 train_visits = train_visits,
                                 starting_constant = list(r0=r0,p0=p0,q0=q0,c0=c0),
                                 penalty_params=list(rho = rho, eps = eps, k = k),
                                 tree_params=list(shrinkage=shrinkage,validation=validation,n_trees=n_trees,sample=sample,min_split=min_split,max_depth=max_depth,cp=cp),
                                 optim_params=list(tol=1e-16,max_iter=500),
                                 return_models = return_models,
                                 return_clusters=return_clusters)
  }
  return(estim_res)
}

estimate_bass_ml = function(train_clusters,
                            train_visits,
                            starting_constant = list(r0 = 4/120,p0 = 0.001,q0 = 0.1,c0 = 100),
                            starting_constraints = c(-7/120,1/120,0,0,0),
                            penalty_params=list(rho = 15000, eps = 0.15, k = 0.2),
                            tree_params=list(shrinkage=0.1,validation=0.2,n_trees=4995,sample=0.8,min_split=20,max_depth=16,cp=0.01),
                            optim_params=list(tol=1e-16,max_iter=500),
                            return_models = FALSE,
                            return_clusters = NULL){
  
  assert_colnames(train_clusters, c("cluster_id","m_factor","pop_size_b0","pop_size_b1","pop_size_b2","pop_size_b3"), only_colnames = FALSE, quiet = FALSE)
  assert_colnames(train_visits, c("cluster_id","cluster_visit_id","demand"), only_colnames = FALSE, quiet = FALSE)
  if(!is.null(return_clusters)){
    assert_colnames(train_clusters, c("cluster_id","m_factor","pop_size_b0","pop_size_b1","pop_size_b2","pop_size_b3"), only_colnames = FALSE, quiet = FALSE)
  }
  if(!is.list(starting_constant) | sum(c("r0","q0","p0","c0") %in% names(starting_constant)) < 4){
    stop("List \"starting_constant\" needs to contain variables r0, p0, q0, and c0")
  }
  if(length(starting_constraints) != 5){
    stop("Vector \"starting_constraints\" needs to contain five elements")
  }
  if(!is.list(penalty_params) | sum(c("rho","eps","k") %in% names(penalty_params)) < 3){
    stop("List \"penalty_params\" needs to contain variables rho, eps, and k")
  }
  if(!is.list(tree_params) | sum(c("shrinkage","validation","n_trees","sample","min_split","max_depth","cp") %in% names(tree_params)) < 7){
    stop("List \"tree_params\" needs to contain variables shrinkage, validation, n_trees, sample, min_split, max_depth, and cp")
  }
  if(!is.list(optim_params) | sum(c("tol","max_iter") %in% names(optim_params)) < 2){
    stop("List \"optim_params\" needs to contain variables tol, and max_iter")
  }
  
  
  start_time = Sys.time()
  mse = rep(Inf,tree_params$n_trees)
  number_models = length(starting_constant)
  
  # Split into training and validation sets
  pick_train = sample(1:nrow(train_visits),size = ceiling((1-tree_params$validation)*nrow(train_visits)))
  training_visits = train_visits[pick_train,]
  training_visits = arrange(training_visits, cluster_id, cluster_visit_id)
  training_clusters = train_clusters[train_clusters$cluster_id %in% unique(training_visits$cluster_id),]
  training_clusters = arrange(training_clusters,cluster_id)
  if(nrow(training_visits) < nrow(train_visits)){
    validation_visits = train_visits[-pick_train,]
    validation_visits = arrange(validation_visits, cluster_id, cluster_visit_id)
    validation_clusters = train_clusters[train_clusters$cluster_id %in% unique(validation_visits$cluster_id),]
    validation_clusters = arrange(validation_clusters,cluster_id)
  }
  
  # Separate regression data and data used to run Bass model
  reg_data = training_clusters[,!colnames(training_clusters) %in% c("cluster_id","visits_to_cluster","training_visits","m_factor","beta","pop_size_b0","pop_size_b1","pop_size_b2","pop_size_b3")]
  importance_df = as.data.frame(matrix(0,nrow=number_models,ncol=length(colnames(reg_data))))
  colnames(importance_df) = colnames(reg_data)
  mf = training_clusters$m_factor
  beta = cbind(training_clusters$pop_size_b0,training_clusters$pop_size_b1,training_clusters$pop_size_b2,training_clusters$pop_size_b3)
  if(nrow(training_visits) < nrow(train_visits)){
    mf_validation = validation_clusters$m_factor
    beta_validation = cbind(validation_clusters$pop_size_b0,validation_clusters$pop_size_b1,validation_clusters$pop_size_b2,validation_clusters$pop_size_b3)
  }
  
  # Generate list of demands and visits (each entry corresponds to a cluster and contains a vector of demands for each visit or visit numbers)
  demand_list = aggregate(demand ~ cluster_id, training_visits, c)
  demand_list = arrange(demand_list, cluster_id)
  demand_list = demand_list$demand
  visit_list = aggregate(cluster_visit_id ~ cluster_id, training_visits, c)
  visit_list = arrange(visit_list, cluster_id)
  visit_list = visit_list$cluster_visit_id
  
  if(nrow(training_visits) < nrow(train_visits)){
    demand_validation = aggregate(demand ~ cluster_id, validation_visits, c)
    demand_validation = arrange(demand_validation, cluster_id)
    demand_validation = demand_validation$demand
    visit_validation = aggregate(cluster_visit_id ~ cluster_id, validation_visits, c)
    visit_validation = arrange(visit_validation, cluster_id)
    visit_validation = visit_validation$cluster_visit_id
  }
  
  # Generate F and h matrices for gradient boosting updates
  number_clusters = nrow(training_clusters)
  F_mat = matrix(0,nrow=number_clusters,ncol=number_models)
  h_mat = matrix(1,nrow=number_clusters,ncol=number_models)
  penalty_m = (length(starting_constraints) + 1) * number_clusters
  
  # Generate constraint matrix for optimization of initial constants
  ui = c(-1,rep(0,number_models-1))
  ui = rbind(ui,diag(number_models))
  ci = starting_constraints
  opt_res = constrOptim(theta = c(starting_constant$r0,starting_constant$p0,starting_constant$q0,starting_constant$c0),
                        f = min_fun, grad = NULL,ui = ui, ci = ci,
                        mf=mf,beta=beta,demand_list=demand_list,visit_list=visit_list,
                        F_mat=F_mat,h_mat=h_mat,
                        penalty_params=NULL,
                        control = list(reltol=optim_params$tol, maxit = optim_params$max_iter))
  constant = opt_res$par
  for(i in 1:number_models){
    F_mat[,i] = constant[i]
  }
  
  # Update prediction on validation set to constant
  if(nrow(training_visits) < nrow(train_visits)){
    pred_validation = matrix(0,nrow=nrow(validation_clusters),ncol= number_models)
    for(i in 1:number_models){
      pred_validation[,i] = constant[i]
    }
  }
  
  # Update prediction on test set to constant
  if(!is.null(return_clusters)){
    pred_test = matrix(0,nrow=nrow(return_clusters),ncol= number_models)
    for(i in 1:number_models){
      pred_test[,i] = constant[i]
    }
    best_pred_test = pred_test
    best_mse = Inf
  }
  
  # Iterate to improve model
  models = list()
  source("Pseudo_Residuals.R")
  
  for(m_it in 1:tree_params$n_trees){
    ## Generate pseudo-residuals of the four functions and predict them
    ## refer to: https://uc-r.github.io/regression_trees#prereq
    pred_F = list()
    for(i in 1:number_models){
      residual_fun = match.fun(paste0("residual_",i))
      e = mapply(residual_fun,mf,beta[,1],beta[,2],beta[,3],beta[,4],
                 F_mat[,1],F_mat[,2],F_mat[,3],F_mat[,4],
                 demand_list,visit_list,
                 MoreArgs = list(penalty_params = NULL, m = penalty_m))
      ## Predict the pseudo-residuals (get both function and resulting values) for a sample - SVR? Simultaneous prediction?
      pred_data = reg_data
      pred_data$e = e
      pred_data = pred_data[sample(nrow(pred_data),round(tree_params$sample*nrow(pred_data))),]
      m = rpart(formula = e ~ .,
                data = pred_data,
                method  = "anova",
                control = list(minsplit=tree_params$min_split,maxdepth=tree_params$max_depth,cp=tree_params$cp))
      h_mat[,i] = predict(m,newdata=reg_data)
      pred_F[[i]] = m
    }
    
    ## Estimate model weights
    opt_res = optim(par = rep(0,number_models),fn = min_fun, gr = NULL,
                    mf=mf,beta=beta,demand_list=demand_list,visit_list=visit_list,
                    F_mat=F_mat,h_mat=h_mat,
                    penalty_params=penalty_params,
                    control = list(reltol=optim_params$tol, maxit = optim_params$max_iter))
    gamma = opt_res$par
    
    ## Update estimated values and model list
    F_mat = F_mat + (h_mat %*% diag(gamma)) * tree_params$shrinkage
    if(return_models){
      models[[m_it]] = list(weights = gamma*tree_params$shrinkage, indiv_models = pred_F)
    }
    
    ## Update validation data prediction and correpsonding MSE computation
    if(nrow(training_visits) < nrow(train_visits)){
      for(i in 1:number_models){
        pred_validation[,i] = pred_validation[,i] + gamma[i] * tree_params$shrinkage * predict(pred_F[[i]],newdata=validation_clusters)
      }
      mse[m_it] = sqrt(min_fun(gamma, mf = mf_validation, beta = beta_validation,
                               demand_list = demand_validation, visit_list = visit_validation,
                               F_mat = pred_validation, h_mat = matrix(0,nrow=nrow(pred_validation),ncol=number_models),
                               penalty_params = penalty_params)) / 
        nrow(pred_validation)
    } else{
      mse[m_it] = sqrt(min_fun(gamma, mf = mf, beta = beta, demand_list = demand_list, visit_list = visit_list,
                               F_mat = F_mat, h_mat = matrix(0,nrow=number_clusters,ncol=number_models),
                               penalty_params = penalty_params)) /
        nrow(F_mat)
    }
    
    ## Update test prediction
    if(!is.null(return_clusters)){
      for(i in 1:number_models){
        pred_test[,i] = pred_test[,i] + gamma[i] * tree_params$shrinkage * predict(pred_F[[i]],newdata=return_clusters)
      }
      if(mse[m_it] < best_mse){
        best_pred_test = pred_test
        best_mse = mse[m_it]
      }
    }
    
    ## Update matrix giving the importance of estimates in the individual residual-DTs
    for(i in 1:number_models){
      importance = pred_F[[i]]$variable.importance
      importance_df[i,names(importance)] = importance_df[i,names(importance)] + importance * gamma[i] * tree_params$shrinkage
    }
    
    ## print result
    if(nrow(training_visits) < nrow(train_visits)){
      print(paste0("Decision tree iteration ", m_it, " completed with average loss function on validation set of ", mse[m_it]))
    } else{
      print(paste0("Decision tree iteration ", m_it, " completed with average loss function on training set of ", mse[m_it]))
    }
    
    ## if the mse explodes (e.g., because of negative r), stop
    if(mse[m_it] > 1000 * mse[max(m_it-1,1)]){
      break
    }
  }
  end_time = Sys.time()
  
  var_importance = importance_df
  result_list = list(optimal_trees = which.min(mse), min_error = min(mse), final_error = mse[length(mse)], runtime = as.numeric(difftime(end_time,start_time,units="min")),
                     penalty_params=penalty_params,tree_params=tree_params,optim_params=optim_params,var_importance = var_importance)
  if(return_models){
    result_list = c(result_list, list(constant = constant,models = models))
  }
  if(!is.null(return_clusters)){
    pred_test = as.data.frame(best_pred_test)
    colnames(pred_test) = c("r","p","q","c")
    pred_test$cluster_id = return_clusters$cluster_id
    pred_test$impact = return_clusters$impact
    result_list = c(result_list, list(test_cluster_predictions = pred_test))
  }
  return(result_list)
}

min_fun = function(gamma,mf,beta,demand_list,visit_list,F_mat,h_mat,penalty_params){
  
  F_mat = F_mat + h_mat %*% diag(gamma)
  main = sum(mapply(sq_error,mf,beta[,1],beta[,2],beta[,3],beta[,4],F_mat[,1],F_mat[,2],F_mat[,3],F_mat[,4],demand_list,visit_list))

  ## Penalty function: Xinsheng Xu, Zhiqing Meng, Jianwu Sun, Rui Shen,
  ## A penalty function method based on smoothing lower order penalty function,
  ## Journal of Computational and Applied Mathematics,Volume 235, Issue 14, 2011
  if(is.null(penalty_params)){
    penalty = 0
  } else{
    t = 1/120 - F_mat[,1]
    t = c(t, F_mat[,1] - 7/120)
    t = c(t, F_mat[,4] - mf*(beta[,1] + beta[,2] * F_mat[,1] + beta[,3] * F_mat[,1]^2 + beta[,4] * F_mat[,1]^3))
    t = c(t,0 - F_mat[,2])
    t = c(t,0 - F_mat[,3])
    t = c(t,0 - F_mat[,4])
    penalty = penalty_fun(t, penalty_params$rho, penalty_params$eps, penalty_params$k)
  }
  main + penalty
}

sq_error = function(mf,b0,b1,b2,b3,F1,F2,F3,F4,y,v){
  mean = g(mf,b0,b1,b2,b3,F1,F2,F3,F4,v)
  sum((mean - y)^2)
}

g = function(mf,b0,b1,b2,b3,F1,F2,F3,F4,v){
  r = F1
  p = F2
  q = F3
  c = F4
  m = mf*(b0 + b1 * r + b2 * r^2 + b3 * r^3)
  theta = (m - c)*p/(m*p+c*q)
  mean = (1-theta*exp(-(p+q)*v))/(1+q/p*theta*exp(-(p+q)*v)) - (1-theta*exp(-(p+q)*(v-1)))/(1+q/p*theta*exp(-(p+q)*(v-1)))
  mean = mean * m
}

penalty_fun = function(t,rho,eps,k){
  m = length(t)
  sum(rho*ifelse(t<0,
                 0,ifelse(t < (eps/(m*rho))^(1/k),
                          (m^2*rho^2*t^(3*k))/(eps^2) - (m^3*rho^3*t^(4*k))/(2*eps^3),t^k - eps/(2*m*rho))))
}

predict_p = function(test_clusters,
                     test_visits,
                     constraints = c(-7/120,1/120,0,0,0),
                     test_cluster_predictions = NULL,
                     constant = NULL,
                     models = NULL,
                     train_visits = NULL){
  
  assert_colnames(test_clusters, c("cluster_id","beta","m_factor","pop_size_b0","pop_size_b1","pop_size_b2","pop_size_b3"), only_colnames = FALSE, quiet = FALSE)
  assert_colnames(test_visits, c("cluster_id","cluster_visit_id","demand"), only_colnames = FALSE, quiet = FALSE)
  if(!is.null(train_visits)){
    assert_colnames(train_visits, c("cluster_id","cluster_visit_id","demand"), only_colnames = FALSE, quiet = FALSE)
  }
  if(length(constraints) != 5){
    stop("Vector \"constraints\" needs to contain five elements")
  }
  if(!is.null(test_cluster_predictions)){
    assert_colnames(test_cluster_predictions, c("cluster_id","r","p","q","c"), only_colnames = FALSE, quiet = FALSE)
  } else{
    if(!is.list(models) | !is.list(models[[1]]) | sum(c("weights","indiv_models") %in% names(models[[1]])) < 2){
      stop("List \"models\" needs to contain lists that in turn contain the vectors weights and indiv_models")
    }
    if(is.null(constant) | length(constant) != length(models[[1]]$weights)){
      stop("If \"test_cluster_predictions\" is not given, a vector of constants, one per model, is needed")
    }
    number_models = length(models[[1]]$weights)
    test_cluster_predictions = matrix(0,nrow=nrow(test_clusters),ncol = number_models)
    for(i in 1:number_models){
      test_cluster_predictions[,i] = constant[i]
      for(j in 1:length(models)){
        test_cluster_predictions[,i] = test_cluster_predictions[,i] + models[[j]]$weights[i] * predict(models[[j]]$indiv_models[[i]],newdata=test_clusters)
      }
    }
    test_cluster_predictions = as.data.frame(test_cluster_predictions)
    colnames(test_cluster_predictions) = c("r","p","q","c")
    test_cluster_predictions$cluster_id = test_clusters$cluster_id
  }
  test_cluster_predictions = merge(test_cluster_predictions,test_clusters[,c("cluster_id","beta","m_factor","pop_size_b0","pop_size_b1","pop_size_b2","pop_size_b3")],by="cluster_id")
  
  r95 = as.numeric(quantile(test_cluster_predictions$r[test_cluster_predictions$r <= -constraints[1]],0.95))
  r5 = as.numeric(quantile(test_cluster_predictions$r[test_cluster_predictions$r >= constraints[2]],0.05))
  p5 = as.numeric(quantile(test_cluster_predictions$p[test_cluster_predictions$p >= constraints[3]],0.05))
  q5 = as.numeric(quantile(test_cluster_predictions$q[test_cluster_predictions$q >= constraints[4]],0.05))
  c5 = as.numeric(quantile(test_cluster_predictions$c[test_cluster_predictions$c >= constraints[5]],0.05))
  
  test_cluster_predictions$r = ifelse(test_cluster_predictions$r > -constraints[1],r95,test_cluster_predictions$r)
  test_cluster_predictions$r = ifelse(test_cluster_predictions$r < constraints[2],r5,test_cluster_predictions$r)
  test_cluster_predictions$p = ifelse(test_cluster_predictions$p < constraints[3],p5,test_cluster_predictions$p)
  test_cluster_predictions$q = ifelse(test_cluster_predictions$q < constraints[4],q5,test_cluster_predictions$q)
  test_cluster_predictions$c = ifelse(test_cluster_predictions$c < constraints[5],c5,test_cluster_predictions$c)
  test_cluster_predictions$m = test_cluster_predictions$m_factor*(test_cluster_predictions$pop_size_b0 + test_cluster_predictions$pop_size_b1 * test_cluster_predictions$r + test_cluster_predictions$pop_size_b2 * test_cluster_predictions$r^2 + test_cluster_predictions$pop_size_b3 * test_cluster_predictions$r^3)
  
  # Correct visits/c based on what was observed in the training data
  visit_list = aggregate(cluster_visit_id ~ cluster_id, test_visits, c)
  cluster_predictions_reduced = test_cluster_predictions[test_cluster_predictions$cluster_id %in% unique(visit_list$cluster_id),]
  cluster_predictions_reduced = arrange(cluster_predictions_reduced,cluster_id)
  visit_list = arrange(visit_list, cluster_id)
  visit_list = visit_list$cluster_visit_id
  if(!is.null(train_visits)){
    previous_totals = train_visits %>% group_by(cluster_id) %>% summarize(total_visits = sum(demand))
    cluster_predictions_reduced = merge(cluster_predictions_reduced,previous_totals,by="cluster_id",all.x = TRUE, all.y = FALSE)
    cluster_predictions_reduced$total_visits[is.na(cluster_predictions_reduced$total_visits)] = 0
    cluster_predictions_reduced$c = cluster_predictions_reduced$c + cluster_predictions_reduced$total_visits
    for(i in 1:length(visit_list)){
      visit_list[[i]] = 1:length(visit_list[[i]])
    }
  }
  
  # Predict with (corrected) Bass-regression parameters
  predicted_demand_list = mapply(g,cluster_predictions_reduced$m_factor,cluster_predictions_reduced$pop_size_b0,cluster_predictions_reduced$pop_size_b1,cluster_predictions_reduced$pop_size_b2,cluster_predictions_reduced$pop_size_b3,
                                 cluster_predictions_reduced$r,cluster_predictions_reduced$p,cluster_predictions_reduced$q,cluster_predictions_reduced$c,
                                 visit_list)
  pred_vec = unlist(predicted_demand_list)
  
  # Correction: negative numbers to zero
  pred_vec[pred_vec < 0] = 0
  
  list(visit_prediction = pred_vec, cluster_prediction = test_cluster_predictions)
}
