library(tidyverse)
library(parallel)
library(assertable)

run_and_predict_NLS = function(){
  load("Data_Ex_1.RData")
  train_visits = arrange(visits[visits$training == 1,],cluster_id,cluster_visit_id)
  train_clusters = arrange(clusters[clusters$cluster_id %in% train_visits$cluster_id,],cluster_id)
  test_visits = arrange(visits[visits$training == 0,],cluster_id,cluster_visit_id)
  test_clusters = arrange(clusters[clusters$cluster_id %in% test_visits$cluster_id,],cluster_id)
  seed = 39218
  
  estim_res = estimate_nls(train_clusters = train_clusters,train_visits = train_visits,
                           validation = 0,  seed=seed)
  
  prediction = predict_nls(test_clusters = test_clusters,test_visits=test_visits,
                         train_visits = train_visits,
                         r_coefs = estim_res$r_coefs, p_coefs = estim_res$p_coefs, q_coefs = estim_res$q_coefs, c_coefs = estim_res$c_coefs)
  hist(prediction - test_visits$demand,breaks=100)
  temp = true_parameters[true_parameters$cluster_id %in% test_clusters$cluster_id,]
  reg_matrix = as.matrix(test_clusters[,!colnames(test_clusters) %in% c("cluster_id","visits_to_cluster","training_visits","m_factor","beta","pop_size_b0","pop_size_b1","pop_size_b2","pop_size_b3")])
  reg_matrix = cbind(rep(1,nrow(reg_matrix)),reg_matrix)
  reg_matrix[is.na(reg_matrix)] = 0
  param_pred_r = (reg_matrix %*% estim_res$r_coefs)[,1]
  param_pred_m = test_clusters$m_factor * (test_clusters$pop_size_b0 + test_clusters$pop_size_b1*param_pred_r + test_clusters$pop_size_b2*param_pred_r^2 + test_clusters$pop_size_b3*param_pred_r^3)
  param_pred_p = (reg_matrix %*% estim_res$p_coefs)[,1]
  param_pred_q = (reg_matrix %*% estim_res$q_coefs)[,1]
  param_pred_c = (reg_matrix %*% estim_res$c_coefs)[,1]
  hist(param_pred_r*120 - temp$r*120,breaks=100)
  hist(param_pred_m - temp$m,breaks=100)
  hist(param_pred_p - temp$p,breaks=100)
  hist(param_pred_q - temp$q,breaks=100)
  hist(param_pred_c - temp$c,breaks=100)
  save(estim_res, prediction, file = "Results_NLS_Ex_1.RData")
}

test_parameters = function(train_clusters, train_visits,seed = 123){
  
  # For reproducibility
  set.seed(seed)
  
  # Generate grid
  hyper_grid = expand.grid(
    seed = sample(1:1000, 1),
    validation = c(0.2,0.4),
    r_0 = c(3/120,5/120),
    p_0 = c(0.001,0.1),
    q_0 = c(0.001,0.1),
    c_0 = c(250,500),
    train_MSE = 0,                    # to dump results
    validation_MSE = 0                # to dump results
  )
  
  # Run NLS for all parameter combinations with partial training set
  calc_cluster = makeCluster(detectCores(),type = 'FORK')
  clusterSetRNGStream(cl = calc_cluster, seed)
  estim_res_list = parApply(calc_cluster,hyper_grid,1,function(x) estimate_nls(train_clusters = train_clusters, train_visits = train_visits,
                                                                           starting_coefs = list(r0 = x["r_0"], p0 = x["p_0"], q0 = x["q_0"], c0 = x["c_0"]),
                                                                           seed = x["seed"], validation = x["validation"]))
  stopCluster(calc_cluster)
  
  # Gather results
  for(i in 1:nrow(hyper_grid)) {
    hyper_grid$train_MSE[i] = estim_res[[i]]$train_MSE
    hyper_grid$validation_MSE[i] = estim_res[[i]]$validation_MSE
  }
  
  # see the results for the tests
  hyper_grid %>% arrange(validation_MSE) %>% head(10)
}

# Estimation procedure
estimate_nls = function(train_clusters,
                        train_visits,
                        starting_coefs = list(r0 = 5/120,p0 = 0.001,q0 = 0.1,c0=250),
                        starting_correction = 0.02,
                        starting_constraints = c(-7/120,1/120,0,0,0),
                        optim_params=list(tol=1e-16,max_iter=5000),
                        validation = 0,
                        seed = 123){
    
  assert_colnames(train_clusters, c("cluster_id","m_factor","pop_size_b0","pop_size_b1","pop_size_b2","pop_size_b3"), only_colnames = FALSE, quiet = FALSE)
  assert_colnames(train_visits, c("cluster_id","cluster_visit_id","demand"), only_colnames = FALSE, quiet = FALSE)
  if(!is.list(starting_coefs) | sum(c("r0","q0","p0","c0") %in% names(starting_coefs)) < 4){
    stop("List \"starting_coefs\" needs to contain variables r0, p0, and q0")
  }
  if(length(starting_constraints) != 5){
    stop("Vector \"starting_constraints\" needs to contain five elements")
  }
  if(!is.list(optim_params) | sum(c("tol","max_iter") %in% names(optim_params)) < 2){
    stop("List \"optim_params\" needs to contain variables tol, and max_iter")
  }
  
  set.seed(seed)
  train_clusters = arrange(train_clusters,cluster_id)
  train_clusters$cluster_id_rebased = 1:nrow(train_clusters)
  train_visits = merge(train_visits,train_clusters[c("cluster_id","cluster_id_rebased")],by="cluster_id")
  train_visits$cluster_id = train_visits$cluster_id_rebased
  train_visits$cluster_id_rebased = NULL
  train_clusters$cluster_id = train_clusters$cluster_id_rebased
  train_clusters$cluster_id_rebased = NULL

  pick_train = sample(1:nrow(train_visits),size = ceiling((1-validation)*nrow(train_visits)))
  training_visits = train_visits[pick_train,]
  training_visits = arrange(training_visits, cluster_id, cluster_visit_id)
  if(nrow(training_visits) < nrow(train_visits)){
    validation_visits = train_visits[-pick_train,]
    validation_visits = arrange(validation_visits, cluster_id, cluster_visit_id)
  }
  
  # Generate X-matrix for the different estimations
  reg_matrix = as.matrix(train_clusters[,!colnames(train_clusters) %in% c("cluster_id","visits_to_cluster","training_visits","m_factor","beta","pop_size_b0","pop_size_b1","pop_size_b2","pop_size_b3")])
  reg_matrix = cbind(rep(1,nrow(reg_matrix)),reg_matrix)
  reg_matrix[is.na(reg_matrix)] = 0
  
  # Generate reduced cluster data frame for remaining tasks
  train_clusters = train_clusters[,c("pop_size_b0","pop_size_b1","pop_size_b2","pop_size_b3","m_factor")]
  
  # Optimization constraints (keep p and q above 0, radius between 1 and 7 kilometers)
  ui = cbind(-reg_matrix,
             matrix(0,nrow=nrow(reg_matrix),ncol=2*ncol(reg_matrix)))
  ui = rbind(ui, cbind(reg_matrix,
                       matrix(0,nrow=nrow(reg_matrix),ncol=2*ncol(reg_matrix))))
  ui = rbind(ui, cbind(matrix(0,nrow=nrow(reg_matrix),ncol=ncol(reg_matrix)),
                       reg_matrix,
                       matrix(0,nrow=nrow(reg_matrix),ncol=ncol(reg_matrix))))
  ui = rbind(ui, cbind(matrix(0,nrow=nrow(reg_matrix),ncol=2*ncol(reg_matrix)),
                       reg_matrix))
  ci = c(rep(starting_constraints[1],nrow(reg_matrix)),rep(starting_constraints[2],nrow(reg_matrix)),rep(starting_constraints[3],nrow(reg_matrix)),rep(starting_constraints[4],nrow(reg_matrix)))
  # Additionally, keep c above 0
  ui_add_1 = matrix(0,nrow=nrow(ui),ncol=ncol(reg_matrix))
  ui_add_2 = matrix(0,nrow=nrow(reg_matrix),ncol=ncol(ui))
  ui_add_3 = reg_matrix
  ui = rbind(cbind(ui,ui_add_1),cbind(ui_add_2,ui_add_3))
  ci = c(ci,rep(starting_constraints[5],nrow(reg_matrix)))
  
  # Correction of starting parameters, in case they violate constraints
  r_coefs_start = c(starting_coefs$r0,rep(0,length.out=ncol(reg_matrix)-1))
  r_coefs = r_coefs_start
  weight = 1-starting_correction
  while(sum(reg_matrix %*% r_coefs <= starting_constraints[2]) > 0 | sum(reg_matrix %*% r_coefs >= -starting_constraints[1]) > 0){
    r_coefs = (1-weight)*c(5/120,rep(0,length(r_coefs_start)-1)) + weight*r_coefs_start
    weight = weight - starting_correction
  }
  p_coefs_start = c(starting_coefs$p0,rep(0,length.out=ncol(reg_matrix)-1))
  p_coefs = p_coefs_start
  weight = 1-starting_correction
  while(sum(reg_matrix %*% p_coefs <= starting_constraints[3]) > 0){
    p_coefs = (1-weight)*c(0.0001,rep(0,lengh(p_coefs_start)-1)) + weight*p_coefs_start
    weight = weight - starting_correction
  }
  q_coefs_start = c(starting_coefs$q0,rep(0,length.out=ncol(reg_matrix)-1))
  q_coefs = q_coefs_start
  weight = 1-starting_correction
  while(sum(reg_matrix %*% q_coefs <= starting_constraints[4]) > 0){
    q_coefs = (1-weight)*c(0.01,rep(0,lengh(q_coefs_start)-1)) + weight*q_coefs_start
    weight = weight - starting_correction
  }
  c_coefs_start = c(starting_coefs$c0,rep(0,length.out=ncol(reg_matrix)-1))
  c_coefs = c_coefs_start
  weight = 1-starting_correction
  while(sum(reg_matrix %*% c_coefs < starting_constraints[5]) > 0){
    c_coefs = (1-weight)*c(100,rep(0,lengh(c_coefs_start)-1)) + weight*c_coefs_start
    weight = weight - starting_correction
  }
  
  # Run estimation
  start = c(r_coefs,p_coefs,q_coefs,c_coefs)
  start_time = Sys.time()
  output = constrOptim(theta = start,f = obj_fun,grad = NULL,
                       ui = ui, ci = ci,
                       input = list(reg_matrix,train_clusters,training_visits),
                       control = list(reltol=optim_params$tol, maxit = optim_params$max_iter))
  
  min = 1
  max = length(r_coefs)
  r_coefs = output$par[min:max]
  
  min = max + 1
  max = max + length(p_coefs)
  p_coefs = output$par[min:max]
  
  min = max + 1
  max = max + length(q_coefs)
  q_coefs = output$par[min:max]
  
  min = max + 1
  max = max + length(c_coefs)
  c_coefs = output$par[min:max]
  
  # Check the MSE in the validation set
  if(nrow(training_visits) < nrow(train_visits)){
    validation_MSE = sqrt(obj_fun(pars = c(r_coefs,p_coefs,q_coefs,c_coefs),
                                  input = list(reg_matrix,train_clusters,validation_visits)))/nrow(validation_visits)
  } else{
    validation_MSE = 0
  }
  
  result = list(train_MSE = sqrt(output$val)/nrow(training_visits), validation_MSE = validation_MSE, runtime = as.numeric(difftime(Sys.time(),start_time,units="min")),
                r_coefs = r_coefs, p_coefs = p_coefs, q_coefs = q_coefs, c_coefs = c_coefs)
  
  return(result)
}

# Minimize the sum of squares
obj_fun = function(pars,input){
  reg_matrix = input[[1]]
  clusters = input[[2]]
  visits = input[[3]]
  
  min = 1
  max = ncol(reg_matrix)
  r_coefs = pars[min:max]
  
  min = max + 1
  max = max + ncol(reg_matrix)
  p_coefs = pars[min:max]
  
  min = max + 1
  max = max + ncol(reg_matrix)
  q_coefs = pars[min:max]
  
  min = max + 1
  max = max + ncol(reg_matrix)
  c_coefs = pars[min:max]
  
  reg_matrix = reg_matrix[visits$cluster_id,]
  clusters = clusters[visits$cluster_id,]
  observations = visits$demand
  visit_no = visits$cluster_visit_id
  
  ssq(r_coefs,p_coefs,q_coefs,c_coefs,
      reg_matrix,clusters,visit_no,observations)
}

# Compute the sum of squares
ssq = function(r_coefs,p_coefs,q_coefs,c_coefs,
               reg_matrix,clusters,visit_no,observations){
  pred = predict_clients(r_coefs,p_coefs,q_coefs,c_coefs,
                         reg_matrix,clusters,visit_no)
  ll =  sum((observations-pred)^2)
  return(ll)
}

# Predict client visits for a given cluster and visit number
predict_clients = function(r_coefs,p_coefs,q_coefs,c_coefs,
                           reg_matrix,clusters,visit_no,
                           c_start = NULL){
  if(is.null(c_start)){
    c_start = (reg_matrix %*% c_coefs)[,1]
  }
  r = (reg_matrix %*% r_coefs)[,1]
  p = (reg_matrix %*% p_coefs)[,1]
  q = (reg_matrix %*% q_coefs)[,1]
  m = 0
  cols = c("pop_size_b0","pop_size_b1","pop_size_b2","pop_size_b3")
  for(i in 1:length(cols)){
    m = m + clusters[[cols[i]]] * r^(i-1)
  }
  m_low = as.numeric(quantile(m[m > 0],0.1,na.rm=TRUE))
  m = ifelse(m<=0,m_low,m)
  m = m * clusters$m_factor
  theta = (m - c_start)*p/(m*p+c_start*q)
  pred = (1-theta*exp(-(p+q)*visit_no))/(1+q/p*theta*exp(-(p+q)*visit_no)) - (1-theta*exp(-(p+q)*(visit_no-1)))/(1+q/p*theta*exp(-(p+q)*(visit_no-1)))
  pred = pred * m
  return(pred)
}

predict_nls = function(test_clusters, test_visits, train_visits = NULL,
                        r_coefs, p_coefs, q_coefs,c_coefs){
  
  test_visits = merge(test_visits,test_clusters,by="cluster_id")
  test_visits = arrange(test_visits,cluster_id,cluster_visit_id)
  
  # Generate matrix of prediction data for p and q
  reg_matrix = as.matrix(test_visits[,colnames(test_visits) %in% colnames(test_clusters) & !colnames(test_visits) %in% c("cluster_id","visits_to_cluster","training_visits","m_factor","beta","pop_size_b0","pop_size_b1","pop_size_b2","pop_size_b3")])
  reg_matrix = cbind(rep(1,nrow(reg_matrix)),reg_matrix)
  reg_matrix[is.na(reg_matrix)] = 0
  
  # Correct visits/c based on what was observed in the training data
  if(!is.null(train_visits)){
    test_visits$c_start = (reg_matrix %*% c_coefs)[,1]
    previous_totals = train_visits %>% group_by(cluster_id) %>% summarize(total_visits = sum(demand))
    test_visits = merge(test_visits,previous_totals,by = "cluster_id")
    test_visits$c_start = test_visits$c_start + test_visits$total_visits
    temp = test_visits %>% group_by(cluster_id) %>% summarize(no_visits = n())
    test_visits = merge(test_visits, temp, by="cluster_id")
    test_visits$cluster_visit_id = test_visits$cluster_visit_id - test_visits$no_visits + 1
  }
  
  predict_clients(r_coefs,p_coefs,q_coefs,c_coefs,
                  reg_matrix,test_visits,test_visits$cluster_visit_id,
                  c_start=test_visits$c_start)
}
