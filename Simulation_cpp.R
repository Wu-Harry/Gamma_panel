
#n = 100              # n= 200, 400;
#case = 1         # case = 1, 2;
#beta2 = 1          # beta2 = 1, -1;
#theta = 0.5          # theta = 0.5, 1.0, 1.5;
#a = 2

args <- as.numeric(commandArgs(TRUE))
#args = c(100, 1, 1, 1, 1)

#rm(list=ls(all=TRUE)); #Args = c(100, 1, 1, 1, 1);

n = args[1]              # n= 200, 400;
case = args[2]         # case = 1, 2;
beta2 = args[3]          # beta2 = 1, -1;
theta = args[4]          # theta = 0.5, 1.0, 1.5;
a = args[5]                # a = 500;


library(statmod)
library(Rcpp)
# library(ggplot2)

#sourceCpp("/Users/wuhanyu/Desktop/Research/LBS/Gamma_frailty/code/E_step.cpp")
#sourceCpp("/Users/wuhanyu/Desktop/Research/LBS/Gamma_frailty/code/M_step_1.cpp")
#sourceCpp("/Users/wuhanyu/Desktop/Research/LBS/Gamma_frailty/code/M_step_2.cpp")
sourceCpp("/root/Simu/code/E_step.cpp")
sourceCpp("/root/Simu/code/M_step_1.cpp")
sourceCpp("/root/Simu/code/M_step_2.cpp")

h_n0 = 1; beta = c(1,beta2); tau = 4; t0 = 0.6;

if(case == 1){ lambda0=function(t){0.5*t^0}; Lambda0=function(t){0.5*t} }else{
  lambda0=function(t){t};   Lambda0=function(t){t^2/2} }

Grid = seq(t0,tau,length.out=200) # 用于计算\Lambda和画图

set.seed(123+2023*a)

##########################
##    Generate Data     ##
##########################

K = NULL; T = matrix(0,n,7); T.vec = NULL;

for(i in 1:n){
  T0 = runif(1,t0,2);
  K_i = (tau - T0) %/% t0 + 1;
  K = c(K,K_i);
  T[i,2:(K_i+1)] = seq(T0, T0+0.6*(K_i -1), 0.6);
  T.vec = c(T.vec, seq(T0, T0+0.6*(K_i -1), 0.6));
}

C = apply(T,1,max)  # stopping time
grid = unique(sort(T.vec))
m = length(grid) # length of t

X = cbind( rbinom(n,1,0.5), rnorm(n,0,1) ); X_len = dim(X)[2];

Lambda <- function(t,i){ 
  if(t!=0){ return( Lambda0(t)*exp( X[i,] %*% beta ) ) } 
  if(t==0){ return(0) }
}

delt = matrix(0,n,6);  delt.vec = NULL; N = matrix(0,n,6);
Xi = rgamma(n, shape=1/theta, scale = theta); # rate = 1/theta

for(i in 1:n){
  Delt = NULL;   
  time = cbind( T[i, 2:(K[i]+1)], T[i, 1:K[i] ] );
  for(j in 1: K[i] ){
    Num = rpois(1, Xi[i]*( Lambda(time[j,1],i) - Lambda(time[j,2],i) ));
    if(Num == 0){ Delt = c( Delt, 0 ) };
    if(Num >= 1){ Delt = c( Delt, 1 ) };
    if(j >= 2){N[i,j] = N[i,j-1] + Num} else{N[i,j] = Num};
  }
  delt[i, 1: K[i]] = Delt;
  delt.vec = c(delt.vec, Delt);
}


##########################
##    Initial Value     ##
##########################

#max(N) / (X[which.max(N)%%dim(N)[1], ] %*% beta) / m
# 这里初始值设定也有问题
initial = c(1/theta, rep(1,X_len), rep(max(N) / m, m)) # c(1/theta, beta, lambda)

#########################################
#######  E-M Algorithm Iteratoion  ######
#########################################

maxiter=1000; ep=1e-3; iter=1; error=1; Converge_Index1=0;
theta_h = initial[1]; beta_h = initial[2:(X_len+1)]; lam_h = initial[(X_len+2):length(initial)];
Q = 50; gq = gauss.quad(Q, kind="laguerre"); nodes = gq$nodes;  whts = gq$weights;
print("EM iteration start...")

theta_0 = theta_h; beta_0 = beta_h; lam_0 = lam_h;
#theta_0 = 1/theta; beta_0 = beta; lam_0 = lam_h;
max_lam = vector()
up_beta = matrix(0,maxiter,2)
up_theta = vector()

start_time <- Sys.time()
while( iter <= maxiter & error > ep  ){
  
  ##########################
  ##        E_step        ##
  ##########################  
  
  Exp_Xbeta = c( exp( X %*% beta_0) );
  E_xi = rep(0,n); E_log_xi = rep(0,n); E_W = matrix(0, nrow=n, m);  # posterior mean
  
  E_all = E_step(T, delt, K, lam_0, c(1/theta, beta_0), grid, Exp_Xbeta, nodes, whts, n, 0)
  E_xi = E_all[,m+1]
  E_log_xi = E_all[,m+2]
  E_W = E_all[,1:m]
  
  #cat(max(E_W))
  #E_W[is.na(E_W)] = 0
  #E_log_xi[is.na(E_log_xi)] = 0
  #E_xi[is.na(E_xi)] = 0
  print("E-step done...")
  
  ##########################
  ##        M_step        ##
  ##########################
  
  ###################
  ##  update theta ##
  ###################
  
  #theta_new = theta_0 - ( mean(E_log_xi) - mean(E_xi) - digamma(theta_0) + 1 - log(1/theta_0) )/( 1/theta_0 - trigamma(theta_0) );
  theta_new = theta_0 + ( mean(E_log_xi) - mean(E_xi) - digamma(theta_0) + 1 + log(theta_0) );
  #theta_new = exp(log(theta_0) - (mean(E_log_xi) - mean(E_xi) + log(theta_0) + 1 - digamma(theta_0)) / (1 - trigamma(theta_0)*theta_0));
  #theta_new = exp(log(theta_0) - (mean(E_log_xi) - mean(E_xi) + log(theta_0) + 1 - digamma(theta_0)) / (mean(E_log_xi) - mean(E_xi) + log(theta_0) + 2 - (trigamma(theta_0)*theta_0 + digamma(theta_0))));
  #THETA = c(THETA, theta_new); 
  error1  = abs(theta_0-theta_new); 
  if(theta_new<=0) break;
  theta_0 = theta_new;
  up_theta[iter] = theta_0;
  # theta_0 = ifelse(theta_new>=0.005,theta_new,0.005) # 防止theta负
  cat("theta =", 1/theta_0,"\n");
  
  ###################
  ##  update beta  ##
  ###################
  
  out = M_step_1(E_W, E_xi, X, grid, C, beta_0, 1/theta, m, n) # theta先设置为真值
  
  H = out[,1:X_len];
  grad = out[,X_len+1];
  if( sum( is.na(grad) ) >= 1 ) break;
  
  if(det(H) > 0.0001 | det(H) < -0.0001) {
    beta_new = beta_0 - solve(H)%*%grad} 
  else{ print("H is Invertible");
    beta_new = beta_0 - solve(H - diag(0.00001, X_len)) %*% grad;}
  cat("M_step done...\n"); cat("beta =",c(beta_new),"\n")
  up_beta[iter,] = beta_new
  error = sqrt( sum((beta_0 - beta_new)^2)); 
  beta_0 = c(min(max(beta_new[1],-3),3),min(max(beta_new[2],-3),3)); 
  
  ###################
  ## update lambda ##
  ###################
  
  theta_0 = out[1,X_len+2];
  lam_new = M_step_2(E_W, E_xi, X, grid, C, beta_0, m, n);
  lam_0 = lam_new;
  #lam_0[which(lam_0==0)] = 1e-8 #in case min(lam_0)=0
  if(sum(is.na(lam_0))>=1){cat("NA in lam_0"); break}; if(max(lam_0)>10000){cat("max(lam) too big"); break}
  cat("max(lam) =",max(lam_0),"\n");
  max_lam[iter] = max(lam_0);
  lam_0[which(lam_0>0.2)] = 0.2; # 防止lam过大，做一个截断
  cat("beta error =",error,", theta error =",error1,", interation",iter,"done\n");
  
  if( (sum( is.na(beta_new) ) >= 1) | sum(abs(beta_new) >= 300 ) >=1) break;
  
  if( error < ep & error1 < ep ){ Converge_Index1 = 1; Conv_knot = 1; break } 
  #print(iter); # print( c(theta_new, eta.hat) );
  #cat("min lam =",min(lam_0))
  iter = iter + 1;
}
end_time <- Sys.time() 
cat("time =", end_time - start_time)


Lam_y = vector(mode = "numeric", length = 200)
for (l in 1:(m-1)) {
  l_index = which((Grid >= grid[l])&(Grid < grid[l+1]))
  Lam_y[l_index] = sum(lam_0[1:l])
}
Lam_y[which(Grid >= grid[m])] = sum(lam_0[1:m]);
lam_plot = data.frame(grid = Grid, est_lam = Lam_y, real_lam = Lambda0(Grid));
# ggplot(data = as.data.frame(lam_plot)) + geom_line(mapping = aes(x=grid,y=real_lam),color="red") + geom_line(mapping = aes(x=grid,y=est_lam),color="blue") + theme_bw() + theme( panel.grid.minor = element_blank(), panel.border = element_blank(), plot.title = element_text(hjust = 0.5)) + ylab("Cumulative baseline hazard function") + labs(title = "Risk Curve Comparison")



##########################
##          SE          ##
########################## 

cat("SE begin...")

logPl_func = function(Y){
  theta_pl = Y[1]; beta_pl = Y[2:(X_len+1)]; lam_pl = initial[(X_len+2):length(initial)];
  
  maxiter=200; ep=1e-2; iter=1; error=1; Converge_Index1=0; 
  Q = 60; gq = gauss.quad(Q, kind="laguerre"); nodes = gq$nodes;  whts = gq$weights;
  
  while( iter <= maxiter & error > ep  ){
    
    ##########################
    ##        E_step        ##
    ##########################  
    
    Exp_Xbeta = c( exp( X %*% beta_pl) );
    E_xi = rep(0,n); E_log_xi = rep(0,n); 
    E_W = matrix(0, nrow=n, m);  # posterior mean
    Loglike = 0;
    
    E_all = E_step(T, delt, K, lam_pl, c(1/theta, beta_pl), grid, Exp_Xbeta, nodes, whts, n, 1);
    E_xi = E_all[,m+1];
    E_log_xi = E_all[,m+2];
    E_W = E_all[,1:m];
    Loglike = E_all[1,m+3];
    cat("E_step done ");
    
    ##########################
    ##        M_step        ##
    ##########################
    
    ###################
    ## update lambda ##
    ###################
    
    lam_new = M_step_2(E_W, E_xi, X, grid, C, beta_pl, m, n);
    error = sqrt( sum((lam_pl - lam_new)^2));
    lam_pl = lam_new
    lam_pl[which(lam_pl==0)] = 1e-8 # 防止lam=0 
    # lam_pl[which(lam_pl>0.2)] = 0.2; # 防止lam过大，做一个截断
    
    cat("lambda error =",error,", interation",iter,"done\n")
    if( error < ep ){ Converge_Index1 = 1; Conv_knot = 1; break } 
    #print(iter); # print( c(theta_new, eta.hat) );
    
    iter = iter + 1;
  }
  
  return(Loglike)
}


################################
h_n = h_n0/sqrt(n);  I_matr = matrix(0, nrow=X_len+1, ncol=X_len+1);
para_hat = c(theta_0, beta_0[1:X_len]);  

start_time <- Sys.time();
logPl_hat = logPl_func(para_hat);

for(i in 1:(X_len+1) ){
  for(j in i:(X_len+1) ){
    v1=v2=rep(0, X_len+1); v1[i]=1; v2[j]=1;
    I_matr[j,i]=I_matr[i,j]=-(logPl_func(para_hat+h_n*v1+h_n*v2)-logPl_func(para_hat+h_n*v1)
                              -logPl_func(para_hat+h_n*v2)+logPl_hat )/(n*(h_n)^2);
    print(c(i,j))
  }
}
end_time <- Sys.time();
cat("time =", end_time - start_time);

if(sum(is.na(I_matr)) >0 ) {print("NA values happend in SE.prof")}

if( abs( det(I_matr)) > 0.00001 ){
  SE_prof = sqrt( diag( solve(I_matr) ) )/sqrt(n); 
}else{
  SE_prof = sqrt( diag( solve(I_matr + diag(0.0001, X_len+1)) ) )/sqrt(n);
}

if(sum(is.na(SE_prof)) >0 ) {print("NA values happend in SE.prof")}


###################################
##    95% Confidence Interval    ##
##  lower bound  &  upper bound  ##
###################################
lower1 = para_hat - qnorm(0.975,0,1)*SE_prof;
upper1 = para_hat + qnorm(0.975,0,1)*SE_prof;
Cpi = 1*( c(1/theta,beta) >= lower1 )*( c(1/theta,beta) <= upper1 );

###########################
##  Simulation Summary   ##
###########################

para_hat = c(1/theta_0, beta_0[1:X_len]);
Output = c(para_hat,SE_prof,Cpi)

if(file.exists("./result.csv")){
  result = read.csv("./result.csv")
  result = rbind(result, Output)
  write.csv(result, "./result.csv", row.names = F, quote = F)
} else{ write.csv(t(Output), "./result.csv", row.names = F, quote = F) }

if(file.exists("./lambda.csv")){
  res_lam = read.csv("./lambda.csv")
  res_lam = rbind(res_lam, Lam_y)
  write.csv(res_lam, "./lambda.csv", row.names = F, quote = F)
} else{ write.csv(t(Lam_y), "./lambda.csv", row.names = F, quote = F) }




