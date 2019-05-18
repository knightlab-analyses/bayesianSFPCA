library(parallel)
library(rstan)
library(loo)
library(Matrix)

options(mc.cores = parallel::detectCores())
source('../../sfpca.R')

#### test on simulated code
load("sim_N50_M80.RData", dat <- new.env())
#ls.str(dat) 

## true values of parameters
J = 2
N = length(dat$Y_SPARSE[[1]])
nknots = dat$params[[1]][[J]] 
Q = nknots + 4 # number of basis
K = dat$params[[2]][[J]] # number of PC
sigma_eps = sqrt(dat$params[[5]][[J]]) # SIGMA_OMEGA_true is error variance

Q1 = dat$params[[1]][[1]] + 4
K1 = dat$params[[2]][[1]]

Theta = dat$params[[6]][(1+Q1):(Q+Q1), (1+K1):(K+K1)] # THETA_true
theta_mu = dat$params[[7]][(1+Q1):(Q+Q1)] # MU_true
ALPHA = t(dat$ALPHA[[1]][(1+K1):(K+K1), ]) # ALPHA
Group = ifelse(ALPHA[, 1] > -0.1, 'G1', 'G2') # 49 in G1 and 51 in G2

# basis
phi_t = dat$phi_t[[J]]
phi_t_stacked=NULL
for(i in 1:N){
	phi_t_stacked = rbind(phi_t_stacked, t(phi_t[[i]]))
}


## simulated data
time = Y = MU = F = Omega = Alpha = list() 
ids=rep(1:N,each=1)

for (i in 1:N){
	time[[i]] = dat$TIME_SPARSE[[1]][[i]][[J]]
	Y[[i]] = dat$Y_SPARSE[[1]][[i]][[J]]
	MU[[i]] = dat$MU_SPARSE[[1]][[i]][[J]]
	F[[i]] = dat$F_SPARSE[[1]][[i]][[J]]
	Omega[[i]] = dat$OMEGA_SPARSE[[1]][[i]][[J]]
}

## convert data to data frame
nrows = length(unlist(time))
nvisits = sapply(time, length)
df = data.frame(matrix(rep(0, nrows * 3), nrow=nrows))
colnames(df) = c('id', 'time', 'response')
for (j in 1:nrows){
	df$id = rep(1:N, nvisits)
	df$time = unlist(time)
	df$response = unlist(Y)
}

# add group assigment
id.list = unique(df$id)
df$group = rep(0, nrows)
for (i in 1:length(id.list)){
	df$group[df$id == id.list[i]] = Group[i]
}
save(df, file='simulated_data_N50_M80_wide.rData')

## data preparation
prepared_data = prepare_data(data=df, unique_subject_id = 'id', time_var='time', response='response')
results_basis = basis_setup_sparse(prepared_data=prepared_data, nknots=nknots, orth=TRUE)

## check whether my generated basis is the same as the simulated one: should always be equal
sum(results_basis$orth_spline_basis_sparse_stacked) == sum(phi_t_stacked) 


##################### run sfpca stan model ##############
Nsamples = 1000
Nchains = 3

model_file = "../../sfpca2.stan"
smod = stan_model(model_file)

pca_data <- list(N = prepared_data$num_subjects, K = K, Q = Q, Y = prepared_data$response.list,
		             V = prepared_data$visits.vector, subject_starts = prepared_data$visits.start,
			   	     subject_stops = prepared_data$visits.stop, cov_starts = prepared_data$cov.start,
			   	     cov_stops = prepared_data$cov.stop,
			   	     cov_size = prepared_data$cov.size, B = phi_t_stacked)
set.seed(31)
sa.list = sampling(smod, data= pca_data, iter=Nsamples, chains=Nchains,init="random")
log_lik.list <- extract(sa.list,"log_lik_marg")[[1]]
looic.list = loo(log_lik.list)
save(sa.list,log_lik.list,looic.list,file="sfpca_results_N50_M80.RData")

##################### sfpca results ##############
load("sfpca_results_N50_M80.RData")
sa = sa.list
Sigma = extract(sa,"Sigma",permuted=FALSE)
W = extract(sa,"W",permuted=FALSE)
sigma_eps = extract(sa,"sigma_eps",permuted=FALSE)
theta_mu = extract(sa,"theta_mu",permuted=FALSE)
alpha = extract(sa,"alpha",permuted=FALSE)
Theta = extract(sa,"Theta",permuted=FALSE)

## Reshape parameters and reorient loadings with PCA rotation 
theta_mu_new = array(0, dim=c(Q, Nchains*Nsamples/2))
alpha_old = alpha_new = array(0, dim=c(K, N, Nchains*Nsamples/2)) 
Theta_old = Theta_new = array(0, dim=c(Q, K, Nchains*Nsamples/2))
W_old = array(0, dim=c(Q, Q, Nchains*Nsamples/2)) 

ind = 0
for(i in 1:dim(W)[1]){
	for(j in 1:dim(W)[2]){
		ind = ind + 1
		theta_mu_new[,ind] = array(theta_mu[i,j,])
		alpha_old[,,ind] = t(array(alpha[i,j,],dim=c(N, K)))
		Theta_old[,,ind] = array(Theta[i,j,],dim=c(Q, K))
		W_old[,,ind] = array(W[i,j,],dim=c(Q,Q)) 

		eigen_temp_sigma=eigen(W_old[,,ind])
		v_temp=eigen_temp_sigma$vectors
		d_temp=eigen_temp_sigma$values 

		for(com in 1:length(d_temp)){
			if(!(d_temp[com]-Re(d_temp[com])==0)){
				d_temp[com]=-1*10^5
			}
		}
		pos_temp=array(0,dim=c(K,1))
		for(pos in 1:K){
			pos_temp[pos]=(1:length(d_temp))[max(d_temp)==d_temp]
			d_temp[pos_temp[pos]]=-1e+5
		}

		Theta_new[,,ind]=v_temp[,pos_temp]
		for(k in 1:K){
			Theta_new[, k, ind]=sign(Theta_new[1,k,ind]) * Theta_new[,k,ind]
		}

		alpha_new[,, ind] = t(Theta_new[,,ind]) %*% Theta_old[,,ind] %*% alpha_old[,,ind]
	}
}

### need to check: Theta_new[,,ind] %*% alpha_new[,, ind] is identital with Theta_old[,,ind] %*% alpha_old[,, ind]
new = Theta_new[,,ind] %*% alpha_new[,, ind]; old = Theta_old[,,ind] %*% alpha_old[,, ind]
which (which (new == old) == FALSE) # should return integer(0), and it did!

save(alpha_new, theta_mu_new, Theta_new, results_basis, prepared_data, 
	 file="estimates_N50_M80.RData")

######### plot sfpca results ######
load('estimates_N50_M80.RData')	

ALPHA_array = alpha_new
MU_array = theta_mu_new
THETA_array = Theta_new
phi_t_cont = results_basis$orth_spline_basis_cont
time_cont = results_basis$time_cont
N = prepared_data$num_subjects
MU_true = dat$params[[7]][(1+Q1):(Q+Q1)]
FPC_true = dat$params[[8]][, (1+K1):(K+K1)]

time = Y = list() 
ids=rep(1:N,each=1)

for (i in 1:N){
	time[[i]] = dat$TIME_SPARSE[[1]][[i]][[J]]
	Y[[i]] = dat$Y_SPARSE[[1]][[i]][[J]]
}

time_sparse = time # coming from wide format: sampling time point for each subject
Y_sparse = Y # coming from wide format: response for each subject
time_unique = dat$time
nloop=dim(ALPHA_array)[3]
first=1
last=nloop

## compare estimated vs. true mean curve ###
MU_mean=MU_array[, 1]
for(iter in 2:nloop){
	MU_mean = MU_mean + MU_array[, iter]
}
MU_mean=cbind(MU_mean/(last-first+1))
Mu_true_functions=t(bdiag(cbind(phi_t_cont)))%*%MU_true
Mu_functions=t(bdiag(cbind(phi_t_cont)))%*%MU_mean
Mu1_true=Mu_true_functions[1:length(time_cont)] # referring to 1st block
Mu1=Mu_functions[1:length(time_cont)]

pdf('Mean_trueVsestimated.pdf', width = 4, height = 4)
plot(time_cont,Mu1_true,type="l",ylim=c(-10,10), xlim=c(0,1),lwd=2,col=4, 
	 xlab='time', ylab='Response', font.lab=2, cex.lab=1.2)
for(i in 1:N){
	lines(time_sparse[[i]],Y_sparse[[i]],type="l",lwd=.25)
}
lines(time_cont,Mu1,type="l",col=2,lwd=2)
#title('True vs. Estimated Mean')
legend('topright', c('True', 'Estimated'), col=c(4, 2), lty=c(1,1), bty='n')
dev.off()


### compare estimated vs. true FPC functions
THETA_mean = THETA_array[,,first]
for(iter in 2:nloop){
	THETA_mean = THETA_mean + THETA_array[,,iter]
}
THETA_mean=cbind(THETA_mean/(last-first+1))

FPC1_mean=t(phi_t_cont)%*%THETA_mean # for 1st block

pdf('FPCs_trueVsestimated.pdf', width = 4, height = 4)
plot(time_unique,FPC_true[,1],type="l",lwd=2,ylim=c(min(FPC_true),max(FPC_true)), 
	 xlab='time', ylab='PC curve', font.lab=2, cex.lab=1.2)
lines(time_unique,FPC_true[,2],type="l",lwd=1)
lines(time_cont,FPC1_mean[,1],type="l",lwd=2,col=2)
lines(time_cont,FPC1_mean[,2],type="l",lwd=1,col=2)
#title(main="True vs. Estimated PFCs")
legend('topright', c('True', 'Estimated'), col=c(1, 2), lty=c(1,1), bty='n')
dev.off()

### plot estimated PC scores ###
ALPHA_mean = ALPHA_array[,,first] 
for(iter in 2:nloop){
	ALPHA_mean = ALPHA_mean + ALPHA_array[,,iter]
}
ALPHA_mean=cbind(ALPHA_mean/(last-first+1))

df = prepared_data$data
Y_sparse = list()
time_sparse = list()
scores = data.frame(t(ALPHA_mean)) 
names(scores)=c("fpc1","fpc2")
df$fpc1=0 # principle component scores
df$fpc2=0

i = 0
for (pid in unique(df$ID)){
	i = i + 1
	Y_sparse[[i]] = df$response[df$ID == pid]
	time_sparse[[i]] = df$time[df$ID == pid]
	df$fpc1[df$ID == pid] = scores[i, 1]
	df$fpc2[df$ID == pid] = scores[i, 2]
}

Fits_sparse=list()
for(i in 1:N){
	Fits_sparse[[i]] = t(phi_t[[i]]) %*% MU_mean + t(phi_t[[i]]) %*% THETA_mean %*% ALPHA_mean[, i]
}

pdf('fpcScores_scatterplots.pdf', width=4, height=4)
colors = rep('black', length(df$group))   
colors[df$group == 'G1'] = 'pink'
colors[df$group == 'G2'] = 'green'
plot(df$fpc1, df$fpc2, pch=20, cex=1,
     xlab='PC1 score', ylab='PC2 score', col=colors,
     font.lab=2, cex.lab=1.2)
legend(x='topright', legend=c("G1", 'G2'), 
       col=c('pink', 'green'), pch=c(20,20), cex=0.8)
dev.off()
