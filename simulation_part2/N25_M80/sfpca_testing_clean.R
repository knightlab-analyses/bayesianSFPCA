library(parallel)
library(rstan)
library(loo)
library(Matrix)

options(mc.cores = parallel::detectCores())
source('../../../../sfpca.R')

#### test on simulated code
load("sim_N100_Orth.RData", dat <- new.env())
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

# compare to Knights method
library(splinectomeR)
result = permuspliner(data=df, xvar='time', yvar='response', cases='id',
	                  category='group', perms = 999, retain_perm = T, quiet = T)
result$pval # 0.001
#permuspliner.plot.permdistance(result, xlabel='time')
permuspliner.plot.permsplines(result, xvar = 'time', yvar = 'response')
ggsave('SplineR_GroupDiff.pdf')

# strange two-group plot: export data to check
write.table(df, 'dat_N25C20_long.txt', sep='\t', row.names=FALSE)

## plot ground truth scatterplot and spaghetti plot
# reference: https://www.r-bloggers.com/my-commonly-done-ggplot2-graphs/
library(ggplot2)
ggplot(df, aes(x=time, y=response)) + geom_line() + guides(colour=FALSE) + 
xlab("Observation Time Point") + ylab("Response") + aes(colour = factor(id)) +
facet_wrap(~ group) + # by different groups
geom_smooth(se=FALSE, colour="black", size=1) # add mean smooth line
ggsave('Observed_spaghetti_groups.pdf')


scores_true = cbind.data.frame(as.vector(ALPHA[,1]), as.vector(ALPHA[,2]), Group, stringsAsFactors = FALSE)
colnames(scores_true ) = c('fpc1', 'fpc2', 'group')
ggplot(scores_true , aes(x=fpc1 , y=fpc2)) + geom_point() + aes(colour=factor(group)) 
ggsave('scatterplot_scores_truth.pdf')

pdf('fpcScores_true_scatterplots.pdf')
colors = rep('black', length(df$group))   
colors[scores_true$group == 'G1'] = 'pink'
colors[scores_true$group == 'G2'] = 'green'
plot(scores_true$fpc1, scores_true$fpc2, pch=20, cex=2,
     xlab='FPC1 Score', ylab='FPC2 Score', col=colors,
     main='True scores')
legend(x='topright', legend=c("G1", 'G2'), 
       col=c('pink', 'green'), pch=c(16,16), cex=0.8)
dev.off()

## data preparation
prepared_data = prepare_data(data=df, unique_subject_id = 'id', time_var='time', response='response')
results_basis = basis_setup_sparse(prepared_data=prepared_data, nknots=nknots, orth=TRUE)

## check whether my generated basis is the same as the simulated one: should always be equal
sum(results_basis$orth_spline_basis_sparse_stacked) == sum(phi_t_stacked) 


##################### run sfpca stan model ##############
Nsamples = 1000
Nchains = 3

model_file = "../../../../sfpca2.stan"
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
save(sa.list,log_lik.list,looic.list,file="sfpca_results_N100_Orth_G2.RData")

################## predictive model checking (need to re-run Stan as it was updated with predictive distribution)
sa = sa.list
Ynew = extract(sa,"Ynew",permuted=FALSE)
V = prepared_data$visits.vector
Ynew_transform = matrix(rep(0, Nsamples/2 * Nchains * sum(V)), ncol=sum(V))
ind = 0
for (i in 1:(Nsamples/2)){
	for (j in 1:Nchains){
		ind = ind + 1
		Ynew_transform[ind, ] = Ynew[i,j,]
	}
}

library("bayesplot")
library("ggplot2")
pdf('predictive_cheking.pdf')
ppc_dens_overlay(df$response, Ynew_transform)
dev.off()
####################

load("sfpca_results_N100_Orth_G2.RData")

sa = sa.list
# traceplot(sa, pars=c("sigma_eps")) # mixed well
# traceplot(sa, pars=c("theta_mu")) # mixed well
# traceplot(sa, pars=c("alpha")) # not, due to unidentifiability
# traceplot(sa, pars=c("Theta")) # not, due to unidentifiability
# traceplot(sa, pars=c("W")) # mixed well: still determined depsite of undetermined Theta
# traceplot(sa, pars=c("Sigma")) # mixed well: still determined depsite of undetermined Theta


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
	 file="estimates_N100_Orth_G2.RData")


##### check estimation of Theta ######

Theta_table= array(0,dim=c(Q,K,3))
for(q in 1:Q){
	for(k in 1:K){
		Theta_table[q,k,] = quantile(Theta_new[q,k,],probs=c(.025,.5,.975))
	}
}

print(Theta_table)

, , 1

           [,1]        [,2]
[1,] 0.03252941  0.87544116
[2,] 0.19101316  0.33637887
[3,] 0.52224252 -0.22195727
[4,] 0.53408304 -0.23336275
[5,] 0.47883126 -0.07217950
[6,] 0.27547901 -0.07599339

, , 2

          [,1]        [,2]
[1,] 0.1345290  0.89595494
[2,] 0.2345997  0.37048392
[3,] 0.5437442 -0.16166684
[4,] 0.5572063 -0.17372493
[5,] 0.4879849 -0.02101163
[6,] 0.2850650 -0.04394199

, , 3

          [,1]        [,2]
[1,] 0.2266213  0.90872296
[2,] 0.2722366  0.40022899
[3,] 0.5587670 -0.09817273
[4,] 0.5736565 -0.10809415
[5,] 0.4965154  0.03433805
[6,] 0.2928808 -0.01054144



### True value of Theta 
# Theta = dat$params[[6]][(1+Q1):(Q+Q1), (1+K1):(K+K1)]

          [,1]        [,2]
[1,] 0.1734486  0.79569747
[2,] 0.2471184  0.31964112
[3,] 0.5260293 -0.16640001
[4,] 0.5316424 -0.17366537
[5,] 0.4836753 -0.04200097
[6,] 0.2758460 -0.05197922


alpha_table= array(0,dim=c(K,N,3))
for(k in 1:K){
	for(n in 1:N){
		alpha_table[k,n,] = quantile(alpha_new[k,n,],probs=c(.025,.5,.975))
	}
}

print(alpha_table[, 1:10,]) # look at only the first 10 subjects
, , 1

           [,1]       [,2]       [,3]       [,4]       [,5]       [,6]
[1,] -1.1812193 -0.3466107 -1.6591680 -1.1333168 -0.1143988 -2.4015149
[2,] -0.3745052 -0.7553360 -0.4661914  0.8155803 -0.2880318  0.4650177
           [,7]      [,8]       [,9]     [,10]
[1,] -4.5390154 -2.886346 -0.8091383  0.986499
[2,] -0.8986483 -0.146596 -1.3904728 -1.045742

, , 2

            [,1]        [,2]        [,3]       [,4]      [,5]       [,6]
[1,] -0.76823567  0.08525775 -1.24714329 -0.6851386 0.3013175 -1.9785236
[2,]  0.02710019 -0.56082946 -0.07485567  1.0109616 0.1740594  0.7520609
           [,7]       [,8]       [,9]      [,10]
[1,] -4.0823438 -2.4766740 -0.3647094  1.4407063
[2,] -0.4087697  0.1787385 -0.6164185 -0.8016778

, , 3

           [,1]       [,2]       [,3]       [,4]      [,5]      [,6]
[1,] -0.3508846  0.5142787 -0.8350015 -0.2218252 0.7101789 -1.510091
[2,]  0.4278610 -0.3655328  0.3014285  1.2082740 0.6799973  1.048718
            [,7]       [,8]       [,9]      [,10]
[1,] -3.62559482 -2.0576273 0.08151103  1.8617919
[2,]  0.05551879  0.5047831 0.18504730 -0.5405619

## true value
dat$ALPHA[[1]][3:4, 1:10] ########

          [,1]        [,2]        [,3]       [,4]       [,5]      [,6]
[1,] -0.802574  0.07270469 -1.27195980 -0.5432805  0.3010709 -2.057059
[2,] -0.335781 -0.58760567 -0.05514026  1.1458458 -0.2241192  1.041724
           [,7]       [,8]       [,9]      [,10]
[1,] -4.1710520 -2.5832875 -0.5900854  1.5061697
[2,] -0.2560843  0.3270311 -1.5527667 -0.8738524

######### plotting ######
load('estimates_N100_Orth_G2.RData')	

ALPHA_array = alpha_new
MU_array = theta_mu_new
THETA_array = Theta_new
phi_t_cont = results_basis$orth_spline_basis_cont
time_cont = results_basis$time_cont
N = prepared_data$num_subjects
MU_true = dat$params[[7]][(1+Q1):(Q+Q1)]
FPC_true = dat$params[[8]][, (1+K1):(K+K1)]

# re-run the following to make sure
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

# pdf('plot_stan.pdf')
# par(mfcol=c(2,2))
#plot mean functions
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


#plot PC functions
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


####### additional plots #########
load('estimates_N100_Orth_G2.RData')	
K = npcs = 2  # revised manually; think about saving its value, also Q, earlier

ALPHA_array = alpha_new
MU_array = theta_mu_new
THETA_array = Theta_new
phi_t_cont = results_basis$orth_spline_basis_cont
phi_t = results_basis$orth_spline_basis_sparse
time_cont = results_basis$time_cont

nloop=dim(ALPHA_array)[3]
first=1
last=nloop
N = prepared_data$num_subjects


MU_mean = MU_array[, first] #mean function across sampling sessions
ALPHA_mean = ALPHA_array[,,first] # mean factor scores
THETA_mean = THETA_array[,,first] # mean factor loading

for(iter in 2:nloop){
	MU_mean = MU_mean + MU_array[, iter]
	ALPHA_mean = ALPHA_mean + ALPHA_array[,,iter]
	THETA_mean = THETA_mean + THETA_array[,,iter]
}

MU_mean=cbind(MU_mean/(last-first+1))
ALPHA_mean=cbind(ALPHA_mean/(last-first+1))
THETA_mean=cbind(THETA_mean/(last-first+1))

Mu_functions = t(bdiag(cbind(phi_t_cont)))%*%MU_mean
FPC_mean=t(phi_t_cont)%*%THETA_mean

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

pdf('mean_spaghetti.pdf')
plot(time_cont, Mu_functions, type="l", ylim=c(min(unlist(Y_sparse)), 
	 max(unlist(Y_sparse))), xlim=c(0,1), ylab='Standardized Response', lwd=5, col=4)
for(i in 1:N){
	lines(time_sparse[[i]],Y_sparse[[i]],type="l",lwd=.25)
}
title(main='Response')
dev.off()


pdf('spaghetti_smoothed.pdf')
par(mfrow=c(1,2))
plot(time_cont, Mu_functions,type="n", ylim=c(min(unlist(Y_sparse)), max(unlist(Y_sparse))),
	 lwd=5, col=4, ylab='Standardized Response', xlab='Time')
for(i in 1:N){
	lines(time_sparse[[i]], Y_sparse[[i]], type="l", lwd=.25)
}
title(main='Observed Response')

plot(time_cont, Mu_functions, type="n", ylim=c(min(unlist(Y_sparse)), max(unlist(Y_sparse))),
	 lwd=2,col=4, ylab='Standardized Response', xlab='Time')
for(i in 1:N){
		lines(time_sparse[[i]],Fits_sparse[[i]], type="l", lwd=.25)
}
title(main='Smoothed Response')
dev.off()

### spaghetti plots with mean curves for each group ###
df_more = df
df_more$Y_sparse = unlist(Y_sparse) # check: sum(df_more$Y_sparse != df_more$response) == 0
df_more$Fits_sparse = unlist(Fits_sparse)

library(ggplot2)
ggplot(df_more, aes(x=time, y=Fits_sparse)) + geom_line() + guides(colour=FALSE) + 
xlab("Observation Time Point") + ylab("Smoothed Response") + aes(colour = factor(id)) +
facet_wrap(~ group) + # by different groups
geom_smooth(se=FALSE, colour="black", size=1) # add mean smooth line
ggsave('Smoothed_spaghetti_groups.pdf')

########## FPC curves seem strange ########
pdf('FPCs.pdf', width=5, height=8)
plot(time_cont, FPC_mean[, 1], type="l", lwd=2, ylim=c(min(unlist(Y_sparse)), max(unlist(Y_sparse))),
     xlab='Time', ylab='Value of PC Curve', cex.lab=1.5)
lines(time_cont, FPC_mean[, 2],type="l",lwd=1, col='blue')
title(main=paste("FPCs for", "Response"))
legend('topright', c('PC1', 'PC2'), lwd=c(3, 3), col=c('black', 'blue'), bty='n')
dev.off()

pdf('FPCs_mean.pdf', width=5, height=8)
par(mfrow=c(2,1))
for (k in 1:K){
	plot(time_cont, Mu_functions, type="l", ylim=c(min(unlist(Y_sparse)), max(unlist(Y_sparse))), 
	lwd=2,col=1, xlab='Time', ylab='Standardized Response', cex.lab=1.5)
	lines(time_cont, Mu_functions + FPC_mean[,k],type="l",lwd=3,lty=2,col=2) # red
	lines(time_cont, Mu_functions - FPC_mean[,k],type="l",lwd=3,lty=2,col=3) # green
	title(main=paste("Effect of FPC", k , "for", "Response"))
}
legend('topright', c('+ pc', '- pc'), lty=c(2,2), lwd=c(3,3), col=c(2, 3), bty='n', cex=0.5)
dev.off()


#### predicted trajectories do not seem as good as before
pdf('predictedTrajectories.pdf')
par(mfrow=c(2,1))
i = 0
for (pid in unique(df$ID)){
	i = i + 1
  	Y_i = Y_sparse[[i]]
  	times_i = time_sparse[[i]]
  	Fitted_i = Fits_sparse[[i]]
  	plot(time_cont, Mu_functions, type="n", ylim=c(min(unlist(Y_sparse)), max(unlist(Y_sparse))), 
  	   	 lwd=2, xlim=c(0, 1), col=4, xlab="Time", ylab="Standardized Response")
  	lines(times_i, Y_i, type="l", col=4, lty=2)
  	lines(times_i, Fitted_i, type="l", col=4)
  	title(main=paste("Fitted", "Response", "for subj", pid))
}
dev.off()

########## analysis with principle component scores ###### 
#### (warning: df is long formatted; convert to wide format for lm? or take unique ID only?)
pdf('fpcScores.pdf')
par(mfrow=c(1, 2))
boxplot(fpc1 ~ group, data=df, main='Response (PC1)', ylab='FPC1 Score')
boxplot(fpc2 ~ group, data=df, main='Response (PC2)', ylab='FPC2 Score')
dev.off()

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

# simple linear regression
fit.1 = lm(df$fpc1 ~ df$group)
summary(fit.1)
# library(xtable)
# obj.1 = xtable(fit.1)
# print.xtable(obj.1, type="latex", file="fit1.tex")

# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.63473    0.09816   16.65   <2e-16 ***
# df$groupG2  -3.14595    0.13612  -23.11   <2e-16 ***

fit.2 = lm(df$fpc2 ~ df$group)
summary(fit.2)

# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)  0.06988    0.06654   1.050    0.294
# df$groupG2  -0.12582    0.09227  -1.364    0.174

### compare predicted and true PC scores
# fpc scores are the same for same subject at different visits
N = length(unique(df$id))
df_scores = data.frame(matrix(rep(0, N*4), nrow=N))
colnames(df_scores) = c('id', 'fpc1', 'fpc2', 'group')
df_scores$id = unique(df$id)
pid.list = df_scores$id
for (i in pid.list){
	df_scores$fpc1[df_scores$id == i] = df$fpc1[df$id == i][1]
	df_scores$fpc2[df_scores$id == i] = df$fpc2[df$id == i][1]
	df_scores$group[df_scores$id == i] = df$group[df$id == i][1]
}


b1 = ggplot(scores_true, aes(x = group, y = fpc1, col=group)) + geom_boxplot() + xlab('true') + 
	 geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)
b2 = ggplot(scores_true, aes(x = group, y = fpc2, col=group)) + geom_boxplot() + xlab('true') + 
	 geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)
b3 = ggplot(df_scores, aes(x = group, y = fpc1, col=group)) + geom_boxplot() + xlab('fit') + 
	 geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)
b4 = ggplot(df_scores, aes(x = group, y = fpc2, col=group)) + geom_boxplot() + xlab('fit') + 
	 geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)

require(gridExtra) # grid for ggplot: https://stackoverflow.com/questions/1249548/side-by-side-plots-with-ggplot2
pdf('boxplots_scores_compare.pdf')
grid.arrange(b1, b2, b3, b4, ncol=2, as.table = FALSE)
dev.off()

##### in fact: the version using repeted fpc scores for the same subject works fine
b1 = ggplot(scores_true, aes(x = group, y = fpc1, col=group)) + geom_boxplot() + xlab('true') 
b2 = ggplot(scores_true, aes(x = group, y = fpc2, col=group)) + geom_boxplot() + xlab('true') 
b3 = ggplot(df, aes(x = group, y = fpc1, col=group)) + geom_boxplot() + xlab('fit')
b4 = ggplot(df, aes(x = group, y = fpc2, col=group)) + geom_boxplot() + xlab('fit') 
pdf('boxplots_scores_compare_v2.pdf')
grid.arrange(b1, b2, b3, b4, ncol=2, as.table = FALSE)
dev.off()










