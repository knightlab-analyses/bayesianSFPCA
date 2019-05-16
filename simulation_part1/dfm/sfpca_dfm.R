library(parallel)
library(rstan)
library(loo) # for waic()
library(Matrix) # for bdiag(): construct a block diagnoal matrix
options(mc.cores = parallel::detectCores())


df=read.csv("dfm_bound.txt",header=TRUE, sep='\t')

library(splinectomeR) 
# p-value = 0.336 
result = permuspliner(data=df, xvar='time', yvar='response', cases='patient',
	                  category='Version', groups = c('baseline','shift_1x'), perms = 999, retain_perm = T)
permuspliner.plot.permsplines(result, xvar = 'time', yvar = 'response')
ggsave('SplineR_GroupDiff_base1x_df.pdf')

# p-value = 0.002 
result = permuspliner(data=df, xvar='time', yvar='response', cases='patient',
	                  category='Version', groups = c('baseline','shift_2x'), perms = 999, retain_perm = T)
permuspliner.plot.permsplines(result, xvar = 'time', yvar = 'response')
ggsave('SplineR_GroupDiff_base2x_df.pdf')

# p-value = 0.001 
result = permuspliner(data=df, xvar='time', yvar='response', cases='patient',
	                  category='Version', groups = c('baseline','shift_4x'), perms = 999, retain_perm = T)
permuspliner.plot.permsplines(result, xvar = 'time', yvar = 'response')
ggsave('SplineR_GroupDiff_base4x_df.pdf')

library(ggplot2)
ggplot(df, aes(x=time, y=response)) + geom_line() + guides(colour=FALSE) + 
xlab("Observed Time") + ylab("response") + aes(colour = factor(patient)) +
geom_smooth(se=FALSE, colour="black", size=1) # add mean smooth line
ggsave('Observed_spaghetti.pdf')

library(ggplot2)
colpal <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")
ggplot(df, aes(x=time, y=response, group=patient, colour=Version)) + 
geom_line(alpha=0.3) + xlab("Observed Time") + ylab("response")  +
theme_classic() + theme(axis.text = element_text(color='black')) +
geom_smooth(se=FALSE, size=1.5, aes(group=Version)) +
scale_color_manual(values=colpal)
ggsave('Observed_spaghetti_group.pdf', width = 4, height = 3, dpi=600)

ggplot(df, aes(x=time, y=response)) + geom_line() + guides(colour=FALSE) + 
xlab("Observed Time") + ylab("response") + aes(colour = factor(patient)) +
facet_wrap(~ Version) + # by different groups
geom_smooth(se=FALSE, colour="black", size=1) # add mean smooth line
ggsave('Observed_spaghetti_each.pdf')

# check each individual's change
pdf('Observed_spaghetti_ind.pdf')
colpal <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")
tmp.1 = df[df$patient %in% c('1', '1_1x', '1_2x', '1_4x'), ]
ggplot(tmp.1, aes(x=time, y=response)) + geom_line() + guides(colour=FALSE) + 
xlab("Observed Time") + ylab("response") + aes(colour = factor(patient)) + 
ggtitle('patient1') + scale_color_manual(values=colpal)

tmp.2 = df[df$patient %in% c('2', '2_1x', '2_2x', '2_4x'), ]
ggplot(tmp.2, aes(x=time, y=response)) + geom_line() + guides(colour=FALSE) + 
xlab("Observed Time") + ylab("response") + aes(colour = factor(patient)) + 
ggtitle('patient2') + scale_color_manual(values=colpal)

tmp.3 = df[df$patient %in% c('3', '3_1x', '3_2x', '3_4x'), ]
ggplot(tmp.3, aes(x=time, y=response)) + geom_line() + guides(colour=FALSE) + 
xlab("Observed Time") + ylab("response") + aes(colour = factor(patient)) + 
ggtitle('patient3') + scale_color_manual(values=colpal)

tmp.4 = df[df$patient %in% c('4', '4_1x', '4_2x', '4_4x'), ]
ggplot(tmp.4, aes(x=time, y=response)) + geom_line() + guides(colour=FALSE) + 
xlab("Observed Time") + ylab("response") + aes(colour = factor(patient)) +
ggtitle('patient4') + scale_color_manual(values=colpal)

tmp.5 = df[df$patient %in% c('5', '5_1x', '5_2x', '5_4x'), ]
ggplot(tmp.5, aes(x=time, y=response)) + geom_line() + guides(colour=FALSE) + 
xlab("Observed Time") + ylab("response") + aes(colour = factor(patient)) + 
ggtitle('patient5') + scale_color_manual(values=colpal)

tmp.6 = df[df$patient %in% c('6', '6_1x', '6_2x', '6_4x'), ]
ggplot(tmp.6, aes(x=time, y=response)) + geom_line() + guides(colour=FALSE) + 
xlab("Observed Time") + ylab("response") + aes(colour = factor(patient)) + 
ggtitle('patient6') + scale_color_manual(values=colpal)

tmp.7 = df[df$patient %in% c('7', '7_1x', '7_2x', '7_4x'), ]
ggplot(tmp.7, aes(x=time, y=response)) + geom_line() + guides(colour=FALSE) + 
xlab("Observed Time") + ylab("response") + aes(colour = factor(patient)) + 
ggtitle('patient7') + scale_color_manual(values=colpal)

tmp.8 = df[df$patient %in% c('8', '8_1x', '8_2x', '8_4x'), ]
ggplot(tmp.8, aes(x=time, y=response)) + geom_line() + guides(colour=FALSE) + 
xlab("Observed Time") + ylab("response") + aes(colour = factor(patient)) + 
ggtitle('patient8') + scale_color_manual(values=colpal)

tmp.9 = df[df$patient %in% c('9', '9_1x', '9_2x', '9_4x'), ]
ggplot(tmp.9, aes(x=time, y=response)) + geom_line() + guides(colour=FALSE) + 
xlab("Observed Time") + ylab("response") + aes(colour = factor(patient)) + 
ggtitle('patient9') + scale_color_manual(values=colpal)

tmp.10 = df[df$patient %in% c('10', '10_1x', '10_2x', '10_4x'), ]
ggplot(tmp.10, aes(x=time, y=response)) + geom_line() + guides(colour=FALSE) + 
xlab("Observed Time") + ylab("response") + aes(colour = factor(patient)) + 
ggtitle('patient10') + scale_color_manual(values=colpal)
dev.off()

# selected patients
ggplot(tmp.1, aes(x=time, y=response)) + geom_line(size=1.2) + guides(colour=FALSE) + 
xlab("time") + ylab("response") + aes(colour = factor(patient)) + 
ggtitle('patient1') + scale_color_manual(values=colpal) +
theme(plot.title = element_text(hjust = 0.5)) + ylim(84, 131) +
theme(axis.title.x = element_text(face="bold"),
      axis.title.y = element_text(face="bold"))
ggsave('Observed_patient1.pdf', width = 4, height = 3, dpi=600)

ggplot(tmp.1, aes(x=time, y=response)) + geom_point() + guides(colour=FALSE) + 
xlab("time") + ylab("response") + aes(colour = factor(patient)) + 
ggtitle('patient1') + scale_color_manual(values=colpal) +
theme(plot.title = element_text(hjust = 0.5)) + ylim(84, 131) +
theme(axis.title.x = element_text(face="bold"),
      axis.title.y = element_text(face="bold"))
ggsave('Observed_patient1_point.pdf', width = 4, height = 3, dpi=600)

ggplot(tmp.2, aes(x=time, y=response)) + geom_line(size=1.2) + guides(colour=FALSE) + 
xlab("time") + ylab("response") + aes(colour = factor(patient)) + 
ggtitle('patient2') + scale_color_manual(values=colpal) +
theme(plot.title = element_text(hjust = 0.5)) + ylim(84, 131) +
theme(axis.title.x = element_text(face="bold"),
      axis.title.y = element_text(face="bold"))
ggsave('Observed_patient2.pdf', width = 4, height = 3, dpi=600)

ggplot(tmp.2, aes(x=time, y=response)) + geom_point() + guides(colour=FALSE) + 
xlab("time") + ylab("response") + aes(colour = factor(patient)) + 
ggtitle('patient2') + scale_color_manual(values=colpal) +
theme(plot.title = element_text(hjust = 0.5)) + ylim(84, 131) +
theme(axis.title.x = element_text(face="bold"),
      axis.title.y = element_text(face="bold"))
ggsave('Observed_patient2_point.pdf', width = 4, height = 3, dpi=600)

ggplot(tmp.3, aes(x=time, y=response)) + geom_line(size=1.2) + guides(colour=FALSE) + 
xlab("time") + ylab("response") + aes(colour = factor(patient)) + 
ggtitle('patient3') + scale_color_manual(values=colpal) +
theme(plot.title = element_text(hjust = 0.5)) + ylim(84, 131) +
theme(axis.title.x = element_text(face="bold"),
      axis.title.y = element_text(face="bold")) 
ggsave('Observed_patient3.pdf', width = 4, height = 3, dpi=600)

ggplot(tmp.3, aes(x=time, y=response)) + geom_point() + guides(colour=FALSE) + 
xlab("time") + ylab("response") + aes(colour = factor(patient)) + 
ggtitle('patient3') + scale_color_manual(values=colpal) +
theme(plot.title = element_text(hjust = 0.5)) + ylim(84, 131) +
theme(axis.title.x = element_text(face="bold"),
      axis.title.y = element_text(face="bold")) 
ggsave('Observed_patient3_point.pdf', width = 4, height = 3, dpi=600)

ggplot(tmp.5, aes(x=time, y=response)) + geom_line(size=1.2) + guides(colour=FALSE) + 
xlab("time") + ylab("response") + aes(colour = factor(patient)) + 
ggtitle('patient5') + scale_color_manual(values=colpal) +
theme(plot.title = element_text(hjust = 0.5)) + ylim(84, 131) +
theme(axis.title.x = element_text(face="bold"),
      axis.title.y = element_text(face="bold"))
ggsave('Observed_patient5.pdf', width = 4, height = 3, dpi=600)

ggplot(tmp.5, aes(x=time, y=response)) + geom_point() + guides(colour=FALSE) + 
xlab("time") + ylab("response") + aes(colour = factor(patient)) + 
ggtitle('patient5') + scale_color_manual(values=colpal) +
theme(plot.title = element_text(hjust = 0.5)) + ylim(84, 131) +
theme(axis.title.x = element_text(face="bold"),
      axis.title.y = element_text(face="bold"))
ggsave('Observed_patient5_point.pdf', width = 4, height = 3, dpi=600)

# outcomes are standardized (problematic if not standardized; why?)
source('../../sfpca.R')
prepared_data = prepare_data(data=df, unique_subject_id = 'patient', time_var='time', 
	                         response='response', standardize.y=TRUE, scale.time=TRUE)
save(prepared_data, file='prepared_data_dfm.RData')

##################### run sfpca stan model ##############
Nsamples = 1000
Nchains = 3

model_file = "../../sfpca2.stan"
smod = stan_model(model_file)


PC_max = 2 # number of PCs
D_max = 2 # number of knots
sa.list = list()
log_lik.list = list()
log_lik_con.list = list()
looic.list = list()
waic.list = list()
waic_con.list = list()
i = 0
for (k in 1:PC_max){
	print(paste('number of PC:', k, sep=''))

	for (d in 1:D_max){
	print(paste('number of knots:', d, sep=''))

	i = i + 1
	print(paste('index i is:', i, sep=''))

	results_basis = basis_setup_sparse(prepared_data=prepared_data, nknots=d, orth=TRUE)
	pca_data <- list(N = prepared_data$num_subjects, K = k, Q = d + 4, Y = prepared_data$response.list,
		             V = prepared_data$visits.vector, subject_starts = prepared_data$visits.start,
			   	     subject_stops = prepared_data$visits.stop, cov_starts = prepared_data$cov.start,
			   	     cov_stops = prepared_data$cov.stop, cov_size = prepared_data$cov.size,
			   	     B = results_basis$orth_spline_basis_sparse_stacked)
	set.seed(31)
	sa.list[[i]] = sampling(smod, data= pca_data, iter=Nsamples, chains=Nchains, init="random")
	log_lik.list[[i]] <- extract(sa.list[[i]],"log_lik_marg")[[1]]
	looic.list[[i]] = loo(log_lik.list[[i]])
	waic.list[[i]] = waic(log_lik.list[[i]])

	log_lik_con.list[[i]] <- extract(sa.list[[i]],"log_lik")[[1]]
	waic_con.list[[i]] = waic(log_lik_con.list[[i]])

	save(Nsamples, Nchains, sa.list, log_lik.list, looic.list, 
		 waic.list,log_lik_con.list, waic_con.list, file="stan_results.RData")
	print("###############################")
	}
}

load('stan_results.RData')
looic.obj = compare(looic.list[[1]],looic.list[[2]],looic.list[[3]], looic.list[[4]])
print(looic.obj) 
#                 elpd_diff se_diff elpd_loo p_loo  looic 
# looic.list[[3]]    0.0       0.0  -646.1     19.7 1292.2
# looic.list[[4]]    0.0       6.5  -646.1     24.5 1292.2
# looic.list[[1]]  -27.6      13.7  -673.7     20.6 1347.4
# looic.list[[2]]  -41.0      15.1  -687.1     49.1 1374.2


#### predictive checking ########
sa = sa.list[[3]]
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
Ynew_mean = colMeans(Ynew_transform)

pdf('predictive_checking.pdf')
library("bayesplot")
library("ggplot2")
color_scheme_set("brightblue")
ppc_dens_overlay(prepared_data$data$response, Ynew_transform)
dev.off()

#######################################
## Extract parameters from the model ##
#######################################
load('prepared_data_dfm.RData')
Nsamples = 1000
Nchains = 3

load('stan_results.RData')
i = 3; k = npcs =2; d = nknots =1; 
sa = sa.list[[i]]
source('../../sfpca.R')
results_basis = basis_setup_sparse(prepared_data=prepared_data, nknots=d, orth=TRUE)

Sigma = extract(sa,"Sigma",permuted=FALSE)
W = extract(sa,"W",permuted=FALSE)
sigma_eps = extract(sa,"sigma_eps",permuted=FALSE)
theta_mu = extract(sa,"theta_mu",permuted=FALSE)
alpha = extract(sa,"alpha",permuted=FALSE)
Theta = extract(sa,"Theta",permuted=FALSE)


## Reshape parameters and reorient loadings with PCA rotation 
N = prepared_data$num_subjects
K = npcs
Q = nknots + 4

theta_mu_new = array(0, dim=c(Q, Nchains*Nsamples/2))
alpha_old = alpha_new = array(0, dim=c(K, N, Nchains*Nsamples/2)) 
Theta_old = Theta_new = array(0, dim=c(Q, K, Nchains*Nsamples/2))
W_old = array(0, dim=c(Q, Q, Nchains*Nsamples/2)) 

ind = 0
prop_var = NULL
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
		prop_var = rbind(prop_var, d_temp/sum(d_temp)) # proportion of variance explained by each PC

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

# proportion of variance explained by each PC (average over all draws)
prop_var_avg_origin = colMeans(prop_var)
prop_var_avg = paste(round(colMeans(prop_var)*100, 2), '%', sep='')

save(N, K, Q, alpha_new, theta_mu_new, Theta_new, prop_var_avg_origin, prop_var_avg, results_basis, prepared_data, 
	 file="post_rotation_results.RData")

pdf('screeplot.pdf')
plot(prop_var_avg_origin, xlab='Principal Component Functions', 
	ylab='Proportion of Variance Explained', ylim=c(0,1),type="b", main='Scree Plot')
dev.off()


######### plotting ######
load("post_rotation_results.RData")
vars_select = c('ID', 'patient', 'time', 'response', 'Version') # selected variables from prepared_data$data

ALPHA_array = alpha_new
MU_array = theta_mu_new
THETA_array = Theta_new
phi_t_cont = results_basis$orth_spline_basis_cont
phi_t = results_basis$orth_spline_basis_sparse
time_cont = results_basis$time_cont

nloop=dim(ALPHA_array)[3]
first=1
last=nloop
# N = prepared_data$num_subjects

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


### create data frame containing needed information ####
df = prepared_data$data[, vars_select]
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

df$Y_sparse = unlist(Y_sparse) # check: sum(df$Y_sparse != df$response) == 0
df$Fits_sparse = unlist(Fits_sparse)
df$residuals = df$Y_sparse - df$Fits_sparse
df$residuals = df$Y_sparse - df$Fits_sparse

# export resulted data #
write.table(df, 'data_sfpca.txt', sep='\t', row.names=F)

## residual analysis ##
pdf('residual_plots.pdf')
par(mfrow=c(1,2))
plot(df$residuals)
library(car); qqPlot(df$residuals)
dev.off()


###### start plotting ######
pdf('mean_spaghetti.pdf')
plot(time_cont, Mu_functions, ylim=c(min(unlist(Y_sparse)), max(unlist(Y_sparse))), ylab='pH', xlab='time', lwd=5, col=4)
for(i in 1:N){
	lines(time_sparse[[i]],Y_sparse[[i]],type="l",lwd=.25)
}
title(main='pH')
dev.off()

pdf('spaghetti_smoothed.pdf')
par(mfrow=c(1,2))
plot(time_cont, Mu_functions,type="n", ylim=c(min(unlist(Y_sparse)), max(unlist(Y_sparse))),
	 lwd=5, col=4, ylab='pH', xlab='Time')
for(i in 1:N){
	lines(time_sparse[[i]], Y_sparse[[i]], type="l", lwd=.25)
}
title(main='Observed pH')

plot(time_cont, Mu_functions, type="n", ylim=c(min(unlist(Y_sparse)), max(unlist(Y_sparse))),
	 lwd=2,col=4, ylab='pH', xlab='Time')
for(i in 1:N){
		lines(time_sparse[[i]],Fits_sparse[[i]], type="l", lwd=.25)
}
title(main='Smoothed pH')
dev.off()


library(ggplot2)
colpal <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")
ggplot(df, aes(x=time*22, y=Fits_sparse, color=Version)) + geom_point(alpha=0) + # use alpha=0 
xlab("Observation Time Point") + ylab("Smoothed Response") + scale_color_manual(values=colpal)+
theme_classic() + theme(axis.text = element_text(color='black')) + 
geom_smooth(se=FALSE, size=1, aes(group=Version)) +
geom_point(data=df, aes(x=time*22, y=response, color=Version), alpha=0.5, size=1.2) 
# theme(legend.position='top') +
ggsave('Smoothed_spaghetti_groups.pdf', width = 4, height = 3, dpi=600)



########## FPC curves ###### (need to automatic ACT different number of PCs)
pdf('FPCs.pdf', width=5, height=8)
plot(time_cont*22, FPC_mean[, 1], type="l", lwd=2, ylim=c(min(unlist(Y_sparse)), max(unlist(Y_sparse))),
     xlab='time', ylab='PC Curve Values', cex.lab=1.2, col='red')
lines(time_cont*22, FPC_mean[, 2],type="l",lwd=2, col='blue')
title(main=paste("FPCs for", "response"))
legend('topright', c('PC1', 'PC2'), lwd=c(3, 3), col=c('red', 'blue'), bty='n')
dev.off()

###### try transform back
data = read.csv("dfm_bound.txt",header=TRUE, sep='\t')
sigma_y = sd(data$response)
mu_y = mean(data$response)

pdf('FPCs_transformed.pdf', width = 4, height = 4)
plot(time_cont*22, FPC_mean[, 1]*sigma_y + mu_y, type="l", lwd=2, ylim=c(84, 131),
     xlab='time', ylab='PC Curve Values', col='red')
lines(time_cont*22, FPC_mean[, 2]*sigma_y + mu_y,type="l",lwd=2, col='blue')
#title(main=paste("FPCs for", "response"))
legend('topright', c('PC1', 'PC2'), lwd=c(2, 2), col=c('red', 'blue'), bty='n')
dev.off()

pdf('FPCs_mean.pdf', width=5, height=8)
par(mfrow=c(2,1))
for (k in 1:K){
	plot(time_cont, Mu_functions, type="l", ylim=c(min(unlist(Y_sparse)), max(unlist(Y_sparse))),
	 lwd=2,col=1, xlab='Time', ylab='pH', cex.lab=1.5)
	lines(time_cont, Mu_functions + FPC_mean[,k],type="l",lwd=3,lty=2,col=2) # red
	lines(time_cont, Mu_functions - FPC_mean[,k],type="l",lwd=3,lty=2,col=3) # green
	title(main=paste("Effect of FPC", k , "for", "pH"))
}
legend('topright', c('+ pc', '- pc'), lty=c(2,2), lwd=c(3,3), col=c(2, 3), bty='n', cex=0.5)
dev.off()

## plot each PC separately with proportion of variation explained
for (k in 1:K){
	pdf(paste(paste('FPCs_mean_PC', k, sep=''), 'pdf', sep='.'), width = 4, height = 4)
	plot(time_cont*22, Mu_functions*sigma_y + mu_y, type="l", ylim=c(84, 121),
	 lwd=2,col=1, xlab='time', ylab='response function', font.lab=2, cex.lab=1.2)
	lines(time_cont*22, (Mu_functions + FPC_mean[,k])*sigma_y + mu_y,type="l",lwd=3,lty=2,col=2) # red
	lines(time_cont*22, (Mu_functions - FPC_mean[,k])*sigma_y + mu_y,type="l",lwd=3,lty=2,col=3) # green
	title(main=paste(paste('PC', k, sep=' '), ' (', prop_var_avg[k], ' )', sep=''))
	#axis(1, font=2) # make x-axis ticks label bold
	legend('topright', c('+ pc', '- pc'), lty=c(2,2), lwd=c(3,3), col=c(2, 3), bty='n', cex=0.5)
	dev.off()
}

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
  	   	 lwd=2, col=4, xlab="Time", ylab="pH")
  	lines(times_i, Y_i, type="l", col=4, lty=2)
  	lines(times_i, Fitted_i, type="l", col=4)
  	title(main=paste("Fitted", "pH", "for subj", pid))
}
dev.off()

########## analysis with principle component scores ###### 

library(ggpubr)
library(gridExtra)
pdf('fpcScores_ggBoxplot.pdf', width = 7, height = 3.8)
### p-values based on T test (same as from simple linear regression)
compare_means(fpc1 ~ Version,  data = df, method = "t.test")
my_comparisons <- list(c("baseline", "shift_1x"), c("baseline", "shift_2x"),
	                   c("baseline", "shift_4x"))
colpal <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")
plot1 = ggboxplot(df, x = "Version", y = "fpc1",
          color = "Version", palette = colpal)+ 
stat_compare_means(comparisons = my_comparisons, method = "t.test") + # Add pairwise comparisons p-value
theme(text = element_text(size=12), 
	  axis.title.x = element_text(size=18, face="bold"),
      axis.title.y = element_text(size=18, face="bold"),
      axis.text.x = element_text(face="bold", size=10),
	  axis.text.y = element_text(face="bold", size=10)) +
theme(legend.position="none") # remove legend

compare_means(fpc2 ~ Version,  data = df, method = "t.test")
my_comparisons <- list(c("baseline", "shift_1x"), c("baseline", "shift_2x"),
	                   c("baseline", "shift_4x"))
colpal <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")
plot2 = ggboxplot(df, x = "Version", y = "fpc2",
          color = "Version", palette = colpal)+ 
stat_compare_means(comparisons = my_comparisons, method = "t.test") +
theme(text = element_text(size=12), 
	  axis.title.x = element_text(size=18, face="bold"),
      axis.title.y = element_text(size=18, face="bold"),
      axis.text.x = element_text(face="bold", size=10),
	  axis.text.y = element_text(face="bold", size=10)) +
theme(legend.position="none") # remove legend
grid.arrange(plot1, plot2, ncol=2)
dev.off()

# scatterplot of PC scores
ggplot(df, aes(x=fpc1 , y=fpc2)) + geom_point() + 
aes(colour=factor(Version)) + theme(legend.position="top") +
scale_color_manual(values=colpal) +
theme(text = element_text(size=12), 
	  axis.title.x = element_text(size=18, face="bold"),
      axis.title.y = element_text(size=18, face="bold"))
ggsave('scatterplot_scores.pdf', width=6, height=5)


### each PC scores separately
library(ggpubr)
compare_means(fpc1 ~ Version,  data = df, method = "t.test")
my_comparisons <- list(c("baseline", "shift_1x"), c("baseline", "shift_2x"),
	                   c("baseline", "shift_4x"))
colpal <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")
plot1 = ggboxplot(df, x = "Version", y = "fpc1",
          color = "Version", palette = colpal)+ 
stat_compare_means(comparisons = my_comparisons, method = "t.test") + # Add pairwise comparisons p-value
theme(text = element_text(size=12), 
	  axis.title.x = element_text(size=18, face="bold"),
      axis.title.y = element_text(size=18, face="bold"),
      axis.text.x = element_text(face="bold", size=10),
	  axis.text.y = element_text(face="bold", size=10)) +
theme(legend.position="right") + ylab('PC 1 scores') 
ggsave('boxplot_PC1_scores.pdf', width=4, height=3.2)

compare_means(fpc2 ~ Version,  data = df, method = "t.test")
my_comparisons <- list(c("baseline", "shift_1x"), c("baseline", "shift_2x"),
	                   c("baseline", "shift_4x"))
colpal <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")
plot2 = ggboxplot(df, x = "Version", y = "fpc2",
          color = "Version", palette = colpal)+ 
stat_compare_means(comparisons = my_comparisons, method = "t.test") +
theme(text = element_text(size=12), 
	  axis.title.x = element_text(size=18, face="bold"),
      axis.title.y = element_text(size=18, face="bold"),
      axis.text.x = element_text(face="bold", size=10),
	  axis.text.y = element_text(face="bold", size=10)) +
theme(legend.position="right") + ylab('PC 2 scores') # remove legend
ggsave('boxplot_PC2_scores.pdf', width=4, height=3.2)
