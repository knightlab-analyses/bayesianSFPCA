library(parallel)
library(rstan)
library(loo) # for waic()
library(Matrix) # for bdiag(): construct a block diagnoal matrix
library(ggplot2)
options(mc.cores = parallel::detectCores())


df=read.csv("df_bound.txt",header=TRUE, sep='\t')

###### plottig the data ############
# compare different versions
colpal <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")
ggplot(df, aes(x=time, y=response, group=patient, colour=Version)) + 
geom_line(alpha=0.3) + xlab("Observed Time") + ylab("response")  +
theme_classic() + theme(axis.text = element_text(color='black')) +
geom_smooth(se=FALSE, size=1.5, aes(group=Version)) +
scale_color_manual(values=colpal)
ggsave('Observed_spaghetti_group.pdf', width = 4, height = 3, dpi=600)


# selected patients
tmp.4 = df[df$patient %in% c('4', '4_1x', '4_2x', '4_4x'), ]
ggplot(tmp.4, aes(x=time, y=response)) + geom_line(size=1.2) + guides(colour=FALSE) + 
xlab("time") + ylab("response") + aes(colour = factor(patient)) + 
ggtitle('patient4') + scale_color_manual(values=colpal) +
theme(plot.title = element_text(hjust = 0.5)) + ylim(84, 121) +
theme(axis.title.x = element_text(face="bold"),
      axis.title.y = element_text(face="bold"))
ggsave('Observed_patient4.pdf', width = 4, height = 3, dpi=600)

tmp.8 = df[df$patient %in% c('8', '8_1x', '8_2x', '8_4x'), ]
ggplot(tmp.8, aes(x=time, y=response)) + geom_line(size=1.2) + guides(colour=FALSE) + 
xlab("time") + ylab("response") + aes(colour = factor(patient)) + 
ggtitle('patient8') + scale_color_manual(values=colpal) +
theme(plot.title = element_text(hjust = 0.5)) + ylim(84, 121) +
theme(axis.title.x = element_text(face="bold"),
      axis.title.y = element_text(face="bold"))
ggsave('Observed_patient8.pdf', width = 4, height = 3, dpi=600)


# apply bayesian SFPCA
source('../../sfpca.R')
prepared_data = prepare_data(data=df, unique_subject_id = 'patient', time_var='time', 
	                         response='response', standardize.y=TRUE, scale.time=TRUE)
save(prepared_data, file='prepared_data_df.RData')

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
# looic.list[[4]]    0.0       0.0  -617.8     12.0 1235.7
# looic.list[[3]]  -19.5       9.0  -637.3      9.5 1274.7
# looic.list[[2]]  -21.4       7.7  -639.3     11.2 1278.5
# looic.list[[1]]  -31.8      11.0  -649.6      8.3 1299.3

looic.list[[4]] ## all good

# Computed from 1500 by 40 log-likelihood matrix

#          Estimate   SE
# elpd_loo   -617.8 17.9
# p_loo        12.0  1.4
# looic      1235.7 35.9
# ------
# Monte Carlo SE of elpd_loo is 0.1.

# All Pareto k estimates are good (k < 0.5).
# See help('pareto-k-diagnostic') for details.

#### posterior predictive checking ########
load('prepared_data_df.RData')
sa = sa.list[[4]]
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

library("bayesplot")
library("ggplot2")
color_scheme_set("brightblue")
ppc_dens_overlay(prepared_data$data$response, Ynew_transform)


#######################################
## Extract parameters from the model ##
#######################################
#load('prepared_data_df.RData')
Nsamples = 1000
Nchains = 3


#load('stan_results.RData')
i = 4; k = npcs =2; d = nknots =2; 
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

# screeplot
plot(prop_var_avg_origin, xlab='Principal Component Functions', 
	ylab='Proportion of Variance Explained', ylim=c(0,1),type="b", main='Scree Plot')


######### plotting ######
#load("post_rotation_results.RData")
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


###### plot results on original data sclae 
data = read.csv("df_bound.txt",header=TRUE, sep='\t')
sigma_y = sd(data$response)
mu_y = mean(data$response)

### plot each PC separately with proportion of variation explained
for (k in 1:K){
	#pdf(paste(paste('FPCs_mean_PC', k, sep=''), 'pdf', sep='.'), width = 4, height = 4)
	plot(time_cont*22, Mu_functions*sigma_y + mu_y, type="l", ylim=c(84, 121),
	 lwd=2,col=1, xlab='time', ylab='response function', font.lab=2, cex.lab=1.2)
	lines(time_cont*22, (Mu_functions + FPC_mean[,k])*sigma_y + mu_y,type="l",lwd=3,lty=2,col=2) # red
	lines(time_cont*22, (Mu_functions - FPC_mean[,k])*sigma_y + mu_y,type="l",lwd=3,lty=2,col=3) # green
	title(main=paste(paste('PC', k, sep=' '), ' (', prop_var_avg[k], ' )', sep=''))
	legend('topright', c('+ pc', '- pc'), lty=c(2,2), lwd=c(3,3), col=c(2, 3), bty='n', cex=0.5)
	#dev.off()
}


########## analysis with principle component scores ###### 
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
theme(legend.position="right") + ylab('PC 2 scores') 
ggsave('boxplot_PC2_scores.pdf', width=4, height=3.2)
