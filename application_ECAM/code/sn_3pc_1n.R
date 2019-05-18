library(parallel)
library(rstan)
library(loo) # for waic()
library(Matrix) # for bdiag(): construct a block diagnoal matrix
options(mc.cores = parallel::detectCores())

#################### looking at results ###########
load('prepared_data_shannon.RData')
load("stan_modelSelecton_shannon.RData")
i = 11; k = npcs =3; d = nknots =1; 
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
	 file="../results_3pc_1n/post_rotation_results.RData")


######### plotting ######
load("../results_3pc_1n/post_rotation_results.RData")
vars_select = c('ID', 'time', 'response', 'abx', 'delivery', 'diet', 'sex') # selected variables from prepared_data$data

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
df$fpc3=0

i = 0
for (pid in unique(df$ID)){
	i = i + 1
	Y_sparse[[i]] = df$response[df$ID == pid]
	time_sparse[[i]] = df$time[df$ID == pid]
	df$fpc1[df$ID == pid] = scores[i, 1]
	df$fpc2[df$ID == pid] = scores[i, 2]
	df$fpc3[df$ID == pid] = scores[i, 3]
}

Fits_sparse=list()
for(i in 1:N){
	Fits_sparse[[i]] = t(phi_t[[i]]) %*% MU_mean + t(phi_t[[i]]) %*% THETA_mean %*% ALPHA_mean[, i]
}

df$Y_sparse = unlist(Y_sparse) 
df$Fits_sparse = unlist(Fits_sparse)
df$residuals = df$Y_sparse - df$Fits_sparse

### plot estimated mean curve on original scale
dat = read.csv('../data/metadata_shannon.txt', sep='\t')
dat = dat[, c('studyid', 'month_of_life', 'shannon', 'antiexposedall', 'delivery', 'diet', 'sex')]
colnames(dat) = c('ID_unique', 'Time', 'response', 'abx', 'delivery', 'diet','sex')
sigma_y = sd(dat$response)
mu_y = mean(dat$response)

pdf('../results_3pc_1n/mean_spaghetti_transform.pdf', width = 4, height = 4)
plot(time_cont, Mu_functions*sigma_y + mu_y, type="l", ylim=c(-0.1, 6.5), 
	 xlab='Month of life', ylab='Shannon diversity', lwd=5, col=4, font.lab=2, cex.lab=1.2)
for(i in 1:N){
	lines(time_sparse[[i]],Y_sparse[[i]]*sigma_y + mu_y,type="l",lwd=.25)
}
#title(main='Shannon')
dev.off()

### plot each PC separately with proportion of variation explained
for (k in 1:K){
	pdf(paste(paste('../results_3pc_1n/FPCs_mean_PC', k, sep=''), 'pdf', sep='.'), width = 4, height = 4)
	plot(time_cont, Mu_functions*sigma_y + mu_y, type="l", ylim=c(-0.1, 6.5),
	 lwd=2,col=1, xlab='Month of life', ylab='Shannon diversity', font.lab=2, cex.lab=1.2)
	lines(time_cont, (Mu_functions + FPC_mean[,k])*sigma_y + mu_y,type="l",lwd=3,lty=2,col=2) # red
	lines(time_cont, (Mu_functions - FPC_mean[,k])*sigma_y + mu_y,type="l",lwd=3,lty=2,col=3) # green
	title(main=paste(paste('PC', k, sep=' '), ' (', prop_var_avg[k], ' )', sep=''))
	#axis(1, font=2) # make x-axis ticks label bold
	legend('topright', c('+ pc', '- pc'), lty=c(2,2), lwd=c(3,3), col=c(2, 3), bty='n', cex=0.5)
	dev.off()
}

######### plot boxplots of PC scores ########
library(ggpubr)
####### antibiotics #########
abx_col = c('#2c7fb8', '#dd1c77')
ggboxplot(df, x = "abx", y = "fpc1", color = "abx") + 
stat_compare_means(method = "t.test")  +
theme(axis.title.x = element_text(face="bold"),
	  axis.title.y = element_text(face="bold")) +
theme(legend.position="right") + ylab('PC 1 scores') + xlab('Antibiotics') +
scale_x_discrete(breaks=c("n", "y"), labels=c("No", "Yes")) +
labs(color='Antibiotics') + scale_color_manual(labels = c("No", "Yes"), values=abx_col) 
ggsave('../results_3pc_1n/boxplot_abx_PC1_scores.pdf', width=4, height=3.2)

ggboxplot(df, x = "abx", y = "fpc2", color = "abx") + 
stat_compare_means(method = "t.test")  +
theme(axis.title.x = element_text(face="bold"),
	  axis.title.y = element_text(face="bold")) +
theme(legend.position="right") + ylab('PC 2 scores') + xlab('Antibiotics') +
scale_x_discrete(breaks=c("n", "y"), labels=c("No", "Yes")) +
labs(color='Antibiotics') + scale_color_manual(labels = c("No", "Yes"), values=abx_col) 
ggsave('../results_3pc_1n/boxplot_abx_PC2_scores.pdf', width=4, height=3.2)

ggboxplot(df, x = "abx", y = "fpc3", color = "abx") + 
stat_compare_means(method = "t.test")  +
theme(axis.title.x = element_text(face="bold"),
	  axis.title.y = element_text(face="bold")) +
theme(legend.position="right") + ylab('PC 3 scores') + xlab('Antibiotics') +
scale_x_discrete(breaks=c("n", "y"), labels=c("No", "Yes")) +
labs(color='Antibiotics') + scale_color_manual(labels = c("No", "Yes"), values=abx_col) 
ggsave('../results_3pc_1n/boxplot_abx_PC3_scores.pdf', width=4, height=3.2)

####### delivery #########
dm_col = c('#636363', '#8856a7')
ggboxplot(df, x = "delivery", y = "fpc1", color = "delivery") + 
stat_compare_means(method = "t.test")  +
theme(axis.title.x = element_text(face="bold"),
	  axis.title.y = element_text(face="bold")) +
theme(legend.position="right") + ylab('PC 1 scores') + xlab('delivery') +
labs(color='Delivery Mode') + scale_color_manual(values=dm_col)
ggsave('../results_3pc_1n/boxplot_delivery_PC1_scores.pdf', width=4, height=3.2)

ggboxplot(df, x = "delivery", y = "fpc2", color = "delivery") + 
stat_compare_means(method = "t.test")  +
theme(axis.title.x = element_text(face="bold"),
	  axis.title.y = element_text(face="bold")) +
theme(legend.position="right") + ylab('PC 2 scores') + xlab('delivery') +
labs(color='Delivery Mode') + scale_color_manual(values=dm_col)
ggsave('../results_3pc_1n/boxplot_delivery_PC2_scores.pdf', width=4, height=3.2)

ggboxplot(df, x = "delivery", y = "fpc3", color = "delivery") + 
stat_compare_means(method = "t.test")  +
theme(axis.title.x = element_text(face="bold"),
	  axis.title.y = element_text(face="bold")) +
theme(legend.position="right") + ylab('PC 3 scores') + xlab('delivery') +
labs(color='Delivery Mode') + scale_color_manual(values=dm_col)
ggsave('../results_3pc_1n/boxplot_delivery_PC3_scores.pdf', width=4, height=3.2)

####### diet #########
diet_col = c('#006837', '#e6550d')
ggboxplot(df, x = "diet", y = "fpc1", color = "diet") + 
stat_compare_means(method = "t.test")  +
theme(axis.title.x = element_text(face="bold"),
	  axis.title.y = element_text(face="bold")) +
theme(legend.position="right") + ylab('PC 1 scores') + xlab('diet') +
scale_x_discrete(breaks=c("bd", "fd"), labels=c('breastfed', 'formula')) +
labs(color='diet') + scale_color_manual(labels=c('breastfed', 'formula'), values=diet_col) 
ggsave('../results_3pc_1n/boxplot_diet_PC1_scores.pdf', width=4, height=3.2)

ggboxplot(df, x = "diet", y = "fpc2", color = "diet") + 
stat_compare_means(method = "t.test")  +
theme(axis.title.x = element_text(face="bold"),
	  axis.title.y = element_text(face="bold")) +
theme(legend.position="right") + ylab('PC 2 scores') + xlab('diet') +
scale_x_discrete(breaks=c("bd", "fd"), labels=c('breastfed', 'formula')) +
labs(color='diet') + scale_color_manual(labels=c('breastfed', 'formula'), values=diet_col) 
ggsave('../results_3pc_1n/boxplot_diet_PC2_scores.pdf', width=4, height=3.2)

ggboxplot(df, x = "diet", y = "fpc3", color = "diet") + 
stat_compare_means(method = "t.test")  +
theme(axis.title.x = element_text(face="bold"),
	  axis.title.y = element_text(face="bold")) +
theme(legend.position="right") + ylab('PC 3 scores') + xlab('diet') +
scale_x_discrete(breaks=c("bd", "fd"), labels=c('breastfed', 'formula')) +
labs(color='diet') + scale_color_manual(labels=c('breastfed', 'formula'), values=diet_col) 
ggsave('../results_3pc_1n/boxplot_diet_PC3_scores.pdf', width=4, height=3.2)

