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

library(splinectomeR)
result = permuspliner(data=df, xvar='time', yvar='response', cases='id',
	                  category='group', perms = 999, retain_perm = T, quiet = T)
result$pval # 0.001

est_g2 = result['v1_interpolated'][[1]] 
est_g2$group = as.character(result['category_1'][[1]])
est_g1 = result['v2_interpolated'][[1]] 
est_g1$group = as.character(result['category_2'][[1]])
colnames(est_g1)[1] = colnames(est_g2)[1] = 'time'
colnames(est_g1)[2] = colnames(est_g2)[2] = 'response' 
est_data = rbind(est_g1, est_g2)

colpal = c('red', 'blue')
ggplot() + geom_line(data=est_data, aes(x=time, y=response, color=factor(est_data$group)), size=1.2) +
scale_color_manual(values=colpal) + labs(colour= "Version") +
theme_classic() + theme(axis.text = element_text(color='black')) + 
geom_point(data=df, aes(x=time, y=response, color=group), alpha=0.5, size=1.2) +
theme(axis.title.x = element_text(face="bold"),
      axis.title.y = element_text(face="bold")) 
ggsave('splinectome_group_curves.pdf', width = 4, height = 3, dpi=600)





