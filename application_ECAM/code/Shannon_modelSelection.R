library(parallel)
library(rstan)
library(loo) # for waic()
library(Matrix) # for bdiag(): construct a block diagnoal matrix
options(mc.cores = parallel::detectCores())


############################################################
############ Data Preparation ##############################
############################################################
dat = read.csv('../data/metadata_shannon.txt', sep='\t')
dat = dat[, c('studyid', 'month_of_life', 'shannon', 'antiexposedall', 'delivery', 'diet', 'sex')]
colnames(dat) = c('ID_unique', 'Time', 'response', 'abx', 'delivery', 'diet','sex')

length(unique(dat$ID_unique)) # 43 babies
summary(dat)

############################################################
############ Check data transformation #####################
############################################################
par(mfrow=c(2,2))
library(car)
qqPlot(dat$response) 
qqPlot(log(dat$response))
qqPlot(sqrt(dat$response)) 
qqPlot((dat$response)^(1/3)) 
## conclusion: original; no transformation

## prepare data for SFPCA
source('../../sfpca.R')
prepared_data = prepare_data(data=dat, unique_subject_id = 'ID_unique', time_var='Time', 
	                         response='response', standardize.y=TRUE, scale.time=FALSE)
save(prepared_data, file='prepared_data_shannon.RData')

############################################################
############ Plot observed trajectories ####################
############################################################
library(ggplot2)
ggplot(dat, aes(x=Time, y=response)) + geom_line(alpha=0.5) + guides(colour=FALSE) + 
xlab("Month of life") + ylab("Shannon diversity") + aes(colour = factor(ID_unique)) +
theme_classic() + theme(axis.text = element_text(color='black')) +
theme(axis.title.x = element_text(face="bold"),
      axis.title.y = element_text(face="bold")) +
scale_fill_distiller(palette = "Spectral") + scale_x_continuous(breaks=seq(0,29,5))
ggsave('Observed_spaghetti.pdf', width = 4, height = 3, dpi=600)


abx_col = c('#2c7fb8', '#dd1c77')
ggplot(dat, aes(x=Time, y=response, color=abx)) + geom_point(alpha=0.2) + 
xlab("Month of life") + ylab("Shannon diversity")  + scale_x_continuous(breaks=seq(0,29,5)) +
geom_smooth(se=FALSE, size=1.5, aes(group=abx)) +
theme_classic() + theme(axis.text = element_text(color='black')) +
theme(axis.title.x = element_text(face="bold"),
      axis.title.y = element_text(face="bold")) +
labs(color='Antibiotics') + scale_color_manual(labels = c("No", "Yes"), values=abx_col)  
ggsave('Observed_shannon_abx.pdf', width=4, height=3)

dm_col = c('#636363', '#8856a7')
ggplot(dat, aes(x=Time, y=response, color=delivery)) + geom_point(alpha=0.2) + 
xlab("Month of life") + ylab("Shannon diversity")  + scale_x_continuous(breaks=seq(0,29,5)) +
geom_smooth(se=FALSE, size=1.5, aes(group=delivery)) +
theme_classic() + theme(axis.text = element_text(color='black')) +
theme(axis.title.x = element_text(face="bold"),
      axis.title.y = element_text(face="bold")) +
labs(color='Delivery Mode') + scale_color_manual(values=dm_col)
ggsave('Observed_shannon_delivery.pdf', width=4, height=3)

diet_col = c('#006837', '#e6550d')
ggplot(dat, aes(x=Time, y=response, color=diet)) + geom_point(alpha=0.2) + 
xlab("Month of life") + ylab("Shannon diversity")  + scale_x_continuous(breaks=seq(0,29,5)) +
geom_smooth(se=FALSE, size=1.5, aes(group=diet)) +
theme_classic() + theme(axis.text = element_text(color='black')) +
theme(axis.title.x = element_text(face="bold"),
      axis.title.y = element_text(face="bold")) +
labs(color='Diet') + scale_color_manual(labels=c('breastfed', 'formula'), values=diet_col)
ggsave('Observed_shannon_diet.pdf', width=4, height=3)


############################################################
############ Results from SFPCA ############################
############################################################
source('../../sfpca.R')
load('prepared_data_shannon.RData')


Nsamples = 1000
Nchains = 3
model_file = "../../sfpca2.stan"
smod = stan_model(model_file)

PC_max = 5 # number of PCs
D_max = 5 # number of knots
sa.list = list()
log_lik.list = list()
log_lik_con.list = list()
looic.list = list()
waic.list = list()
waic_con.list = list()
i = 0
for (k in 1:PC_max){
	for (d in 1:D_max){
	i = i + 1
	print(paste('index i is:', i, 'number of PC:', k, 'number of knots:', d))

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
		 waic.list,log_lik_con.list, waic_con.list, file="stan_modelSelecton_shannon.RData")
	print("###############################")
	}
}

#################################
## Model selection using LOOIC ##
#################################
load("stan_modelSelecton_shannon.RData")
looic.obj = compare(looic.list[[1]],looic.list[[2]],looic.list[[3]],
				    looic.list[[4]],looic.list[[5]],looic.list[[6]],
					looic.list[[7]],looic.list[[8]],looic.list[[9]],
					looic.list[[10]],looic.list[[11]],looic.list[[12]], looic.list[[13]])

print(looic.obj) 
#                  elpd_diff se_diff elpd_loo p_loo   looic  
# looic.list[[11]]     0.0       0.0 -1051.7    132.5  2103.4  3 PC 1 knots
# looic.list[[7]]    -18.4      17.7 -1070.1    113.5  2140.2  2 PC 1 knots
# looic.list[[12]]   -24.6      16.0 -1076.3    157.1  2152.7
# looic.list[[8]]    -50.1      18.3 -1101.8    144.8  2203.7
# looic.list[[9]]    -56.0      18.8 -1107.7    158.4  2215.4
# looic.list[[6]]    -58.0      16.3 -1109.7    148.8  2219.4
# looic.list[[13]]   -58.7      16.2 -1110.4    191.1  2220.7
# looic.list[[10]]   -63.1      19.2 -1114.8    159.8  2229.6
# looic.list[[2]]   -104.6      24.2 -1156.3    145.9  2312.5
# looic.list[[5]]   -134.3      27.7 -1186.0    163.1  2372.1
# looic.list[[4]]   -135.9      28.2 -1187.6    164.9  2375.2
# looic.list[[3]]   -139.1      29.2 -1190.8    161.2  2381.6
# looic.list[[1]]   -339.7      47.0 -1391.4    376.1  2782.7

#################################
## predictive checking ##
#################################
sa = sa.list[[11]]
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
ppc_dens_overlay(prepared_data$data$response, Ynew_transform) + ggtitle("shannon_3pc_1knots")


