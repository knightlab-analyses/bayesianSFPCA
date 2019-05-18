############################################################
############ Data Preparation ##############################
############################################################
dat = read.csv('../data/metadata_shannon.txt', sep='\t')
dat = dat[, c('studyid', 'month_of_life', 'shannon', 'antiexposedall', 'delivery', 'diet', 'sex')]
colnames(dat) = c('ID_unique', 'Time', 'response', 'abx', 'delivery', 'diet','sex')


############################################################
############ Plotting  #####################################
############################################################
library(splinectomeR)
library(ggplot2)

##### antibiotic exposure #######
result_abx= permuspliner(data=dat, xvar='Time', yvar='response', cases='ID_unique',
	                  category='abx', perms = 999, retain_perm = T)
### splinectomeR cannot run, since abx status is time changing


##### delivery mode #######
result_delivery = permuspliner(data=dat, xvar='Time', yvar='response', cases='ID_unique',
	                  category='delivery', perms = 999, retain_perm = T)
# p-value = 0.021

dm_g2 = result_delivery['v1_interpolated'][[1]] 
dm_g2$group = as.character(result_delivery['category_1'][[1]])
dm_g1 = result_delivery['v2_interpolated'][[1]] 
dm_g1$group = as.character(result_delivery['category_2'][[1]])
colnames(dm_g1)[1] = colnames(dm_g2)[1] = 'time'
colnames(dm_g1)[2] = colnames(dm_g2)[2] = 'response' 
dm_data = rbind(dm_g1, dm_g2)

dm_col = c('#636363', '#8856a7')
ggplot() + geom_line(data=dm_data, aes(x=time, y=response, color=factor(dm_data$group)), size=1.2) +
scale_color_manual(values=dm_col) + labs(colour= "delivery") + 
xlab("Month of life") + ylab("Shannon diversity") + ylim(-0.1, 7) + scale_x_continuous(breaks=seq(0,29,5)) +
theme_classic() + theme(axis.text = element_text(color='black')) + 
geom_point(data=dat, aes(x=Time, y=response, color=delivery), alpha=0.2) +
theme(axis.title.x = element_text(face="bold"),
      axis.title.y = element_text(face="bold")) +
labs(color='Delivery Mode') + scale_color_manual(values=dm_col)
ggsave('splinectome_delivery_curves.pdf', width = 4, height = 3, dpi=600)


######### diet ########
result_diet = permuspliner(data=dat, xvar='Time', yvar='response', cases='ID_unique',
	                  category='diet', perms = 999, retain_perm = T)
# p-value = 0.112

diet_g2 = result_diet['v1_interpolated'][[1]] 
diet_g2$group = as.character(result_diet['category_1'][[1]])
diet_g1 = result_diet['v2_interpolated'][[1]] 
diet_g1$group = as.character(result_diet['category_2'][[1]])
colnames(diet_g1)[1] = colnames(diet_g2)[1] = 'time'
colnames(diet_g1)[2] = colnames(diet_g2)[2] = 'response' 
diet_data = rbind(diet_g1, diet_g2)

diet_col = c('#006837', '#e6550d')
ggplot() + geom_line(data=diet_data, aes(x=time, y=response, color=factor(diet_data$group)), size=1.2) +
scale_color_manual(values=diet_col) + labs(colour= "diet") + scale_x_continuous(breaks=seq(0,29,5)) +
xlab("Month of life") + ylab("Shannon diversity") + ylim(-0.1, 7) +
theme_classic() + theme(axis.text = element_text(color='black')) + 
geom_point(data=dat, aes(x=Time, y=response, color=diet), alpha=0.2) +
theme(axis.title.x = element_text(face="bold"),
      axis.title.y = element_text(face="bold")) +
labs(color='Diet') + scale_color_manual(labels=c('breastfed', 'formula'), values=diet_col)
ggsave('splinectome_diet_curves.pdf', width = 4, height = 3, dpi=600)

