## plot estimated mean curve from splinectomeR
library(splinectomeR)
library(ggplot2)


load('simulated_data_N25_M80_wide.rData')

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





