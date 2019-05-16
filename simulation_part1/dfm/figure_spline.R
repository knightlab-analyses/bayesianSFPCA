library(splinectomeR) 
require(ggplot2)
require(reshape2)

df=read.csv("dfm_bound.txt",header=TRUE, sep='\t')
result1 = permuspliner(data=df, xvar='time', yvar='response', cases='patient',
                    category='Version', groups = c('baseline','shift_1x'), perms = 999, retain_perm = T)
result2 = permuspliner(data=df, xvar='time', yvar='response', cases='patient',
                    category='Version', groups = c('baseline','shift_2x'), perms = 999, retain_perm = T)
result4 = permuspliner(data=df, xvar='time', yvar='response', cases='patient',
                    category='Version', groups = c('baseline','shift_4x'), perms = 999, retain_perm = T)

est_bs = result4['v1_interpolated'][[1]] # baseline (x = time, var1 = estimated baseline values)
est_bs$group = as.character(result4['category_1'][[1]])
est_s4 = result4['v2_interpolated'][[1]] # shift_x (x = time, var2 = estimated shift_x values)
est_s4$group = as.character(result4['category_2'][[1]])
est_s1 = result1['v2_interpolated'][[1]] 
est_s1$group = as.character(result1['category_2'][[1]])
est_s2 = result2['v2_interpolated'][[1]] 
est_s2$group = as.character(result2['category_2'][[1]])
colnames(est_bs)[1] = colnames(est_s4)[1] = colnames(est_s1)[1] = colnames(est_s2)[1] = 'time'
colnames(est_bs)[2] = colnames(est_s4)[2] = colnames(est_s1)[2] = colnames(est_s2)[2]= 'response' 
est_data = rbind(est_bs, est_s1, est_s2, est_s4)

colpal <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")
#colpal = c('yellow', 'blue', 'green', 'purple')
ggplot() + geom_line(data=est_data, aes(x=time, y=response, color=factor(est_data$group)), size=1.2) +
scale_color_manual(values=colpal) + labs(colour= "Version") + 
theme_classic() + theme(axis.text = element_text(color='black')) + 
geom_point(data=df, aes(x=time, y=response, color=Version), alpha=0.5, size=1.2)
ggsave('splinectome_dfm_curves.pdf', width = 4, height = 3, dpi=600)




