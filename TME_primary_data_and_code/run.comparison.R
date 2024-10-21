
library("tidyverse")

#feature importance comparison
pal_imp <- read_tsv("../human_pal-data_20230907/signature.RF.feat_importance.txt")
wu_imp <- read_tsv("../human_wu-data_20230911/signature.age_cutoff=50.RF.feat_importance.txt")

#
pal_imp_prep <- pal_imp %>% 
  select(Variable, Importance) %>% 
  dplyr::rename(Importance_Pal = Importance)
wu_imp_prep <- wu_imp %>% 
  dplyr::rename(Importance_Wu = Importance)

#
merge <- wu_imp_prep %>% 
  inner_join(pal_imp_prep, by = "Variable")
#replace 0s
Importance_Pal_min <- min(abs(merge$Importance_Pal[merge$Importance_Pal != 0]))
merge$Importance_Pal[merge$Importance_Pal == 0] <- Importance_Pal_min
Importance_Wu_min <- min(abs(merge$Importance_Wu[merge$Importance_Wu != 0]))
merge$Importance_Wu[merge$Importance_Wu == 0] <- Importance_Wu_min

#label on the top in both
merge <- merge %>% 
  mutate(label_top_both = ifelse(Importance_Wu >= 0.001 & Importance_Pal >= 0.001, Variable, ""))

#save data
write_tsv(x = merge, file = "merge.RF_importance.txt")

#SCC in the title
#linear and log10 plot
#highlight the top features

cor_res <- cor.test(merge$Importance_Pal, merge$Importance_Wu, method = "spearman")
tit <- paste0("SCC = ", formatC(cor_res$estimate), "; p = ", formatC(cor_res$p.value))

scatter_plt <- merge %>% ggplot(aes(x=Importance_Wu, y=Importance_Pal, col=label_top_both != "")) + 
  geom_point() + 
  geom_text(aes(label=label_top_both), size=3, vjust = -1, col="red") +
  scale_colour_manual(values = c("black", "red")) +
  theme_bw() +
  NoLegend() +
  ggtitle(tit)

#save plot
pdf("merge.RF_importance.scatter.pdf", width = 16, height = 8)
plot(scatter_plt | scatter_plt + scale_x_log10() + scale_y_log10())
dev.off()

#Fisher test at increasingly stringent threshold on the feature importance
#selected by Wu [y, n] vs selected by Pal [y, n]

ft_res <- c()

ts <- seq(0.0001, 0.01, 0.0001)

for (thresh in ts) {
  
  both <- sum(merge$Importance_Wu >= thresh & merge$Importance_Pal >= thresh)
  only_wu <- sum(merge$Importance_Wu >= thresh & merge$Importance_Pal < thresh)
  only_pal <- sum(merge$Importance_Wu < thresh & merge$Importance_Pal >= thresh)
  none <- sum(merge$Importance_Wu < thresh & merge$Importance_Pal < thresh)
  
  if (both > 0) {
  
    ft <- fisher.test(matrix(c(both, only_wu, only_pal, none), ncol = 2))
  
    ft_res <- rbind(ft_res, c(thresh, both, as.numeric(ft$estimate), ft$p.value))
    
  }

}

ft_res <- tibble(thrshold_feat_imp = ft_res[,1], 
                 feat_common_n = ft_res[,2], 
                 odds_ratio = ft_res[,3],
                 p_value = ft_res[,4])

#save results
write_tsv(x = ft_res, file = "merge.RF_importance.fisher_tests.txt")

#plot results
plt_odds <- ft_res %>% ggplot(aes(x = thrshold_feat_imp, y = odds_ratio)) + geom_point() + theme_bw() + geom_vline(xintercept=0.001)
plt_feats <- ft_res %>% ggplot(aes(x = thrshold_feat_imp, y = feat_common_n)) + geom_point() + theme_bw() + scale_y_log10() + geom_vline(xintercept=0.001)
plt_pval <- ft_res %>% ggplot(aes(x = thrshold_feat_imp, y = -log10(p_value))) + geom_point() + theme_bw() + geom_vline(xintercept=0.001)

#save plots
pdf("merge.RF_importance.fisher_tests.stats.pdf", width = 3, height = 7)
plot(plt_odds / plt_feats / plt_pval)
dev.off()

