library(readxl)
df <- read_xlsx('X://volume2//noam//gfp_1664//results2.xlsx')

library(ggpubr)

symnum.args <- list(
  cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
  symbols = c("xxxx", "***", "**", "*", "ns"))

my_comparisons <- list( c("pUC57-A1664G", "pUC57-WT"), c("pUC57-A1664G", "no plasmid"), c("pUC57-WT", "no plasmid") )

ggdotplot(df, x = "sample", y = "fluorescence", fill='sample', palette='jco', ylim = c(0, 12.5), size=4) +
  stat_compare_means(comparisons = my_comparisons, method='t.test',symnum.args = symnum.args, label.y=c(9.3,10.5,11.7)) + # Add pairwise comparisons p-value
  #stat_compare_means(comparisons = my_comparisons, method='t.test', label.y=c(9.1,10.3,11.5)) + # Add pairwise comparisons p-value
  #stat_compare_means(label.y = 13, method='anovafill     # Add global p-value
  labs(y= "GFP (A.U.)", x = "Sample")+
  theme(legend.position = "none", panel.background = element_rect(colour = "black", size=1.5))+
  scale_fill_manual(values=c('black', '#F49D09','white'), labels=c("pUC57-WT", "pUC57-A1664G", 'no plasmid'))
ggsave('X://volume2//noam//gfp_1664//gfp_final.png', dpi = 800)
