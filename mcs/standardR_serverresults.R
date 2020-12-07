
save <- rbind(
  readRDS("./results/revision/standard_ser1.Rda"),
  readRDS("./results/revision/standard_ser2.Rda")
)

source("./R/plots.R")

# ## show results
plot_on_T(nice.names = TRUE,revision.levels = T,
          subset(save,
                 test %in% c("spline","quantiles_inst","quantile")),size.adjusted = F,alpha.level = .1)+theme_classic(20)
ggsave("size_revision3.pdf", height = 5, width=8)


plot_on_T(nice.names = TRUE,revision.levels = T,
          subset(save, test %in% c("spline","expectiles_inst","expectile")),size.adjusted = F,alpha.level = .1)+theme_classic(20)
ggsave("size_revision_supp3.pdf", height = 5, width=8)


plot_on_T(nice.names = TRUE,revision.levels = 2,
          subset(save, test %in% c("spline","quantiles","expectiles")),size.adjusted = F,alpha.level = .1)+theme_classic(20)
ggsave("size_revision_supp_inst3.pdf", height = 5, width=8)



