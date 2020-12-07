save <- plyr::rbind.fill(
  readRDS("./results/revision/friction_ser1.Rda"),
  readRDS("./results/revision/friction_ser2.Rda"),
  readRDS("./results/revision/friction_ser3.Rda"),
  readRDS("./results/revision/friction_ser4.Rda")
)

source("./R/plots.R")

save$test <- as.character(save$test)

plot_on_T(nice.names = TRUE,revision.levels = T,
          subset(save, test %in% c("splineP","splineP_lag","quantiles","quantiles_lag","quantile","quantile_lag"#,"median_lag","median"
                                   )),
          size.adjusted = TRUE,alpha.level = .1)+theme_classic(20)
ggsave("power_revision2000.pdf", height = 5, width=8)


plot_on_T(nice.names = TRUE,revision.levels = T,
          subset(save, test %in% c("splineP","splineP_lag","expectiles","expectiles_lag","expectile","expectile_lag"#,"mean","mean_lag"
                                   )),
          size.adjusted = TRUE,alpha.level = .1)+theme_classic(20)
ggsave("power_revision_supp2000.pdf", height = 5, width=8)


plot_on_T(nice.names = TRUE,revision.levels = 2,
          subset(save, test %in% c("splineP","expectilesP","quantilesP","splineP_lag","expectilesP_lag","quantilesP_lag")),
          size.adjusted = TRUE,alpha.level = .1)+theme_classic(20)
ggsave("power_revision_supp_inst2000.pdf", height = 5, width=8)



