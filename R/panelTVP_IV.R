panelTVP_IV <- function(formula_stage1 = NULL,
                        formula_stage2 = NULL,
                        data = NULL,
                        id = NULL,
                        t = NULL,
                        prior.reg_stage1 = list(
                          a.tau = 0.1, a.xi = 0.1,
                          learn.a.tau = TRUE, learn.a.xi = TRUE,
                          alpha.a.tau = 2, alpha.a.xi = 2,
                          beta.a.tau = 1, beta.a.xi = 1,
                          iota.a.tau = 1, iota.a.xi = 1,
                          target.rate.a.tau = 0.44, target.rate.a.xi = 0.44,
                          c.tau = 0.1, c.xi = 0.1,
                          learn.c.tau = TRUE, learn.c.xi = TRUE,
                          alpha.c.tau = 2, alpha.c.xi = 2,
                          beta.c.tau = 1, beta.c.xi = 1,
                          iota.c.tau = 1, iota.c.xi = 1,
                          target.rate.c.tau = 0.44, target.rate.c.xi = 0.44,
                          kappa.tau = 10, kappa.xi = 10,
                          learn.kappa.tau = TRUE, learn.kappa.xi = TRUE,
                          d.tau = 0.001, d.xi = 0.001,
                          e.tau = 0.001, e.xi = 0.001,
                          type = "rw-t1", c = 1, B0 = 1, TG = FALSE, TG.alternative = FALSE
                          ),
                        prior.load_stage1 = list(
                          a.phi = 0.1, a.zeta = 0.1,
                          kappa.phi = 10, kappa.zeta = 10,
                          learn.kappa.phi = TRUE, learn.kappa.zeta = TRUE,
                          d.phi = 0.001, d.zeta = 0.001,
                          e.phi = 0.001, e.zeta = 0.001,
                          type = "rw-t1", c = 1, L0 = 1
                          ),
                        prior.reg_stage2 = list(
                          a.tau = 0.1, a.xi = 0.1,
                          learn.a.tau = TRUE, learn.a.xi = TRUE,
                          alpha.a.tau = 2, alpha.a.xi = 2,
                          beta.a.tau = 1, beta.a.xi = 1,
                          iota.a.tau = 1, iota.a.xi = 1,
                          target.rate.a.tau = 0.44, target.rate.a.xi = 0.44,
                          c.tau = 0.1, c.xi = 0.1,
                          learn.c.tau = TRUE, learn.c.xi = TRUE,
                          alpha.c.tau = 2, alpha.c.xi = 2,
                          beta.c.tau = 1, beta.c.xi = 1,
                          iota.c.tau = 1, iota.c.xi = 1,
                          target.rate.c.tau = 0.44, target.rate.c.xi = 0.44,
                          kappa.tau = 10, kappa.xi = 10,
                          learn.kappa.tau = TRUE, learn.kappa.xi = TRUE,
                          d.tau = 0.001, d.xi = 0.001,
                          e.tau = 0.001, e.xi = 0.001,
                          type = "rw-t1", c = 1, B0 = 1, TG = FALSE, TG.alternative = FALSE
                          ),
                        prior.load_stage2 = list(
                          a.phi = 0.1, a.zeta = 0.1,
                          kappa.phi = 10, kappa.zeta = 10,
                          learn.kappa.phi = TRUE, learn.kappa.zeta = TRUE,
                          d.phi = 0.001, d.zeta = 0.001,
                          e.phi = 0.001, e.zeta = 0.001,
                          type = "rw-t1", c = 1, L0 = 1
                          ),
                        prior.var_stage2 = list(
                          learn.C0.hyp = list(g0 = 5, G0 = 3.333333), c0 = 2.5
                        ),
                        mcmc.opt = list(
                          chain.length = 12000, burnin = 2000, thin = 10, asis = TRUE
                          ),
                        HPD.coverage = 0.95,
                        random.effects_stage1 = FALSE,
                        random.effects_stage2 = TRUE,
                        progress.bar = FALSE
){

  # HERE INPUT CHECKS ARE NEEDED !!!

  # ordering dataset and removing gaps in id
  data$id <- as.numeric(factor(id))
  data$t <- t
  data <- data[order(data$t, data$id),]

  if(prior.reg_stage1$type == "rw-t0") prior.reg_stage1$type <- "rw1"
  if(prior.reg_stage1$type == "rw-t1") prior.reg_stage1$type <- "rw2"
  if(prior.load_stage1$type == "rw-t0") prior.load_stage1$type <- "rw1"
  if(prior.load_stage1$type == "rw-t1") prior.load_stage1$type <- "rw2"

  if(prior.reg_stage2$type == "rw-t0") prior.reg_stage2$type <- "rw1"
  if(prior.reg_stage2$type == "rw-t1") prior.reg_stage2$type <- "rw2"
  if(prior.load_stage2$type == "rw-t0") prior.load_stage2$type <- "rw1"
  if(prior.load_stage2$type == "rw-t1") prior.load_stage2$type <- "rw2"

  result <- fit_panelTVP_IV(formula_stage1 = formula_stage1,
                            formula_stage2 = formula_stage2,
                            data = data,
                            prior.reg_stage1 = prior.reg_stage1,
                            prior.load_stage1 = prior.load_stage1,
                            prior.reg_stage2 = prior.reg_stage2,
                            prior.load_stage2 = prior.load_stage2,
                            prior.var_stage2 = prior.var_stage2,
                            mcmc.opt = mcmc.opt,
                            HPD.coverage = HPD.coverage,
                            random.effects_stage1 = random.effects_stage1,
                            random.effects_stage2 = random.effects_stage2,
                            progress.bar = progress.bar)

}
