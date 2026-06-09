library(dplyr)
library(here)
library(ggplot2)
get_AR <- function(x, level = 0.95) {
  clevs <- c((1 - level) / 2, (1 + level) / 2)
  theta <- getME(x, "theta")
  vv <- vcov(x, full = TRUE)
  cpos <- which(colnames(vv) == "theta_YEAR3+0|group.2")
  s2 <- sqrt(vv[cpos, cpos])
  t2 <- theta[2]
  
  list(
    est = theta2phi(t2),
    ci = theta2phi(qnorm(clevs, t2, s2))
  )
}


biotime=fread(here("data/data_filtered_19_nov.csv"))

insects_biotime_models = biotime %>% filter(class=="Insecta") %>% 
  group_by(MODEL_ID) %>% 
  summarise (duration=length(YEAR)) 

model_1 = biotime %>% filter(MODEL_ID=="788-NA-Eristalis tenax")

plot(model_1$ABUNDANCE ~ model_1$YEAR, type = "b")

model_1$year = model_1$YEAR
model_1$YEAR3 =  factor(model_1$YEAR2)
model_1$abund =  model_1$ABUNDANCE
model_1$repi = 1

df=model_1
i=1

    m <- glmmTMB(abund ~ YEAR2 + ar1(YEAR3 + 0 | group), data = df, family = poisson)
    m_aut_null <- glmmTMB(abund ~ YEAR2, data = df, family = poisson)
    m_null <- glmmTMB(abund ~ 1, data = df, family = poisson)
    
    df_model=data.frame(
      repi=i,
      var = var(m$frame$abund),
      duration = length(unique(df$YEAR2)),
      mean = mean(m$frame$abund),
      slope_pval = summary(m)$coefficients$cond[2, "Pr(>|z|)"],
      slope_estimated = fixef(m)$cond[[2]],
      r_estimated = exp(fixef(m)$cond[[2]]) - 1,
      sdev_estimated = as.vector(exp(m$sdr$par.fixed[3])), 
      phi_estimated = attr(VarCorr(m)$cond$group, "correlation")[2], 
      aic = AIC(m), 
      aic_null = AIC(m_null), 
      r2_conditional  = 1 - logLik(m)[1] / logLik(m_null)[1], 
      r2_marginal     = 1 - logLik(m_aut_null)[1] / logLik(m_null)[1], 
      convergence = m$fit$convergence == 0, 
      message = m$fit$message
    )
    

df_model

df_model$sdev_estimated

res_sim=simulate_ts (r=0, phi=df_model$phi_estimated, sdev=df_model$sdev_estimated, mu_zero=mu_zero, years=1:6, n_rep=100)
df_model_sim=run_model(res_sim)


check_results <- function (model_empiric, model_simulated) {
  
  quantiles <- quantile (model_simulated$r_estimated,probs=c(0.025,0.975))
  compatible_with_null_model <- quantiles[1] < model_empiric$r_estimated & model_empiric$r_estimated < quantiles[2]
  cat("Is it compatible with stability?", compatible_with_null_model, "\n")
  z_score <- (model_empiric$r_estimated - mean(model_simulated$r_estimated)) / sd(model_simulated$r_estimated)
  p_value_z <- 2 * (1 - pnorm(abs(z_score)))
  cat("p value null model = ", p_value_z, "\n")
  cat("p value glm = ", round(model_empiric$slope_pval,7), "\n")

}

check_results(model_empiric = df_model, model_simulated = df_model_sim)

plot_results <- function (model_empiric, model_simulated) {
  
  q <- quantile(
    model_simulated$r_estimated,
    c(0.025, 0.975)
  )

ggplot(model_simulated, aes(r_estimated)) +
  geom_text(aes(x=0, y=13,label=paste0("model p-value= ",round(model_empiric$slope_pval,3), "\n",
                                       "null model p-value= ",round(2 * (1 - pnorm(abs(z_score))),3)))) +
  geom_histogram() +
  geom_vline(xintercept = q,
             linetype = 2,
             colour = "blue") +
  geom_vline(
    xintercept = model_empiric$r_estimated,
    colour = "red") + theme_minimal()+
  labs(x = "Estimated trend (r)", y = "Frequency")


}

plot_results(model_empiric = df_model, model_simulated = df_model_sim)

