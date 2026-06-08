#### POISSON log-normal SIMULATIONS
library(dplyr)
library(mvtnorm)
library(glmmTMB)
library(MASS)
library(data.table)

# ------------------------------------------------------------
# 1. Simulate n Poisson-lognormal AR1 time series
# ------------------------------------------------------------

n_years=10
years=1:n_years
r_obs=-0.02
phi=0.2
sdev=0.20
mu_zero = runif(1, 0.2, 90)

simulate_ts <- function(r, phi, sdev, mu_zero, years, n_rep) {
  
  df_all=NULL
  for (i in 1:n_rep) {
    n_years <- length(years)
    
    trend <- mu_zero * (1 + r) ^ (years - 1)
    
    R <- phi ^ as.matrix(dist(years))
    
    lognormal_noise <- MASS::mvrnorm(
      n = 1,
      mu = rep(0, n_years),
      Sigma = (sdev^2) * R
    )
    
    mu <- exp(log(trend) + lognormal_noise)
    
    df=data.frame(
      repi = i,
      year = years,
      YEAR2 = years - min(years),
      YEAR3 = factor(years - min(years)),
      group = factor(1),
      trend = trend,
      mu = mu,
      abund = rpois(n_years, mu),
      r_sim = r,
      phi_sim = phi,
      sdev_sim = sdev,
      mu_zero = mu_zero
    )
    df_all=rbind(df_all,df)
    
  }
  
  df_all
}

res_obs=simulate_ts (r=r_obs, phi=phi, sdev=sdev, mu_zero=mu_zero, years=years, n_rep=1)
res_sim=simulate_ts (r=0, phi=phi, sdev=sdev, mu_zero=mu_zero, years=years, n_rep=100)

plot(res_obs$abund ~ res_obs$year, type = "b", main = paste0("phi=", phi, ", sdev=", sdev, ", r=", r_obs))

# ------------------------------------------------------------
# 2. Run model and extract time-series characteristics
# ------------------------------------------------------------
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

run_model <- function (df_full) {
  models=NULL
  for(i in unique(df_full$repi)) {
  df=df_full %>% filter(repi == i)
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
    r_sim = unique(df$r_sim),
    sdev_estimated = as.vector(exp(m$sdr$par.fixed[3])), 
    sdev_sim = unique(df$sdev),
    phi_estimated = attr(VarCorr(m)$cond$group, "correlation")[2], 
    phi_sim = unique(df$phi),
    aic = AIC(m), 
    aic_null = AIC(m_null), 
    r2_conditional  = 1 - logLik(m)[1] / logLik(m_null)[1], 
    r2_marginal     = 1 - logLik(m_aut_null)[1] / logLik(m_null)[1], 
    convergence = m$fit$convergence == 0, 
    message = m$fit$message
  )
    
  models=rbind(models, df_model)
  }
  models
}

df_model_sim=run_model(res_sim)
df_model_obs=run_model(res_obs)

quantiles=quantile(df_model_sim$r_estimated,probs=c(0.025,0.975))
hist(df_model_sim$r_estimated)
compatible_with_null_model <- quantiles[1] < df_model_obs$r_estimated & df_model_obs$r_estimated < quantiles[2]
compatible_with_null_model

mean(df_model_sim$r_estimated)
sd(df_model_sim$r_estimated)
df_model_obs$r_estimated

z_score <- (df_model_obs$r_estimated - mean(df_model_sim$r_estimated)) / sd(df_model_sim$r_estimated)
df_model_obs$slope_pval
2 * (1 - pnorm(abs(z_score)))

q <- quantile(
  df_model_sim$r_estimated,
  c(0.025, 0.975)
)

ggplot(df_model_sim, aes(r_estimated)) +
  geom_text(aes(x=0, y=13,label=paste0("model p-value= ",round(df_model_obs$slope_pval,3), "\n",
                                       "null model p-value= ",round(2 * (1 - pnorm(abs(z_score))),3)))) +
  geom_histogram() +
  geom_vline(xintercept = q,
             linetype = 2,
             colour = "blue") +
  geom_vline(
    xintercept = df_model_obs$r_estimated,
    colour = "red") +
  labs(x = "Estimated trend (r)", y = "Frequency")


