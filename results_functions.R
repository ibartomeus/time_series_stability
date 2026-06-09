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


model_1 = biotime %>% filter(MODEL_ID=="788-NA-Eristalis intricaria")

model_1$year = model_1$YEAR
model_1$YEAR3 =  factor(model_1$YEAR2)
model_1$abund =  model_1$ABUNDANCE
model_1$repi = 1
model_1$phi_sim=rep(2, nrow(model_1))
model_1$sdev_sim=rep(2, nrow(model_1))
model_1$r_sim =rep(2, nrow(model_1))

plot(model_1$ABUNDANCE ~ model_1$YEAR, type = "b")

window=1:6
n_models=nrow(model_1)-length(window)+1
model_chunks <- vector("list", n_models)

for(i in 1:n_models){
  model <- model_1[window + i - 1, ]
  model$repi <- i
  model_chunks[[i]] <- model
}

models_window_df=do.call(rbind,model_chunks)

LTS = run_model(model_1)
STS = run_model(models_window_df)
STS_sim_df = simulate_ts (r=0, phi=LTS$phi_estimated, sdev=LTS$sdev_estimated, mu_zero=LTS$mean, years=1:6, n_rep=100)
STS_sim = run_model(STS_sim_df)


check_results <- function (model_empiric, model_simulated) {
  
  quantiles <- quantile (model_simulated$r_estimated,probs=c(0.025,0.975))
  compatible_with_null_model <- quantiles[1] < model_empiric$r_estimated & model_empiric$r_estimated < quantiles[2]
  quantile (model_simulated$r_estimated)
  model_empiric$r_estimated
  #cat("Is it compatible with stability?", compatible_with_null_model, "\n")
  perc=length(model_simulated$r_estimated[which(model_simulated$r_estimated < model_empiric$r_estimated)])/100
  
  # z_score <- (model_empiric$r_estimated - mean(model_simulated$r_estimated)) / sd(model_simulated$r_estimated)
  # p_value_z <- 2 * (1 - pnorm(abs(z_score)))
  # cat("p value null model = ", p_value_z, "\n")
  # cat("p value glm = ", round(model_empiric$slope_pval,7), "\n")
  return(perc)

}


table_full=NULL

for(i in 1:length(STS$repi)){
  quantiles <- quantile (STS_sim$r_estimated,probs=c(0.025,0.975))
  
table=data.frame(
  percentil=check_results(model_empiric = STS[i,], model_simulated = STS_sim),
  quantil_5_null =quantiles[1][1],
  quantil_95_null =quantiles[2][1],
  pval_glm_STS=STS[i,]$slope_pval,
  pvalue_LTS = LTS$slope_pval, 
  r_LTS = LTS$r_estimated, 
  sdev_LTS = LTS$sdev_estimated, 
  phi_LTS = LTS$phi_estimated,
  slope_glm_STS = STS[i,]$slope_estimated,
  GLM_significant = STS[i,]$slope_pval<0.05,
  LTS_significant =  LTS$slope_pval<0.05,
  compatibility = quantiles[1] < STS[i,]$r_estimated & STS[i,]$r_estimated < quantiles[2]
  )

table_full=rbind(table_full, table)

}

table_full

plot_results <- function (model_empiric, model_simulated) {

  q <- quantile(model_simulated$r_estimated, c(0.025, 0.975))
  z_score <- (model_empiric$r_estimated - mean(model_simulated$r_estimated)) / sd(model_simulated$r_estimated)
  p_value_z <- 2 * (1 - pnorm(abs(z_score)))
  compatible_with_null_model <- q[1] < model_empiric$r_estimated & model_empiric$r_estimated < quantiles[2]

ggplot(model_simulated, aes(r_estimated)) +
  geom_text(aes(x=0, y=13,label=paste0("model p-value= ",round(model_empiric$slope_pval,3), "\n",
                                       "null model p-value= ",round(p_value_z,3), "\n",
                                       "is it compatible?", compatible_with_null_model))) +
  geom_histogram() +
  geom_vline(xintercept = q,
             linetype = 2,
             colour = "blue") +
  geom_vline(
    xintercept = model_empiric$r_estimated,
    colour = "red") + theme_minimal()+
  labs(x = "Estimated trend (r)", y = "Frequency")


}

plot_results(model_empiric = STS[6,], model_simulated = STS_sim)




