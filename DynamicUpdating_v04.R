################################################################################
#################### Load libraries  ###########################################
################################################################################
library(here)
library(readxl)
library(dplyr)
library(ggplot2)
library(pROC)
library(furniture)
library(tidyverse)
library(rms)
library(gtsummary)
library(ggalt)
library(ggsci)
library(dcurves)
library(ggpubr)
library(dma)
library(mice)
library(Hmisc)

################################################################################
#################### Functions  ################################################
################################################################################

yearly.validation <- function(glm_pred, method, offset, val_alldata){
  # Loop generating validation results for each year separately 
  glm_val_results <- c()
  for(i in 2016:2018){
    dat <- glm_pred %>% filter(year==i)
    mod1 <- glm(outcome ~ prob, family="binomial",x=TRUE,y=TRUE, data= dat)
    c_slope <- mod1$coefficients[2]
    mod1_ci <- confint(mod1)
    mod2 <- glm(outcome~offset(prob),family="binomial", data = dat)
    citl <- mod2$coefficients
    mod2_ci <- confint(mod2)
    g <- roc(outcome ~ risk, data = dat, direction = "<")
    disc <- g$auc
    gci <- ci(roc(outcome ~ risk, data = dat, direction = "<"))
    O <- prop.table(table(dat$outcome))[2]*length(dat$outcome)
    E <- mean(dat$risk)*length(dat$outcome)
    o_e <- O/E
    o_e_u <- ((sqrt(O)+1.96*0.5)^2)/E
    o_e_l <- ((sqrt(O)-1.96*0.5)^2)/E
    bs <- 1/length(dat[,1]) * sum((dat$risk - dat$outcome)^2) 
    glm_val_results <- rbind(glm_val_results, c(i,citl, mod2_ci[1], mod2_ci[2],
                                                c_slope, mod1_ci[2,1], 
                                                mod1_ci[2,2],
                                                o_e, o_e_u, o_e_l, bs, disc, 
                                                gci[1],gci[3], 0))
  }
  
  glm_val_results <- rbind(glm_val_results, val_alldata)
  glm_val_results <- as.data.frame(glm_val_results)
  names(glm_val_results)[9] <- "2"
  names(glm_val_results)[10] <- "3"
  old_names <- names(glm_val_results)
  new_names <- c("year", "citl", "citl_li", "citl_ui", "c_slope", "c_slope_li", 
                 "c_slope_ui", "o_e", "o_e_ui", "o_e_li", "brier", "disc", 
                 "disc_li", "disc_ui", "logLik")
  glm_val_results <- glm_val_results %>% data.table::setnames(old = old_names, 
                                                              new = new_names)
  # view(glm_val_results)
  
  
  vcdf <- glm_val_results[-9,]
  vcdf[] <- lapply(vcdf, function(x) as.numeric(as.character(x)))
  vcdf$year <- vcdf$year + offset
  vcdf$mod <- method
  return(vcdf)
}

period.validation <- function(glm_pred){
  
  # c-slope
  mod1 <- glm(outcome ~ prob, family="binomial",x=TRUE,y=TRUE, data= glm_pred)
  #mod1 # 
  mod1_ci <- confint(mod1)
  
  # CITL
  mod2 <- glm(outcome~offset(prob),family="binomial", data = glm_pred)
  mod2_ci <- confint(mod2)
  
  # Discrimination
  g <- roc(outcome ~ risk, data = glm_pred, direction = "<")
  gci <- ci(roc(outcome ~ risk, data = glm_pred, direction = "<"))
  
  # E/O
  prop.table(table(glm_pred$outcome))
  O <- prop.table(table(glm_pred$outcome))[2]*length(glm_pred$outcome)
  E <- mean(glm_pred$risk)*length(glm_pred$outcome)
  o_e_u <- ((sqrt(O)+1.96*0.5)^2)/E
  o_e_l <- ((sqrt(O)-1.96*0.5)^2)/E
  
  # Log/likelihood
  loglik <- sum(dbinom(glm_pred$outcome, glm_pred$risk, size=1, log=TRUE))
  
  bs <- 1/length(glm_pred[,1]) * sum((glm_pred$risk - glm_pred$outcome)^2) 
  
  glm_val_alldat <- c("16-18", mod2$coefficients, mod2_ci[1],
                      mod2_ci[2], mod1$coefficients[2],
                      mod1_ci[2,1], mod1_ci[2,2],
                      O/E, o_e_u, o_e_l, bs, g$auc, gci[1],gci[3], loglik)
  
  
  return(glm_val_alldat)
  
}

logLik_test <- function(logLik1, logLik2){
  
  teststat <- -2 * (as.numeric(logLik1)-as.numeric(logLik2))
  p.val <- pchisq(teststat, df = 1, lower.tail = FALSE)
  return(p.val)
}
  

################################################################################
#################### Load UHR 1000 tables  #####################################
################################################################################

uhr1000 <- read.csv("U:/PostDoc/Research/UHR1000+/UHR1000_MasterFile_V02_18032024.csv", 
                    header = TRUE, sep=",")

################################################################################
#################### Convert and prepare variables  ############################
################################################################################

# Convert GAF scores to SOFAS scores
gaf <- read.csv(paste(dirname(rstudioapi::getSourceEditorContext()$path),"/GAF_SOFAS.csv", sep = ""), header = TRUE, sep=",") # Load conversion table
names(gaf) <- c("GAF","SOFAS")
gaf$SOFAS <- round(gaf$SOFAS)

uhr1000$gaf_sofas <- uhr1000$sofas_currscore_0

for(index in c(1:nrow(gaf))){
  uhr1000[uhr1000$gafrat_0 == gaf[index,"GAF"] & !is.na(uhr1000$gafrat_0 == gaf[index,"GAF"]) & is.na(uhr1000$sofas_currscore_0), "gaf_sofas"] <- gaf[index,"SOFAS"]
}

# Calculate UHR category based on dominating UHR group
uhr1000$uhr_cat <- uhr1000$uhr_group

uhr1000[(uhr1000$uhr_group == 3 | uhr1000$uhr_group == 5  | uhr1000$uhr_group == 7) & !is.na(uhr1000$uhr_group), c('uhr_cat')] <- 1
uhr1000[uhr1000$uhr_group == 6 & !is.na(uhr1000$uhr_group), c('uhr_cat')] <- 2
uhr1000[uhr1000$uhr_group == 4 & !is.na(uhr1000$uhr_group), c('uhr_cat')] <- 3

uhr1000$uhr_cat <- as.factor(uhr1000$uhr_cat)

# Extract just Melbourne sample
uhr1000 <- subset(uhr1000, site == 1)

uhr1000$year <- format(as.Date(uhr1000$assessment_date_0, 
                               format="%d/%m/%Y"),"%Y")
uhr1000$year <- as.numeric(uhr1000$year)

# Exclude invdividuals with no follow-up time information
uhr1000 <- uhr1000[!is.na(uhr1000$transdays),]

# Exclude invdividuals who did not transition and drop out of study 
# prior to 1 year
uhr1000 <- uhr1000[!(uhr1000$transdays < 365 & uhr1000$transtat == 0),]

# Set 1-year as follow-up, anyone who transitions after 1 year is marked as 
# non-transition
uhr1000$transtat <- uhr1000$transtat * (uhr1000$transday < 1*365.25)
uhr1000[uhr1000$transdays > 1*365.25 & !is.na(uhr1000$transdays), 
        'transdays'] <- 1*365.25

datOrig <- uhr1000 %>% filter(year<2015)
datval <- uhr1000 %>% filter(year>2014)

################################################################################
#################### Prepare data ##############################################
################################################################################

uhr1000.org <- uhr1000

# Exclude individuals with missing data
#uhr1000 <- uhr1000[!is.na(uhr1000$caarms_DS_sev_0),]
#uhr1000 <- uhr1000[!is.na(uhr1000$caarms_UTC_sev_0),]
#uhr1000 <- uhr1000[!is.na(uhr1000$gaf_sofas),]
#uhr1000 <- uhr1000[!is.na(uhr1000$uhr_cat),]

vars <- c("caarms_DS_sev_0", "caarms_UTC_sev_0", "gaf_sofas", "sans_tot_0", "transtat",
          "year")

uhr1000 <- uhr1000[,vars]

#Calculating imputed values with aregImpute
#impute_arg <- aregImpute(~ caarms_DS_sev_0 + caarms_UTC_sev_0 + gaf_sofas + sans_tot_0 + transtat + year, data = uhr1000, n.impute = 5)

#imputed <- impute.transcan(impute_arg, imputation=1, data=uhr1000, list.out=TRUE,pr=FALSE, check=FALSE) 
# convert the list to the database
#imputed.data <- as.data.frame(do.call(cbind,imputed))

# arrange the columns accordingly
#uhr1000.imputed <- imputed.data[, colnames(uhr1000), drop = FALSE]

uhr1000.imputed <- mice(uhr1000, m = 1, maxit = 30, 
                    method = "norm", print =  FALSE, seed = 1235)

uhr1000.imputed <- complete(uhr1000.imputed, action = "long", include = FALSE)



dat_f <- uhr1000.imputed

dat_f <- dat_f[order(dat_f$year),]
dat_f <- dat_f[dat_f$year < 2019,]

df_val <- dat_f %>% filter(year>2014)
df <- dat_f %>% filter(year<2015)

################################################################################
#################### Method 1: No updating #####################################
################################################################################

# Run model
fit <- lrm(transtat ~ caarms_DS_sev_0 + caarms_UTC_sev_0 + gaf_sofas + sans_tot_0, 
           data = df, x=TRUE, y=TRUE)

dat_orig <- df[,c('transtat','year')]
# Calculate risk
dat_orig$risk <- exp(fit$linear.predictors)/(1+exp(fit$linear.predictors))

# Calibration plot in development set
glmOld_cal = lowess(dat_orig$risk, dat_orig$transtat, iter=0)

# dca plot in development set
tmp <- dca(transtat ~ risk, dat_orig, thresholds = seq(0.05, 0.5, 0.01)) %>%
  plot(smooth = TRUE) 
glmOld_dca <- tmp$plot_env$.[,c(1,4,9)]

# Calculate probabilities in validation data
glm_pred <- predict(fit, newdata = df_val[vars])
glm_pred <- as.data.frame(glm_pred)
names(glm_pred)[1] <- "prob"

# Transfor to risk
glm_pred$risk <- exp(glm_pred$prob)/(1+exp(glm_pred$prob))

glm_pred <- cbind(glm_pred,df_val$transtat,df_val$year)

names(glm_pred)[3] <- "outcome"
names(glm_pred)[4] <- "year"

# Remove 2015 data from validation
glm_pred <- subset(glm_pred, year != 2015)

####################
#### Validation ####
####################

# Calculate performance across entire period
val_alldata <- period.validation(glm_pred)

# Calibration plot in validation period
glm_cal = lowess(glm_pred$risk, glm_pred$outcome, iter=0)

# dca plot in validation period
tmp <- dca(outcome ~ risk, glm_pred, thresholds = seq(0.05, 0.5, 0.01)) %>%
  plot(smooth = TRUE) 
glm_dca <- tmp$plot_env$.[,c(1,4,9)]

# Calculate performance across years
df_original <- yearly.validation(glm_pred, "a_orig", -0.3, val_alldata)

################################################################################
#################### Method 2: Yearly recalibration ############################
################################################################################

# Run model
fit <- lrm(transtat ~ caarms_DS_sev_0 + caarms_UTC_sev_0 + gaf_sofas + sans_tot_0, 
           data = df, x=TRUE, y=TRUE)

# Calculate probabilities in validation data
lp <- predict(fit, newdata = df_val[vars])

df_val$lp14 <- lp

df_val15 <- df_val %>% filter(year==2015)

fit <- lrm(transtat ~ lp14, 
           data = df_val15, x=TRUE, y=TRUE)

lp <- predict(fit, newdata = df_val[,c('lp14')])

df_val$lp15 <- lp

# now update using 2016 data and predict LP in all val data
df_val16 <- df_val %>% filter(year==2016)

fit <- lrm(transtat ~ lp15, 
           data = df_val16, x=TRUE, y=TRUE)

lp <- predict(fit, newdata = df_val[,c('lp15')])

df_val$lp16 <- lp

# now update using 2017 data and predict LP in all val data
df_val17 <- df_val %>% filter(year==2017)

fit <- lrm(transtat ~ lp16, 
           data = df_val17, x=TRUE, y=TRUE)

lp <- predict(fit, newdata = df_val[,c('lp16')])

df_val$lp17 <- lp


# now update using 2018 data and predict LP in all val data
df_val18 <- df_val %>% filter(year==2018)

fit <- lrm(transtat ~ lp17, 
           data = df_val18, x=TRUE, y=TRUE)

lp <- predict(fit, newdata = df_val[,c('lp17')])

df_val$lp18 <- lp

# Store probability estimate for each year in data frame
df_val$prob <- ifelse(df_val$year==2015, df_val$lp14, 99)
df_val$prob <- ifelse(df_val$year==2016, df_val$lp15, df_val$prob)
df_val$prob <- ifelse(df_val$year==2017, df_val$lp16, df_val$prob)
df_val$prob <- ifelse(df_val$year==2018, df_val$lp17, df_val$prob)


glm_pred <- as.data.frame(df_val$prob)
names(glm_pred)[1] <- "prob"

# Transfor to risk
glm_pred$risk <- exp(glm_pred$prob)/(1+exp(glm_pred$prob))

glm_pred <- cbind(glm_pred,df_val$transtat,df_val$year)

names(glm_pred)[3] <- "outcome"
names(glm_pred)[4] <- "year"

# Remove 2015 data from validation
glm_pred <- subset(glm_pred, year != 2015)

####################
#### Validation ####
####################

val_alldata <- period.validation(glm_pred)

# Calibration plot
glmRecal_cal = lowess(glm_pred$risk, glm_pred$outcome, iter=0)

# dca plot
tmp <- dca(outcome ~ risk, glm_pred, thresholds = seq(0.05, 0.5, 0.01)) %>%
  plot(smooth = TRUE) 

glmRecal_dca <- tmp$plot_env$.[,c(1,4,9)]

df_recal <- yearly.validation(glm_pred, "b_recal", -0.1, val_alldata)

################################################################################
#################### Method 3: Continuous refitting ############################
################################################################################

df_val <- dat_f %>% filter(year>2014)

# Run model
fit <- lrm(transtat ~ caarms_DS_sev_0 + caarms_UTC_sev_0 + gaf_sofas + sans_tot_0, 
           data = df, x=TRUE, y=TRUE)

# Calculate probabilities in validation data
p_glm <- predict(fit, newdata = df_val[vars])
p_glm <- as.data.frame(p_glm)

dat_train <- na.omit(df[vars])
dat_val <- na.omit(df_val[vars])

for(i in c(1:(nrow(df_val)-2))){
  
  dat <- rbind(dat_train[i:nrow(dat_train),], dat_val[1:i,])
  
  fit <- lrm(transtat ~ caarms_DS_sev_0 + caarms_UTC_sev_0 + gaf_sofas + sans_tot_0, 
             data = dat, x=TRUE, y=TRUE)

  # Calculate probabilities in validation data
  p_glm[(i+1):nrow(p_glm),] <- predict(fit, newdata = df_val[(i+1):nrow(dat_val),
                                                             vars])
}

glm_pred <- p_glm

# Sum over rows
names(glm_pred)[1] <- "prob"

# Transfor to risk
glm_pred$risk <- exp(glm_pred$prob)/(1+exp(glm_pred$prob))

glm_pred <- cbind(glm_pred,dat_val$transtat,dat_val$year)

names(glm_pred)[3] <- "outcome"
names(glm_pred)[4] <- "year"

# Remove 2015 data from validation
glm_pred <- subset(glm_pred, year != 2015)

####################
#### Validation ####
####################

val_alldata <- period.validation(glm_pred)

# Calibration plot
glmRefit_cal = lowess(glm_pred$risk, glm_pred$outcome, iter=0)

# dca plot
tmp <- dca(outcome ~ risk, glm_pred, thresholds = seq(0.05, 0.5, 0.01)) %>%
  plot(smooth = TRUE) 

glmRefit_dca <- tmp$plot_env$.[,c(1,4,9)]

df_refit <- yearly.validation(glm_pred, "c_refit", 0.1, val_alldata)

################################################################################
#################### Method 4: Dynamic updating ################################
################################################################################
### Start with DM with forgetting

x <- dat_f[c("caarms_DS_sev_0", "caarms_UTC_sev_0", "gaf_sofas", "sans_tot_0")]
x <- as.matrix(x)
y <- dat_f$transtat
mmat<-matrix(c(1,1,1,1),1,4,byrow=TRUE)

dma.test<-logistic.dma(x,y,mmat,lambda=.99,alpha=.99, autotune = TRUE, 
                       initialsamp=524)

betas_dm <- t(as.data.frame(dma.test$theta))

dm <- matrix(betas_dm,nrow = dim(dat_f)[1],ncol = 5)

# drop last coef obs
dm[525,] <- dm[526,]

dm_val <- dm[(dim(dm)[1]-dim(df_val)[1]+1):(dim(dm)[1]),1:5]

# generate val data with variables in DM and add intercept term to val data
x <- df_val[vars]
x <- x %>% replace(is.na(.), 0)
dat <- na.omit(x)
y <- dat$transtat
x1 <- as.data.frame(sapply(dat[,1:(dim(dat)[2]-2)], as.numeric))
x1 <- cbind(1,x1)

# dm_val and val data now match.
# multiply and then row sum to get LP
dm_val <- dm_val * x1

# Sum over rows
dm_val <- as.data.frame(rowSums(dm_val))
names(dm_val)[1] <- "prob"

# Transfor to risk
dm_val$risk <- exp(dm_val$prob)/(1+exp(dm_val$prob))

dm_val <- cbind(dm_val,dat$transtat,dat$year)

names(dm_val)[3] <- "outcome"
names(dm_val)[4] <- "year"

# Remove 2015 data from validation
dm_val <- subset(dm_val, year != 2015)

dm_val <- dm_val %>% filter(risk!="NaN")

####################
#### Validation ####
####################

val_alldata <- period.validation(dm_val)

# Calibration plot
glmBayes_cal = lowess(dm_val$risk, dm_val$outcome, iter=0)

# dca plot
tmp <- dca(outcome ~ risk, dm_val, thresholds = seq(0.05, 0.5, 0.01)) %>%
  plot(smooth = TRUE) 

glmBayes_dca <- tmp$plot_env$.[,c(1,4,9)]

df_bayes <- yearly.validation(dm_val, "d_bayes", 0.3, val_alldata)

################################################################################
#################### Plot results ##############################################
################################################################################

#### Performance metrics across years ####
df_plot <- rbind(df_original, df_recal, df_refit, df_bayes)

df_long <- gather(df_plot, pm, pm_num ,c(citl, c_slope, disc, o_e), factor_key=TRUE)
lower <- gather(df_plot, pm, lower ,c(citl_li, c_slope_li, disc_li, o_e_li), factor_key=TRUE)
upper <- gather(df_plot, pm, upper ,c(citl_ui, c_slope_ui, disc_ui, o_e_ui), factor_key=TRUE)
df_long$lower <- lower$lower
df_long$upper <- upper$upper

df_long <- df_long[!is.na(df_long$year),]
df_long$intercept <- c(rep(0,12), rep(1,36))

# Plot CITL, c_slope and discrimination
# Plot CITL, c_slope and discrimination
g <- ggplot(df_long, aes(x=year, y=pm_num, color=mod)) + geom_point(shape=16,size=3) +
  geom_segment(aes(x = year, y = lower, xend = year, yend = upper, color=mod), size = 1, data = df_long)

g.labs <- c("Calibration-in-the-large", "Calibration slope", "Discrimination", "Observed-expected ratio")
names(g.labs) <- c("citl","c_slope", "disc", "o_e")

g + facet_wrap(~ pm, labeller = labeller(pm = g.labs), scales = "free")  + 
  scale_color_manual(labels = c("Original model",  "Yearly recalibration", "Refitting", "Dynamic Updating"),
                     values = c("#374e55",  "#00a1d5", "#b24745", "#df8f44")) +
  labs(x = "Year", y = "Performance value", color = "Model") + theme_light() + 
  scale_x_continuous(breaks=seq(2016,2018,1), labels=seq(2016,2018,1)) + theme(legend.position = "bottom") + 
  geom_abline(aes(intercept = intercept, slope = 0), linetype = 2 ) +
  theme(plot.subtitle = element_text(size = 14),
         plot.title = element_text(size = 16, hjust = 0.5, face="bold", margin=margin(b = 10, unit = "pt")),
         axis.title.x = element_text(size = 12),
         axis.title.y = element_text(size = 12),
         legend.text = element_text(size = 12),
         legend.title = element_text(size = 12),
         axis.text.x = element_text(size = 11),
         axis.text.y = element_text(size = 11),
         panel.grid.major.y = element_line(),
         legend.key.width = unit(2, "line"),
         panel.grid.minor.y = element_line(),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill="#36454F"))

#### Calculate logLike test ####
logLik_vec <- c( df_original$logLik[4],  df_recal$logLik[4],  df_refit$logLik[4], df_bayes$logLik[4])

logLik_p <- matrix(0,4,4)

for(i in c(1:4)){
  for(j in c(1:4)){
    logLik_p[i,j] <- logLik_test(logLik_vec[i], logLik_vec[j])
  }
}

logLik_p
#### Save results ####

write.csv(df_plot, "U:/PostDoc/Research/DynamicUpdating/results.csv", col.names = TRUE)

#### Plot calibration ####

glm_cal <- as.data.frame(glm_cal)
glm_cal$mod <- 'glm'
#glmOld_cal <- as.data.frame(glmOld_cal)
#glmOld_cal$mod <- 'glmOld'
glmBayes_cal <- as.data.frame(glmBayes_cal)
glmBayes_cal$mod <- 'xdm'
glmRecal_cal <- as.data.frame(glmRecal_cal)
glmRecal_cal$mod <- 'Recal'
glmRefit_cal <- as.data.frame(glmRefit_cal)
glmRefit_cal$mod <- 'Refit'

#cal.all <- rbind(glm_cal, glmOld_cal, glmRecal_cal, glmRefit_cal, glmBayes_cal)
cal.all <- rbind(glm_cal, glmRecal_cal, glmRefit_cal, glmBayes_cal)

cal.plot <- ggplot(cal.all, aes(x = x, y = y, color = mod)) +
  geom_line(linewidth = 1.5) +
  scale_x_continuous(name =paste("Predicted 1-year Probability of Transition to Psychosis",sep = ""), breaks = c(0.1,0.2,0.3,0.4,0.5), limits=c(0, 0.5)) +
  scale_y_continuous(name=paste("Observed 1-year Probability of Transition to Psychosis",sep = ""), breaks = c(0.1,0.2,0.3,0.4,0.5), limits=c(0, 0.5)) +
  labs(color = '') +
  scale_color_manual(labels = c("Original model (2016-2018)", "Yearly recalibration", "Refitting", "Dynamic Updating"),
                     values = c("#374e55", "#00a1d5", "#b24745", "#df8f44")) +
  geom_abline(aes(intercept = 0, slope = 1), guide = "none")	+
  theme_classic() +
  theme(plot.subtitle = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5, face="bold", margin=margin(b = 20, unit = "pt")),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        panel.grid.major.y = element_line(),
        panel.grid.minor.y = element_line()) + 
        guides(color = guide_legend(nrow = 1)) +
  theme(legend.position="bottom") 


#### Plot decision curve analysis ####
dca.all <- subset(glm_dca, variable != 'risk')
glm_dca2 <- subset(glm_dca, variable == 'risk')
glm_dca2$variable <- 'glm'
#glmOld_dca <- subset(glmOld_dca, variable == 'risk')
#glmOld_dca$variable <- 'glmOld'
glmRecal_dca2 <- subset(glmRecal_dca, variable == 'risk')
glmRecal_dca2$variable <- 'Recal'
glmRefit_dca2 <- subset(glmRefit_dca, variable == 'risk')
glmRefit_dca2$variable <- 'Refit'
glmBayes_dca2 <- subset(glmBayes_dca, variable == 'risk')
glmBayes_dca2$variable <- 'xdm'

#dca.all <- rbind(dca.all, glmOld_dca, glm_dca, glmRecal_dca, glmRefit_dca, glmBayes_dca)
dca.all <- rbind(dca.all, glm_dca2, glmRecal_dca2, glmRefit_dca2, glmBayes_dca2)

dca.plot <- ggplot(dca.all, aes(x = threshold, y = net_benefit, color = variable)) +
  #geom_line(linewidth = 1.5) +
  geom_smooth(method = "loess", se = FALSE , linewidth = 1.5) +
  scale_x_continuous(name = "Threshold Probability", breaks = c(0,0.1,0.2,0.3,0.4), limits=c(0.05, 0.4)) +
  scale_y_continuous(name="Net benefit", breaks = c(0.00, 0.025, 0.05, 0.075), limits=c(-0.02, 0.075)) +
  scale_color_manual(labels = c("all" = "Treat All", "none" = "Treat None", "glm" = "Original model (2016-2018)", 
                                "Recal" = "Yearly recalibration",
                                "Refit" = "Refitting", "xdm" = "Dynamic Updating"),
                     values = c("all" = "#AEAEAE", "none" = "black", "glm" = "#374e55", 
                                "Recal" = "#00a1d5", "Refit"="#b24745", "xdm"="#df8f44" )) +
  labs(colour = "") + 
  theme_classic() +
  theme(plot.subtitle = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5, face="bold", margin=margin(b = 20, unit = "pt")),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        panel.grid.major.y = element_line(),
        panel.grid.minor.y = element_line()) +
        guides(color = guide_legend(nrow = 2)) +
  theme(legend.position="bottom") 

ggarrange(cal.plot, dca.plot, ncol= 2)

#### UHR1000 descriptive table ####

uhr1000.org <- uhr1000.org[uhr1000.org$year < 2019,]

uhr1000.org$year2 <- uhr1000.org$year
uhr1000.org[uhr1000.org$year < 2015, "year2" ] <- 2014

uhr1000.org$year2 <- as.factor(uhr1000.org$year2)

table_desc <- uhr1000.org %>% 
  dplyr::select(year2, assessment_age, gender, timeSxService, uhr_cat, transtat, group, caarms_DS_sev_0, caarms_UTC_sev_0, gaf_sofas, sans_tot_0) %>% # keep only columns of interest
  mutate(
    gender = factor(gender, labels = c("Male", "Female")) ,
    uhr_cat = factor(uhr_cat, labels = c("BLIPS","Attenuated Psychosis","Vulnerability")) ,
    group = factor(group, labels = c("Standard intervention treatment", "Non-standard intervention treatment")),
    transtat = factor(transtat, labels = c("Not transitioned", "Transitioned"))
  ) %>% 
  tbl_summary(  
    by = year2,
    statistic = list(all_continuous() ~ "{mean} ({sd})",        # stats and format for continuous columns
                     all_categorical() ~ "{n} / {N} ({p}%)"),   # stats and format for categorical columns
    digits = all_continuous() ~ 1,                              # rounding for continuous columns
    type   = list(c(gender, group, transtat) ~ "dichotomous",
                  c(assessment_age, timeSxService, caarms_DS_sev_0, caarms_UTC_sev_0, sans_tot_0, gaf_sofas) ~ "continuous",
                  c(uhr_cat) ~ "categorical"),
    value = list(gender ~ "Female",
                 transtat ~ "Transitioned",
                 group ~ "Non-standard intervention treatment"),
    label  = list(                                              # display labels for column names
      assessment_age   ~ "Age (years)",                           
      gender ~ "Gender",
      timeSxService  ~ "Time between first symptoms and acceptance/treatment at UHR/PACE service",    
      uhr_cat ~ "UHR syndrome",
      group ~ "Received non-standard intervention treatment as part of trial",
      caarms_DS_sev_0 ~ "CAARMS Disorganized Speech, Severity",
      caarms_UTC_sev_0 ~ "CAARMS Unusual Thought Content, Severity",
      gaf_sofas ~ "SOFAS score",
      sans_tot_0 ~ "SANS total score", 
      transtat ~ "Transition status"
    ),
    missing = "no", # don't list missing data separately;
  ) %>%
  add_n() %>% # add column with total number of non-missing observations
  modify_header(label = "Variable") %>% # update the column header
  bold_labels()

table_desc %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = "U:/PostDoc/Research/DynamicUpdating/descriptive_table.docx")

pairwise.wilcox.test(uhr1000.org$assessment_age, uhr1000.org$year2, p.adjust.method = "bonferroni")
