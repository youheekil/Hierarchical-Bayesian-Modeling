library(metafor)
library(dplyr)
library(ggplot2)
library(brms)

df <- escalc(measure="ZCOR", ri=ri, ni=ni, data=dat.molloy2014)
data<- df %>%
  mutate(study = paste0(authors, " (", year, ")"), sei = sqrt(vi)) %>%
select(study, yi, sei) %>% 
  slice(1:16)


par(mfrow=c(1,1))
ggplot(data, aes(x=yi, y=study)) +
  geom_segment(aes(x=yi-sei*2, xend = yi+sei*2, 
                   y = study, yend= study)) +
  geom_point()



## Two level model (standard random effects model)

mafixed<- rma( yi=yi, sei=sei,  data= data, method="FE")
print(mafixed, digits=3)
ma<- rma( yi=yi, sei=sei, slab=study, data= data)
print(ma, digits=3)
forest(ma)
ml.ma <- rma.mv(yi, sei, random = ~1 | study, data = data)
print(ml.ma, digits =3)

#BHMA model
set.seed(1004)
prior_c <- c(set_prior("normal(0, 1)", class = "Intercept"),
             set_prior("cauchy(0, 0.2)", class = "sd"))
brm1<- brm(
  yi|se(sei) ~ 1 + (1|study), 
  prior = prior_c, 
  data = data, 
  cores = 2,
  file = NULL)

summary(brm1)


#BHMA forest graph
a<-posterior_summary(brm1)
hist(a[,1])
# Study-specific effects are deviations + average
out_r <- spread_draws(brm1, r_study[study,term], b_Intercept) %>% 
  mutate(b_Intercept = r_study + b_Intercept) 
# Average effect
out_f <- spread_draws(brm1, b_Intercept) %>% 
  mutate(study = "Average")
# Combine average and study-specific effects' data frames
out_all <- bind_rows(out_r, out_f) %>% 
  ungroup() %>%
  # Ensure that Average effect is on the bottom of the forest plot
  mutate(study = fct_relevel(study, "Average"))
# Data frame of summary numbers
out_all_sum <- group_by(out_all, study) %>% 
  mean_qi(b_Intercept)
#> Warning: unnest() has a new interface. See ?unnest for details.
#> Try `cols = c(.lower, .upper)`, with `mutate()` needed
# Draw plot
out_all %>%   
  ggplot(aes(b_Intercept, study)) +
  geom_density_ridges(
    rel_min_height = 0.01, 
    col = NA,
    scale = 1
  ) +
  geom_pointintervalh(
    data = out_all_sum, size = 1
  ) +
  geom_text(
    data = mutate_if(out_all_sum, is.numeric, round, 2),
    # Use glue package to combine strings
    aes(label = glue::glue("{b_Intercept} [{.lower}, {.upper}]"), x = Inf),
    hjust = "inward"
  ) + 
  geom_vline(xintercept = 0, linetype = "longdash") +
  theme_classic()

#probabilty check 
brm1 %>%
  plot(
    combo = c("hist", "trace"), widths = c(1, 1.5),
    theme = theme_classic(base_size = 10))
# investigate model fit
mcmc_plot(brm1)


#application 2

#install.packages("MetaStan")
library("MetaStan")
# Loading required package: Rcpp
data("dat.Berkey1995", package = "MetaStan")
data<-(dat.Berkey1995)
library(ggplot2)
# Calculating log odds ratios and variances from data
logodds <- function(x) log((x[1] * (x[4] - x[3]))/((x[2] - x[1]) * x[3]))
stdes   <- function(x) sqrt(1/x[1] + 1/(x[2] - x[1]) + 1/x[3] + 1/(x[4] - x[3]))
r_ind   <- apply(cbind(dat.Berkey1995$rt, dat.Berkey1995$nt, 
                       dat.Berkey1995$rc, dat.Berkey1995$nc), 1, logodds)
se_ind  <- apply(cbind(dat.Berkey1995$rt, dat.Berkey1995$nt, 
                       dat.Berkey1995$rc, dat.Berkey1995$nc), 1, stdes)
lower95_ind <- r_ind + qnorm(.025) * se_ind
upper95_ind <- r_ind + qnorm(.975) * se_ind
# Comparison of the results
trials  <- c("1", "2" ,"3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13")
trials <- ordered(trials, levels = trials)
table(dat.Berkey1995)
d <- data.frame(x = trials,
                y = r_ind,
                sei = se_ind,
                ylo = lower95_ind,
                yhi = upper95_ind)
forest.plot <- ggplot(d, aes(x = x, y = y, ymin = ylo, ymax = yhi)) +
  geom_pointrange() +
  coord_flip() +
  geom_hline(aes(yintercept=0), lty = 2) +
  xlab("Studies") +
  ggtitle("Forest Plot (BCG vaccines)") +
  theme_classic()

plot(forest.plot)


library(rstan)
library(MetaStan)
library("shinystan")
library(shiny)

bnhm  <- meta_stan(ntrt = nt, 
                             nctrl = nc, 
                             rtrt = rt,
                             rctrl = rc,
                             data = dat.Berkey1995,
                             tau_prior_dist = "half-normal",
                             tau_prior = 0.5,
                             theta_prior = c(0, 2.82),
                             model = "BNHM1",
                             chains = 4,
                             iter = 2000,
                             warmup = 1000)
library("shiny")
bnhm.shinystan = as.shinystan(bnhm$fit)
launch_shinystan(bnhm.shinystan)
# print(bnhm.shinystan)


