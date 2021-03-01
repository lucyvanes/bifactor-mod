# bifactor modelling on ePrime 4-year outcome data
# June 2020
# Lucy Vanes

library(nFactors)
library(psych)
library(GPArotation)
library(semTools)

# call functions
source("C:/Users/k1456336/Documents/KCL BMEIS/ePrime/scripts/EBFA_CBFA_functions.R")

# =======a======== = Import Data = ===============
setwd("C:/Users/k1456336/Documents/KCL BMEIS/ePrime/data")
dat <- read.csv("ePrime_merged_data2.csv", header=T)
dat$id <- factor(dat$id)
dat$sex <- factor(dat$sex)



composite_psych_vars <- c("ecbq_neg_affect", "ecbq_surgency",  "ecbq_effort_control","emque_emo_contagion","emque_attent_others","emque_prosocial", "sdq_emo", "sdq_conduct", "sdq_adhd", "sdq_peer", "sdq_prosocial", "adhd_inattention", "adhd_hyper", "srs_soc_awareness", "srs_soc_cog", "srs_soc_comm", "srs_soc_mot", "srs_rrb")

composite_psych_vars_excl <- c("ecbq_neg_affect", "emque_emo_contagion","emque_attent_others","emque_prosocial", "sdq_emo", "sdq_conduct", "sdq_adhd", "sdq_peer", "sdq_prosocial", "adhd_inattention", "adhd_hyper", "srs_soc_awareness", "srs_soc_cog", "srs_soc_comm", "srs_soc_mot", "srs_rrb")

item_vars <- names(dat)[25:188]

# decide which variables to include 
which_psych <- "composite" # "composite", "composite_excl", or "item"

if (which_psych=="composite"){
	bi_dat <- dat[,c(composite_psych_vars)]
	bi_dat_ids <- dat[,c("id",composite_psych_vars)]
	bi_dat <- bi_dat[complete.cases(bi_dat),]
	bi_dat_ids <- bi_dat_ids[complete.cases(bi_dat_ids),]
	bi_dat_cor <- cor(bi_dat)	
} else if (which_psych=="composite_excl"){
	bi_dat <- dat[,c(composite_psych_vars_excl)]
	bi_dat_ids <- dat[,c("id",composite_psych_vars_excl)]
	bi_dat <- bi_dat[complete.cases(bi_dat),]
	bi_dat_ids <- bi_dat_ids[complete.cases(bi_dat_ids),]
	bi_dat_cor <- cor(bi_dat)
} else if (which_psych=="item"){
	bi_dat <- dat[,c(item_vars)]
	bi_dat_ids <- dat[,c("id",item_vars)]
	bi_dat <- bi_dat[complete.cases(bi_dat),]
	bi_dat_ids <- bi_dat_ids[complete.cases(bi_dat_ids),]
	bi_dat_cor <- cor(bi_dat)
}

# General procedure:
#=======================================
# apply exploratory bifactor analysis to get an idea of factor structure
# then run confirmatory bifactor analysis using this factor structure
# then test things like measurement invariance
# since measurement invariance across sex was not met here, I repeated
# those steps after excluding two variables that seemed 
# responsible (ecbq surgency and ecbq effortful control)


#==================================================================
# 				Exploratory bifactor analysis
#==================================================================

# exploratory bifactor analysis (EBFA) largely following this:
# https://www.tandfonline.com/doi/full/10.1080/10705511.2019.1622421?casa_token=CB6g7O-U4wkAAAAA%3AqkJZVL7Bkf-rS859WFyzEjPzOvTFHd8vEMyzAKlcz2VS_KHzaicdvdpQmO_YvhGZe49Fb7aaEDSm

# and corresponding example syntax here:
# https://osf.io/ne795/

# ===================== = Number of Factors = =====================
# determine number of factors
ev <- eigen(bi_dat_cor)
ap <- parallel(subject=nrow(bi_dat),var=ncol(bi_dat), rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS)


# MAP and BIC
n <- length(bi_dat[,1])
nfactors(bi_dat_cor, n.obs = n)


# parallel analysis--principal axis
dat.pa.mean <- fa.parallel(
  bi_dat_cor,
  n.obs = n,
  fa = "both",
  fm = "pa")

dat.pa.95 <- fa.parallel(
  bi_dat_cor,
  n.obs = n,
  fa = "both",
  fm = "pa",
  quan = 0.95)



# ===================== = Factor Extraction = =====================

# looks like 3-5 factors might be good; played around with all options. If 3 factors seems good, need to extract one extra factor which will later become the general factor (after rotation), so total factors = 4

n_factors <- 4


# extraction of 4 factors (via principal axis) with different communality
# starting values.  (see link above)
# here starting values correspond to the SMCs of the data

f.pa.smc <- fa(
  bi_dat_cor,
  nfactors = n_factors,
  n.obs = n,
  rotate = "none",
  fm = "pa",
  max.iter = 100,
  SMC = smc(bi_dat_cor))

# here starting values range from 0.3 to 0.9, in increments of 0.1. 
# Will compare communalities later and see which starting values did the best job.

starting_values <- seq(0.3, 0.9, 0.1)
communalities <- data.frame(var=dimnames(bi_dat_cor)[[1]])
communalities[,1:length(starting_values)+1] <- NA
names(communalities) <- c("var",starting_values)

eval <- data.frame(starting_value=starting_values, fit=rep(NA,length(starting_values)), objective=rep(NA,length(starting_values)))
	
for (s in 1:length(starting_values)){
	 f.pa.s <- fa(
				   bi_dat_cor, 
				   nfactors = n_factors,
				   n.obs = n,
				   rotate = "none",
				   fm = "pa",
				   max.iter = 100,
				   SMC = rep(starting_values[s], length(bi_dat_cor[1,]))  )
	 communalities[s+1] <- f.pa.s$communality
	 eval$fit[s] <- f.pa.s$fit
	 eval$objective[s] <- f.pa.s$objective 
  }

View(communalities)	# see communalities for each starting value
View(eval)			# see fit and objective function

f.pa.smc$communality	# see communalities when SMC used as starting value
f.pa.smc$fit			# see fit
f.pa.smc$objective		# see objective function

# communalities look consistent across all starting values; i.e. doesn't seem to make a big difference.
# I will proceed with the option using SMC since this looks nearly identical.

# extract factor loadings
f.pa.smc.load <- f.pa.smc$loadings


# ============ = Rotation = ============

# here we apply the bifactor rotation. There are two options, geomin and quartimin, which I will compare

# bi-quartimin, 1000 random starts, return up to the 10 best solutions
# this calles the FindBifactor.orth function which was called at the top of the script

bifac_quartimin_rot <- FindBifactor.orth(f.pa.smc.load, reps = 1000, rotation = "bifactor", solutions = 10)

# bi-geomin, 1000 random starts, return up to the 10 best solutions

geomin_rot <- FindBifactor.orth(f.pa.smc.load, reps = 1000, rotation = "geomin", solutions = 10)

# display results
bifac_quartimin_rot
geomin_rot

# geomin gives 2 solutions, quartimin 1. Will compare first geomin solution and quartimin solution

loadings_quart <- as.data.frame(bifac_quartimin_rot$loadings[[1]])
loadings_geo <- as.data.frame(geomin_rot$loadings[[1]])

# write.csv(loadings_geo, "loadings_f4_composite_geomin.csv", row.names=T, quote=F)
# write.csv(loadings_quart, "loadings_f4_composite_quartimin.csv", row.names=T, quote=F)

# I write out the loadings into csv files so I can highlight loadings greater .3(ish) (both pos and neg) to see
# what the factors represent. In the first attempt I continue with the geomin solution because
# it makes most sense to me conceptually. Later when I have removed problematic variables I end up using the 
# quartimin solution

#==================================================================
# 				Confirmatory bifactor analysis
#==================================================================
# Inspiration for confirmatory analysis taken from here:
# https://maelliott1010.github.io/Confirmatory_Factor_Analysis_ADHD/

# here I specify the confirmatory model, using the factors I identified in the previous
# step and visualised for myself in the csv file. 
# the adhd factor seems to reoccur quite robustly with lots of versions, the remaining ones are more mysterious.

# 4 factor (1 p 4 specific), geomin rotation
bimod_4_geomin <- " p =~ ecbq_neg_affect + ecbq_surgency + ecbq_effort_control + emque_emo_contagion + emque_attent_others + emque_prosocial + sdq_emo + sdq_conduct + sdq_adhd + sdq_peer + sdq_prosocial + adhd_inattention + adhd_hyper + srs_soc_awareness + srs_soc_cog + srs_soc_comm + srs_soc_mot + srs_rrb

		s1 =~ emque_emo_contagion + emque_attent_others + sdq_emo + sdq_peer	# emotion sensitivity?
		s2 =~ sdq_conduct + sdq_adhd + adhd_inattention + adhd_hyper			# adhd
		s3 =~ ecbq_neg_affect + ecbq_surgency + ecbq_effort_control + emque_prosocial	# ecbq mystery

		p ~~ 0*s1
		p ~~ 0*s2
		p ~~ 0*s3
		s1 ~~ 0*s2
		s1 ~~ 0*s3
		s2 ~~ 0*s3
"


fit_bi_4_geomin <- cfa(bimod_4_geomin, data=bi_dat)

summary(fit_bi_4_geomin, fit.measures=T,rsquare=TRUE) # TLI=0.869; CFI=.895; RMSEA=.093; SRMR=.068
#  value close to or greater than .95 for the TLI and the CFI, and a value close to or less than .06 for the RMSEA (Hu & Bentler, 1998) and SRMR less than .08 (Kline, 2016) indicate a good fit between the model and the observed data

reliability(fit_bi_4_geomin) 
# Recommendations by Reise et al. suggest a minimual acceptable value of omega-h of .50, though a value of .75 is more acceptable
# Comparison of omega and omega-h values can indicate how much reliable variance could be attributed to general vs. specific factors.

HancockMueller(fit_bi_4_geomin)
# High H values (>.70) indicate the factor is stable and is less likely to fluctuate between various samples.


inspect(fit_bi_4_geomin,what="std")$lambda


# measurement invariance
#==========================
# very useful links:
# https://rstudio-pubs-static.s3.amazonaws.com/194879_192b64ad567743d392b559d650b95a3b.html
# https://users.ugent.be/~yrosseel/lavaan/multiplegroup6Dec2012.pdf

# this is to test whether the factor structure is invariant across different subsamples; 
# e.g. is the factor structure the same in boys and girls? Only if this is the case can we
# use the model across both boys and girls and subsequently compare factor scores between these groups
# weak invariance: set loadings to be equal across both groups
# strong invariance: set loadings AND intercept to be equal across both groups

# first I test for invariance across two random split halves of the sample; don't know if this is 
# a sensible thing to do but I'm doing it anyway:
# create random subsamples (50/50)
samp_rows <- sample(nrow(bi_dat), round(length(bi_dat[,1])/2))
bi_dat$split_half <- 0
bi_dat$split_half[samp_rows] <- 1

# random split half
#=====================
configural <- cfa(bimod_4_geomin, data=bi_dat,estimator="MLR", std.lv=TRUE, missing="FIML", group="split_half")
weak.invariance <- cfa(bimod_4_geomin, data=bi_dat,estimator="MLR", std.lv=TRUE, missing="FIML", group="split_half", group.equal=c("loadings"))
strong.invariance <- cfa(bimod_4_geomin, data=bi_dat, group="split_half", group.equal=c("loadings", "intercepts"), missing="FIML", estimator="MLR", std.lv=TRUE)
lavTestLRT(configural,weak.invariance) # test weak invariance / metric invariance
lavTestLRT(weak.invariance,strong.invariance) # test strong invariance / scalar invariance

# then I test for invariance across boys and girls
# sex
#=================================================
bi_dat_ids$sex <- NA
for (i in levels(bi_dat_ids$id)){
	bi_dat_ids$sex[bi_dat_ids$id==i] <- dat$sex[dat$id==i]
}
# configural invariance
configural <- cfa(bimod_4_geomin, data=bi_dat_ids,estimator="MLR", group="sex")
weak.invariance <- cfa(bimod_4_geomin, data=bi_dat_ids, group="sex", group.equal=c("loadings"))
strong.invariance <- cfa(bimod_4_geomin, data=bi_dat_ids, group="sex", group.equal=c("loadings", "intercepts"))

lavTestLRT(configural,weak.invariance)		# test weak/metric invariance
lavTestLRT(weak.invariance,strong.invariance) # does not pass scalar invariance for sex

lavTestScore(strong.invariance)
# look for greatest X2 value - here p61==p141, X2=8.724
# also .p60. == .p140. 6.879 

parTable(strong.invariance)
# find it in table: p61==p141 refers to ecbq_effort_control
# 					p60==p140 refers to ecbq_surgency

# test partial invariance
strong.invariance.effort_control <- cfa(bimod_4_geomin, data=bi_dat_ids, group = "sex", group.equal = c("loadings", "intercepts"), group.partial = c("ecbq_effort_control ~ 1"))
lavTestScore(strong.invariance.effort_control)

strong.invariance.effort_control_surgency <- cfa(bimod_4_geomin, data=bi_dat_ids, group = "sex", group.equal = c("loadings", "intercepts"), group.partial = c("ecbq_effort_control ~ 1", "ecbq_surgency ~ 1"))
lavTestScore(strong.invariance.effort_control_surgency)

anova(strong.invariance.effort_control_surgency, weak.invariance)


# Rerun, excluding ECBQ surgency and ECBQ effortful control
#==========================================================
# go back up and repeat EBFA steps, using which_psych <- "composite_excl" this time.
# Once again write out the geomin and quartimin loadings as csv files and highlight
# loadings > .30 or .25 to get an idea of the structure
# here I go with the quartimin solution. Geomin has a bunch of variables (mostly SRS) that 
# do not load on to a specific factor, suggesting that the p-factor has a strong SRS component
# Quartimin has specific factors that cover all variables.
# have here also included variables with loadings closer to .2 on specific factors

bimod_4_quartimin <- " p =~ ecbq_neg_affect +  emque_emo_contagion + emque_attent_others + emque_prosocial + sdq_emo + sdq_conduct + sdq_adhd + sdq_peer + sdq_prosocial + adhd_inattention + adhd_hyper + srs_soc_awareness + srs_soc_cog + srs_soc_comm + srs_soc_mot + srs_rrb

		s1 =~ sdq_conduct + sdq_adhd + adhd_inattention + adhd_hyper + srs_soc_mot			# adhd
		s2 =~ ecbq_neg_affect + emque_emo_contagion +emque_attent_others + emque_prosocial	+ sdq_emo + sdq_peer + sdq_prosocial
		s3 =~  srs_soc_awareness + srs_soc_cog + srs_soc_comm + srs_rrb
		

		p ~~ 0*s1
		p ~~ 0*s2
		p ~~ 0*s3
		s1 ~~ 0*s2
		s1 ~~ 0*s3
		s2 ~~ 0*s3
"
# Rerun the above invariance tests across sex using this model
# configural model does not converge, presumably due to sample size in males and females separately (?). However, strong invariance holds (is this valid?). 

fit_bi_4 <- cfa(bimod_4_quartimin, data=bi_dat)

summary(fit_bi_4, fit.measures=T,rsquare=TRUE) # TLI=0.869; CFI=.895; RMSEA=.093; SRMR=.068
reliability(fit_bi_4) 
HancockMueller(fit_bi_4)
inspect(fit_bi_4,what="std")$lambda



# Extract bifactor scores
#=====================================
bifac_scores_4fac <- as.data.frame(predict(fit_bi_4))
bi_dat$bifac4_p <- bifac_scores_4fac$p
bi_dat$bifac4_s1 <- bifac_scores_4fac$s1
bi_dat$bifac4_s2 <- bifac_scores_4fac$s2
bi_dat$bifac4_s3 <- bifac_scores_4fac$s3

dat1 <- merge(dat, bi_dat, all.x=T, all.y=T)

# s1 = ADHD factor
# s2 = social/emotional factor
# s3 = ASD factor




