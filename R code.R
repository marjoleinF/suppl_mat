###################################################################
###
###
### R code for reproducing analyses from:
###
###   One Model May Not Fit All: Subgroup Detection Using
###   Model-Based Recursive Partitioning
###




#########################
##
## Introduction
##
##

set.seed(3)
x1 <- runif(250, min = 5, max = 11)
x2 <- runif(250, min = 2, max = 8)
x3 <- runif(250, min = 4, max = 10)
error <- rnorm(250, sd = 1.5)
y <- 7 + ifelse(x1 > 8 & x2 > 5, yes = 2.5, no = -2.5) + error
toy_data <- round(data.frame(x1, x2, x3, y), digits = 2L)
toy_data$x2 <- factor(ifelse(x2 > 5, yes = ifelse(x2 > 6.5, "A", "C"), no = "B"))
boxplot(toy_data$y, cex = .7)
library("partykit")
tree <- lmtree(y ~ 1 | x1 + x2 + x3, data = toy_data)
plot(tree, gp = gpar(cex = .7), ylim = c(1, 14))



######################################################################
##
## Application Example 1: Subgroup Detection in Mixed-Effects Models 
##
##

## Load data
load("HS_dat.Rda")
head(HS_dat, 3)

## Fit traditional GLMM
library("lme4")
lmm <- lmer(PPVT ~ Program*Age + (1|MotherID/ChildID), data = HS_dat)
beta <- fixef(lmm)
plot(jitter(HS_dat$Age), HS_dat$PPVT, col = HS_dat$Program,
     cex = .5,
     xlab = "Age (years)", ylab = "PPVT score", xaxt = "n",
     main = paste0("Global Model (", length(unique(HS_dat$ChildID)), " children)"))
abline(a = beta["(Intercept)"], b = beta["Age"], col = "black")
abline(a = beta["(Intercept)"]+beta["ProgramHS"], 
       b = beta["Age"]+beta["ProgramHS:Age"], col = "red")
axis(1, cex.axis = .8)

## Fit GLMM tree
library("glmertree")
HS_tree <- lmertree(PPVT ~ Program*Age | (1|MotherID/ChildID) | AFTQ + Race + 
                      Income + Mom_edu_yrs + Mom_height, 
                    data = HS_dat, cluster = MotherID, minsize = 250)
plot(HS_tree, type = "simple", which = "tree", gp = gpar(cex = .5), 
     nodesize_level = 2)

## Extract parameter estimates and predictions (not shown in main paper)
coef(HS_tree)
fixef(HS_tree)
VarCorr(HS_tree)
ranef(HS_tree)
predict(HS_tree, newdata = HS_dat[1:10, ]) ## random effects included
predict(HS_tree, newdata = HS_dat[1:10, ], re.form = NA) ## ranef excluded

## Create node-specific plots
nodes <- predict(HS_tree, type = "node")
beta <- coef(HS_tree)
par(mfrow = c(1, 4))
for (node in sort(unique(nodes))) {
  plot(HS_dat$Age[nodes == node], HS_dat$PPVT[nodes == node],
       col = HS_dat$Program[nodes == node],
       cex = .5, cex.lab = .7, cex.axis=.7, cex.main = .7,
       xlab = "Age", ylab = "PPVT", #xlim = c(sqrt(3), sqrt(20)), 
       ylim = c(0, 150))
  abline(a = beta[as.character(node), "(Intercept)"], b = beta[as.character(node), "Age"], col = "black")
  abline(a = beta[as.character(node), "(Intercept)"] + beta[as.character(node), "ProgramHS"], 
         b = beta[as.character(node), "Age"] + beta[as.character(node), "ProgramHS:Age"], col = "red")
}

## Create node-specific predictions
coefs <- round(coef(HS_tree), digits = 2)
design_df <- data.frame(.tree = rep(sort(unique(HS_tree$data$.tree)), each = 2),
                        Program = factor(rep(c("None", "HS"))))
design_df$`PPVT at age 6` <- predict(HS_tree$lmer, newdata = cbind(design_df, Age = 6), 
                                     re.form = ~0)
design_df$`PPVT at age 12` <- predict(HS_tree$lmer, newdata = cbind(design_df, Age = 12), 
                                      re.form = ~0)
names(design_df)[1] <- "Node"
knitr::kable(design_df, format = "latex", digits = 2, 
             caption = "Node-specific average predicted PPVT scores at different ages.", 
             label = "predictions", centering = FALSE, linesep = "",
             row.names = FALSE, escape=TRUE, align = c("cccc"), booktabs = TRUE)






######################################################################
##
## Application Example 2: Subgroup Detection in Rasch models 
##
##

## Load data
load("dat_SPISA.Rda")
head(dat_SPISA, 3)
levels(dat_SPISA$Area) <- c("Lang & Cult", "\nLaw and Econ",
                            "\nMed & Health", "\nEngin",
                            "\nSci, Phar, Geo",
                            "\nAgri & Nutri",
                            "\nSports", "Arts",
                            "\nno Info")

## Fit Rasch tree
library("psychotree")
Raschtree_culture <- raschtree(culture ~  Gender + Age + Area,
                               data = dat_SPISA, maxdepth = 4)
plot(Raschtree_culture)

## Install raschtreeMH
library("devtools")
install_github("mirka-henninger/raschtreeMH")
library("raschtreeMH")

## Fit Raschtree using MH effect-size stopping criterion
Raschtree_MH_culture <- raschtree(culture ~  Gender + Age + Area, 
                                  data = dat_SPISA,
                                  stopfun= stopfun_mantelhaenszel(purification = "iterative"))
Raschtree_MH_culture <- add_mantelhaenszel(Raschtree_MH_culture,
                                           purification = "iterative")
plot(Raschtree_MH_culture)
plot(Raschtree_MH_culture, color_by_node = 1)
plot(Raschtree_MH_culture, color_by_node = 4)
Raschtree_MH_culture$info$mantelhaenszel