library(lme4)
library(rstatix)
library(rsq)


# Search performance per condition (linear model)

dir_behav = "/rds/projects/j/jenseno-visual-search-rft/visual_search_rift/results/behavior"


# change name of .csv file to check analyse based on split during search

file_name = "fig1_performance_condi.csv"

file_path = file.path(dir_behav,file_name)
data <- read.csv(file_path)

data$id <- factor(data$id)

# recode condition into set size and guided/unguided

data$gui.ung[data$condi==1] <- 0
data$gui.ung[data$condi==2] <- 1
data$gui.ung[data$condi==3] <- 0
data$gui.ung[data$condi==4] <- 1

data$gui.ung[data$condi==5] <- 0
data$gui.ung[data$condi==6] <- 1
data$gui.ung[data$condi==7] <- 0
data$gui.ung[data$condi==8] <- 1


data$set.size[data$condi==1] <- 1
data$set.size[data$condi==3] <- 2
data$set.size[data$condi==2] <- 1
data$set.size[data$condi==4] <- 2

data$set.size[data$condi==5] <- 1
data$set.size[data$condi==7] <- 2
data$set.size[data$condi==6] <- 1
data$set.size[data$condi==8] <- 2


## RT #################################

# models
bsl_model <- glmer(RT~ 1 +(1|id),data=data, family=Gamma,control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                                                                optCtrl = list(method = "L-BFGS-B", starttests = FALSE, kkt = FALSE)))


model_setsize <- glmer(RT~1 + set.size + (1 +set.size|id),data=data, family=Gamma(identity),
                       control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                              optCtrl = list(method = "L-BFGS-B", starttests = FALSE, kkt = FALSE)))


model_setsize_gui <- glmer(RT~ 1+ set.size + gui.ung + (1 +set.size|id) + (1 +gui.ung|id),data=data, family=Gamma(identity),
                           control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                                  optCtrl = list(method = "L-BFGS-B", starttests = FALSE, kkt = FALSE)))

model_setsize_gui_int <- glmer(RT~ 1+ set.size * gui.ung + (1 +set.size|id) + (1 +gui.ung|id),data=data, family=Gamma(identity),
                           control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                                  optCtrl = list(method = "L-BFGS-B", starttests = FALSE, kkt = FALSE)))

anova(bsl_model,model_setsize,model_setsize_gui,model_setsize_gui_int)

rsq(bsl_model)

rsq(model_setsize)
rsq(model_setsize_gui)

rsq(model_setsize_gui_int)


summary(model_setsize_gui)


# DF for fixed effects
n <- nrow(model_setsize_gui@frame)
p <- length(coef(summary(model_setsize_gui)))
df <- n-p

pwc = pairwise_wilcox_test(data, RT~condi,p.adjust.method = "bonferroni",paired=TRUE)
pwc
# SE of difference
diff_se <- matrix(NA,4,4)

for (i in 1:4) {
  for (j in 1:4){
    if (i != j){
      subdata1 <- data$RT[data$condi == i]
      subdata2 <- data$RT[data$condi == j]
      
      diff_se[i,j] <- sd(subdata1 - subdata2)/sqrt(31)
      
    }
  }
}


ung16 = data$RT[data$condi == 1]
gui32 = data$RT[data$condi == 4]

pwc= wilcox.test(ung16,gui32,paired=TRUE)

wilcox_effsize(ung16,gui32,paired=TRUE)

eff = wilcox_effsize(data,RT~condi, paired = TRUE, hedges.correction = TRUE)



Z = eff$effsize*sqrt(31)



## d' ################################

shapiro_test(data$d)


# models
bsl_model <- glmer(d~1 +(1|id),data=data)
model_setsize <- glmer(d~1 + set.size + (1 +set.size|id),data=data)
model_setsize_gui <- glmer(d~1+ set.size + gui.ung + (1 +set.size|id) + (1 +gui.ung|id),data=data)
model_setsize_gui_int <- glmer(d~1+ set.size * gui.ung + (1 +set.size|id) + (1 +gui.ung|id),data=data)


anova(bsl_model,model_setsize,model_setsize_gui,model_setsize_gui_int)

summary(model_setsize_gui)

rsq(bsl_model)
rsq(model_setsize)
rsq(model_setsize_gui)
rsq(model_setsize_gui_int)

pwc = pairwise_t_test(data, d~condi,p.adjust.method = "bonferroni",paired=TRUE,detailed=TRUE)
pwc
pwc_cohd = cohens_d(data, d~condi,hedges.correction = TRUE)
pwc_cohd



# SE of difference
diff_se <- matrix(NA,4,4)

for (i in 1:4) {
  for (j in 1:4){
    if (i != j){
      subdata1 <- data$d[data$condi == i]
      subdata2 <- data$d[data$condi == j]
      
      diff_se[i,j] <- sd(subdata1 - subdata2)/sqrt(31)
      
    }
  }
}diff_se
