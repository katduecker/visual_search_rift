library(lme4)
library(rstatix)
library(rsq)
library(dplyr)
dir_behav = "Z:\\Visual Search RFT\\results\\eyelink"

# Compare ocular artefacts for fast vs slow trials


## Blinks ##########################

file_name = "blinks_fast_slow.csv"

file_path = file.path(dir_behav,file_name)
data <- read.csv(file_path)

data$id <- factor(data$id)

# decode condition into set size and guided/unguided

data$gui.ung[data$condi==1] <- 0
data$gui.ung[data$condi==2] <- 1
data$gui.ung[data$condi==3] <- 0
data$gui.ung[data$condi==4] <- 1

data$set.size[data$condi==1] <- 1
data$set.size[data$condi==3] <- 2
data$set.size[data$condi==2] <- 1
data$set.size[data$condi==4] <- 2

# models
bsl_model <- lmer(num_bl~ 1 +(1|id),data=data)
model_fast_slow<- lmer(num_bl~fast + (1|id),data=data)
anova(bsl_model,model_fast_slow)

model_setsize <- lmer(num_bl~set.size + (1|id),data=data)
model_setsize_gui <- lmer(num_bl~set.size * gui.ung + (1|id),data=data)

anova(bsl_model,model_setsize,model_setsize_gui)


## Saccades ##########################

file_name = "saccades_fast_slow.csv"

file_path = file.path(dir_behav,file_name)
data <- read.csv(file_path)

data$id <- factor(data$id)

# decode condition into set size and guided/unguided

data$gui.ung[data$condi==1] <- 0
data$gui.ung[data$condi==2] <- 1
data$gui.ung[data$condi==3] <- 0
data$gui.ung[data$condi==4] <- 1

data$set.size[data$condi==1] <- 1
data$set.size[data$condi==3] <- 2
data$set.size[data$condi==2] <- 1
data$set.size[data$condi==4] <- 2

# models
bsl_model <- lmer(num_sac~ 1 +(1|id),data=data)
model_fast_slow<- lmer(num_sac~fast + (1|id),data=data)

anova(bsl_model,model_fast_slow)


rsq(bsl_model)
rsq(model_fast_slow)

# t-test
pwc <- data %>%
  group_by(gui.ung,set.size) %>%
  pairwise_t_test(num_sac ~ fast, paired=TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
pwc

## Gaze bias #########################

file_name = "bias_fast_slow.csv"

file_path = file.path(dir_behav,file_name)
data <- read.csv(file_path)

data$id <- factor(data$id)

# decode condition into set size and guided/unguided

data$gui.ung[data$condi==1] <- 0
data$gui.ung[data$condi==2] <- 1
data$gui.ung[data$condi==3] <- 0
data$gui.ung[data$condi==4] <- 1

data$set.size[data$condi==1] <- 1
data$set.size[data$condi==3] <- 2
data$set.size[data$condi==2] <- 1
data$set.size[data$condi==4] <- 2


bsl_model <- lmer(bias~ 1 +(1|id),data=data)
model_fast_slow<- lmer(bias~fast + (1|id),data=data)

anova(bsl_model,model_fast_slow)
rsq(bsl_model)
rsq(model_fast_slow)