
#COLLIDER BIAS in MAIHDA analyses demonstration

#code adapted from https://mixtape.scunning.com/03-directed_acyclical_graphs by Dr Andrew Bell and Dr Gareth J. Griffith

library(tidyverse)
library(lme4)
library(sjPlot)
library(merTools)
library(ggpubr)

set.seed(33)

tb <- tibble(
  minoritised = ifelse(runif(100000)>=0.8,1,0),
  helpseek = ifelse(rnorm(100000)>0,1,0),
  education = -1*minoritised + 2*helpseek + rnorm(100000),
  degree = ifelse(education>0,1,0),
  sex = ifelse(runif(100000)>=0.5,1,0),
  age = runif(100000, min=18, max=80),
  agecat = cut(age, breaks=c(0,40,70,100))
) %>%
  mutate(
    
    health = 6 + 2*as.numeric(helpseek) + 1*degree - 1*as.numeric(agecat) + 1*sex + rnorm(100000, 0, 3),
    
    # Numeric IDs for modeling
    strata_id     = group_indices(., as.numeric(agecat), sex, degree, minoritised),
    stratanoed_id = group_indices(., as.numeric(agecat), sex, minoritised),
    stratahelp_id = group_indices(., helpseek, as.numeric(agecat), sex, degree, minoritised),
    
    # Labels
    minoritised_label = ifelse(minoritised == 1, "Minority", "Majority"),
    degree_label      = ifelse(degree == 1, "Degree", "NoDegree"),
    sex_label         = ifelse(sex == 1, "Male", "Female"),
    agecat_label = case_when(
      as.numeric(agecat) == 1 ~ "Under 41",
      as.numeric(agecat) == 2 ~ "41-70",
      as.numeric(agecat) == 3 ~ "Over 71"
    ),
    helpseek_label = ifelse(helpseek == 1, "HelpSeeking", "NoHelpSeek"),
    
    # NEW: stratum label WITHOUT age category or minoritised label (for graphics)
    strata_label_nage = paste(sex_label, degree_label, sep = "_"),
    
    # Preserve labels in case needed elsewhere
    strata_label     = paste(agecat_label, sex_label, degree_label, minoritised_label, sep = "_"),
    stratanoed_label = paste(agecat_label , sex_label, minoritised_label, sep = "_"),
    stratahelp_label = paste(helpseek_label, agecat_label, sex_label, degree_label, minoritised_label, sep = "_")
  )


# Lookup includes agecat_label AND the new label without age/minority status
strata_lookup <- tb %>%
  distinct(strata_id, agecat_label, strata_label_nage)


#########################
### Figure 1 analyses ###
#########################

# Fit models (typical MAIHDA models 1 & 2, with and without help-seeking 
# stratification) to show how the inclusion of the collider (help-seeking) 
# changes the results - and also show how the results change if we stratify 
# excluding education

MAIHDA1     <- lmer(health ~ (1|strata_id), tb)
MAIHDA2     <- lmer(health ~ minoritised  + sex + agecat + degree + (1|strata_id), tb)
MAIHDA1NOED <- lmer(health ~ (1|stratanoed_id), tb)
MAIHDA2NOED <- lmer(health ~ minoritised + sex + agecat + (1|stratanoed_id), tb)
MAIHDA1whelp <- lmer(health ~ (1|stratahelp_id), tb)
MAIHDA2whelp <- lmer(health ~ minoritised + sex + agecat + degree + helpseek + (1|stratahelp_id), tb)

tab_model(MAIHDA1, MAIHDA2, MAIHDA1NOED, MAIHDA2NOED, MAIHDA1whelp, MAIHDA2whelp,
          p.style="stars", show.se=TRUE, show.ci=FALSE)

# Predictions from models 1 and 2
m1m <- predictInterval(MAIHDA1, level=0.95, include.resid.var=FALSE) %>% mutate(id=row_number())
m2m <- predictInterval(MAIHDA2, level=0.95, include.resid.var=FALSE) %>% mutate(id=row_number())

# Merge predictions
tb$id <- seq.int(nrow(tb))
tb2 <- tb %>%
  left_join(m1m, by="id") %>%
  rename(m1mfit=fit, m1mupr=upr, m1mlwr=lwr) %>%
  left_join(m2m, by="id") %>%
  rename(m2mfit=fit, m2mupr=upr, m2mlwr=lwr)

# Stratum-level summary for plotting
stratum_level <- aggregate(
  x = tb2[c("health")],
  by = tb2[c("sex","minoritised","degree","agecat","strata_id",
             "m1mfit","m1mupr","m1mlwr","m2mfit","m2mupr","m2mlwr")],
  FUN = mean
)

# Ensure minority is a factor so ggplot2 recognises it as a variable to "color" by
stratum_level$minoritised <- as.factor(stratum_level$minoritised)


# Bring in generated agecat_label and map new label WITHOUT age/minority for plotting
stratum_level <- stratum_level %>%
  left_join(strata_lookup, by = "strata_id") %>%
  arrange(strata_id)

# Gen interpretable labels for graphics
stratum_level <- stratum_level %>%
  mutate(agecat_label = fct_relevel(agecat_label, "Under 41", "41-70", "Over 71"))

# Build a named label map for x-axis ticks
lab_map <- stratum_level %>%
  distinct(strata_id, strata_label_nage) %>%
  { setNames(.$strata_label_nage, .$strata_id) }




# Pivot the m1m* and m2m* columns into long format for ggplot2 faceting and make
# model names interpretable factors
stratum_long <- stratum_level %>%
  pivot_longer(
    cols = matches("^m[12]m(fit|upr|lwr)$"),  
    names_to = c("model", ".value"),          
    names_pattern = "^(m[12]m)(fit|upr|lwr)$" 
  )
stratum_long$model <- fct_recode(stratum_long$model,"Model 1A" = "m1m", "Model 2A" = "m2m")

dodge <- position_dodge(width = 0.6)

### Plot: Model 1 (faceting by age cat and model no.)
pred1 <- ggplot(
  stratum_long,
  aes(y = fit, x = as.factor(strata_id), color = minoritised, group = strata_label_nage)) +
  geom_point(size=1) +
  geom_errorbar(aes(ymax = upr, ymin = lwr),
                position = dodge, width = 0.3)+
  ylab("Predicted Health Score") +
  xlab("") + ylim(2,10) +
  theme_bw() +
  scale_x_discrete(labels = lab_map) +
  facet_grid(vars(model), vars(agecat_label), scales = "free_x") +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
# For checking
pred1

ggsave("maihda_fig3.png", pred1, dpi = 500, width = 6, height = 6)

#############################
### End Figure 1 analyses ###
#############################


#################################
### MNAR Example for Figure 2 ###
#################################
set.seed(1)

tb <- tibble(
  minoritised = ifelse(runif(100000)>=0.8,1,0), #minoritised ethnic, binary variable
  degree = ifelse(runif(100000)>=0.6,1,0), #binary version of education
  sex = ifelse(runif(100000)>=0.5,1,0), #sex is a binary variable
  age = runif(100000, min=18, max=80), #age is a uniform variable...
  agecat = cut(age, breaks=c(0,40,70,100)), #...which we cut into 3 categories
  health = 6 + 1*degree - 1*as.numeric(agecat) + 1*sex + rnorm(100000,0,3), #outcome variable generated as function of age, sex and education (unobserved latent value underpinning degree)
  strata = minoritised + 10*degree + 100*sex + 1000*as.numeric(agecat),
  unhealthy = 2 - as.numeric(cut(health, breaks=c(-10,4,100))), #create a binary measure of unhealthiness
  random = rbinom(n=100000, size=1, prob=0.7), #create a random 0/1 identifier to help define missingness
  keepMAR =  1-(random*sex*minoritised), #this will define missingness that is MAR - ie remove 70% of all non-degree having men
  keepMNAR =  1-(random*unhealthy*sex*minoritised)) %>%   #this will define MNAR missingness- ie remove 70% of all unhealthy non-degree having men
  
  # Numeric IDs for modeling
  mutate(strata_id = group_indices(., as.numeric(agecat), sex, degree, minoritised),
         stratanoed_id = group_indices(., as.numeric(agecat), sex, minoritised),
         stratahealth_id = group_indices(., unhealthy, as.numeric(agecat), sex, degree, minoritised),
         
         # Labels
         minoritised_label = ifelse(minoritised == 1, "Minority", "Majority"),
         degree_label      = ifelse(degree == 1, "Degree", "NoDegree"),
         sex_label         = ifelse(sex == 1, "Male", "Female"),
         unhealth_label    = ifelse(unhealthy ==1, "Unhealthy", "Healthy"),
         agecat_label = case_when(
           as.numeric(agecat) == 1 ~ "Under 41",
           as.numeric(agecat) == 2 ~ "41-70",
           as.numeric(agecat) == 3 ~ "Over 71"
         ),
         
         # stratum label WITHOUT age category or minoritised label (to use on x-axis)
         strata_label_nage = paste(sex_label, degree_label, sep = "_"),
         
         # (Keep the others if needed elsewhere)
         strata_label     = paste(agecat_label, sex_label, degree_label, minoritised_label, sep = "_"),
         stratanoed_label = paste(agecat_label , sex_label, minoritised_label, sep = "_"),
         strata_unhealth_label = paste(agecat_label, sex_label, unhealth_label, sep="_")
  )

tb <- tb[order(tb$strata),]

MAIHDA1full <- lmer(health ~ (1|strata), tb)
MAIHDA2full <- lmer(health ~ minoritised  + sex + agecat + degree + (1|strata), tb)
summary(MAIHDA2full)

#MAR: cut 70% of the minoritised female individuals (MAR with interaction)
#sense check bias on re effects?

tbMAR <- subset(tb, keepMAR==1) 

MAIHDA1mar <- lmer(health ~ (1|strata), tbMAR)
MAIHDA2mar <- lmer(health ~ minoritised  + sex + agecat + degree + (1|strata), tbMAR)
summary(MAIHDA2mar)

#MNAR: cut 70% of the *least healthy* (ie health score below X) minoritised female individuals 
#should be biased.

tbMNAR <- subset(tb, keepMNAR==1) 

MAIHDA1mnar <- lmer(health ~ (1|strata), tbMNAR)
MAIHDA2mnar <- lmer(health ~ minoritised  + sex + agecat + degree + (1|strata), tbMNAR)
summary(MAIHDA2mnar)

tab_model(MAIHDA1full, MAIHDA2full, MAIHDA1mar, MAIHDA2mar, MAIHDA1mnar, MAIHDA2mnar, 
          p.style="stars", show.se=T, show.ci=F)


#make predictions from MAIHDA2

m2mfull <- predictInterval(MAIHDA2full, level=0.95, include.resid.var=FALSE)
m2mfull <- mutate(m2mfull, id=row_number())

m2mMAR <- predictInterval(MAIHDA2mar, level=0.95, include.resid.var=FALSE)
m2mMAR <- mutate(m2mMAR, id=row_number())

m2mMNAR <- predictInterval(MAIHDA2mnar, level=0.95, include.resid.var=FALSE)
m2mMNAR <- mutate(m2mMNAR, id=row_number())

#merge predictions into main dataframe
tb$id <- seq.int(nrow(tb))
tb2full <- merge(tb, m2mfull, by="id")
tb2full <- tb2full %>%
  rename(
    m2mfit=fit,
    m2mupr= upr,
    m2mlwr=lwr
  )

tbMAR$id <- seq.int(nrow(tbMAR))
tb2mar <- merge(tbMAR, m2mMAR, by="id")
tb2mar <- tb2mar %>%
  rename(
    m2mfit=fit,
    m2mupr= upr,
    m2mlwr=lwr
  )

tbMNAR$id <- seq.int(nrow(tbMNAR))
tb2mnar <- merge(tbMNAR, m2mMNAR, by="id")
tb2mnar <- tb2mnar %>%
  rename(
    m2mfit=fit,
    m2mupr= upr,
    m2mlwr=lwr
  )



#create stratum level dataframe
stratum_levelfull <- aggregate(x=tb2full[c("health")], 
                               by=tb2full[c("sex", "minoritised", "degree", "agecat", "strata_id",
                                            "m2mfit", "m2mupr", "m2mlwr","strata_label_nage" 
                               )],
                               FUN=mean)

stratum_levelmar <- aggregate(x=tb2mar[c("health")], 
                              by=tb2mar[c("sex", "minoritised", "degree", "agecat", "strata_id",
                                          "m2mfit", "m2mupr", "m2mlwr", "strata_label_nage"
                              )],
                              FUN=mean)

stratum_levelmnar <- aggregate(x=tb2mnar[c("health")], 
                               by=tb2mnar[c("sex", "minoritised", "degree", "agecat", "strata_id",
                                            "m2mfit", "m2mupr", "m2mlwr", "strata_label_nage" 
                               )],
                               FUN=mean)


# Build the lookup again before joins, preserve needed columns
strata_lookup <- tb %>%
  distinct(strata_id, agecat_label, strata_label_nage)

# Join lookup into each strata df
stratum_levelfull <- stratum_levelfull %>%
  left_join(strata_lookup, by = "strata_id") %>%
  mutate(agecat_label = fct_relevel(agecat_label, "Under 41", "41-70", "Over 71"))

stratum_levelmar <- stratum_levelmar %>%
  left_join(strata_lookup, by = "strata_id") %>%
  mutate(agecat_label = fct_relevel(agecat_label, "Under 41", "41-70", "Over 71"))

stratum_levelmnar <- stratum_levelmnar %>%
  left_join(strata_lookup, by = "strata_id") %>%
  mutate(agecat_label = fct_relevel(agecat_label, "Under 41", "41-70", "Over 71"))

# Long form with model indicator — preserve strata_label_nage across sims
stratum_long <- bind_rows(
  stratum_levelfull %>% mutate(model = "Full"),
  stratum_levelmar  %>% mutate(model = "MAR"),
  stratum_levelmnar %>% mutate(model = "MNAR")
)

# Prepare x-axis levels and gen factor labels (not IDs)
x_levels <- strata_lookup %>%
  arrange(strata_id) %>%
  pull(strata_label_nage) %>%
  unique()

#combine using suffixed labels
stratum_long <- stratum_long %>%
  mutate(
    x_lab = factor(strata_label_nage.y, levels = x_levels),         
    minoritised = factor(minoritised, levels = c(0, 1),
                         labels = c("Majority", "Minority"))       
  )

# Plot w/ consistent dodge value and combined X ticks

dodge <- position_dodge(width = 0.6)

pred2 <- ggplot(
  stratum_long %>% filter(model != "MAR"),
  aes(y = m2mfit, x = x_lab, color = model, group = model)
) +
  geom_point(size = 1.8, position = dodge) +
  geom_errorbar(aes(ymax = m2mupr, ymin = m2mlwr),
                position = dodge, width = 0.3) +
  ylab("Predicted Health (MNAR Selection)") +
  xlab("") + ylim(2.7, 7.9) +
  theme_bw() +
  scale_x_discrete(drop = TRUE) + 
  facet_grid(
    rows = vars(minoritised),     
    cols  = vars(agecat_label),
    scales = "free_x"
  ) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )
# Checking
pred2

ggsave("MAIHDA_MNAR_fig4.png", pred2, dpi=500, height = 6, width = 6)
