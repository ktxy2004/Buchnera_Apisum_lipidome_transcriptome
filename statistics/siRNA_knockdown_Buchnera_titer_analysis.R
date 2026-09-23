
############################################################

library(dplyr)
library(ggplot2)
library(smplot2)
library(DescTools)

############################################################
# Load data

titer <- read.delim("data/siRNA_titer_results.txt", sep = "\t", header = TRUE)

############################################################
# 12 µM siRNA (30 h)
# (ApB744, ApW584)

titer_12 <- titer %>%
  filter(siRNA_conc_uM == 12,
         Duration == 30,
         Treatment %in% c("Universal", "ApB744", "ApW584")) %>%
  mutate(Treatment = factor(Treatment, levels = c("Universal", "ApB744", "ApW584")))

# Summary statistics
titer_12 %>%
  group_by(Treatment) %>%
  summarise(
    n = n(),
    mean = mean(LOG_Buchnera_titer, na.rm = TRUE),
    sd = sd(LOG_Buchnera_titer, na.rm = TRUE)
  )

# ANOVA + Dunnett’s post hoc
anova_12 <- aov(LOG_Buchnera_titer ~ Treatment, data = titer_12)
summary(anova_12)
DunnettTest(LOG_Buchnera_titer ~ Treatment, data = titer_12, control = "Universal")

# Visualization
ggplot(titer_12, aes(x = Treatment, y = LOG_Buchnera_titer, fill = Treatment)) +
  sm_raincloud() +
  scale_fill_manual(values = sm_color("orange", "blue", "pink")) +
  labs(y = "log(Buchnera titer)", x = NULL)

############################################################
# 20 µM siRNA (50 h)
# knockdown (ApB744, ApW584)

titer_20 <- titer %>%
  filter(siRNA_conc_uM == 20,
         Duration == 50,
         Treatment %in% c("Universal", "ApB744", "ApW584")) %>%
  mutate(Treatment = factor(Treatment, levels = c("Universal", "ApB744", "ApW584")))

# Summary statistics
titer_20 %>%
  group_by(Treatment) %>%
  summarise(
    n = n(),
    mean = mean(LOG_Buchnera_titer, na.rm = TRUE),
    sd = sd(LOG_Buchnera_titer, na.rm = TRUE)
  )

# ANOVA + Dunnett’s post hoc
anova_20 <- aov(LOG_Buchnera_titer ~ Treatment, data = titer_20)
summary(anova_20)
DunnettTest(LOG_Buchnera_titer ~ Treatment, data = titer_20, control = "Universal")

# Visualization
ggplot(titer_20, aes(x = Treatment, y = LOG_Buchnera_titer, fill = Treatment)) +
  sm_raincloud() +
  scale_fill_manual(values = sm_color("orange", "blue", "pink")) +
  labs(y = "log(Buchnera titer)", x = NULL)
