############################################################
# Cholesterol supplementation: Buchnera titer analysis


library(dplyr)
library(ggplot2)
library(smplot2)
library(DescTools)

############################################################
# Load data

cholesterol <- read.delim("data/cholesterol_injection_titer.txt", sep = "\t", header = TRUE)

############################################################
# 1 µg/µl cholesterol – 24 h post injection

chol_24 <- cholesterol %>%
  filter(conc_ug_ul == 1,
         Duration == 24,
         Treatment %in% c("Cholesterol", "NegControl")) %>%
  mutate(Treatment = factor(Treatment, levels = c("NegControl", "Cholesterol")))

# Summary statistics
chol_24 %>%
  group_by(Treatment) %>%
  summarise(
    n = n(),
    mean = mean(LOG_Buchnera_titer, na.rm = TRUE),
    sd = sd(LOG_Buchnera_titer, na.rm = TRUE)
  )

# Two-tailed Student’s t-test (matches Fig. 5I)
t.test(LOG_Buchnera_titer ~ Treatment,
       data = chol_24,
       var.equal = TRUE,
       alternative = "two.sided")

# Visualization
ggplot(chol_24, aes(x = Treatment, y = LOG_Buchnera_titer, fill = Treatment)) +
  sm_raincloud() +
  scale_fill_manual(values = sm_color("orange", "blue")) +
  labs(y = "log(Buchnera titer, 24 h)", x = NULL)

############################################################
# 1 µg/µl cholesterol – 48 h post injection

chol_48 <- cholesterol %>%
  filter(conc_ug_ul == 1,
         Duration == 48,
         Treatment %in% c("Cholesterol", "NegControl", "Ptx")) %>%
  mutate(Treatment = factor(Treatment, levels = c("Ptx", "NegControl", "Cholesterol")))

# Summary statistics
chol_48 %>%
  group_by(Treatment) %>%
  summarise(
    n = n(),
    mean = mean(LOG_Buchnera_titer, na.rm = TRUE),
    sd = sd(LOG_Buchnera_titer, na.rm = TRUE)
  )

# One-way ANOVA + Dunnett’s post hoc vs. control (matches Fig. 5J)
anova_48 <- aov(LOG_Buchnera_titer ~ Treatment, data = chol_48)
summary(anova_48)
DunnettTest(LOG_Buchnera_titer ~ Treatment, data = chol_48, control = "Ptx")

# Visualization
ggplot(chol_48, aes(x = Treatment, y = LOG_Buchnera_titer, fill = Treatment)) +
  sm_raincloud() +
  scale_fill_manual(values = sm_color("red", "orange", "blue")) +
  labs(y = "log(Buchnera titer, 48 h)", x = NULL)
