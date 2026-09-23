############################################################
############################################################

library(dplyr)
library(ggplot2)
library(smplot2)

# Load data -------------------------------------------------------------
siRNA <- read.delim("data/siRNA_rtqPCR_results.txt", sep = "\t", header = TRUE)

############################################################
# 12 µM siRNA treatments


# --- ApB744 ------------------------------------------------------------
b744_12 <- siRNA %>%
  filter(Incubation == 24,
         siRNA_conc_uM == 12,
         Treatment %in% c("Universal", "ApB744")) %>%
  mutate(Treatment = factor(Treatment, levels = c("Universal", "ApB744")))

# Summary statistics
b744_12 %>%
  group_by(Treatment) %>%
  summarise(
    n = n(),
    mean = mean(LOG_B744_expression, na.rm = TRUE),
    sd = sd(LOG_B744_expression, na.rm = TRUE)
  )

# Welch’s one-tailed t-test
t.test(LOG_B744_expression ~ Treatment,
       data = b744_12,
       var.equal = FALSE,
       alternative = "less")

# --- ApW584 ------------------------------------------------------------
w584_12 <- siRNA %>%
  filter(Incubation == 24,
         siRNA_conc_uM == 12,
         Treatment %in% c("Universal", "ApW584")) %>%
  mutate(Treatment = factor(Treatment, levels = c("Universal", "ApW584")))

w584_12 %>%
  group_by(Treatment) %>%
  summarise(
    n = n(),
    mean = mean(LOG_W584_expression, na.rm = TRUE),
    sd = sd(LOG_W584_expression, na.rm = TRUE)
  )

# One-tailed Student’s t-test
t.test(LOG_W584_expression ~ Treatment,
       data = w584_12,
       var.equal = TRUE,
       alternative = "less")

############################################################
# 20 µM siRNA treatments


# --- ApB744 ------------------------------------------------------------
b744_20 <- siRNA %>%
  filter(Incubation == 24,
         siRNA_conc_uM == 20,
         Treatment %in% c("Universal", "ApB744")) %>%
  mutate(Treatment = factor(Treatment, levels = c("Universal", "ApB744")))

b744_20 %>%
  group_by(Treatment) %>%
  summarise(
    n = n(),
    mean = mean(LOG_B744_expression, na.rm = TRUE),
    sd = sd(LOG_B744_expression, na.rm = TRUE)
  )

# One-tailed Student’s t-test
t.test(LOG_B744_expression ~ Treatment,
       data = b744_20,
       var.equal = TRUE,
       alternative = "less")

# --- ApW584 ------------------------------------------------------------
w584_20 <- siRNA %>%
  filter(Incubation == 24,
         siRNA_conc_uM == 20,
         Treatment %in% c("Universal", "ApW584")) %>%
  mutate(Treatment = factor(Treatment, levels = c("Universal", "ApW584")))

w584_20 %>%
  group_by(Treatment) %>%
  summarise(
    n = n(),
    mean = mean(LOG_W584_expression, na.rm = TRUE),
    sd = sd(LOG_W584_expression, na.rm = TRUE)
  )

# One-tailed Student’s t-test
t.test(LOG_W584_expression ~ Treatment,
       data = w584_20,
       var.equal = TRUE,
       alternative = "less")

############################################################
# Visualization (example for ApB744 12 µM)

ggplot(b744_12, aes(x = Treatment, y = LOG_B744_expression, fill = Treatment)) +
  sm_raincloud() +
  scale_fill_manual(values = sm_color("orange", "blue")) +
  labs(y = "log(Expression, ApB744)", x = NULL)
