library(tidyverse)
library(scales)

# data1 <- read.csv("Comparisons.csv")
# 
# data1 %>%
#   ggplot(aes(x = 10^Lowest_aPOD, y = 10^Lowest_tPOD)) +
#   geom_point(size = 3.5, aes(shape = Type)) +
#   scale_x_log10(breaks = c(10^-4, 10^-2, 10^0, 10^2, 10^4),
#                 limits = c(10^-5, 10^5),
#                 labels = trans_format("log10", math_format(10^.x))) +
#   scale_y_log10(breaks = c(10^-4, 10^-2, 10^0, 10^2, 10^4),
#                 limits = c(10^-5, 10^5),
#                 labels = trans_format("log10", math_format(10^.x))) +
#   geom_abline(intercept = 0, slope = 1) +
#   # geom_abline(intercept = 1.5, slope = 1) +
#   labs(x = "Most Sensitive Acute aPOD", y = "Most Sensitive tPOD") +
#   theme_classic()
  
data2 <- read.csv("Comparisons_Ratio.csv")

data2 %>% 
  ggplot(aes(x = Chemical, y = Ratio, label = Chemical)) +
  geom_point(size = 3, aes(color = factor(Type), shape = factor(Type)), show.legend = F) +
  geom_text(aes(label = round(Ratio, digits = 2)), vjust = -1) +
  facet_wrap(~Type, scales = "free_x") +
  geom_hline(yintercept = 0, color = "red") +
  scale_color_manual(values = (c("black", "blue"))) +
  scale_y_continuous(breaks = c(-6, -4, -2, 0, 2, 4),
                     limits = c(-6, 4.5)) +
  labs(y = "log10 ( tPOD/aPOD )") +
  theme_classic()
