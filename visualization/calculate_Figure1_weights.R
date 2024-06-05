library(tidyverse)

exponential_weight <- function(theta, delta) {
  alpha <- log(theta)
  weight <- exp(-alpha * delta)
  return(weight)
}

normalized_exponential_weight <- function(theta, delta, generation) {
  alpha <- log(theta) / 2
  weight <- exp((-alpha * delta) / 2^(-(generation+1)/2))
  return(weight)
}

normal_distribution_weight <- function(theta, parent1_ancestry0,
                                       parent2_ancestry0) {
  calculate_weight <- function(theta, parent1_ancestry0, parent2_ancestry0) {
    denom <- 2 * log(1 / theta)
    sigma <- sqrt(-1 / denom)
    weight <- dnorm(parent2_ancestry0, mean = parent1_ancestry0, sd = sigma)
    return(weight)
  }

  weight <- suppressWarnings(case_when(theta == 1 ~ 1,
                      .default = calculate_weight(theta, parent1_ancestry0,
                                                  parent2_ancestry0)))
  return(weight)
}

ancestry_psi <- expand_grid(x1 = c(0, 0.5),
                              xi = seq(0, 1, by = 0.05),
                              theta = c(1, 4, 10)) %>%
  mutate(delta = abs(x1 - xi),
         exponential_decay = exponential_weight(theta, delta),
         exponential2 = normalized_exponential_weight(theta, delta, 2),
         exponential3 = normalized_exponential_weight(theta, delta, 3),
         normal_distribution = normal_distribution_weight(theta, x1, xi)) %>%
  group_by(theta) %>%
  mutate(across(exponential_decay:normal_distribution, ~ .x / sum(.x)),
         n = 1 / n()) %>%
  ungroup() %>%
  mutate(across(exponential_decay:normal_distribution, ~ .x / n)) %>%
  select(-n) %>%
  pivot_longer(cols = exponential_decay:normal_distribution, names_to = "model",
               values_to = "psi") %>%
  mutate(theta = factor(theta, levels = c(1, 4, 10)),
         generation = factor(case_when(model == "exponential2" ~ 2,
                                       model == "exponential3" ~ 3,
                                       .default = 1), levels = c(1, 2, 3)))

normalized_exp <- ancestry_psi %>%
  filter(grepl("exponential", model), theta %in% c(1, 4)) %>%
  mutate(model = "exponential_decay_normalized")

ancestry_psi <- ancestry_psi %>%
  filter(!(model %in% c("exponential2", "exponential3"))) %>%
  rbind(., normalized_exp)


social_psi <- tibble(theta = c(1, 4, 10),
                                  different = 1 / (theta + 1),
                                  same = 1 - different) %>%
  pivot_longer(cols = c(same, different), names_to = "social_group",
               values_to = "psi") %>%
  mutate(psi = psi / 0.5)

parent_dir <- "data/"

write_csv(ancestry_psi, paste0(parent_dir, "/example_weights_ancestry.csv"))

write_csv(social_psi, paste0(parent_dir, "/example_weights_social_group.csv"))