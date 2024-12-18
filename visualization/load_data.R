library(magrittr)
library(dplyr)

load_csv <- function(filename, models, thetas, migration_rates) {
  df <- readr::read_csv(filename, col_types = "d") %>%
    mutate(simulation = factor(simulation),
           model = factor(model, levels = models),
           theta = factor(theta, levels = c(seq(1,10), 25, 50, 75, 100)),
           seed = factor(seed),
           migration = factor(migration, levels = migration_rates),
           initial = factor(initial, levels = seq(0.5, 0.9, by = 0.1)))
  
  return(df)
  
}

symbol_theta <- "<span style='font-family: Georgia'>\u03b1</span>"
symbol_psi <- "<span style='font-family: Georgia'>\u03c8</span>"
symbol_psi_i <- "<span style='font-family: Georgia'>\u03c8</span>~i~,~j~"
symbol_sigma <- paste0("<span style='font-family: Georgia'>",
                       "\u03c3</span><sup>2</sup>~x~(t)")

values <- list()
labels <- list()

thetas <- c(seq(1, 10, by = 1), seq(25, 100, by = 25))
thetas_shown <- c(1, seq(2, 10, by = 2), seq(25, 100, by = 25))
values$theta <- set_names(c("#737373", "#c6dbef", "#9ecae1", "#6baed6",
                            "#3182bd", "#08519c", "#a1d99b", "#74c476",
                            "#31a354", "#006d2c"), thetas_shown)

labels$theta <- lapply(thetas, function(x) {
  if(x == 1) { return("Random") } else {
    return(paste0(symbol_theta, " = ", x)) }
})

labels$theta <- set_names(labels$theta, thetas)

models <- c("random", "exponential_decay", "exponential_decay_normalized",
            "normal_distribution", "subpop")
values$model <- set_names(c("#737373", "#4477AA", "#EE6677", "#228833",
                            "#CCBB44"), models)

labels$model <- set_names(c("Random mating",
                            "Stationary-preference model",
                            "Increasing-preference model", 
                            "Broad-preference model",
                            "Social group model"), models)

labels$short_model <- set_names(c("Random mating",
                                  "Stationary-preference",
                                  "Increasing-preference", 
                                  "Broad-preference",
                                  "Social group"), models)

generations <- c(1, 2, 3)
values$generation <- set_names(c("solid", "dashed", "dotted"), generations)

migration_rates <- c(0, 0.01)
values$migration <- set_names(c("solid", "dashed"), migration_rates)
labels$migration <- set_names(c("Single-pulse admixture",
                                "Continuous migration"),
                              migration_rates)

subpops <- c(0, 1)
values$subpop <- set_names(c("0" = "#FC8D62", "1" = "#8DA0CB"), subpops)
labels$subpop <- set_names(c("0" = "Social group A", "1" = "Social group B"),
                           subpops)

labels$mu <- c("observed" = bquote(bar(Delta)[obs]),
               "permuted" = bquote(bar(Delta)[perm]))
values$mu <- c("observed" = "#2C7FB8", "permuted" = "#636363")

parent_dir <- "data/"

binned_ancestry0 <- load_csv(paste0(parent_dir, "/binned_ancestry0.csv"),
                             models, thetas, migration_rates)

summary_data <- load_csv(paste0(parent_dir, "/summary.csv"),
                         models, thetas, migration_rates)

paired_examples <- readr::read_csv(paste0(parent_dir,
                                          "/example_parameters_paired.csv"),
                                   col_names = c("seed", "model", "theta",
                                                 "initial", "migration", "N",
                                                 "generations"),
                                   col_types = "ffffff") %>%
  select("seed", "model", "theta", "migration", "initial") %>%
  mutate(example = factor(c(1, 2, 2, 3, 3, 4, 4)),
         model = case_when(model == "exponential_decay" ~ "random",
                           .default = model),
         model = factor(model, levels = models))

summary_data <- summary_data %>%
  left_join(., paired_examples,
            by = c("seed", "model", "theta", "migration", "initial"))

summary_data <- summary_data %>%
  mutate(tracts_exp_min = tracts_exp_min * 100,
         tracts_exp_q25 = tracts_exp_q25 * 100,
         tracts_exp_q50 = tracts_exp_q50 * 100,
         tracts_exp_q75 = tracts_exp_q75 * 100,
         tracts_exp_max = tracts_exp_max * 100,
         p_outliers = p_outliers * 100)

binned_ancestry0 <- binned_ancestry0 %>%
  left_join(., paired_examples, by = c("seed", "model", "theta", "migration",
                                       "initial"))

rm(paired_examples)

example_pedigree <- load_csv(paste0(parent_dir, "/pedigree.csv"),
                             models, thetas, migration_rates) %>%
  mutate(delta = abs(parent1_ancestry0 - parent2_ancestry0),
         offspring = offspring1 + offspring2,
         example = as.factor(example))

initial_proportion <- load_csv(
  paste0(parent_dir, "/pedigree_example_initial_proportion.csv"),
  models, thetas, migration_rates) %>%
  mutate(offspring = offspring1 + offspring2)

example_psi_ratio <- load_csv(paste0(parent_dir, "/pedigree_example_psi.csv"),
                              models, thetas, migration_rates)

shuffled_rho <- load_csv(paste0(parent_dir, "/shuffled_rho.csv"),
                         models, thetas, migration_rates)

ancestry_psi <- readr::read_csv(paste0(parent_dir,
                                       "/example_weights_ancestry.csv"),
                                col_types = "ddfdfdf")

social_psi <- readr::read_csv(paste0(parent_dir,
                                     "/example_weights_social_group.csv"),
                              col_types = "ffd")

titles <- c("ancestry0" = "x",
            "proportion" = "Frequency",
            "generation" = "Generations, t",
            "parent_corr" = bquote(r(x[i],x[j])),
            "ancestry0_var" = paste0("Variance in x, ", symbol_sigma),
            "obs_delta" = bquote(bar(Delta)[obs]),
            "parent1_ancestry0" = bquote(Parent~1*","~x[i]),
            "parent2_ancestry0" = bquote(Parent~2*","~x[j]),
            "expressed_mating_bias" = "B",
            "timing_discrepancy" = "Discrepancy in estimated t",
            "psi" = paste0("Scaled ", symbol_psi_i),
            "xi" = bquote(x[j]),
            "social_group" = "",
            "tract_median" = "Median local-ancestry tract length, cM",
            "timing_0" = "Estimated t",
            "p_migrant_pairs" = "% matings between migrants",
            "theta" = symbol_theta)
