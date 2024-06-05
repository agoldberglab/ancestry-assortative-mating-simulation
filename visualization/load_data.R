library(magrittr)

load_csv <- function(filename, models, thetas, migration_rates) {
  df <- readr::read_csv(filename, col_types = "d") %>%
    dplyr::mutate(simulation = factor(simulation),
           model = factor(model, levels = models),
           theta = factor(theta, levels = c(seq(1,10), 25, 50, 75, 100)),
           seed = factor(seed),
           migration = factor(migration, levels = migration_rates))
  
  return(df)
  
}

symbol_theta <- "<span style='font-family: Georgia'>\u03b1</span>"
symbol_psi <- "<span style='font-family: Georgia'>\u03c8</span>"
symbol_psi_i <- "<span style='font-family: Georgia'>\u03c8</span>~i~"

values <- list()
labels <- list()

thetas <- c(seq(1, 10, by = 1), seq(25, 100, by = 25))
thetas_shown <- c(1, seq(2, 10, by = 2), seq(25, 100, by = 25))
values$theta <- set_names(c("#737373", "#C6DBEF", "#9ECAE1", "#6BAED6",
                            "#3182BD", "#08519C", "#A1D99B", "#74C476",
                            "#31A354", "#006D2C"), thetas_shown)

labels$theta <- lapply(thetas, function(x) {
  if(x == 1) { return("Random") } else {
    return(paste0(symbol_theta, " = ", x)) }
})

labels$theta <- set_names(labels$theta, thetas)

models <- c("random", "exponential_decay", "exponential_decay_normalized",
            "normal_distribution", "subpop")
values$model <- set_names(c("#737373", "#4477AA", "#EE6677", "#228833",
                            "#CCBB44"), models)

labels$model <- set_names(c("Random mating", "Stationary-preference model",
                            "Increasing-preference model", 
                            "Broad-preference model",
                            "Social group model"), models)

labels$short_model <- set_names(c("Random mating", "Stationary-preference",
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
  select("seed", "model", "theta", "migration") %>%
  mutate(example = factor(c(1, 2, 2, 3, 3, 4, 4)),
         model = case_when(model == "exponential_decay" ~ "random",
                           .default = model),
         model = factor(model, levels = models))

summary_data <- summary_data %>%
  left_join(., paired_examples, by = c("seed", "model", "theta", "migration"))

binned_ancestry0 <- binned_ancestry0 %>%
  left_join(., paired_examples, by = c("seed", "model", "theta", "migration"))

rm(paired_examples)

example_pedigree <- load_csv(paste0(parent_dir, "/pedigree.csv"),
                     models, thetas, migration_rates) %>%
  dplyr::filter(generation == 20) %>%
  dplyr::mutate(delta = abs(parent1_ancestry0 - parent2_ancestry0))

example_psi_ratio <- load_csv(paste0(parent_dir, "/pedigree_example_psi.csv"),
                             models, thetas, migration_rates)

shuffled_rho <- load_csv(paste0(parent_dir, "/shuffled_rho.csv"),
                              models, thetas, migration_rates)

shuffled_delta <- load_csv(paste0(parent_dir, "/shuffled_delta.csv"),
                         models, thetas, migration_rates)

ancestry_psi <- read_csv(paste0(parent_dir, "/example_weights_ancestry.csv"),
                         col_types = "ddfdfdf")

social_psi <- read_csv(paste0(parent_dir, "/example_weights_social_group.csv"),
                       col_types = "ffd")

titles <- c("ancestry0" = "Ancestry proportion",
            "proportion" = "Proportion",
            "generation" = "Generations post-admixture",
            "parent_corr" = "Correlation between mates",
            "ancestry0_var" = "Variance in ancestry",
            "obs_delta" = bquote(bar(Delta)[obs]),
            "parent1_ancestry0" = "Parent 1 ancestry proportion",
            "parent2_ancestry0" = "Parent 2 ancestry proportion",
            "expressed_mating_bias" = "B",
            "timing_discrepancy" = "Discrepancy in est. time",
            "psi" = symbol_psi_i,
            "xi" = bquote(x[i]),
            "social_group" = "",
            "tract_median" = "Median tract length, cM",
            "timing_0" = "Est. time post-admixture",
            "p_migrant_pairs" = "% matings between migrants",
            "theta" = symbol_theta)