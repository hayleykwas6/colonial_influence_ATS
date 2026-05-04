## ---------------------------
##
## Script name: ATS_code.R
##
## Purpose of script: Load the data. Characterize State patterns by ATS participation.
##                    Examine drivers of participation patterns. Visualize patterns.
##
## Authors: Hayley Kwasniewski, Dr. Peyton Thomas
##
## Date Created: 04-08-2026
##
## Copyright (c) Hayley Kwasniewski, 2026
## Email: hayley.kwasniewski@colorado.edu
##
## ---------------------------

##### ------------------------------------------------------
##### (1) Load and Prep Data for Analysis
##### ------------------------------------------------------

## ===== Load Packages =====
library(vegan)
library(car)
library(dplyr)
library(ggplot2)
library(factoextra)
library(cluster)
library(patchwork)
library(countrycode)
library(maps)
library(mice)

## ===== Read data =====
# Update this path to point to your local copy of the data file
df <- read.csv("data/....csv")

## ===== Factorize Categorical Variables and Scale Numerical =====
df$AT                 <- as.factor(df$AT)
df$CCAMLR             <- as.factor(df$CCAMLR)
df$Colonization.Status <- as.factor(df$Colonization.Status)
df$InNATO             <- as.factor(df$InNATO)
df$Continent          <- as.factor(df$Continent)
df$BRICS              <- as.factor(df$InBRICS)

df_use <- df %>%
  mutate(across(c(GDP, GDPperCap, Mil_Total, Mil_Per, Edu_Total, Edu_Per, Age, NatResRents),
                scale))

##### ------------------------------------------------------
##### (2) Create and Analyze Gower Distance Matrix
##### ------------------------------------------------------

## ===== Create distance matrix  =====
dummy_AT     <- model.matrix(~ AT - 1,     data = df_use)
dummy_CCAMLR <- model.matrix(~ CCAMLR - 1, data = df_use)

response      <- cbind(dummy_AT, dummy_CCAMLR)
dist_vector   <- vegdist(response, method = "gower")
dist_matrix_1 <- as.matrix(dist_vector)

## ===== Run NMDS =====
set.seed(123) # same random starting point for reproducibility
nmds_2d <- metaMDS(dist_vector, k = 2, trymax = 100)
set.seed(123)
nmds_3d <- metaMDS(dist_vector, k = 3, trymax = 100)

# check stress for valid/useable NMDS
cat("2D Stress:", nmds_2d$stress, "\n")
cat("3D Stress:", nmds_3d$stress, "\n")

nmds_result <- nmds_2d

## ===== Extract NMDS scores =====
nmds_scores                    <- as.data.frame(scores(nmds_result))
nmds_scores$Country            <- df_use$Countries
nmds_scores$AT                 <- df_use$AT
nmds_scores$CCAMLR             <- df_use$CCAMLR
nmds_scores$Continent          <- df_use$Continent
nmds_scores$InNATO             <- df_use$InNATO
nmds_scores$Colonization.Status <- df_use$Colonization.Status
nmds_scores$InBRICS            <- df_use$InBRICS

##### ------------------------------------------------------
##### (3) Run PERMANOVAs with Three Different Approaches
##### ------------------------------------------------------

## ===== Check Collinearity =====
nmds1_score <- scores(nmds_result)[, 1]

lm_model <- lm(nmds1_score ~ GDP + GDPperCap + Mil_Total + Mil_Per + Edu_Total + Edu_Per + Age + NatResRents +
                 Colonization.Status + InNATO + InBRICS + Continent,
               data = df_use, na.action = na.omit)

vif_results <- vif(lm_model)
print(vif_results)

## ===== Approach 1: Complete Cases =====
vars_used <- c("GDPperCap", "Edu_Per", "Mil_Per", "NatResRents",
               "Colonization.Status", "InNATO", "InBRICS", "Age", "Continent")

complete_rows <- complete.cases(df_use[, vars_used])
cat("Approach 1 — Complete cases:", sum(complete_rows), "of", nrow(df_use), "countries\n")
cat("Dropped:\n"); print(df_use$Countries[!complete_rows])

df_complete   <- df_use[complete_rows, ]
dist_complete <- vegdist(response[complete_rows, ], method = "gower") # subset response matrix

adonis_complete <- adonis2(dist_complete ~ GDPperCap + Edu_Per + Mil_Per + NatResRents +
                             Colonization.Status + InNATO + InBRICS + Age + Continent,
                           data = df_complete,
                           by = "margin", permutations = 9999)
print(adonis_complete)

## ===== Approach 2: Imputation =====

# define ALL variables needed — both imputed and analysis variables
vars_impute  <- c("GDPperCap", "Edu_Per", "Mil_Per", "NatResRents")
vars_analysis <- c("Colonization.Status", "InNATO", "InBRICS", "Age", "Continent")
vars_all      <- c(vars_impute, vars_analysis)

# prepare the data — include all analysis variables in imputation model
df_imputed <- df_use
df_imputed[, vars_impute] <- lapply(df_imputed[, vars_impute], as.numeric)
df_imputed$Age <- as.numeric(df_imputed$Age)

# check missing rates to determine number of imputations
missing_rates <- colSums(is.na(df_imputed[, vars_impute])) / nrow(df_imputed) * 100
cat("Missing rates (%):\n")
print(round(missing_rates, 1))

# number of imputations should be >= max % missing
m_imputations <- max(25, ceiling(max(missing_rates)))
cat("\nRunning", m_imputations, "imputations\n")

# run MICE with all analysis variables included in imputation model
imp <- mice(
  df_imputed[, vars_all],
  m          = m_imputations,
  method     = "pmm",
  maxit      = 20,
  seed       = 123,
  printFlag  = FALSE
)

# plot to assess convergence
plot(imp)

# analyse each imputed dataset and pool results
# Note: adonis2 is not directly compatible with mice's with() and pool() workflow,
# so we run the model on each completed dataset and extract results manually
adonis_results <- lapply(1:m_imputations, function(i) {
  df_complete_i <- complete(imp, i)  # extract the i-th completed dataset
  
  adonis2(
    dist_matrix_1 ~ GDPperCap + Edu_Per + Mil_Per + NatResRents +
      Colonization.Status + InNATO + InBRICS + Age + Continent,
    data         = df_complete_i,
    by           = "margin",
    permutations = 9999
  )
})

# pool R2 and F statistics across imputations (Rubin's rules)
r2_values <- do.call(rbind, lapply(adonis_results, function(x) x$R2))
f_values  <- do.call(rbind, lapply(adonis_results, function(x) x$F))

cat("\nPooled R2 (mean across imputations):\n")
print(round(colMeans(r2_values, na.rm = TRUE), 4))

cat("\nPooled F statistics (mean across imputations):\n")
print(round(colMeans(f_values, na.rm = TRUE), 4))

# derive pooled summary statistics before building results_table
predictor_names <- colnames(r2_values)
pooled_r2       <- colMeans(r2_values, na.rm = TRUE)
pooled_f        <- colMeans(f_values,  na.rm = TRUE)
pooled_sd       <- apply(r2_values, 2, sd, na.rm = TRUE)

results_table <- data.frame(
  Predictor = predictor_names,
  Pooled_R2 = round(pooled_r2, 4),
  Pooled_F  = round(pooled_f,  4),
  SD_R2     = round(pooled_sd, 4)
)

print(results_table)

## ===== Approach 3: Binning with Missing as a Category =====

# define function
bin_quantile <- function(x, n_bins = 4, var_name = "") {
  probs  <- seq(0, 1, length.out = n_bins + 1)
  breaks <- quantile(x, probs = probs, na.rm = TRUE)
  breaks <- unique(breaks)  # remove duplicate breaks if distribution is very discrete
  cat(var_name, "quantile breaks:", round(breaks, 3), "\n")
  binned <- cut(x, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  binned <- as.character(binned)
  binned[is.na(binned)] <- "Missing"
  as.factor(binned)
}

# plot histograms and inspect
dev.off()
par(mfrow = c(2, 2), mar = c(3, 3, 2, 1))
hist(as.numeric(df_use$GDPperCap),   breaks = 100, main = "GDPperCap",   xlab = "z-score")
hist(as.numeric(df_use$Edu_Per),     breaks = 100, main = "Edu_Per",     xlab = "z-score")
hist(as.numeric(df_use$Mil_Per),     breaks = 100, main = "Mil_Per",     xlab = "z-score")
hist(as.numeric(df_use$NatResRents), breaks = 100, main = "NatResRents", xlab = "z-score")
par(mfrow = c(1, 1))

# test different n_bins and inspect counts before deciding
for (n in 3:5) {
  cat("\nGDPperCap n_bins =", n, "\n")
  print(table(bin_quantile(as.numeric(df_use$GDPperCap), n_bins = n, var_name = "")))
}
for (n in 3:5) {
  cat("\nEdu_Per n_bins =", n, "\n")
  print(table(bin_quantile(as.numeric(df_use$Edu_Per), n_bins = n, var_name = "")))
}
for (n in 3:5) {
  cat("\nMil_Per n_bins =", n, "\n")
  print(table(bin_quantile(as.numeric(df_use$Mil_Per), n_bins = n, var_name = "")))
}
for (n in 3:5) {
  cat("\nNatResRents n_bins =", n, "\n")
  print(table(bin_quantile(as.numeric(df_use$NatResRents), n_bins = n, var_name = "")))
}

# after inspecting histograms and n States in bins, set n_bins per variable and bin
# n_bins chosen to balance bin counts while accounting for missingness:
# Mil_Per uses n_bins = 3 due to 50 missing countries (~25% of dataset)
# all others use n_bins = 4 (~50 countries per bin)
df_binned <- df_use
df_binned$GDPperCap_bin   <- bin_quantile(as.numeric(df_use$GDPperCap),   n_bins = 4, var_name = "GDPperCap")
df_binned$Edu_Per_bin     <- bin_quantile(as.numeric(df_use$Edu_Per),     n_bins = 4, var_name = "Edu_Per")
df_binned$Mil_Per_bin     <- bin_quantile(as.numeric(df_use$Mil_Per),     n_bins = 3, var_name = "Mil_Per")
df_binned$NatResRents_bin <- bin_quantile(as.numeric(df_use$NatResRents), n_bins = 4, var_name = "NatResRents")

# use dist_matrix_1 from section 2
# use binned variable names and print the correct object
adonis_binned <- adonis2(dist_matrix_1 ~ GDPperCap_bin + Edu_Per_bin + Mil_Per_bin + NatResRents_bin +
                           Colonization.Status + InNATO + InBRICS + Age + Continent,
                         data         = df_binned,
                         by           = "margin", permutations = 9999)
print(adonis_binned)

# remove clearly insignificant variables
adonis_binned_sig <- adonis2(dist_matrix_1 ~ Colonization.Status + InNATO + InBRICS + Continent,
                             data         = df_binned,
                             by           = "margin", permutations = 9999)
print(adonis_binned_sig)

##### ------------------------------------------------------
##### (4) Conduct Cluster Analysis
##### ------------------------------------------------------

## ===== Run Hierarchical Clustering =====
attr(dist_vector,   "Labels") <- df_use$Countries
attr(dist_matrix_1, "Labels") <- df_use$Countries
hclust_result <- hclust(dist_vector, method = "ward.D2")

## ===== Determine optimal clusters: Elbow Method =====
fviz_nbclust(as.matrix(dist_vector),
             FUN    = hcut,
             method = "wss",
             k.max  = 10) +
  labs(title = "Elbow Method for Optimal Clusters")

clusters_3 <- cutree(hclust_result, k = 3)

nmds_scores$Cluster3 <- as.factor(clusters_3)
df_use$Cluster3      <- as.factor(clusters_3)
df$Cluster3          <- as.factor(clusters_3)

## ===== Visualize NMDS =====

p1 <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = AT)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_bw() +
  labs(title = "AT Status") +
  theme(
    plot.title       = element_text(size = 12),
    legend.position  = c(0.8, 0.85),
    legend.text      = element_text(size = 8),
    legend.title     = element_blank()
  )

p2 <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = CCAMLR)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_bw() +
  labs(title = "CCAMLR Status") +
  theme(
    plot.title       = element_text(size = 12),
    axis.title.y     = element_blank(),
    legend.position  = c(0.8, 0.85),
    legend.text      = element_text(size = 8),
    legend.title     = element_blank()
  )

p3 <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = Cluster3)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_bw() +
  labs(title = "Cluster") +
  theme(
    plot.title       = element_text(size = 12),
    axis.title.y     = element_blank(),
    legend.position  = c(0.85, 0.85),
    legend.text      = element_text(size = 8),
    legend.title     = element_blank()
  )

p1 + p2 + p3

## ===== Summarize Cluster Characterization =====
cluster_summary_3 <- df %>%
  group_by(Cluster3) %>%
  summarise(
    n               = n(),
    AT_mode         = names(sort(table(AT),                  decreasing = TRUE)[1]),
    CCAMLR_mode     = names(sort(table(CCAMLR),              decreasing = TRUE)[1]),
    Continent_mode  = names(sort(table(Continent),           decreasing = TRUE)[1]),
    NATO_mode       = names(sort(table(InNATO),              decreasing = TRUE)[1]),
    BRICS_mode      = names(sort(table(InBRICS),             decreasing = TRUE)[1]),
    Colonial_mode   = names(sort(table(Colonization.Status), decreasing = TRUE)[1]),
    GDPperCap_mean  = mean(GDPperCap,   na.rm = TRUE),
    Edu_Per_mean    = mean(Edu_Per,     na.rm = TRUE),
    Mil_Per_mean    = mean(Mil_Per,     na.rm = TRUE),
    Age_mean        = mean(Age,         na.rm = TRUE),
    NatResRents_mean = mean(NatResRents, na.rm = TRUE)
  )

print(cluster_summary_3)

## ===== Check Cluster Validity =====
# determine differences between clusters using PERMANOVA
cluster_permanova <- adonis2(dist_matrix_1 ~ Cluster3,
                             data         = df_use,
                             permutations = 9999)
print(cluster_permanova)

# get silhouette values
sil_3   <- silhouette(clusters_3, dist_vector)
avg_sil <- mean(sil_3[, 3])
cat("Average Silhouette Width:", round(avg_sil, 3), "\n")

summary(sil_3)

##### ------------------------------------------------------
##### (5) Identify Drivers of Clusters
##### ------------------------------------------------------

## ===== Determine Drivers of Clusters =====
var_tests <- list()
var_tests$Colonization <- adonis2(dist_matrix_1 ~ Colonization.Status, data = df_use, permutations = 9999, na.action = na.omit)
var_tests$NATO         <- adonis2(dist_matrix_1 ~ InNATO,              data = df_use, permutations = 9999, na.action = na.omit)
var_tests$BRICS        <- adonis2(dist_matrix_1 ~ InBRICS,             data = df_use, permutations = 9999, na.action = na.omit)
var_tests$Continent    <- adonis2(dist_matrix_1 ~ Continent,           data = df_use, permutations = 9999, na.action = na.omit)
var_tests$GDPperCap    <- adonis2(dist_matrix_1 ~ GDPperCap,           data = df_use, permutations = 9999, na.action = na.omit)
var_tests$Mil_Per      <- adonis2(dist_matrix_1 ~ Mil_Per,             data = df_use, permutations = 9999, na.action = na.omit)
var_tests$Edu_Per      <- adonis2(dist_matrix_1 ~ Edu_Per,             data = df_use, permutations = 9999, na.action = na.omit)
var_tests$Age          <- adonis2(dist_matrix_1 ~ Age,                 data = df_use, permutations = 9999, na.action = na.omit)
var_tests$NatRes       <- adonis2(dist_matrix_1 ~ NatResRents,         data = df_use, permutations = 9999, na.action = na.omit)

r2_values <- sapply(var_tests, function(x) x$R2[1])
p_values  <- sapply(var_tests, function(x) x$`Pr(>F)`[1])

var_importance <- data.frame(
  Variable    = names(r2_values),
  R2          = round(r2_values, 4),
  P_value     = p_values,
  Significant = ifelse(p_values < 0.001, "***",
                       ifelse(p_values < 0.01, "**",
                              ifelse(p_values < 0.05, "*", "ns")))
)

var_importance <- var_importance[order(-var_importance$R2), ]
print(var_importance)

# plot drivers of clustering
var_labels <- c(
  "Colonization" = "Colonial Status",
  "NATO"         = "NATO Status",
  "BRICS"        = "BRICS Status",
  "Continent"    = "Continent",
  "Age"          = "Age",
  "GDPperCap"    = "GDP per Capita",
  "NatRes"       = "Natural Resources Rents",
  "Mil_Per"      = "Military Spending",
  "Edu_Per"      = "Education Spending"
)

ggplot(var_importance, aes(x = reorder(Variable, R2), y = R2, fill = Significant)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_x_discrete(labels = var_labels) +
  labs(y = "R² (Proportion of Variance Explained)") +
  theme_minimal() +
  scale_fill_manual(values = c("***" = "darkgreen", "**" = "green",
                               "*"   = "lightgreen", "ns" = "grey"))

##### ------------------------------------------------------
##### (6) Visualize Global Clusters
##### ------------------------------------------------------

## ===== Mapping AT/CCAMLR Clusters =====
# load ISO3 codes
df$iso3 <- countrycode(df$Countries,
                       origin      = "country.name",
                       destination = "iso3c")
# Note: Kosovo and Channel Islands cannot be matched to ISO3 codes and will
# appear as NA; these countries are too small to be visible at map resolution

# get world map with ISO codes
world_map      <- map_data("world")
world_map$iso3 <- countrycode(world_map$region,
                              origin      = "country.name",
                              destination = "iso3c")
# Note: some regions in the map data (e.g. island territories) cannot be matched
# to ISO3 codes; these are too small to be visible at map resolution

# combine
map_clusters <- world_map %>%
  left_join(dplyr::select(df, iso3, Cluster3), by = "iso3", relationship = "many-to-many")

# plot
ggplot(map_clusters, aes(x = long, y = lat, group = group, fill = Cluster3)) +
  geom_polygon(color = "white", linewidth = 0.1) +
  scale_fill_manual(values   = c("1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A"),
                    na.value = "grey90",
                    name     = "Cluster") +
  theme_void() +
  labs(title = "Global Distribution of AT/CCAMLR Participation Clusters") +
  coord_fixed(1.3)
