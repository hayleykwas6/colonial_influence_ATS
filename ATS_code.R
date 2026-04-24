## ---------------------------
##
## Script name: ATS_code.R
##
## Purpose of script: Load the data. Characterize State patterns by ATS participation.
##                    Examine drivers of participation patterns. Visualize patterns.
##
## Author: Hayley Kwasniewski, Dr. Peyton Thomas
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

## ===== Read data =====
df <- read.csv("~/Documents/Data/ATS/ATSData.csv") #rename path

## ===== Factorize Categorical Variables and Scale Numerical =====
df$AT <- as.factor(df$AT)
df$CCAMLR <- as.factor(df$CCAMLR)
df$Colonization.Status <- as.factor(df$Colonization.Status)
df$In.NATO <- as.factor(df$In.NATO)
df$Continent <- as.factor(df$Continent)

df_use <- df %>%
  mutate(across(c(GDP, GDPperCap, Mil_Total, Mil_Per, Edu_Total, Edu_Per, Age, NatResRents), 
                scale))

##### ------------------------------------------------------
##### (2) Create and Analyze Gower Distance Matrix
##### ------------------------------------------------------

## ===== Create distance matrix  =====
dummy_AT <- model.matrix(~ AT - 1, data = df_use)
dummy_CCAMLR <- model.matrix(~ CCAMLR - 1, data = df_use)

response <- cbind(dummy_AT, dummy_CCAMLR)
dist_vector <- vegdist(response, method = "gower")
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

## ===== Check for Co-linearity =====
nmds1_score <- scores(nmds_result)[, 1]

lm_model <- lm(nmds1_score ~ GDP + GDPperCap + Mil_Total + Mil_Per + Edu_Total + Edu_Per + Age + NatResRents + 
                 Colonization.Status + In.NATO + Continent, 
               data = df_use,
               na.action = na.omit)  # make the omission explicit

vif_results <- vif(lm_model)
print(vif_results) # omit relevant variables from PERMANOVA

## ===== Run PERMANOVA =====
adonis_1 <- adonis2(dist_matrix_1 ~ GDPperCap + Edu_Per + Mil_Per + NatResRents + 
                      Colonization.Status + In.NATO + Age + Continent,
                    data = df_use,
                    by = "margin",
                    permutations = 9999)
summary(adonis_1)

## ===== Extract NMDS scores =====
nmds_scores <- as.data.frame(scores(nmds_result))
nmds_scores$Country <- df_use$Countries
nmds_scores$AT <- df_use$AT
nmds_scores$CCAMLR <- df_use$CCAMLR
nmds_scores$Continent <- df_use$Continent
nmds_scores$In.NATO <- df_use$In.NATO
nmds_scores$Colonization.Status <- df_use$Colonization.Status

##### ------------------------------------------------------
##### (3) Conduct Cluster Analysis
##### ------------------------------------------------------

## ===== Run Hierarchical Clustering =====
attr(dist_vector, "Labels") <- df_use$Countries
attr(dist_matrix_1, "Labels") <- df_use$Countries
hclust_result <- hclust(dist_vector, method = "ward.D2")

## ===== Determine optimal clusters: Elbow Method =====
fviz_nbclust(as.matrix(dist_vector), 
             FUN = hcut,
             method = "wss",
             k.max = 10) +
  labs(title = "Elbow Method for Optimal Clusters") 

clusters_3 <- cutree(hclust_result, k = 3)

nmds_scores$Cluster3 <- as.factor(clusters_3)
df_use$Cluster3 <- as.factor(clusters_3)
df$Cluster3 <- as.factor(clusters_3)

## ===== Visualize NMDS =====

p1 <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = AT)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_bw() +
  labs(title = "AT Status") +
  theme(
    plot.title = element_text(size = 12),
    legend.position = c(0.8, 0.85),
    legend.text = element_text(size = 8),
    legend.title = element_blank()
  )

p2 <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = CCAMLR)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_bw() +
  labs(title = "CCAMLR Status") +
  theme(
    plot.title = element_text(size = 12),
    axis.title.y = element_blank(),
    legend.position = c(0.8, 0.85),    
    legend.text = element_text(size = 8),
    legend.title = element_blank()
  )

p3 <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = Cluster3)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_bw() +
  labs(title = "Cluster") +
  theme(
    plot.title = element_text(size = 12),
    axis.title.y = element_blank(),
    legend.position = c(0.85, 0.85),
    legend.text = element_text(size = 8),
    legend.title = element_blank()
  )

p1 + p2 + p3


## ===== Summarize Cluster Characterization =====
cluster_summary_3 <- df %>%
  group_by(Cluster3) %>%
  summarise(
    n = n(),
    AT_mode = names(sort(table(AT), decreasing = TRUE)[1]),
    CCAMLR_mode = names(sort(table(CCAMLR), decreasing = TRUE)[1]),
    Continent_mode = names(sort(table(Continent), decreasing = TRUE)[1]),
    NATO_mode = names(sort(table(In.NATO), decreasing = TRUE)[1]),
    Colonial_mode = names(sort(table(Colonization.Status), decreasing = TRUE)[1]),
    GDPperCap_mean = mean(GDPperCap, na.rm = TRUE),
    Edu_Per_mean = mean(Edu_Per, na.rm = TRUE),
    Mil_Per_mean = mean(Mil_Per, na.rm = TRUE),
    Age_mean = mean(Age, na.rm = TRUE),
    NatResRents_mean = mean(NatResRents, na.rm = TRUE)
  )

print(cluster_summary_3)

## ===== Check Cluster Validity =====
# determine differences between clusters using PERMANOVA
cluster_permanova <- adonis2(dist_matrix_1 ~ Cluster3,
                             data = df_use,
                             permutations = 9999)
print(cluster_permanova)

# get silhouette values
sil_3 <- silhouette(clusters_3, dist_vector)
avg_sil <- mean(sil_3[, 3])
cat("Average Silhouette Width:", round(avg_sil, 3), "\n")

summary(sil_3)

##### ------------------------------------------------------
##### (4) Identify Drivers of Clusters
##### ------------------------------------------------------

## ===== Determine Drivers of Clusters =====
var_tests <- list()

var_tests$Colonization <- adonis2(dist_matrix_1 ~ Colonization.Status, data = df_use, permutations = 9999)
var_tests$NATO <- adonis2(dist_matrix_1 ~ In.NATO, data = df_use, permutations = 9999)
var_tests$Continent <- adonis2(dist_matrix_1 ~ Continent, data = df_use, permutations = 9999)
var_tests$GDPperCap <- adonis2(dist_matrix_1 ~ GDPperCap, data = df_use, permutations = 9999)
var_tests$Mil_Per <- adonis2(dist_matrix_1 ~ Mil_Per, data = df_use, permutations = 9999)
var_tests$Edu_Per <- adonis2(dist_matrix_1 ~ Edu_Per, data = df_use, permutations = 9999)
var_tests$Age <- adonis2(dist_matrix_1 ~ Age, data = df_use, permutations = 9999)
var_tests$NatRes <- adonis2(dist_matrix_1 ~ NatResRents, data = df_use, permutations = 9999)

r2_values <- sapply(var_tests, function(x) x$R2[1])
p_values <- sapply(var_tests, function(x) x$`Pr(>F)`[1])

var_importance <- data.frame(
  Variable = names(r2_values),
  R2 = round(r2_values, 4),
  P_value = p_values,
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
##### (5) Visualize Global Clusters
##### ------------------------------------------------------

## ===== Mapping AT/CCAMLR Clusters =====
# load ISO3 codes
df$iso3 <- countrycode(df$Countries, 
                       origin = "country.name", 
                       destination = "iso3c")
# Note: Kosovo and Channel Islands cannot be matched to ISO3 codes and will
# appear as NA; these countries are too small to be visible at map resolution

# get world map with ISO codes
world_map <- map_data("world")
world_map$iso3 <- countrycode(world_map$region, 
                              origin = "country.name", 
                              destination = "iso3c")
# Note: some regions in the map data (e.g. island territories) cannot be matched
# to ISO3 codes; these are too small to be visible at map resolution

# combine
map_clusters <- world_map %>%
  left_join(dplyr::select(df, iso3, Cluster3), by = "iso3", relationship = "many-to-many")

# plot
ggplot(map_clusters, aes(x = long, y = lat, group = group, fill = Cluster3)) +
  geom_polygon(color = "white", linewidth = 0.1) +
  scale_fill_manual(values = c("1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A"),
                    na.value = "grey90",
                    name = "Cluster") +
  theme_void() +
  labs(title = "Global Distribution of AT/CCAMLR Participation Clusters") +
  coord_fixed(1.3)
