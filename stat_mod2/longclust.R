

# Installs & Load ---------------------------------------------------------

install.packages("tidyverse")
install.packages("patchwork")
install.packages("mvtnorm")
install.packages("nlme")
install.packages("TSdist")
install.packages("dtwclust")
install.packages("stat_mod2/longclust_1.2.3.tar.gz") # need to be installed from the binary version

library(tidyverse)
library(patchwork)
library(mvtnorm)
library(dtwclust)
library(longclust)



# Custom functions --------------------------------------------------------

get_clusters <- function(cluster_score, threshold = .5) {
	
	n_clus <- ncol(cluster_score)
	n_obs <- nrow(cluster_score)
	res <- tibble::tibble(id = 1:n_obs, mb_clus = vector("numeric", n_obs))
	for (j in 1:n_clus) {
		cluster_score[, j] <- ifelse(cluster_score[, j] >= threshold, j, 0)
	}
	res$mb_clus <- rowSums(cluster_score)
	return(res)
	
}

log_standardize_vec <- function(x) {
	y <- log1p(x)
	y <- (y - mean(y)) / sd(y)
	return(y)
}

standardize_vec <- function(x) {
	y <- (x - mean(x)) / sd(x)
	return(y)
}



# Parameters --------------------------------------------------------------

# what method to use for choosing the optimal number of clusters in k-means and 
# hierarchical clustering
# possible values from cvi - type 
?cvi
method <- "D"



# Simulated Data ----------------------------------------------------------

# * Simulation ------------------------------------------------------------

# simulate time series with 6 observations each
nsim <- 5

m1 <- c(23, 34, 39, 45, 51, 56)
S1 <- matrix(c(1.00, -0.90, 0.18, -0.13, 0.10, -0.05, -0.90, 
							 1.31, -0.26, 0.18, -0.15, 0.07, 0.18, -0.26, 4.05, -2.84, 
							 2.27, -1.13, -0.13, 0.18, -2.84, 2.29, -1.83, 0.91, 0.10, 
							 -0.15, 2.27, -1.83, 3.46, -1.73, -0.05, 0.07, -1.13, 0.91, 
							 -1.73, 1.57), 6, 6)

m2 <- c(16, 18, 15, 17, 21, 17)
S2 <- matrix(c(1.00, 0.00, -0.50, -0.20, -0.20, 0.19, 0.00, 
							 2.00, 0.00, -1.20, -0.80, -0.36,-0.50, 0.00, 1.25, 0.10, 
							 -0.10, -0.39, -0.20, -1.20, 0.10, 2.76, 0.52, -1.22,-0.20, 
							 -0.80, -0.10, 0.52, 1.40, 0.17, 0.19, -0.36, -0.39, -1.22, 
							 0.17, 3.17), 6, 6)

m3 <- c(8, 11, 16, 22, 25, 28)
S3 <- matrix(c(1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 
							 1.00, -0.20, -0.64, 0.26, 0.00, 0.00, -0.20, 1.04, -0.17, 
							 -0.10, 0.00, 0.00, -0.64, -0.17, 1.50, -0.65, 0.00, 0.00, 
							 0.26, -0.10, -0.65, 1.32, 0.00, 0.00, 0.00, 0.00, 0.00, 
							 0.00, 1.00), 6, 6)

m4 <- c(12, 9, 8, 5, 4 ,2)
S4 <- diag(c(1, 1, 1, 1, 1, 1))

set.seed(123456)
data <- rbind(
	rmvnorm(nsim, m1, S1), rmvnorm(nsim, m2, S2),	
	rmvnorm(nsim, m3, S3), rmvnorm(nsim, m4, S4)
)


# * Visualization ---------------------------------------------------------

data_tbl <- data |> 
	t() |> 
	as_tibble() |> 
	set_names(1:ncol(t(data))) |> 
	mutate(time = 1:n(), .before = 1) |> 
	pivot_longer(-time) |> 
	mutate(name = as.integer(name)) |> 
	arrange(name, time)

data_tbl |> 
	ggplot(aes(x = time, y = value, group = name, col = name)) +
	geom_line(alpha = 0.8, show.legend = FALSE) +
	# facet_wrap(. ~ name) + 
	labs(title = "", y = "") +
	theme_bw()


# * Modelling -------------------------------------------------------------

data_mat <- apply(t(data), 2, standardize_vec) |> t()

# Model-based
mb_clus <- longclustEM(data_mat, 2, 6, gaussian = TRUE, criteria = "BIC")
mb_clus$Gbest
mb_clus$bicres # EVI (-6.540292)
summary(mb_clus)
get_clusters(mb_clus$zbest)

# K-medoids with DTW
km_dtw <- tsclust(data_mat, type = "partitional", k = 2L:6L, distance = "dtw", centroid = "pam")
best_km_dtw <- sapply(km_dtw, cvi, type = method) |> which.max() + 1
km_dtw <- tsclust(data_mat, type = "partitional", k = best_km_dtw, distance = "dtw", centroid = "pam")

# Hierarchical with DTW
hie_dtw <- tsclust(data_mat, type = "hierarchical", k = 2L:6L, distance = "dtw")
best_hie_dtw <- sapply(hie_dtw, cvi, type = method) |> which.max() + 1
hie_dtw <- tsclust(data_mat, type = "hierarchical", k = best_hie_dtw, distance = "dtw")

# Imposing the number of clusters
mb_clus <- longclustEM(data_mat, 4, 4, gaussian = TRUE, criteria = "BIC")
km_dtw <- tsclust(data_mat, type = "partitional", k = 4, distance = "dtw", centroid = "pam")
hie_dtw <- tsclust(data_mat, type = "hierarchical", k = 4, distance = "dtw")

clusters <- get_clusters(mb_clus$zbest) |> 
	mutate(km_clus = km_dtw@cluster, hie_clus = hie_dtw@cluster) |>
	mutate(
		mb_clus = as.factor(mb_clus),
		km_clus = as.factor(km_clus),
		hie_clus = as.factor(hie_clus)
	)
	
data_clus <- data_tbl |> 
	left_join(clusters, by = c("name" = "id"))

g_mb <- data_clus |> 
	ggplot(aes(x = time, y = value, group = name, col = mb_clus)) +
	geom_line(alpha = 0.5, show.legend = FALSE) +
	facet_wrap(. ~ mb_clus) + 
	labs(title = "Model-based (EVI)", y = "") +
	theme_bw()
g_km <- data_clus |> 
	ggplot(aes(x = time, y = value, group = name, col = km_clus)) +
	geom_line(alpha = 0.5, show.legend = FALSE) +
	facet_wrap(. ~ km_clus) + 
	labs(title = "K-medoids with DTW", y = "") +
	theme_bw()
g_hie <- data_clus |> 
	ggplot(aes(x = time, y = value, group = name, col = hie_clus)) +
	geom_line(alpha = 0.5, show.legend = FALSE) +
	facet_wrap(. ~ hie_clus) + 
	labs(title = "Hierarchical with DTW", y = "") +
	theme_bw()
g_mb / g_km / g_hie


# Replication Study - Rats Data -------------------------------------------

# * Data ------------------------------------------------------------------

# time series of rats' body weight for 4 different diets  
data_rats <- nlme::BodyWeight |> 
	as_tibble() |> 
	set_names(c("value", "time", "name", "diet")) |> 
	mutate(name = as.integer(name)) |> 
	select(name, time, value, diet) |> 
	arrange(name, time)


# * Visualization ---------------------------------------------------------

data_rats |> 
	ggplot(aes(x = time, y = value, group = name, col = diet)) + 
	geom_line() + 
	facet_wrap(. ~ diet) +
	labs(title = "Rats body weight over time by diet", y = "") +
	theme_bw()


# * Modelling -------------------------------------------------------------

# data in matrix form by row
data_rats_mat <- data_rats |> 
	select(-diet) |> 
	pivot_wider(names_from = name, values_from = value) |> 
	select(-time) |> 
	mutate(across(everything(), standardize_vec)) |> 
	as.matrix() |> 
	t()

# Model-based
mb_clus_rats <- longclustEM(data_rats_mat, 2, 6, gaussian = TRUE, criteria = "BIC")
mb_clus_rats$Gbest
mb_clus_rats$bicres # EVI (-172.6680)
summary(mb_clus_rats)
get_clusters(mb_clus_rats$zbest)

# K-medoids with DTW
km_dtw_rats <- tsclust(data_rats_mat, type = "partitional", k = 2L:6L, distance = "dtw", centroid = "pam")
best_km_dtw_rats <- sapply(km_dtw_rats, cvi, type = method) |> which.max() + 1
km_dtw_rats <- tsclust(data_rats_mat, type = "partitional", k = best_km_dtw_rats, distance = "dtw", centroid = "pam")

# Hierarchical with DTW
hie_dtw_rats <- tsclust(data_rats_mat, type = "hierarchical", k = 2L:6L, distance = "dtw")
best_hie_dtw_rats <- sapply(hie_dtw_rats, cvi, type = method) |> which.max() + 1
hie_dtw_rats <- tsclust(data_rats_mat, type = "hierarchical", k = best_hie_dtw_rats, distance = "dtw")

# Imposing the number of clusters
mb_clus_rats <- longclustEM(data_rats_mat, 3, 3, gaussian = TRUE, criteria = "BIC")
km_dtw_rats <- tsclust(data_rats_mat, type = "partitional", k = 3, distance = "dtw", centroid = "pam")
hie_dtw_rats <- tsclust(data_rats_mat, type = "hierarchical", k = 3, distance = "dtw")

clusters_rats <- get_clusters(mb_clus_rats$zbest) |> 
	mutate(km_clus = km_dtw_rats@cluster, hie_clus = hie_dtw_rats@cluster) |> 
	mutate(
		mb_clus = as.factor(mb_clus),
		km_clus = as.factor(km_clus),
		hie_clus = as.factor(hie_clus)
	)

data_clus_rats <- data_rats |> 
	left_join(clusters_rats, by = c("name" = "id"))

g_rats <- data_clus_rats |> 
	ggplot(aes(x = time, y = value, group = name, col = diet)) +
	geom_line(alpha = 0.5, show.legend = FALSE) +
	facet_wrap(. ~ diet) + 
	labs(title = "Real Groups - Diets", y = "") +
	theme_bw()
g_mb_rats <- data_clus_rats |> 
	ggplot(aes(x = time, y = value, group = name, col = mb_clus)) +
	geom_line(alpha = 0.5, show.legend = FALSE) +
	facet_wrap(. ~ mb_clus) + 
	labs(title = "Model-based (EVI)", y = "") +
	theme_bw()
g_km_rats <- data_clus_rats |> 
	ggplot(aes(x = time, y = value, group = name, col = km_clus)) +
	geom_line(alpha = 0.5, show.legend = FALSE) +
	facet_wrap(. ~ km_clus) + 
	labs(title = "K-medoids with DTW", y = "") +
	theme_bw()
g_hie_rats <- data_clus_rats |> 
	ggplot(aes(x = time, y = value, group = name, col = hie_clus)) +
	geom_line(alpha = 0.5, show.legend = FALSE) +
	facet_wrap(. ~ hie_clus) + 
	labs(title = "Hierarchical with DTW", y = "") +
	theme_bw()
g_rats / g_mb_rats 
g_rats / g_km_rats
g_rats / g_hie_rats



# Real Application - TSdist Data ------------------------------------------

# * Data ------------------------------------------------------------------

# dataset 1: 100 clustered time series from TSdist package
data("example.database2", package = "TSdist")

labels1 <- tibble(
	name = 1:length(example.database2$classes), 
	group = example.database2$classes
)

data1_tbl <- t(example.database2$data) |> 
	as_tibble() |> 
	set_names(labels1$name) |> 
	mutate(time = 1:n(), .before = 1) |> 
	slice(1:40) |> # first 40 observations otherwise the algorithm does not work
	pivot_longer(-time) |> 
	mutate(name = as.integer(name)) |> 
	left_join(labels1, by = "name") |> 
	filter(group != 6) |> 
	mutate(group = as.factor(group)) |> 
	arrange(name, time)

rm(example.database2)


# * Visualization ---------------------------------------------------------

idx1 <- labels1 |> group_by(group) |> slice(1) |> pull(name)
data1_tbl |> 
	filter(name %in% idx1) |> # 1 per cluster
	ggplot(aes(x = time, y = value, color = group)) +
	geom_line(alpha = 0.5, show.legend = FALSE) +
	facet_wrap(. ~ group) + 
	labs(title = "One sampled time series by cluster", y = "") +
	theme_bw()

data1_tbl |> 
	group_by(name) |> 
	ggplot(aes(x = time, y = value, group = name, col = group)) +
	geom_line(alpha = 0.5, show.legend = FALSE) +
	facet_wrap(. ~ group) + 
	labs(title = "85 time series by clusters", y = "") +
	theme_bw()


# * Modelling -------------------------------------------------------------

# in matrix and rowise format for clustering modelling
data1_mat <- data1_tbl |> 
	select(-group) |> 
	pivot_wider(names_from = name, values_from = value) |> 
	select(-time) |> 
	mutate(across(everything(), standardize_vec)) |> # standardization 
	as.matrix() |> 
	t()

# Model-based
mb_clus1 <- longclustEM(data1_mat, 2, 8, gaussian = TRUE, criteria = "BIC")
mb_clus1$Gbest
mb_clus1$bicres # EVI (-7242.576)
summary(mb_clus1)
get_clusters(mb_clus1$zbest)

# K-medoids with DTW
km_dtw1 <- tsclust(data1_mat, type = "partitional", k = 2L:8L, distance = "dtw", centroid = "pam")
best_km_dtw1 <- sapply(km_dtw1, cvi, type = method) |> which.max() + 1
km_dtw1 <- tsclust(data1_mat, type = "partitional", k = best_km_dtw1, distance = "dtw", centroid = "pam")

# Hierarchical with DTW
hie_dtw1 <- tsclust(data1_mat, type = "hierarchical", k = 2L:8L, distance = "dtw")
best_hie_dtw1 <- sapply(hie_dtw1, cvi, type = method) |> which.max() + 1
hie_dtw1 <- tsclust(data1_mat, type = "hierarchical", k = best_hie_dtw1, distance = "dtw")

# Imposing the number of clusters
mb_clus1 <- longclustEM(data1_mat, 5, 5, gaussian = TRUE, criteria = "BIC")
km_dtw1 <- tsclust(data1_mat, type = "partitional", k = 5, distance = "dtw", centroid = "pam")
hie_dtw1 <- tsclust(data1_mat, type = "hierarchical", k = 5, distance = "dtw")

clusters1 <- get_clusters(mb_clus1$zbest) |> 
	mutate(km_clus = km_dtw1@cluster, hie_clus = hie_dtw1@cluster) |> 
	mutate(
		mb_clus = as.factor(mb_clus),
		km_clus = as.factor(km_clus),
		hie_clus = as.factor(hie_clus)
	)

data_clus1 <- data1_tbl |> 
	left_join(clusters1, by = c("name" = "id"))

g1 <- data_clus1 |> 
	ggplot(aes(x = time, y = value, group = name, col = group)) +
	geom_line(alpha = 0.5, show.legend = FALSE) +
	facet_wrap(. ~ group) + 
	labs(title = "Real Groups", y = "") +
	theme_bw()
g_mb1 <- data_clus1 |> 
	ggplot(aes(x = time, y = value, group = name, col = mb_clus)) +
	geom_line(alpha = 0.5, show.legend = FALSE) +
	facet_wrap(. ~ mb_clus) + 
	labs(title = "Model-based (EVI)", y = "") +
	theme_bw()
g_km1 <- data_clus1 |> 
	ggplot(aes(x = time, y = value, group = name, col = km_clus)) +
	geom_line(alpha = 0.5, show.legend = FALSE) +
	facet_wrap(. ~ km_clus) + 
	labs(title = "K-medoids with DTW", y = "") +
	theme_bw()
g_hie1 <- data_clus1 |> 
	ggplot(aes(x = time, y = value, group = name, col = hie_clus)) +
	geom_line(alpha = 0.5, show.legend = FALSE) +
	facet_wrap(. ~ hie_clus) + 
	labs(title = "Hierarchical with DTW", y = "") +
	theme_bw()
g1 / g_mb1 
g1 / g_km1
g1 / g_hie1

