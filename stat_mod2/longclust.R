

# Installs & Load ---------------------------------------------------------

install.packages("tidyverse")
install.packages("patchwork")
install.packages("mvtnorm")
install.packages("nlme")
install.packages("TSdist")
install.packages("stat_mod2/longclust_1.2.3.tar.gz")

library(tidyverse)
library(patchwork)
library(mvtnorm)
library(longclust)



# Custom functions --------------------------------------------------------

get_clusters <- function(cluster_score, threshold = .5) {
	
	n_clus <- ncol(cluster_score)
	n_obs <- nrow(cluster_score)
	res <- tibble::tibble(id = 1:n_obs, cluster = vector("numeric", n_obs))
	for (j in 1:n_clus) {
		cluster_score[, j] <- ifelse(cluster_score[, j] >= threshold, j, 0)
	}
	res$cluster <- rowSums(cluster_score)
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
	ggplot(aes(x = time, y = value, group = name, color = name)) +
	geom_line(alpha = 0.8, show.legend = FALSE) +
	facet_wrap(. ~ name) + 
	labs(title = "Simulated time series", y = "") +
	theme_bw()


# * Modelling -------------------------------------------------------------

clus <- longclustEM(data, 2, 6, gaussian = TRUE, criteria = "BIC")

clus$Gbest
clus$bicres # EEI

summary(clus)

clus_res <- get_clusters(clus$zbest)

data_clus <- data_tbl |> 
	left_join(clus_res, by = c("name" = "id")) |> 
	mutate(cluster = as.factor(cluster))

data_clus |> 
	ggplot(aes(x = time, y = value, group = name, color = cluster)) +
	geom_line(alpha = 0.5, show.legend = FALSE) +
	facet_wrap(. ~ cluster) + 
	labs(title = "Simulated time series", y = "") +
	theme_bw()



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
	labs(title = "Rats body weight over time by diet", y = "") +
	theme_bw()


# * Modelling -------------------------------------------------------------

# data in matrix form by row
data_rats_mat <- data_rats |> 
	select(-diet) |> 
	pivot_wider(names_from = name, values_from = value) |> 
	select(-time) |> 
	as.matrix() |> 
	t()

clus <- longclustEM(data_rats_mat, 2, 6, gaussian = TRUE)

clus$Gbest
clus$bicres # EEA

clus_res <- get_clusters(clus$zbest)

data_rats_clus <- data_rats |> 
	left_join(clus_res, by = c("name" = "id"))

data_rats_clus |> 
	mutate(diet = as.numeric(diet)) |> 
	pivot_longer(cols = c(diet, cluster), names_to = "type", values_to = "group_value") |> 
	mutate(group_value = as.factor(group_value)) |> 
	ggplot(aes(x = time, y = value, group = name, col = group_value)) + 
	geom_line() + 
	facet_wrap(~ type) +
	labs(y = "") +
	theme_bw()



# Real Application - TSdist Data ------------------------------------------

# * Data ------------------------------------------------------------------

# dataset 1: 100 clustered time series from TSdist package
data("example.database2", package = "TSdist")

labels1 <- tibble(
	name = 1:length(example.database2$classes), 
	group = example.database2$classes
)

data1 <- t(example.database2$data) |> 
	as_tibble() |> 
	set_names(labels1$name) |> 
	mutate(across(everything(), log_standardize_vec)) |> # log + standardization
	mutate(time = 1:n(), .before = 1) |> 
	slice(1:40) # first 40 observations otherwise the algorithm does not work

data1_tbl <- data1 |> 
	pivot_longer(-time) |> 
	mutate(name = as.integer(name)) |> 
	left_join(labels1, by = "name") |> 
	arrange(name, time)

rm(example.database2)


# * Visualization ---------------------------------------------------------

idx1 <- labels1 |> group_by(group) |> slice(1) |> pull(name)
data1_tbl |> 
	filter(name %in% idx1) |> # 1 per cluster
	ggplot(aes(x = time, y = value, color = name)) +
	geom_line(alpha = 0.5, show.legend = FALSE) +
	facet_wrap(. ~ group) + 
	labs(title = "Dataset 1: one sampled time series by cluster", y = "") +
	theme_bw()

data1_tbl |> 
	group_by(name) |> 
	ggplot(aes(x = time, y = value, group = name, col = name)) +
	geom_line(alpha = 0.5, show.legend = FALSE) +
	facet_wrap(. ~ group) + 
	labs(title = "Dataset 1: 100 time series by cluster", y = "") +
	theme_bw()


# * Modelling -------------------------------------------------------------

# in matrix and rowise format for clustering modelling
data1_mat <- data1 |> select(-time) |> as.matrix() |> t()

clus1 <- longclustEM(data1_mat, 2, 8, gaussian = TRUE)

clus1$Gbest
clus1$bicres # EVI

clus1_res <- get_clusters(clus1$zbest)

data1_clus <- data1_tbl |> 
	left_join(clus1_res, by = c("name" = "id"))

g1 <- data1_clus |> 
	mutate(group = as.factor(group)) |> 
	ggplot(aes(x = time, y = value, group = name, col = group)) + 
	geom_line() + 
	facet_wrap(. ~ group) +
	labs(title = "Real Clusters", y = "") +
	theme_bw()
g2 <- data1_clus |> 
	mutate(cluster = as.factor(cluster)) |> 
	ggplot(aes(x = time, y = value, group = name, col = cluster)) + 
	geom_line() + 
	facet_wrap(. ~ cluster) +
	labs(title = "Estimated Clusters", y = "") +
	theme_bw()
g1 / g2
