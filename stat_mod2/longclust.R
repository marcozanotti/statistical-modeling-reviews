

# Installs & Load ---------------------------------------------------------

install.packages("tidyverse")
install.packages("mvtnorm")
install.packages("TSdist")
install.packages("stat_mod2/longclust_1.2.3.tar.gz")

library(tidyverse)
library(mvtnorm)
library(TSdist)
library(longclust)



# Import & Cleaning -------------------------------------------------------

# dataset 1: 100 clustered time series 
data("example.database2") # get example data from TSdist package
labels1 <- tibble(
	name = 1:length(example.database2$classes), 
	cluster = example.database2$classes
)
data1 <- t(example.database2$data) |> 
	as_tibble() |> 
	set_names(labels1$name) |> 
	mutate(time = 1:n(), .before = 1)

# dataset 2: 50 clustered time series
data("example.database3") # get example data from TSdist package
labels2 <- tibble(
	name = 1:length(example.database3$classes), 
	cluster = example.database3$classes
)
data2 <- t(example.database3$data) |> 
	as_tibble() |> 
	set_names(labels2$name) |> 
	mutate(time = 1:n(), .before = 1)

# in matrix and rowise format for clustering modelling
data1_clus <- example.database2$data
data2_clus <- example.database3$data
rm(example.database2)
rm(example.database3)



# Visualization -----------------------------------------------------------

idx1 <- labels1 |> group_by(cluster) |> slice(1) |> pull(name)
data1 |> 
	pivot_longer(-time) |> 
	mutate(name = as.integer(name)) |> 
	left_join(labels1, by = "name") |> 
	filter(name %in% idx1) |> # 1 per cluster
	ggplot(aes(x = time, y = value, group = cluster, color = name)) +
	geom_line(alpha = 0.5, show.legend = FALSE) +
	facet_wrap(. ~ cluster) + 
	labs(title = "Dataset 1: 1 time series by cluster", y = "") +
	theme_bw()

data1 |> 
	pivot_longer(-time) |> 
	mutate(name = as.integer(name)) |> 
	left_join(labels1, by = "name") |> 
	ggplot(aes(x = time, y = value, group = cluster, color = name)) +
	geom_line(alpha = 0.5, show.legend = FALSE) +
	geom_smooth(color = "red", show.legend = FALSE) +
	facet_wrap(. ~ cluster) + 
	labs(title = "Dataset 1: 100 time series by cluster", y = "") +
	theme_bw()


idx2 <- labels2 |> group_by(cluster) |> slice(1) |> pull(name)
data2 |> 
	pivot_longer(-time) |> 
	mutate(name = as.integer(name)) |> 
	left_join(labels2, by = "name") |> 
	filter(name %in% idx2) |> # 1 per cluster
	ggplot(aes(x = time, y = value, group = cluster, color = name)) +
	geom_line(alpha = 0.5, show.legend = FALSE) +
	facet_wrap(. ~ cluster) + 
	labs(title = "Dataset 2: 1 time series by cluster", y = "") +
	theme_bw()

data2 |> 
	pivot_longer(-time) |> 
	mutate(name = as.integer(name)) |> 
	left_join(labels2, by = "name") |> 
	ggplot(aes(x = time, y = value, group = cluster, color = name)) +
	geom_line(alpha = 0.5, show.legend = FALSE) +
	geom_smooth(color = "red", show.legend = FALSE) +
	facet_wrap(. ~ cluster) + 
	labs(title = "Dataset 2: 50 time series by cluster", y = "") +
	theme_bw()



# Modelling ---------------------------------------------------------------

?longclustEM

clus1 <- longclustEM(t(data1_clus), 2, 8, gaussian = TRUE)
clus1$Gbest
clus1$bicres
clus1$zbest |> round(2)
plot(clus1, t(data1_clus))

clus2 <- longclustEM(t(data2_clus), 2, 8, gaussian = TRUE)
clus2$Gbest
clus2$bicres
clus2$zbest |> round(2)
plot(clus2, t(data2_clus))






# Simulation --------------------------------------------------------------

# mu and sigma of 4 different multivariate normal distributions
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

data <- rbind(
	rmvnorm(25, m1, S1), rmvnorm(25, m2, S2),	
	rmvnorm(25, m3, S3), rmvnorm(25, m4, S4)
)



# Visualization -----------------------------------------------------------

data |> 
	as_tibble() |> 
	set_names(1:6) |> 
	mutate(time = 1:n(), .before = 1) |> 
	pivot_longer(-time) |> 
	mutate(name = as.integer(name)) |> 
	arrange(name, time) |> 
	ggplot(aes(x = time, y = value, group = name, color = name)) +
	geom_line(alpha = 0.5, show.legend = FALSE) +
	facet_wrap(. ~ name) + 
	labs(title = "Simulated time series", y = "") +
	theme_bw()



# Modelling ---------------------------------------------------------------

clus <- longclustEM(data, 3, 5, gaussian = TRUE)

clus$Gbest
clus$bicres
clus$zbest |> round(2)

summary(clus)
plot(clus, data)
print(clus)

