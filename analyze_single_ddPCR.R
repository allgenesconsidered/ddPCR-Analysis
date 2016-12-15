annotate_clusters_for_ASE <- function(cluster_centers){
  cluster_centers_temp = cluster_centers
  cluster_list = vector(mode='integer')
  
  cluster_list["FAM ++"] = cluster_centers_temp[which(cluster_centers_temp$Ch1.Amplitude == max(cluster_centers_temp$Ch1.Amplitude)),3]
  cluster_centers_temp = cluster_centers_temp[which(cluster_centers_temp$Cluster != cluster_list["FAM ++"]),]
  cluster_centers_temp$a = cluster_centers_temp$Ch1.Amplitude + cluster_centers_temp$Ch2.Amplitude
  
  cluster_list["FAM+/VIC+"] = cluster_centers_temp[which(cluster_centers_temp$a == max(cluster_centers_temp$a)),3]
  cluster_centers_temp = cluster_centers_temp[which(cluster_centers_temp$Cluster != cluster_list["FAM+/VIC+"]),]
  
  cluster_list["FAM +"] = cluster_centers_temp[which(cluster_centers_temp$Ch1.Amplitude == max(cluster_centers_temp$Ch1.Amplitude)),3]
  cluster_centers_temp = cluster_centers_temp[which(cluster_centers_temp$Cluster != cluster_list["FAM +"]),]
  
  cluster_list["VIC +"] = cluster_centers_temp[1,3]
  cluster_centers_ID = vector(mode='character', length = nrow(cluster_centers))
  
  for(row in 1:length(cluster_centers_ID)){
    cluster_centers_ID[row] = names(cluster_list)[which(cluster_list == cluster_centers[row,3])]
  }
  cluster_centers$ID = cluster_centers_ID
  return(cluster_centers)
}

annotate_dataset <- function(cluster_centers, dataset){
  dataset_ID = vector(mode='character', length = nrow(dataset))
  dataset_c1_center = vector(mode='numeric', length = nrow(dataset))
  dataset_c2_center = vector(mode='numeric', length = nrow(dataset))
  for(row in 1:length(dataset_ID)){
    dataset_ID[row] = cluster_centers$ID[which(cluster_centers$Cluster == dataset$Cluster[row])]
    dataset_c1_center[row] = cluster_centers$Ch1.Amplitude[which(cluster_centers$Cluster == dataset$Cluster[row])]
    dataset_c2_center[row] = cluster_centers$Ch2.Amplitude[which(cluster_centers$Cluster == dataset$Cluster[row])]
  }
  dataset$ID = dataset_ID
  dataset$Ch1.center = dataset_c1_center
  dataset$Ch2.center = dataset_c2_center
  return(dataset)
}

point_distance <- function(p1,p2,q1,q2){
  return(sqrt((q1 -p1)^2 + (q2 - p2)^2))
}

file = './data/2016.11.29/2016.11.29 Angela_C01_Amplitude.csv'
k = 4


require(ggplot2)
dat <- read.csv(file)

#dat = dat[!dat$Cluster == 1,]

for(col in c(1,2)){
  dat[,col] = (dat[,col] - min(dat[,col]))
  dat[,col] = (dat[,col]/max(dat[,col]))
}  

reclustered = kmeans(dat[,c(1,2)], k, nstart = 1000)
dat$Cluster = reclustered$cluster
#dat$Cluster = as.factor(dat$Cluster)

centers <- as.data.frame(reclustered$centers)
centers$Cluster <- c(1:k)

centers = annotate_clusters_for_ASE(centers)
dat = annotate_dataset(centers, dat)

dat_dist = vector(mode='numeric', length = nrow(dat))
for
for(row in 1:nrow(dat)){
  dat_dist[row] = point_distance(dat[row,1],dat[row,2],dat[row,5],dat[row,6])
}
dat$distance = dat_dist
for(row in 1:nrow(dat)){
  if(dat$distance[row] >= quantile(dat$distance, .95)){
    dat$ID[row] = 'Noise'
  }
}

ggplot() +
  geom_point(data = dat,aes(Ch2.Amplitude,Ch1.Amplitude, color = ID), alpha = 0.4) +
  geom_point(data = centers, aes(Ch2.Amplitude,Ch1.Amplitude), size = 3) +
  geom_point(data = centers, aes(Ch2.Amplitude,Ch1.Amplitude, color = ID)) +
  scale_x_continuous("Scaled VIC") +
  scale_y_continuous("Scaled FAM") +
  scale_fill_brewer(palette="Set2") +
  ggtitle(file) +
  theme_grey() 
#ggsave('ddPCR_test_RAW.png')

