# read_ddpcr_data.R for analysis of Allele-specific Editing using droplet digital PCR
# (c) 2016 Michael Olvera
# micahel.olvera@gladstone.ucsf.edu


#' Collects raw .csv data from each droplet and returns the merged dataframe, including info on
#' clusters. Inputs include:
#' path - string of path to the filder holding all the raw droplet data (in .csv format).
#' remove_empties - Logical whether or not to remove the empty droplets. Default is True.
#' recluster - should the algorithm attempt to recluster the data. Default is True.
#' k - if reclustering, how many clusters are you trying for.
#'
#' The algorithm works by:
#'   * importing all the .csv files from the directory given in the path. This loads
#' in the well ID as well.
#'   For each well:
#'     * empty wells are removed from the analysis. Right now this is done using the
#'     emty wells as determined by the machine, however it can be changed later.
#'     * the axis are normalized (seting the maximum to 1 and minimum to 0).
#'     * with normalized axis, the data is reclustered
#'     * the data is then annotated (done only for current experiment).
#'     
#' Imports come from the external funtion \code{read_ddpcr_data}.
parse_files_from_path <- function(path, remove_empties = T, recluster = T, k){
  files = list.files(path, pattern = '.csv', full.names = F, recursive = F)
  full_plate = data.frame()
  for (file in files) {
    importedFile = parseFileHelper(path, file, remove_empties, recluster, k)
    if(nrow(full_plate)==0){
      full_plate = importedFile
    } else { full_plate = rbind(full_plate,importedFile) }
  }
  cat(paste(length(files), 'files merged together.\n'))
  return(full_plate)
}

parse_files_from_list <- function(list_of_paths, remove_empties = T, recluster = T, k, 
                                  fileNames){
  full_plate = data.frame()
  for (i in 1:length(list_of_paths)) {
    importedFile = parseFileHelper(path = "", list_of_paths[i], remove_empties, recluster, k,
                                   fileName = fileNames[i])
    if(nrow(full_plate)==0){
      full_plate = importedFile
    } else { full_plate = rbind(full_plate,importedFile) }
  }
  
  cat(paste(length(list_of_paths), 'files merged together.\n'))
  return(full_plate)
}

#' Helper funciton for parseFiles. Should not be used by itself.
parseFileHelper <- function(path, file, rmv_emp, recluster, k, fileName = file){
  importedFile = read.csv(paste(path,file, sep = ""), stringsAsFactors = F)
  if(rmv_emp){ # Revome empty droplets
    importedFile = importedFile[!importedFile$Cluster == 1,]
  }
  for(col in c(1,2)){ # Maximum normalize the data
    importedFile[,col] = maximumNormalize(importedFile[,col])
  }
  if(recluster){ # kmeans reclustering
    importedFile_and_centers = recluster(importedFile, k)
    
    centers = annotate_clusters_for_ASE(as.data.frame(importedFile_and_centers[2]))
    importedFile = annotate_dataset(centers, as.data.frame(importedFile_and_centers[1]))
  }
  
  importedFile$Well <- findWellnumber(fileName)
  return(importedFile)
}

#' Extracts the well number from the experimental .csv file name. This .csv should
#' come from the ddPCR BioRad software and should not be editied.  
#' Assumes:
#'     The well number is flanked by '_' character.
#'     The well matches the pattern 'alphabet-number-number' (ex. A01)
#' @param file The string form of the filename.
#' @return The string corresponding to the well number.
findWellnumber <- function(file){
  split = unlist(strsplit(file, '_'))
  return(grep('[A-Z][0-9]{2}',split, value = T)) # Return well number (ex. "A01")
}

#' Handles reclustering of ddPCR-assigned clusters via k-means. The k is defined
#' by the user ahead of time.
#' @TODO Allow for muliple k's for different wells.
#' @param dataset The dataset to recluster. Should be the ddPCR plate
#' @parap k The numer of clusters to find.
#' @return A list containing the reclustered dataset
recluster <- function(dataset, k){
  reclustered = kmeans(dataset[,c(1,2)], k, nstart = 100)
  dataset$Cluster = reclustered$cluster
  
  centers <- as.data.frame(reclustered$centers)
  centers$Cluster <- c(1:k)
  return(list(dataset, centers))
}

#' Maximum normalize each intensity read. The logic being that ddPCR intensity
#' values are inherently arbitrary. Normalizing the data so that the most
#' intense values are 
#' x_normalized = (x/max-min) - min
maximumNormalize <- function(column){
  # Maximum normalize each intensity read. 
  # x_normalized = (x/max-min) - min

  column = (column - min(column))
  column = (column/max(column))
  return(column)
}

### Automatic Annotation ###

#' Function specific for Angela's ddPCR experiments. Annotates each clused based on where you
#' would predic them to fall relative to eachother. For instance, it is assumed that the 
#' 'FAM++' poplation will be the population with the largest Ch1 amplitude. All other groups can
#' be extrapolated based on these properties.
#' 
#' The funciton relies on only the cluster centers returned from the kmeans function. Th
annotate_clusters_for_ASE <- function(cluster_centers){
  cluster_centers_temp = cluster_centers
  cluster_list = vector(mode='integer')
  
  cluster_list["FAM ++"] = cluster_centers_temp[which(cluster_centers_temp$Ch1.Amplitude == max(cluster_centers_temp$Ch1.Amplitude)),3]
  cluster_centers_temp = cluster_centers_temp[which(cluster_centers_temp$Cluster != cluster_list["FAM ++"]),]
  
  cluster_list["VIC +"] = cluster_centers_temp[which(cluster_centers_temp$Ch1.Amplitude == min(cluster_centers_temp$Ch1.Amplitude)),3]
  cluster_centers_temp = cluster_centers_temp[which(cluster_centers_temp$Cluster != cluster_list["VIC +"]),]
  
  cluster_list["FAM +"] = cluster_centers_temp[which(cluster_centers_temp$Ch2.Amplitude == min(cluster_centers_temp$Ch2.Amplitude)),3]
  cluster_centers_temp = cluster_centers_temp[which(cluster_centers_temp$Cluster != cluster_list["FAM +"]),]
  
  cluster_list["FAM+/VIC+"] = cluster_centers_temp[which(cluster_centers_temp$Ch1.Amplitude == min(cluster_centers_temp$Ch1.Amplitude)),3]
  cluster_centers_temp = cluster_centers_temp[which(cluster_centers_temp$Cluster != cluster_list["FAM+/VIC+"]),]
  
  if(nrow(cluster_centers_temp) != 0){
    for( i in nrow(cluster_centers_temp)){
      cluster_list[paste("Unknown", i)] = cluster_centers_temp[i,3]
    }
  }
  
  cluster_centers_ID = vector(mode='character', length = nrow(cluster_centers))
  for(row in 1:length(cluster_centers_ID)){
    cluster_centers_ID[row] = names(cluster_list)[which(cluster_list == cluster_centers[row,3])]
  }
  cluster_centers$ID = cluster_centers_ID
  return(cluster_centers)
}

annotate_dataset <- function(cluster_centers, dataset){
  dataset_ID = vector(mode='character', length = nrow(dataset))
  for(row in 1:length(dataset_ID)){
    dataset_ID[row] = cluster_centers$ID[which(cluster_centers$Cluster == dataset$Cluster[row])]
  }
  dataset$ID = dataset_ID
  return(dataset)
}

###

processPlatemap <- function(platemap_file){
  # Read in the metadata of the experiment from the platemap, including
  # # of replicates and experimental condition.
  platemap <- read.csv(platemap_file, stringsAsFactors = F)
  platemap <- seperateWellNumbers(platemap)
  return(platemap)
}


plotFacet <- function(full_plate, platemap = NULL, path = ""){
  # Plot the data as a facet plot together, seperating the plots by row letter
  # and column number as displayed on the plate.
  full_plate <- seperateWellNumbers(full_plate)
  
  require(ggplot2) 
  full_plate$Cluster <- as.factor(full_plate$Cluster)
  g = ggplot(full_plate, aes(x=Ch2.Amplitude, y=Ch1.Amplitude), environment = environment()) + 
    geom_point(shape=1, size=1, aes(color = ID)) + 
    facet_grid(Letter ~ Number, switch = 'y') +
    scale_x_continuous("VIC Intensity (AU)") +
    scale_y_continuous("FAM Intensity (AU)") +
    ggtitle(path)
  if(!is.null(platemap)){
    g =  g + geom_text(data = platemap, aes(x = 0.5, y =0.9, label = Condition, group = NULL))
  }
  g = g +
    theme_bw() +
    theme(axis.ticks = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_blank(), strip.text.y = element_text(angle = 180))
  return(g)
}

seperateWellNumbers <- function(dataset){
  # With the input dataset, takes the 'Well' column and seperates the number and the letter. 
  # Argument: dataset, with a column names 'Well'.
  # Returns: the dataset, with two new columns.
  # Assumes:
  #     The well matches the pattern 'alphabet-number-number' (ex. A01)
  try(if(!'Well' %in% colnames(dataset)) stop('No column named \"Well\" in dataset.'))
  dataset$Letter <- substr(dataset$Well,1,1)
  dataset$Number <- substr(dataset$Well,2,3)
  return(dataset)
}

cleanPlate <- function(plate, platemap){
  # From the raw data (every point a row), return a data.frame with every row being a
  # well, and the sum of every cluster type (ie # of droplets per condition).
  
  require(reshape2)
  condensed = as.data.frame(table(plate$Well, plate$ID))
  condensed = as.data.frame(acast(condensed, Var1~Var2, value.var="Freq"))
  condensed$Well = rownames(condensed)
  merged = merge(condensed,platemap, by = 'Well')
  
  return(merged[ , -which(names(merged) %in% c('Letter','Number'))])
}

read_ddpcr_data <- function(path_to_data_folder, path_to_platemap, output_dir = '~/Desktop/',
                            remove_empties = T, recluster = T,k = 4, graph_data = T, output = F){
  
  testPlate <- parse_files_from_path(path_to_data_folder, k = k)
  platemap <- processPlatemap(path_to_platemap)
  outputPlate <- cleanPlate(testPlate, platemap)
  if(graph_data){
    plotFacet(testPlate, platemap, path_to_data_folder)
    #ggsave('facet_11_29_16_Angela.png')
  }
  if(output){
    write.csv(outputPlate,paste(output_dir,basename(path_to_data_folder),'_autoclustered.csv', sep = ''), row.names = F)
  }
  cat('Plate read sucessfully!')
  return(outputPlate)
}


#read_ddpcr_data('./data/2016.11.30/','./data/2016_11_30_platemap.csv', k=5)


