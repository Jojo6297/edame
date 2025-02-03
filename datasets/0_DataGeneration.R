# import libraries
library(curatedMetagenomicData)
library(dplyr)


# set output directory
dir <- getwd()


# import all the studies and data related to them
SampleInfo <- data.frame(sampleMetadata)


# check for which body site maximum samples are available -> result: stool samples
table(SampleInfo$body_site, useNA = 'always')


# filter the data for only stool samples
stool_samp <- SampleInfo[which(SampleInfo$body_site == "stool"),]
stool_samp <- Filter(function(x)!all(is.na(x)), stool_samp)


# Find the different types of diseases available in the data
stool_tab1 <- table(stool_samp$disease, useNA = 'always')
stool_tab1 <- data.frame(stool_tab1[order(stool_tab1,       # Decreasing order of table
                                decreasing = TRUE)])

## Based of sample size (high) and relationship of diseases, chosen diseases are: IBD and CRC, and control: Healthy state.
disease_states <- c('healthy', 'CRC', 'IBD')

## selected studies
h_studies <- list("LifeLinesDeep_2016", "AsnicarF_2021", "MehtaRS_2018", "ZeeviD_2015") #healthy
c_studies <- list("YachidaS_2019", "WirbelJ_2018", "ZellerG_2014", "VogtmannE_2016")  #CRC
i_studies <- list("HMP_2019_ibdmdb", "VilaAV_2018", "HallAB_2017", "NielsenHB_2014") #IBD
studies <- list(h_studies, c_studies, i_studies)

# retrieve disease datasets
for (j in 1:3){
  disease_state <- disease_states[j]
  sortedStudies <- studies[[j]]
    
  # For each chosen study export the abundance counts and phylo.trees. 
  for (i in 1:length(sortedStudies)){
    stool.disease.study <-
      filter(SampleInfo) |>
      filter(body_site == "stool") |>
      filter(study_name == sortedStudies[[i]]) |>
      filter(disease == disease_state) |>
      select(where(~ !all(is.na(.x)))) |>
      returnSamples("relative_abundance", rownames = "short")
    
    disease.counts <- stool.disease.study@assays@data@listData[["relative_abundance"]]
    disease.phylo.tree <- data.frame(stool.disease.study@rowLinks@listData[["nodeLab"]])
    #disease.row.names <- data.frame(stool.disease.study@rowLinks@rownames)
    #disease.ColData <- data.frame(colData(stool.disease.study))
    
    write.csv(disease.counts, file = paste0(dir,i, '_', sortedStudies[i], "_", disease_state, ".relative.abundance.csv"), row.names = TRUE)
    write.csv(disease.phylo.tree, file = paste0(dir,i, '_', sortedStudies[i], "_", disease_state, ".phylo.tree.csv"), row.names = TRUE)
    #write.csv(disease.row.names, file = paste0(dir,i, '_', sortedStudies[i], "_disease.row.names.csv"), row.names = TRUE)
    #write.csv(disease.ColData, file = paste0(dir,i, '_', sortedStudies[i], "_disease.ColData.csv"), row.names = TRUE)
  }
}