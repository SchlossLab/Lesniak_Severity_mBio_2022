pack_used <- c('dplyr','Boruta')
for (dep in pack_used){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), repos = 'http://cran.us.r-project.org', 
                     quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}

meta_file   <- '~/Documents/Github/Schubert_humanCdGF_XXXX_2016/data/process/human_CdGF_metadata.txt'
shared_file <- '~/Documents/Github/Schubert_humanCdGF_XXXX_2016/data/process/human_CdGF.an.unique_list.0.03.subsample.shared'

# read in files
meta_file   <- read.table(meta_file, sep = '\t', header = T, row.names = 'sample_id')
shared_file <- read.table(shared_file, sep = '\t', header = T, row.names = 'sample_id')

# select only mice with 431 cdiff
meta_file <- meta_file[meta_file$cdiff_strain == '431', ]

# subset shared
#Create df with relative abundances
stool_shared <- shared_file[rownames(meta_file[meta_file$day %in% 0, ]), ]
# use %in% instead of == because to ignore NA
# rename day 0 community rows by mouse id - since meta rows have day in row names, use mouse id for both
rownames(stool_shared) <- meta_file$mouse_id[rownames(meta_file) %in% 
                                               rownames(stool_shared) ]
rel_abund_stool <- 100 * stool_shared / unique(apply(stool_shared, 1, sum))
#Create vector of OTUs with median abundances >1%
OTUs_stool_0_01 <- apply(rel_abund_stool, 2, max) > 1
#df of OTUs with abundances >1% - by cage and inoculum
rel_abund_stool <- rel_abund_stool[ ,OTUs_stool_0_01]

#random forest regression model to determine Day 0 OTU correlated with increased CFU
seed <- -592166825
n_trees <- 2001
iters <- 100

meta_file$Euth_Early <- factor(ifelse(meta_file[ , 'Early_Euth'] == FALSE, 1, 0)) 

Predict_df <- merge(meta_file[meta_file$day %in% 0, 
                              c('cage_id', 'mouse_id', 'Euth_Early')],
                    rel_abund_stool, by.x = 'mouse_id', by.y = 'row.names')

Predict_early_euth_df <- select(Predict_df, -cage_id, -mouse_id)

#create dataframe for RF-classification
Validate_early_euth_df <- Predict_early_euth_df
set.seed(seed)
Validate_early_euth_df$Euth_Early <- sample(Validate_early_euth_df$Euth_Early)

# boruta feature selection
random_model_features <- c()

for(i in 1:iters){
  set.seed(seed)
  Boruta(Euth_Early ~ ., data = Validate_early_euth_df)
  feature_df <- attStats(random_bor)
  feature_df$OTU <- rownames(feature_df)
  random_model_features[[i]] <- feature_df
}
random_model_features <- do.call('rbind', random_model_features)

write.table(random_model_features, file = '~/Documents/Github/Schubert_humanCdGF_XXXX_2016/data/process/random_rf_features_boruta.txt', sep = '\t')
