# ======================================{ + }====================================== #
#                                                                                   #
#       An assessment about of prognostic immunity markers in breast cancer
#
#                                                                                   #
# ======================================{ + }====================================== #




# Plan when analysisng dataset:
# EDA: Explanatory data analysis
#     --> Graphs (simple ones)
#     --> Build hypothesis
#     --> Test with models




# =============================================================================================


# Libraries
library(tidyverse)    # Data manipulation and visualization
library(rpart)        # To define tree models
library(rpart.plot)   # To plot them
library(randomForest) # Random forest for classification trees
library(corrplot)     # Correlation plot


# =============================================================================================



          # ==== PART I: First observations ==== #
          # ==================================== #


# Time processing study
median(GSE113863_data_analysis.csv$DMFS_time..YEARS., na.rm = TRUE) # Time median of the study

# Distribution of breast cancer subtypes 
table(GSE113863_data_analysis.csv$iRDM_subtype)
barplot(table(GSE113863_data_analysis.csv$iRDM_subtype))

table(GSE113863_data_analysis.csv$PAM50_subtype)
barplot(table(GSE113863_data_analysis.csv$PAM50_subtype))

df_subtypes <- GSE113863_data_analysis.csv %>%
  pivot_longer(
    cols = c(iRDM_subtype, PAM50_subtype),
    names_to = "Evaluation_type",
    values_to = "Subtype") %>%
  select(Evaluation_type, Subtype)

ggplot(data = df_subtypes, aes(x = Subtype, fill = Evaluation_type)) +
  geom_bar(position = 'dodge') +
  theme_light()

# Find most important genes (more and less expressed)
# Correlation between groups


          # ==== PART II: Hypothesis and protocol ==== #
          # ========================================== #


# Hypothesis: Differences in gene expression according to different subtypes of breast cancer



          # ==== PART III: Classification & Prediction ==== #
          # =============================================== #

# ---- Classification model: Predicting subtypes ~ genes ---- #

## Prepare dataframe
# Fuse dataframes
GSE113863_dwd_normalized_std_expression.csv <- GSE113863_dwd_normalized_std_expression.csv %>%  # Rename column to same name
  rename(Sample_name = CLID)
GSE113863 <- GSE113863_data_analysis.csv %>%                                                    # Merge dataframes
  right_join(GSE113863_dwd_normalized_std_expression.csv, by = join_by(Sample_name))

# Cleanse dataframe
str(GSE113863)                          # Observe dataframe
sapply(GSE113863[,14:85], is.numeric)   # Verify that every gene expression variables are numeric
summary(GSE113863[,14:85])              # Description

verif_gene_expression <- GSE113863 %>%  # Select only gene expression
  pivot_longer(
    cols = CD27:TFRC,
    names_to = "Genes",
    values_to = "Gene_expression"
  ) %>%
  select(Genes, Gene_expression)

# I'm not sure that was a good idea, maybe make a heatmap with gene_expression over patient and/or group with other variables
ggplot(data = verif_gene_expression, aes(x = Genes, y = Gene_expression, fill = Genes)) +
  geom_boxplot() +
  theme_light() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    legend.position = "none")





#### Classification model for each evaluation (iRDM // PAM50)

## (1) Classification: response_variable = iRDM


# Select data iRDM against gene expression
Class_iRDM_df <- GSE113863 %>%   # Select only iRDM//genes corresponding
  dplyr::select(iRDM_subtype, 14:85)

# Factorize the response_variable VARIABLE DE FDP QUI M4A BOUSILLE 3H DE MA VIE
Class_iRDM_df$iRDM_subtype <- factor(Class_iRDM_df$iRDM_subtype)


# ///////////////////////////////////////////////////////////////////////////  #
# 1. Distribution des classes (CAUSE #1)
print("=== DISTRIBUTION iRDM_subtype ===")
print(table(Class_iRDM_df$iRDM_subtype, useNA = "always"))

# 2. Taille du dataset
print(paste("Nombre d'échantillons:", nrow(Class_iRDM_df)))

# 3. NA dans la réponse
print(paste("NA dans iRDM_subtype:", sum(is.na(Class_iRDM_df$iRDM_subtype))))
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #


# Separate a dataset in a train and test subsets
sample_iRDM <- sample(c(TRUE, FALSE), nrow(Class_iRDM_df), replace=TRUE, prob=c(0.7,0.3))

train  <- Class_iRDM_df[sample_iRDM, ]
test   <- Class_iRDM_df[!sample_iRDM, ]

# Build Class tree
max_tree_iRDM = rpart(iRDM_subtype ~ ., data = Class_iRDM_df, method = "class", cp = 0) # CP needs to be low to have a complex tree

# Plot Tree
plotcp(max_tree_iRDM)  # Represents the cross-validation error as a function of cp
printcp(max_tree_iRDM) # Table version of cross-validation errors as a function of cp

# Determines the value of cp (iRDM) for the lowest cross-validation error (specific gene)
cpOpt = max_tree_iRDM$cptable[which.min(max_tree_iRDM$cptable[,4]),1]
cpOpt = max_tree_iRDM$cptable[which.min(max_tree_iRDM$cptable[,"xerror"]), "CP"] # Answer from forum and IA (same effect)
print(paste("CP optimal:", round(cpOpt, 6)))

# The best tree
pruned_tree_iRDM = rpart(
  iRDM_subtype ~ ., 
  data = Class_iRDM_df, 
  method = "class",
  cp = cpOpt)

# Plot the tree obtained with the choices made
printcp(pruned_tree_iRDM)
rpart.plot(pruned_tree_iRDM, cex = 0.6, box.palette = "Blues", extra = 1)


# Prediction on a new non-annotated dataset

x = pruned_tree_iRDM$variable.importance
df_importance_iRDM = data.frame(
  Variable = factor(names(x),
  levels = names(x)),
  Importance = as.numeric(x))

# y_test = predict(pruned_tree, newdata = New_CT_iRDM_df, type = 'class') # output: most probable class
# y_test = predict(pruned_tree, newdata = New_CT_iRDM_df, type = "prob") # with all probabilities


ggplot(df_importance_iRDM, aes(x = Importance, y = Variable)) +
  geom_col(fill = "skyblue") +
  labs(x = "Variable importance", size = 0.1) +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 6))



# Random forest
forest_iRDM = randomForest(iRDM_subtype ~ ., data = Class_iRDM_df)
plot(forest_iRDM)


# ROC Curve




#### (2) Classification: response_variable = PAM50


#### (3) Comparison between two evaluation

# Comparison between them
# Visualization

















## Survival model: Predicting DMFS (Year/event) ~ genes and subtypes

# Transform DMFS_event in binary / Cleanse NA
# Build tree: DMFS_event ~ iRDM_subtypes + PAM50_subtypes + Proliferation_group + Immunity_group + 72 gènes
# Does any genes are more important ? 

# Visualization + Visualization clusters using UMAP


## Observation of 10-15 most important genes from classification tree

# Take the 10-15 most important genes from classification tree
# Re-do regression tree for them
# Analysis
# Visualization





