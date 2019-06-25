# ST-Composer

Generates paired synthetic Spatial Transcriptomics and Single Cell Data 

The ST-data uses an image as the "histology" defining regions with
different cell compositions of K-specified celltypes. 

The SC-data will be in the form of a balanced set of single cells
coming from either of the K-celltypes.

The software is CLI based but the main functions can also be imported and used in isolation for other means, though this has not been the main usage kept in mind upon design.

By providing an image equal or larger to the dimensions of 625x625 pixels and the number of regions (corresponding to colors) to be used, these will be assigned unique celltype compositions of the
same K celltypes (specified by the user). By default, all pixels with the rgb-color (128,128,128) or equivalent HEX \#808080 will be considered as 'background' and not included in the generated data.
The user can specify an alternative color to be used as background by using the "--background" command and the desired HEX-color value (including the hashtag).

The process for genereating the synthetic ST-data is in short

1. Generate spot-coordinates and regions based on provided image
2. Specify the number of celltypes present within the dataset
3. Generate poisson parameters for each gene (shared among all cell types) 
4. Pick a number of outlier genes and rescale these (shared among all cell types)
5. Assign DE status to a number of genes and rescale these (specific to cell type)
6. Assign probability of observing a celltype within each region (shared between spots of the same region)
7. Draw cells that contribute to every spot based on the probabilities from (6.)
8. Rescale observations as to make them partial
9. Generate a SC-dataset consisting of the same cell types

This will generate regions where the same cell types are present but at different proportions (which are known to the user)


## Examples
To be included



