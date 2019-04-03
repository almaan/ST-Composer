# ST-Composer

Generates paired synthetic Spatial Transcriptomics and Single Cell Data 

The ST-data uses an image as the "histology" defining regions with
different cell compositions of K-specified celltypes. 

The SC-data will be in the form of a balanced set of single cells
coming from either of the K-celltypes.

The software is CLI based but the main functions can also be imported and used in isolation for other means, though this has not been the main usage kept in mind upon design.

By providing an image equal or larger to the dimensions of 625x625 pixels and the number of regions (corresponding to colors) to be used, these will be assigned unique celltype compositions of the
same K celltypes (specified by the user). By default, all pixels with the rgb-color (128,128,128) or equivalent HEX \#808080 will be considered as 'background' and not included in the generated data.
The user can specify an alternative color to be used as background by using the "--backgroun" command and the desired HEX-color value (including the hashtag).

The process for genereating the synthetic ST-data is in short

1. Generate spot-coordinates and regions based on provided image
2. Specify the number of celltypes present within the dataset
3. make  n\_region draws from Dir(alpha) where alpha is a vector of length K (celltypes) with the same value, default 0.5, of all entries; these are the priors for each region. Denote the
   K-dimensional vector for region r as p\_r
4. For each spot s within region r draw $n_{sr}$ from Unif(1,upper\_limit) representing the number of cells. 
5. For each spot sample $z$ from Mult($n_{sr},p_r)$ 

