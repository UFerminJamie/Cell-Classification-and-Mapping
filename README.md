# Cell-Classification-and-Mapping

This is a **group assignment**.

## Code Implementation & Technical Report

The final deliverables include a 2-page IEEE-format report, code implementation and a detailed GitHub readme file.

## Data

* insert link to the data

## About The Project

EEL5934 (Spring 2023) Final Project 

Title: Single Cell Classification and Mapping

Project Description: Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.

## Usage 

To replicate the results from the report, follow these steps:

1.  Clone the repository and navigate to the project directory.
	
	`git clone https://github.com/UFerminJamie/Cell-Classification-and-Mapping.git`
	 
2.  Ensure that all the necessary dependencies are installed.

3. Within `data_extraction_visualization.ipynb`, run the `mask_to_xml` function to generate an XML file containing all the annotated cell regions. This XML file can be used to visualize the quality of the nuclei segmentation mask by reading the H&E image and its corresponding XML file in ImageScope Aperio. The `mask_to_xml` function uses the following arguments : `xml_path`(directory where to save the XML file), and the `mask`(nuclei mask image).

	Secondly, in the same notebook, run the `cell_extraction` function, which extracts patches where there are nuclei in the 17-channel CODEX image. This function accepts the following arguments: `tif_file`( 17-channel CODEX image), `xml_file`(corresponding XML file generated earlier), `mask_file`(nuclei mask image), `patch_folder` (path to where image patches get saved). The output of this function would be N number of image patches, each containing a single nucleus.

4. To do feature extraction, run the `extract_intensity.m` file after completing the previous steps. Before running, edit line 19 to specify the directory where the image patches are located.  This script calculates the mean and standard deviation for each channel in the CODEX image for all nuclei, producing an N x M feature matrix where N is the number of features and M is the number of nuclei. This outputs two CSV files 1) 17 X M feature matrix (mean intensities only), and 2) 33 x M feature matrix (mean & SD). Please note that the first column of the extracted features CSV files corresponds to the ID of the nucleus and is not considered a feature.

5. Next, execute the `UMAP.R` script to perform unsupervised clustering. This script uses the feature matrix generated in the previous step and produces a CSV file that contains predicted labels for each cell. The Seurat package, which is commonly used for scRNA-seq data analysis, is utilized in this file. Prior to clustering, the D most variable features were selected, followed by data normalization (log normalization), scaling, and PCA. To identify the different clusters, a graph-based approach was utilized, which constructs a k-nearest neighbor (kNN) graph based on the intensity features. In this approach, each cell is connected to its k nearest neighbors.

6. To complete the pipeline, open the `Mapping.ipynb` notebook and load the CSV file containing the predicted cluster labels. Then, use the overlay function with the following arguments: `tif_path` (the path to the CODEX image), `xml_path` (the path to the corresponding XML file), `num_cluster` (the number of clusters), and `dataframe` (the cluster ID column in the CSV file). The `overlay` function generates an RGB overlay image on top of a CODEX image where each cluster is assigned a different color for visualization purposes.

## Authors
Jamie L. Fermin - j.fermin@ufl.edu

Hannah Kempfert - hannahkempfert@ufl.edu

Antonio Hendricks - a.hendricks1@ufl.edu
