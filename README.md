# HTP

---

## Overview

**HTP: Hashtag Pipeline**

This repo was made to streamline execution of hashtagged single cell data, all the way from basic QC to clustering, markers and various plot outputs

- Performs Hashtaq QC (Doublet & Empty Removal) and outputs QC plots
- Performs Standard Single Cell RNA QC such as filtering on features, percent.mt etc.
- Runs Seurat Pipeline including UMAP generation, Harmony Batch Correction, Clustering & various plots
- Outputs RDS files at each stage of analysis for more custom analysis

## Project Structure

```sh
└── hashtag-pipeline/
    └── htp
        ├── raw-data
        ├── routes
        ├── run
        └── scripts
```

### Project Index

<details open>
	<summary><b><code>HASHTAG-PIPELINE/</code></b></summary>
	<!-- htp Submodule -->
	<details>
		<summary><b>htp</b></summary>
		<blockquote>
			<div class='directory-path' style='padding: 8px 0; color: #666;'>
				<code><b>⦿ htp</b></code>
			<table style='width: 100%; border-collapse: collapse;'>
			<thead>
				<tr style='background-color: #f8f9fa;'>
					<th style='width: 30%; text-align: left; padding: 8px;'>File Name</th>
					<th style='text-align: left; padding: 8px;'>Summary</th>
				</tr>
			</thead>
				<tr style='border-bottom: 1px solid #eee;'>
					<td style='padding: 8px;'><b><a href='/Users/aleksandrprystupa/Projects/Alex-BINF-Pipelines/hashtag-pipeline/blob/master/htp/run'>run</a></b></td>
					<td style='padding: 8px;'>- Execute a bash script to orchestrate High-Throughput Processing (HTP) tasks in the project<br>- The script triggers sequential jobs for HTP data quality control and processing based on dependencies, enhancing workflow efficiency.</td>
				</tr>
			</table>
			<!-- scripts Submodule -->
			<details>
				<summary><b>scripts</b></summary>
				<blockquote>
					<div class='directory-path' style='padding: 8px 0; color: #666;'>
						<code><b>⦿ htp.scripts</b></code>
					<table style='width: 100%; border-collapse: collapse;'>
					<thead>
						<tr style='background-color: #f8f9fa;'>
							<th style='width: 30%; text-align: left; padding: 8px;'>File Name</th>
							<th style='text-align: left; padding: 8px;'>Summary</th>
						</tr>
					</thead>
						<tr style='border-bottom: 1px solid #eee;'>
							<td style='padding: 8px;'><b><a href='/Users/aleksandrprystupa/Projects/Alex-BINF-Pipelines/hashtag-pipeline/blob/master/htp/scripts/HTO-QC.R'>HTO-QC.R</a></b></td>
							<td style='padding: 8px;'>- Generate HTO-QC statistics, filter singlets, and plot feature scatter graphs for each sample in the project<br>- Save resulting data and plots in designated directories<br>- This script automates quality control processes for high-throughput single-cell sequencing data, enhancing efficiency and reproducibility.</td>
						</tr>
						<tr style='border-bottom: 1px solid #eee;'>
							<td style='padding: 8px;'><b><a href='/Users/aleksandrprystupa/Projects/Alex-BINF-Pipelines/hashtag-pipeline/blob/master/htp/scripts/HTP-Processing.R'>HTP-Processing.R</a></b></td>
							<td style='padding: 8px;'>- Generate downstream analysis plots, identify markers, and run a multi-sample pipeline on RNA sequencing data using Seurat and other libraries<br>- Organize directories, load data, and save processed objects for further analysis.</td>
						</tr>
						<tr style='border-bottom: 1px solid #eee;'>
							<td style='padding: 8px;'><b><a href='/Users/aleksandrprystupa/Projects/Alex-BINF-Pipelines/hashtag-pipeline/blob/master/htp/scripts/RNA-QC.R'>RNA-QC.R</a></b></td>
							<td style='padding: 8px;'>- Generate RNA quality control (QC) plots, calculate summary stats, and merge Seurat objects<br>- Filter merged object based on RNA features and mitochondrial content<br>- Save QC results and merged objects for downstream analysis<br>- Key steps include QC plotting, stats calculation, object merging, and filtering based on feature thresholds.</td>
						</tr>
					</table>
					<!-- helper-functions Submodule -->
					<details>
						<summary><b>helper-functions</b></summary>
						<blockquote>
							<div class='directory-path' style='padding: 8px 0; color: #666;'>
								<code><b>⦿ htp.scripts.helper-functions</b></code>
							<table style='width: 100%; border-collapse: collapse;'>
							<thead>
								<tr style='background-color: #f8f9fa;'>
									<th style='width: 30%; text-align: left; padding: 8px;'>File Name</th>
									<th style='text-align: left; padding: 8px;'>Summary</th>
								</tr>
							</thead>
								<tr style='border-bottom: 1px solid #eee;'>
									<td style='padding: 8px;'><b><a href='/Users/aleksandrprystupa/Projects/Alex-BINF-Pipelines/hashtag-pipeline/blob/master/htp/scripts/helper-functions/utils.R'>utils.R</a></b></td>
									<td style='padding: 8px;'>Create new save directories and load multiple sheets from Excel files.</td>
								</tr>
								<tr style='border-bottom: 1px solid #eee;'>
									<td style='padding: 8px;'><b><a href='/Users/aleksandrprystupa/Projects/Alex-BINF-Pipelines/hashtag-pipeline/blob/master/htp/scripts/helper-functions/hto-utils.R'>hto-utils.R</a></b></td>
									<td style='padding: 8px;'>- Enable full HTO workflow by loading data, performing downstream analysis with HTO demux, adding metadata, and obtaining HTO demux stats<br>- Implement filtering stats per hashtag and plot feature scatter plots for comprehensive analysis<br>- Use these functions to streamline HTO data processing and analysis effortlessly.</td>
								</tr>
								<tr style='border-bottom: 1px solid #eee;'>
									<td style='padding: 8px;'><b><a href='/Users/aleksandrprystupa/Projects/Alex-BINF-Pipelines/hashtag-pipeline/blob/master/htp/scripts/helper-functions/hto-testing.R'>hto-testing.R</a></b></td>
									<td style='padding: 8px;'>- Load and normalize data from HTO samples, then display metadata and sample groups<br>- This script interacts with Seurat and other libraries to process and analyze HTO data, facilitating data loading and manipulation for downstream analysis in the project architecture.</td>
								</tr>
								<tr style='border-bottom: 1px solid #eee;'>
									<td style='padding: 8px;'><b><a href='/Users/aleksandrprystupa/Projects/Alex-BINF-Pipelines/hashtag-pipeline/blob/master/htp/scripts/helper-functions/hto-main-function.R'>hto-main-function.R</a></b></td>
									<td style='padding: 8px;'>- Define the main function <code>run.hto.workflow</code> to execute the HTO workflow, including loading data, downstream analysis, and HTO demultiplexing<br>- This function integrates helper functions from <code>hto-utils.R</code> and <code>utils.R</code> to process data efficiently within the projects architecture.</td>
								</tr>
								<tr style='border-bottom: 1px solid #eee;'>
									<td style='padding: 8px;'><b><a href='/Users/aleksandrprystupa/Projects/Alex-BINF-Pipelines/hashtag-pipeline/blob/master/htp/scripts/helper-functions/hto-testing.Rmd'>hto-testing.Rmd</a></b></td>
									<td style='padding: 8px;'>- Generate a summary detailing the function of the hto-testing file within the project's architecture<br>- The file processes Seurat objects, calculates various metrics, and generates quality control summaries for single-cell RNA sequencing data<br>- It plays a crucial role in assessing data quality and ensuring downstream analysis accuracy.</td>
								</tr>
								<tr style='border-bottom: 1px solid #eee;'>
									<td style='padding: 8px;'><b><a href='/Users/aleksandrprystupa/Projects/Alex-BINF-Pipelines/hashtag-pipeline/blob/master/htp/scripts/helper-functions/processing-utils.R'>processing-utils.R</a></b></td>
									<td style='padding: 8px;'>- Define a multi-sample pipeline to process data, including normalization, feature selection, scaling, PCA, harmony correction, UMAP, and clustering<br>- Implement functions to find markers using statistical tests and visualize UMAP plots with options for labeling and saving<br>- These functions enhance data analysis capabilities within the projects architecture.</td>
								</tr>
								<tr style='border-bottom: 1px solid #eee;'>
									<td style='padding: 8px;'><b><a href='/Users/aleksandrprystupa/Projects/Alex-BINF-Pipelines/hashtag-pipeline/blob/master/htp/scripts/helper-functions/rna-utils.R'>rna-utils.R</a></b></td>
									<td style='padding: 8px;'>- Generate RNA quality control statistics and plots for Seurat objects, aiding in data analysis and visualization<br>- The functions in this file help assess RNA quality by plotting violin and scatter plots, and summarizing key metrics like feature counts and mitochondrial gene percentages.</td>
								</tr>
							</table>
						</blockquote>
					</details>
				</blockquote>
			</details>
			<!-- routes Submodule -->
			<details>
				<summary><b>routes</b></summary>
				<blockquote>
					<div class='directory-path' style='padding: 8px 0; color: #666;'>
						<code><b>⦿ htp.routes</b></code>
					<table style='width: 100%; border-collapse: collapse;'>
					<thead>
						<tr style='background-color: #f8f9fa;'>
							<th style='width: 30%; text-align: left; padding: 8px;'>File Name</th>
							<th style='text-align: left; padding: 8px;'>Summary</th>
						</tr>
					</thead>
						<tr style='border-bottom: 1px solid #eee;'>
							<td style='padding: 8px;'><b><a href='/Users/aleksandrprystupa/Projects/Alex-BINF-Pipelines/hashtag-pipeline/blob/master/htp/routes/run-htp-rna-qc.sh'>run-htp-rna-qc.sh</a></b></td>
							<td style='padding: 8px;'>- Execute a bash script to run RNA quality control analysis<br>- Set job parameters and load necessary modules<br>- Create log directories and define the R script command<br>- Finally, execute the R script for Step-1 of the process.</td>
						</tr>
						<tr style='border-bottom: 1px solid #eee;'>
							<td style='padding: 8px;'><b><a href='/Users/aleksandrprystupa/Projects/Alex-BINF-Pipelines/hashtag-pipeline/blob/master/htp/routes/run-htp-hto-qc.sh'>run-htp-hto-qc.sh</a></b></td>
							<td style='padding: 8px;'>- Execute a script to run quality control analysis for HTO data<br>- Set up job parameters and load necessary modules<br>- The script runs an R command to perform QC tasks.</td>
						</tr>
						<tr style='border-bottom: 1px solid #eee;'>
							<td style='padding: 8px;'><b><a href='/Users/aleksandrprystupa/Projects/Alex-BINF-Pipelines/hashtag-pipeline/blob/master/htp/routes/run-htp-processing.sh'>run-htp-processing.sh</a></b></td>
							<td style='padding: 8px;'>- Facilitates running HTP processing job with specified configurations<br>- Sets up job parameters and executes R script for processing<br>- Handles job scheduling and resource allocation efficiently<br>- Simplifies running HTP processing tasks within the project architecture.</td>
						</tr>
					</table>
				</blockquote>
			</details>
		</blockquote>
	</details>
</details>

---

## Usage

Run the project with:

1. htp/run (runs all 3 scripts, which depend on each other to continue)
