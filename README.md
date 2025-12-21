# User manual for 16S RNA analysis pipeline
Last update: DD/MM/2025
Created by: Alžbeta Žiarovská, as a part of bachleor's thesis
Supervised by: Mgr. Martin Blažek
Bachleor's thesis ID: FIIT-16768-127323
___
## About the pipeline
This pipeline implements a workflow for 16S RNA data analysis. It is written using Nextflow sotware and implements Anaconda distribution for virtual environments and tools used for the processes.
___
## How to use the pipeline
#### 1. Before you start
To use this pipeline you have to clone the repository into your local storage using command: 
```git clone https://github.com/```

The the cloned folder download a pretrained 16S taxonomic classifier of your choice. For the testing of the pipeline was used: [SILVA 138 99% OTUs full-length sequences](https://library.qiime2.org/data-resources)

Anaconda Distribution is required for running the pipeline. To install Anaconda follow official [Anaconda Installation tutorial](https://www.anaconda.com/docs/getting-started/anaconda/install).
#### 2. Input parameters
In order to change the default behaviour of the pipeline you can modify the input parameters. In this section you can find the explanation for all of the input parameters that can be passed to the workflow along with their explanation and possible values.
##### pwd
This parameter hold a current working directory. **Do not change the parameter.** It is set automatically and changing it could fatally effect the pipeline.

##### download
Based on the value of this parameter, the pipeline either loads the data from a local directory or downloads it from the SRA archive based on SRR ID.

Default value: *true*

Possible values: *true*/*false*

Parameter configuration: ```nextflow run main.nf --download true```

##### data_file
If params.download is set to value true, the data is downloaded based on a SRR IDs written in a text file.

Example data file: [SRR_Acc_list.txt](SRR_Acc_list_70.txt)

Default value: *SRR_Acc_list.txt*

Possible values: *file_name.txt*

Parameter configuration: ```nextflow run main.nf --data_file SRR_Acc_list.txt```

##### data_dir
This parameter is used when the data is loaded from a local storage. 

Needed files: *<SRR_ID>_1.fastq.gz* and *<SRR_ID>_2.fastq.gz*
**Please use this format of file names for correct pipeline function**

Default value: *data*

Possible values: *data_folder_name*

Parameter configuration: ```nextflow run main.nf --data_dir data_folder```

##### out_dir
Parameter to determine the folder name where the results are stored.

Default value: *results*

Possible values: *results_folder_name*

Parameter configuration: ```nextflow run main.nf --out_dir results_folder```

##### name
This parameter names the output plots. 

Default value: *run1*

Possible values: any string of choice

Parameter configuration: ```nextflow run main.nf --name name_of_run```

##### adapter
Adapter sequence that gets cut from the sequences. Use only forward adapter by entering only one adapter sequence or both forward and reverse adapter by diving their sequences with ..., the pipeline automatically formats them afterwards.

Default value: *TACGGRAGGCAGCAG...AGGGTATCTAATCCT*

Example value: *TACGGRAGGCAGCAG...AGGGTATCTAATCCT*

Possible values: *TACGGRAGGCAGCAG* or *TACGGRAGGCAGCAG...AGGGTATCTAATCCT*

Parameter configuration: ```nextflow run main.nf --adpter TACGGRAGGCAGCAG...AGGGTATCTAATCCT```

##### qiime_classifier
The path to a Qiime2 pretraind classifier that will be used for taxonomic classification. This file is needed to be located in the working directory.

Example classifier: [SILVA 138 99% OTUs full-length sequences](https://library.qiime2.org/data-resources)

Defalut value: *${params.pwd}/silva-138-99-nb-classifier.qza*

Possible values: *absolute/path/to/classifier.qza*

Parameter configuration: ```nextflow run main.nf --qiime_classifier absolute/path/to/classifier.qza```

##### top_taxa
This parameter determines a number of top represented taxons that get visualized in some of the plots for higher readibility.

Default value: *20*

Possible values: *number from 1 to number of taxons*

Parameter configuration: ```nextflow run main.nf --top_taxa 20```

```nextflow run main.nf --parameter_name parameter_value```

#### 3. Running the pipeline
1. Open the folder, where ```main.nf``` is located in terminal
2. Make sure you have Anaconda installed
3. Make sure the classifier is located in the working directory
4. Make sure you have either ```data_file``` or ```data_dir``` in the working directory in needed format
5. Run the pipeline by ```nextflow run main.nf -with-conda```

#### 4. Output files
All of the output files can be found in the folder determined by parameter *out_dir*.

In this folder you can find following subfolders and files:
- data - downloaded or locally found .fastq.gz files for pair-end reads
- fastqc_trimmed - fastqc reports of trimmed data files
- fastqc_raw - fastqc reports of raw data
- statistics - plots in .pdf format and two .csv tables with taxonomy classification and abundance data
- taxonomy - feature, taxonomy and merged table in .csv format for each sequence
- trimmed - trimmed .fastq.gz files
- pipeline_run_info - automatically created file by [Nextflow CO₂footprint plugin](https://github.com/nextflow-io/nf-co2footprint)
- multiqc_raw_report.html
- multiqc_trimmed_report.html

___
## Example input and output

Input:

[Example input file](SRR_Acc_list.txt)
[Example input folder](data)

Output:
[Example output folder](results)

___
## References
`J. Carl, N. Volkmann, J. Mir-Pedrol, P. Ewels, S. Nahnsen, S. Krakau nextflow-io/nf-co2footprint v1.0.0. (Jun., 2025). nextflow-io. Available: https://github.com/nextflow-io/nf-co2footprint`

`Bokulich, N. A., Kaehler, B. D., Rideout, J. R., Dillon, M., Bolyen, E., Knight, R., Huttley, G. A., & Gregory Caporaso, J. (2018). Optimizing taxonomic classification of marker-gene amplicon sequences with QIIME 2’s q2-feature-classifier plugin. Microbiome, 6(1).`

`Robeson, M. S., II, O’Rourke, D. R., Kaehler, B. D., Ziemski, M., Dillon, M. R., Foster, J. T., & Bokulich, N. A. (2020). RESCRIPt: Reproducible sequence taxonomy reference database management for the masses. In bioRxiv. bioRxiv.`

`Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P, Peplies J, Glöckner FO (2013) The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucl. Acids Res. 41 (D1): D590-D596.`

`Yilmaz P, Parfrey LW, Yarza P, Gerken J, Pruesse E, Quast C, Schweer T, Peplies J, Ludwig W, Glöckner FO (2014) The SILVA and "All-species Living Tree Project (LTP)" taxonomic frameworks. Nucl. Acids Res. 42:D643-D648`

`Glöckner FO, Yilmaz P, Quast C, Gerken J, Beccati A, Ciuprina A, Bruns G, Yarza P, Peplies J, Westram R, Ludwig W (2017) 25 years of serving the community with ribosomal RNA gene reference databases and tools. J. Biotechnol.`

`Daniel Straub, Nia Blackwell, Adrian Langarica-Fuentes, Alexander Peltzer, Sven Nahnsen, Sara Kleindienst (2020) Interpretations of Environmental Microbial Community Studies Are Biased by the Selected 16S rRNA (Gene) Amplicon Sequencing Pipeline, Frontiers in Microbiology 2020, 11:2652 doi: 10.3389/fmicb.2020.550420.`

`Huang, X., Zeng, J., Li, S., Chen, J., Wang, H., Li, C., & Zhang, S. (2024). 16S rRNA, metagenomics and 2bRAD-M sequencing to decode human thanatomicrobiome. Scientific Data, 11(1), 736.`