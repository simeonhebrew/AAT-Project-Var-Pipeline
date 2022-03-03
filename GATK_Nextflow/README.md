## Reproducing Nextflow GATK Variant Calling Pipeline Uisng Conda or Docker ##

---
### Using a Conda environment 

This following is stepwise guideline on how to reproduce this Nextflow variant calling pipeline for paired-end datasets using a Conda environment.

1. Kindly ensure that ```Conda``` is installed in your computer. You can follow these [installation steps](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)
   if you haven't installed it yet.

2. Create a new conda environment that has all the required packages and tools. You can find the environment file in the ```env``` directory. 

    ```conda env create -n gatkenv -f env/gatkenv.yaml```

3. Activate the conda environment

    ```conda activate gatkenv``` 

4. In the ```nextflow.config``` file, add the paths to your:

   - Paired-end reads

   - Reference sequence

   - Known variants file

   - Output directory (directory containing all your output files)

   - ```Tracedir``` directory (directory containing all the nextflow reports)

5. Run Nextflow using the following command

    ```nextflow run main.nf```

### Using a Docker container

1. To reproduce the workflow using the required environment from Docker, kindly ensure that Docker is installed on your machine. You can follow these 
[installation steps](https://docs.docker.com/engine/install/) if you haven't installed it yet.

2. Pull the container into your machine 
   
   ```docker pull hebrewsimeon/condatest:latest```


---

If you are calling variants from multiple samples, consider using the `split.sh` script provided which splits your multisample vcf file to individual
files based on sample identification, extracts SNPs from the splitted files and converts them to nexus matrices ready for phylogenetic analysis.

---
   
