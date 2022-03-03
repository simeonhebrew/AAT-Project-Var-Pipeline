## Checking for decontamination

This is a detailed documentation of how fastQscreen tool was used to check for suspected decontaminate reads from other genomes.

Documentation on the usage of fastQscreen can be found [here](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/_build/html/index.html) as well as video tutorials.

fastQscreen was downloaded in the respective working directory.

As fastQscreen requires a config file to work, the fastq_screen.conf.example was used in the creation of the working config file.

`cp fastq_screen.conf.example fastq_screen.conf`

This config file `fastq_screen.conf` was edited according to the preferred aligner and relevant databases need for alignment.

### Steps in editing the config file.
1. Selecting the preferred aligner -  Uncommenting the preferred aligner and giving path to the executable file.

`BWA /opt/apps/bwa/0.7.17/bin/bwa`

2. Adding paths to the reference genomes databases.


        ## Glossinapallidipes - sequences available from
        ##https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/688/715/GCA_000688715.1_Glossina_pallidipes-1.0.3/GCA_000688715.1_Glossina_pallidipes-1.0.$
        #DATABASE	Glossinapallidipes   /data/rose/AATproject/Glossina/GCA_000688715.1_Glossina_pallidipes-1.0.3_genomic.fna

        ## Trypanasomacongolense -sequences available from
        ##https://tritrypdb.org/common/downloads/Current_Release/TcongolenseIL3000/fasta/data/TriTrypDB-52_TcongolenseIL3000_Genome.fasta
        #DATABASE        Trypanasomacongolense   /data/rose/AATproject/T.congolense/TriTrypDB-52_TcongolenseIL3000_Genome.fasta

        ## Trypanasomavivax -sequences available from
        ##https://tritrypdb.org/common/downloads/Current_Release/TvivaxY486/fasta/data/TriTrypDB-52_TvivaxY486_Genome.fasta
        #DATABASE	Trypanasomavivax   /data/rose/AATproject/T.vivax/TriTrypDB-52_TvivaxY486_Genome.fasta

        ## Bos_Taurus -sequences available from
        ##https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/263/795/GCF_002263795.1_ARS-UCD1.2/GCF_002263795.1_ARS-UCD1.2_genomic.fna.gz
        #DATABASE	Bostaurus   /data/rose/AATproject/Bos_taurus/GCF_002263795.1_ARS-UCD1.2_genomic.fna

3. Creating index files of the respective reference genomes

`bwa index GCA_000688715.1_Glossina_pallidipes-1.0.3_genomic.fna`

`bwa index TriTrypDB-52_TcongolenseIL3000_Genome.fasta`

`bwa index TriTrypDB-52_TvivaxY486_Genome.fasta`

`bwa index GCF_002263795.1_ARS-UCD1.2_genomic.fna`

### Confirming for the presence of decontaminated reads

After the config file was customized to our sample dataset, contaminated reads were checked using fastQscreen.

1. Checking for decontaminated reads in Tryp-B87_S2 sample.

`FastQ-Screen-0.14.1/fastq_screen Tryp-B87_S2_1.fastq.gz Tryp-B87_S2_2.fastq.gz`

2. Tryp-T70_S4 sample.

`FastQ-Screen-0.14/fastq_screen Tryp-T70_S4_1.fastq.gz Tryp-T70_S4_2.fastq.gz`

The above commands generated a `html` report of each read with a bar graph indicating the level of contamination against the suspected genomes. Additionally, it also generated a `txt` file containing a table showing the number of reads that mapped to each genome.

Using fastQscreen, the reads were also tagged against the reference genomes to get an overall representation as the first two commands only checked a subset of the reads.

 `FastQ-Screen-0.14/fastq_screen --tag Tryp-T70_S4_1.fastq.gz`

 This decontamination analysis confirmed a high percentage of decontaminate reads and using `samtools` reads that only mapped to *Trypanasoma congolense* were extracted and proceeded to downstream analysis.
