# TOGA README: step-by-step instructions
**Authors: Tanya Lama, Bill Thomas**
with contributions from Diana Moreno Santillan, David Ray, Ariadna Morales (thanks!)

## Table of Contents

[TOC]

## Step 1. Get genomes 
### Create a folder for the raw fasta [here](https://drive.google.com/drive/folders/1xD49VJuw6NxfhGmiyQkbVXXYF5kfgvbk?usp=sharing)
You can ask collaborators to deposit their assembly fastas in this folder. 

### Deposit the raw fasta folder (gzipped) in the /project space 
Create a folder with the genSpe name in /project: 
/project/tlama_umass_edu/projects/project_toga/data/genomes/hg38/trackData/
Deposit zipped fasta here using rclone 
`rclone copy sbugoogledrive:unity_backup/project_toga/inprogress_togas/musFur_tpursell/musFur_hifiasm.bp.p_ctg.fna.gz /project/tlama_umass_edu/projects/project_toga/data/genomes/hg38/trackData/musFur/`

Create a folder with the genSpe name in /work: 
/work/tlama_umass_edu/project_shrew_genome/data/genomes/hg38/trackData/

### Calculate relevant stats on your genome assembly
contig N50 and L50; BUSCO

### Loading dependencies 

Load the anaconda module
```
module load anaconda/2022.10

```
You need to be in a conda environment called `toga` with preloaded dependencies to execute tools like faToTwoBit and faSize in the script below.
You only need to do this ONCE. 

### Create the toga environment from toga.yml
```
conda env create -f /project/tlama_umass_edu/projects/project_toga/toga.yml
```
Once the environment has been created, you can activate it using 
```
conda activate toga
```
You should see (toga) to the left of your username on Unity. 
![](https://hackmd.io/_uploads/Hy5ZpcaWT.png)
You are now in the toga environment and you may proceed!

### Get chromsizes 
You can follow the script `mask_prep.sh` or edit the code below and submit with sbatch 
```
#!/bin/bash
#SBATCH --output=prep.out
#SBATCH --error=prep.err
#SBATCH --job-name=prep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=cpu

#this is our working directory
GL=/project/tlama_umass_edu/projects/project_toga/data/genomes/hg38/trackData #do not edit this

#genSpe name is musFur
idx=musFur #change this
G1=musFur_hifiasm.bp.p_ctg.fna.gz #change this
G2=musFur_hifiasm.bp.p_ctg.fna #change this
#mkdir $GL/$idx
#wget  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/829/155/GCA_009829155.1>
gunzip $GL/$idx/$G1
/project/tlama_umass_edu/bin/faToTwoBit $GL/$idx/$G2  $GL/$idx/$idx.2bit
/project/tlama_umass_edu/bin/faSize -detailed -tab $GL/$idx/$G2 > $GL/$idx/$idx.
chromsize
gzip $GL/$idx/$G2 #zip it back up
seff ${SLURM_JOB_ID} >> ./logs/${SLURM_JOB_ID}.log
```
This should have created a .2bit and .chromsize file in your folder.

### Sanity Check!!!
Do you see a bunch of chromosome names and their sizes?
`head musFur_hifiasm.chromsize`
Is your genome.fna file at least 2.0G in size (for mammals)
`ls -sh musFur_hifiasm.fna`

### Document relevant metadata [here](https://docs.google.com/spreadsheets/d/14R93B03xRgjdVJGi53Jds7WhwaxREW9qbgt6ts9eAmo/edit?usp=sharing)

If **masked**, skip Step 2.
If **unmasked**, complete Step 2.

## Step 2. RepeatModel & RepeatMasker
Please review the [RepeatModeler](https://www.pnas.org/doi/10.1073/pnas.1921046117) and [RepeatMasker](https://github.com/rmhubley/RepeatMasker) papers to understand what is happening in that step. 

Note that the RepeatModeler step will run for **3-5 days**. 
RepeatMasker may run for up to **24h**.
This is one of the most time intensive steps in the TOGA pipeline.

### Load the mask_genomes environment 
```
conda activate /project/tlama_umass_edu/bin/anaconda3/envs/mask_genomes
```
You can create the mask_genomes environment from mask_genomes.yml if necessary (this did not work)
```
conda env create -f /project/tlama_umass_edu/projects/project_toga/mask_genomes.yml
```
Once the environment has been created, you can activate it using 
```
conda activate mask_genomes
```
If that doesn't work you can use my mask_genomes environment in our /bin.
`conda activate /project/tlama_umass_edu/bin/anaconda3/envs/mask_genomes`

## Execute RepeatModeler 
Working directory for this step is: /work/tlama_umass_edu/project_toga/genomes/hg38/trackData/musFur
Edit and sbatch `slurm_mask_genome_NEWUNITY.sh` or edit the code below and submit with sbatch
```
#!/bin/bash
#SBATCH --output=masked.out
#SBATCH --error=masked.err
#SBATCH --job-name=mask-stc319-artjam
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=cpu-long
#SBATCH --time=120:00:00 #5 days

#Are you in the mask_genomes environment? Have you loaded the slurm module? Are you executing this job from the /work directory? If NO to any of the above, go back.

result=${PWD##*/}
echo $result
G=stc319-artjam.fna.gz #your genome.fna.gz
G1=stc319-artjam.fna #your genome.fna 
gunzip /project/tlama_umass_edu/projects/project_toga/data/genomes/hg38/trackData/$result/$G

python /project/tlama_umass_edu/projects/project_shrew_genome/scripts/mask_genomes_toga_NEWUNITY.py -qy /project/tlama_umass_edu/projects/project_toga/data/genomes/hg38/trackData/$result/$G1 -c 2 -o $result'_masked'

gzip /project/tlama_umass_edu/projects/project_toga/data/genomes/hg38/trackData/
$result/$G1

seff ${SLURM_JOB_ID} >> ./logs/${SLURM_JOB_ID}.log
```
## Edit the README to keep a log of what has been executed/when

### Sanity Check!!!
RepeatModeler should have created a series of folders: 1_masking2/RepeatModeler/RM_output. 
Check out the .err and .out files in each folder -- see anything worrisome?
`cat *.err`
RM_output includes your final results. Let's go and look at them. 
Is consensi.fa betweek 700-900K in size? 
`ls -sh consensi.fa`
grep a couple of chromosome names from families-classified.stk -- do you see results for each chromosome? 
`grep 'ptg000001l' families-classified.stk` #yes
`grep 'ptg000201l' families-classified.stk` #nope
cd to the /RepeatModeler folder
`cat stc319-artjam_repmod.translation`
Do you see the full list of chromosomes? If yes, then RepeatModeler has scanned the whole genome for repeats. If NO (e.g., I only see chromosomes 1-6) then your .fna file is likely corrupted or truncated, and you need to go to /project folder and take a look at it. 

## Execute RepeatMasker
cd to /1_masking2 and make a new folder called RepeatMasker 
You will need two scripts in this folder: Repeatmasker_slurm.sh and rm2bed.sh
You need to edit both scripts as below, then execute `Repeatmasker_slurm.sh` with sbatch 

Repeatmasker_slurm.sh
```#!/bin/bash
#SBATCH --job-name=artJam_repmask #change your job name
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=cpu-long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

echo makes variables used in this script
gunzip /project/tlama_umass_edu/projects/project_toga/data/genomes/hg38/trackData/artJam/stc319-artjam.fna.gz #unzip your genome

GENOME=/project/tlama_umass_edu/projects/project_toga/data/genomes/hg38/trackData/artJam/stc319-artjam.fna #path to your genome.fna

RUNTYPE=stc319-artjam_RM #your genome prefix_RM

DIR=/work/tlama_umass_edu/project_toga/genomes/hg38/trackData/artJam/1_masking2/
RepeatMasker/$RUNTYPE
mkdir -p $DIR
cd $DIR
ln -s $GENOME
echo "fa to 2bit"
#/project/tlama_umass_edu/bin/faToTwoBit $GENOME stc319-artjam.2bit #we may not 
need to run this because we converted fa to 2bit in the mask_prep step
cp /project/tlama_umass_edu/projects/project_toga/data/genomes/hg38/trackData/artJam/stc319-artjam.2bit ./
echo "Use slurm_clusterrunC8.py to generate all batches needed to run RepeatMasker and the doLift.sh to compile the results."
python3 /work/tlama_umass_edu/project_shrew_genome/scripts/slurm_clusterrunC8.py -i $GENOME -b 50 -lib /work/tlama_umass_edu/project_toga/genomes/hg38/trackData
/artJam/1_masking2/RepeatModeler/RM_output/consensi.fa.classified -dir . -xsmall #insert the path to your RepeatModeler outputs

echo "Submits the RepeatMasker jobs created by slurm_clusterrunC8.py"

bash qsub.sh

#Creates a list of jobIDs to keep track of the RepeatMasker batches being run.
jobIDs=""; for i in `squeue  | grep stc319-artjam | awk '{print $1}'`; do jobIDs=$jobIDs,$i; done; jobIDs="${jobIDs:1}" 

#Submits the doLift script to the queue but holds it until all jobs in the jobIDs list (the RepeatMasker batches) have cleared.
sleep 10m
X=`squeue -u tlama_umass_edu | wc -l` #change to your username  
while [ $X -gt 2 ]; do sleep 1m; X=`squeue -u tlama_umass_edu | wc -l`; done ##change to your username  
sbatch doLift.sh

#Creates a list of job IDs to keep track of the doLift script.
jobIDs=""; for i in `squeue  | grep stc319-artjam | awk '{print $1}'`; do jobIDs=$jobIDs,$i; done; jobIDs="${jobIDs:1}"

#Submits the rm2bed script to the queue but holds it until all jobs in the jobIDs list (the doLift job) have cleared.
sleep 10m
X=`squeue -u tlama_umass_edu | wc -l`  
while [ $X -gt 2 ]; do sleep 1m; X=`squeue -u tlama_umass_edu | wc -l`; done
sbatch /work/tlama_umass_edu/project_toga/genomes/hg38/trackData/artJam/1_masking2/RepeatMasker/rm2bed.sh

seff ${SLURM_JOB_ID} >> ./logs/${SLURM_JOB_ID}.log
```
rm2bed.sh
```
#!/bin/bash
#SBATCH --job-name=stc319-artjam_RM2 #change to genome name_RM2
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=cpu-long
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks=1

echo makes variables used in this script
GENOME=stc319-artjam # your genome prefix
RUNTYPE=stc319-artjam_RM # your genome prefix _RM
DIR=/work/tlama_umass_edu/project_toga/genomes/hg38/trackData/artJam/1_masking2/RepeatMasker/$RUNTYPE
cd $DIR

#Runs a python script RM2bed to generate one complete .bed file and several subfiles subdivided by TE class. Merges overlapping hits based using lower_divergence criterion.
[ ! -f ${GENOME}_rm.bed ] && python /project/tlama_umass_edu/projects/project_shrew_genome/scripts/RM2bed.py -d . -sp class -p ${GENOME} -o lower_divergence ${G
ENOME}.fa.align.gz

module load gcc/9.3.0 bedtools/2
RMPATH=/project/tlama_umass_edu/bin/anaconda3/envs/mask_genomes/bin/

gunzip -c ${GENOME}.fa.out.gz > ${GENOME}.fa.out
perl $RMPATH/rmOutToGFF3.pl ${GENOME}.fa.out > ${GENOME}.gff
bedtools maskfasta -soft -fi /project/tlama_umass_edu/projects/project_toga/data/genomes/hg38/trackData/artJam/stc319-artjam.fna -bed ${GENOME}.gff -fo ${GENOME}.masked.fa

echo "Masking genome done!

seff ${SLURM_JOB_ID} >> ./logs/${SLURM_JOB_ID}.log
```
### Sanity Check!!!
RepeatMasker should have created a folder called your-genome_RM

Check out the .err and .out files -- see anything worrisome?
`cat *.err`
RM_output includes your final result. Let's go and look at it. 
stc319-artjam.masked.fa is your final masked genome. It should be the same size as your original .fna file in /project 
`ls -sh stc319-artjam.masked.fa` #yes
Take a peek at stc319-artjam.masked.fa. Do you see a mix of capital and lowercase (masked) bases? (GAAAGCAACAgaatagttagggaag)
`head -2 stc319-artjam.masked.fa`
There should be a folder within RM_output called RMPart
Do the folder names run from 000 to 047 (or whatever number of scaffolds your genomes has)?
grep a couple of chromosome names from the folders in RMPart
`grep ptg000001l ./*/*.fa` #any results from chrom1? yes

If all of the above looked good, your genome has been successfully masked, copy stc319-artjam.masked.fa to your /project folder
`cp stc319-artjam.masked.fa /project/tlama_umass_edu/project/project_toga/data/genomes/hg38/trackdata/artJam`
Run chain_prep.sh again to generate 2.bit and chromsizes for your masked genome 
```
#!/bin/bash
#SBATCH --output=prep.out
#SBATCH --error=prep.err
#SBATCH --job-name=prep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=cpu

#this is our working directory
GL=/project/tlama_umass_edu/projects/project_toga/data/genomes/hg38/trackData #do not edit this

#genSpe name is artJam
idx=artJam #change this
G1=stc319-artjam.masked.fna.gz #change this
G2=stc319-artjam.masked.fna #change this

/project/tlama_umass_edu/bin/faToTwoBit $GL/$idx/$G2  $GL/$idx/$idx.masked.2bit
/project/tlama_umass_edu/bin/faSize -detailed -tab $GL/$idx/$G2 > $GL/$idx/$idx.masked.chromsize
#gzip $GL/$idx/$G2 #zip it back up
#seff ${SLURM_JOB_ID} >> ./logs/${SLURM_JOB_ID}.log
```
Proceed to Step 3

## Step 3. Chains 

## Add login1 and login2 to your ~/.ssh config file? 
You will need to make a public/private key and add it to your account on Unity so that TOGA can execute arrays from the login nodes. Navigate to `~/.ssh`
`ssh-keygen` ( select no password)
`cat id_rsa.pub` (copy and paste to your account settings on https://unity.rc.umass.edu/panel/account.php)
try `ssh login1`
Repeat this process for login2.

## Edit your DEF file
Navigate to your /work folder if you're not already there
Edit your DEF file in /work
```
# Hg38_vs_SPECIES #here human vs mm10
PATH=/usr/local/bin:/project/tlama_umass_edu/bin:/project/tlama_umass_edu/projects/project_shrew_genome/scripts
BLASTZ=/project/tlama_umass_edu/bin/lastz-1.04.00
BLASTZ_O=400
BLASTZ_E=30
BLASTZ_M=254
BLASTZ_Q=/project/tlama_umass_edu/bin/scripts/diana_scripts/scores/scores_human_mammal.q

# TARGET: Human_HG38 #borrowing these from the toga annotation data folder
SEQ1_DIR=/project/tlama_umass_edu/projects/project_shrew_genome/data/genomes/hg38/hg38.2bit
SEQ1_LEN=/project/tlama_umass_edu/projects/project_shrew_genome/data/genomes/hg38/hg38.chrom.sizes
SEQ1_CHUNK=175000000
SEQ1_LAP=0
SEQ1_LIMIT=18

# QUERY: SPECIES #here we have eptFus3
SEQ2_DIR=/project/tlama_umass_edu/projects/project_shrew_genome/data/genomes/hg38/trackData/eptFus3/eptFus3_masked.2bit
SEQ2_LEN=/project/tlama_umass_edu/projects/project_shrew_genome/data/genomes/hg38/trackData/eptFus3/eptFus3.masked.chromsize
SEQ2_CHUNK=50000000
SEQ2_LIMIT=2000
SEQ2_LAP=10000

BASE=/work/tlama_umass_edu/project_shrew_genome/data/genomes/hg38/trackData/eptFus3
TMPDIR=/work/tlama_umass_edu/project_shrew_genome/data/genomes/hg38/trackData/eptFus3/tmp
```
Make a directory called tmp 
```
mkdir tmp
```
You can copy `blastz_slurm.sh` or edit the code below and submit with sbatch
```
#!/bin/bash
#SBATCH --output=batch_%A_%a.out
#SBATCH --error=batch_%A.%a.err
#SBATCH --job-name=doBlast
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=cpu-long

/project/tlama_umass_edu/projects/project_shrew_genome/scripts/tanya_modified_doBlastzChain_NU.pl `pwd`/DEF -verbose=2 -noDbNameCheck -workhorse=login1 -bigClusterHub=login1 -skipDownload -dbHost=login1 -smallClusterHub=login1 -trackHub -fileServer=login1 -syntenicNet > do.log
```
## Step 4. TOGA 

## Step 5. 

## Step 6. BUSCO annotation completeness

## Step 7. Calculate n intact genes
The **number of intact genes** is a better metric of assembly (and annotation) quality than BUSCO provides. 
<16k = bad
16-17k = ok 
>17k = good
Set threshold(s) for inclusion as you see fit, ours is set at 16k minimum n intact genes.

Using TOGA for annotations, the file "loss_summ_data.tsv" has the classification of each gene/transcript. Everything on the third column that has an an "I" for intact are the genes you want to use. The citation for using this threshold would be Kirilenko et al. 2022

There are a couple of ways of doing this: 
```
awk '{print $1}' loss_summ_data.tsv | grep -n GENE > just_genes.txt
awk '{print}' just_genes.txt | grep -w "I" > intact_genes.txt
wc -l intact_genes.txt
```
OR (fewer lines of code):
```
cat loss_summ_data.tsv | grep "GENE" | grep -v P | grep "I" | wc -l
```

## Step 8. Update relevant metadata [here](https://docs.google.com/spreadsheets/d/14R93B03xRgjdVJGi53Jds7WhwaxREW9qbgt6ts9eAmo/edit?usp=sharing)

## Step 9. Backup

## Description

NOTE FOR BILL: the mask_genomes_toga_copy.py script needs to be finished (everything after line 170 needs to be reviewed). Already run RModeler on sorAra -- pickup where left off. 
Remember to call source activate mask_genomes before starting for necessary dependencies

This is a python script to mask genomes for annotation in TOGA. Briefly, this script will run RepeatModeler to identify De Novo repeats on the genome and use that genome-specific library to mask the genome with RepeatMasker.
:::warning
:bulb: <span style="color:red">**Important**</span>: This script is written to be run on Texas Tech University cluster, using the "nocona" partition, that works in a SLURM environment.
:::

### Software required:
*  [RepeatModeler](https://github.com/Dfam-consortium/RepeatModeler/blob/master/RepeatModeler)
*  [RepeatMasker](https://www.repeatmasker.org/)
*  Python

### Installation 
module load docker/singularity
curl -sSLO https://github.com/Dfam-consortium/TETools/raw/master/dfam-tetools.sh
chmod +x dfam-tetools.sh 
./dfam-tetools.sh --singularity

#### Usage

##call the conda environment that we made for this 
called mask_genomes before running

To run the mask_genomes_toga.py script just type the command as follows

`python /home/tlama_umass_edu/scratch/project_shrew_genome/scripts/mask_genomes_toga_copy.py -qy mSorAra1.pri.asm.20211103.fasta -c 4 -o sorAra_masked`

```
nohup python -u mask_genomes_toga.py -qy genome.fa -c <number of cores> > mask.log &
```


```python=
##########################################################INFO###############################################################################################
#  mask_genomes_toga.py			
#
#  Setup script for mask genomes with Repeat modeler and Repeat masker.
#
#  Author: Diana Moreno
#  v.1.0 November 2019
#  v.2.0 January 2020 --> add xsmall and lib as Flags.
#  v.3.0 April 2021 --> set up for TOGA and add David Ray scripts.
#
#  Description:
#
#	This script will perform the following steps for TOGA functional annotation:
#		I.- Identify de novo repeats with RepeatModeler.
#        The output for this step will be stored at 1_Repeat_modeler directory
#	   II.- Mask  genomes with RepeatMasker, using the RepeatModeler output.
#  
#    syntax: python variant_call_annotation.py
#           -qy <path/to/query/assembly.fa> \
#           -rf <path/to/reference/genome/sample.fa> \
#           -px <prefix for output files> \
#           -c <# number of cores to use> \
#           -o <path to main working directory. all files will be output to subdirectories here> \
#           -ct <# compute processors. should be < nt> \
#           -nt <# data threads. nct x nt = -c>
################################################################################################################################################################



#----------------------------------------------------------------- SETTING ENVIRONMENT -----------------------------------------------------------------

#1)Import modules:

import argparse
import os
import os.path
import subprocess
from subprocess import check_output
import sys 
import time


#2)Set path to software packages:
FATOBIT = ('/lustre/work/daray/software/faToTwoBit')
CLUSTERRUN = ('/home/dmorenos/scripts/slurm_clusterrunC8.py')
REPMODE= ('/lustre/work/daray/software/RepeatModeler') 
REPMASK = ('/lustre/work/daray/software/RepeatMasker-4.1.0') 
AUGUSTUS = ('/lustre/work/daray/software/augustus/bin')
MAINDIR = os.getcwd() 


#3)Define input arguments:

def get_args():
	parser = argparse.ArgumentParser(description='Input needed for variant call annotation', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-qy', '--query', type=str, help='Path/to/query/assembly.fa with variant call', required=True)
	parser.add_argument('-c', '--cores', type=int, help='Define number of cores to use for multi-threaded process', default=1, required=True) 
	parser.add_argument('-o', '--outdir', type=str, help='Define the path/to/output/directory', default='./') 
	parser.add_argument('-sp', '--species', type=str, help='Define the repeat library for RepeatMasker search.', required=False)
	parser.add_argument('-xsmall', help='Select a RepeatMasker masking option as lowercase bases [-xsmall], default is to mask as Ns. ', required=False, action='store_true')
	args = parser.parse_args()
	QRY = args.query 
	CORES = args.cores
	OUTDIR = args.outdir
	SPECIES=args.species
	XSMALL=args.xsmall
	return QRY, CORES, OUTDIR, SPECIES, XSMALL

QRY, CORES, OUTDIR, SPECIES, XSMALL = get_args() 

QUERY = os.path.abspath(QRY) #Set the complete path for query fasta file. Due to during all the scripts we are switching directories, by assigning the full path we avoid errors. 
PREFIX = os.path.basename(QRY).split(".")[0]

#4)Sanity check

print("The query genome is {}\n".format(QRY))
print("The query genome is located in {}\n".format(QUERY))
print("Output files will be generated with {} prefix".format(PREFIX))
print("For this run you are going to use {} cores".format(CORES))
print("Results will be saved on {}".format(OUTDIR))

#5)Create working directories for each task
print ('Creating main working directory')
os.mkdir('1_masking')
os.chdir('1_masking')

#***I. RepeatModeler***

os.mkdir('RepeatModeler')
os.chdir('RepeatModeler')

#  **Ia.- Create qsub for Repeatmodeler

REPMODENAME=(PREFIX + '_RM_slurm.sh')

with open(REPMODENAME, 'w') as f:
        f.write('#!/bin/bash' + '\n')
        f.write('#SBATCH --job-name=' + PREFIX + '_RM' + '\n')
        f.write('#SBATCH --output=%x.%j.out' + '\n')
        f.write('#SBATCH --error=%x.%j.err' + '\n')
        f.write('#SBATCH --partition=nocona' + '\n')
        f.write('#SBATCH --nodes=1' + '\n')
        f.write('#SBATCH --ntasks=' + str(CORES) + '\n')
        f.write('\n' + '\n')
        f.write('module load gcc/9.2.0 ncbi-rmblastn/2.9.0ls' + '\n')
        f.write('echo "Build database using"' + QUERY + 'as a genome reference' + '\n' + '\n')
        f.write((REPMODE) + '/BuildDatabase -name ' + PREFIX + '_repmod ' + '-engine ncbi ' + QUERY + '\n')
        f.write('echo "Running repeat modeler..."' + '\n' + '\n')
        f.write((REPMODE) + '/RepeatModeler -database ' + PREFIX + '_repmod -engine ncbi -pa ' + str(CORES) + '\n')
        f.write('mv RM* RM_output' + '\n')
        f.write('rm -r RM_output/round-*' + '\n')
        f.write('echo "Repeat modeler finished on $ (date)" > repmod_done.txt' + '\n')
        
#  **Ib.-  Excecute sbatch to run RepeatModeler and safety check points.

subprocess.check_call(['sbatch', REPMODENAME])

print ("Running repeat modeler, this might take a while...")

while not os.path.exists('repmod_done.txt'):  #Whit this loop Python is paused until RepeatModeler finish.
        time.sleep(1)
        if os.path.exists('repmod_done.txt'):
                print ('Repeat Modeler finished')

#  **IIIb.-  Check if the output file required for RepeatMasker exists, if not the script will stop.
#    NOTE: This code is going to work only if you are located in the RepeatModeler directory.

try:
        open('./RM_output/consensi.fa.classified')
except IOError:  #file not found
        print ('File not found, Repeat Modeler failed. \n This Python script will die!')
        sys.exit()

REPMODOUT = (MAINDIR + '/1_masking/RepeatModeler/RM_output/consensi.fa.classified') #Assign the output from Repeat Modeler to a path for further analysis.


#II. REPEAT MASKER -----------------------------------------------------------------

#  **IIa.- Run RepeatMasker using the denovo library built with Repeat Modeler.
#        - IMPORTANT: This script uses a python program. You need to make sure python is loaded and ready to go by activating conda.
#        -The following script was developed by David A. Ray

os.chdir(MAINDIR) #Moving from /1_Repeat_modeler directory to /2_Repeat_masker.
os.chdir('1_masking')
os.mkdir('RepeatMasker')
os.chdir('RepeatMasker')

RMSLURM=('Repeatmasker_slurm.sh')

with open(RMSLURM, 'w') as f:
    f.write('#!/bin/bash' + '\n')
    f.write('#SBATCH --job-name=' + PREFIX + '_repmask' + '\n')
    f.write('#SBATCH --output=%x.%j.out' + '\n')
    f.write('#SBATCH --error=%x.%j.err' + '\n')
    f.write('#SBATCH --partition=nocona' + '\n')
    f.write('#SBATCH --nodes=1' + '\n')
    f.write('#SBATCH --ntasks=1' + '\n')
    f.write('\n')
    f.write('source ~/conda/etc/profile.d/conda.sh' + '\n')
    f.write('conda activate' + '\n' + '\n') 
    f.write('echo makes variables used in this script' + '\n')
    f.write('GENOME=' + QRY + '\n' + '\n')
    f.write('RUNTYPE=' + PREFIX + '_RM' + '\n')
    f.write('DIR=' + MAINDIR + '/1_masking/RepeatMasker/$RUNTYPE' + '\n')
    f.write('mkdir -p $DIR' + '\n')
    f.write('cd $DIR' + '\n' + '\n')
    f.write('ln -s ' + QUERY + '\n' + '\n')
    f.write('echo "fa to 2bit"'+ '\n')
    f.write(FATOBIT + ' ' + QRY + ' ' + PREFIX + '.2bit' + '\n' + '\n')
    f.write('echo "Use slurm_clusterrunC8.py to generate all batches needed to run RepeatMasker and the doLift.sh to complile the results."'+ '\n')
    f.write('python3 ' + CLUSTERRUN + ' -i ' + QRY + ' -b 50 -lib ' + REPMODOUT + ' -dir . -xsmall' + '\n')
    f.write('echo "Submits the RepeatMasker jobs created by slurm_clusterrunC8.py"'+ '\n')
    f.write('sh qsub.sh'+ '\n')
    f.write('#Creates a list of jobIDs to keep track of the RepeatMasker batches being run.' + '\n')
    f.write('jobIDs=""; for i in `squeue  | grep ' + PREFIX + """ | awk '{print $1}'`; do jobIDs=$jobIDs,$i; done; jobIDs=\"${jobIDs:1}\" """ + '\n')
    f.write('#Submits the doLift script to the queue but holds it until all jobs in the jobIDs list (the RepeatMasker batches) have cleared.' + '\n')
    f.write('sbatch -d afterok:$jobIDs doLift.sh' + '\n')
    f.write('#Creates a list of job IDs to keep track of the doLift script.' + '\n')
    f.write('jobIDs=""; for i in `squeue  | grep ' + PREFIX + """ | awk '{print $1}'`; do jobIDs=$jobIDs,$i; done; jobIDs=\"${jobIDs:1}\" """ + '\n')
    f.write('#Submits the rm2bed script to the queue but holds it until all jobs in the jobIDs list (the doLift job) have cleared.' + '\n')
    f.write ('sbatch -d afterok:$jobIDs ' + MAINDIR + '/1_masking/RepeatMasker/rm2bed.sh' + '\n')
  
#Generate rm2bed.sh file

RM2BED=('rm2bed.sh')

with open(RM2BED, 'w') as f:
    f.write('#!/bin/bash' + '\n')
    f.write('#SBATCH --job-name=' + PREFIX + '_RM2' + '\n')
    f.write('#SBATCH --output=%x.%j.out' + '\n')
    f.write('#SBATCH --error=%x.%j.err' + '\n')
    f.write('#SBATCH --partition=nocona' + '\n')
    f.write('#SBATCH --mem-per-cpu=100G' + '\n')
    f.write('#SBATCH --nodes=1' + '\n') 
    f.write('#SBATCH --ntasks=1' + '\n')
    f.write('\n')
    f.write('source ~/conda/etc/profile.d/conda.sh' + '\n')
    f.write('conda activate' + '\n') 
    f.write('echo makes variables used in this script' + '\n')
    f.write('GENOME=' + PREFIX + '\n')
    f.write('RUNTYPE=' + PREFIX + '_RM' + '\n')
    f.write('DIR=' + MAINDIR + '/1_masking/RepeatMasker/$RUNTYPE' + '\n')
    f.write('cd $DIR' + '\n')
    f.write('#Runs a python script RM2bed to generate one complete .bed file and several subfiles subdivided by TE class. Merges overlapping hits based using lower_divergence criterion.' + '\n')
    f.write('[ ! -f ${GENOME}_rm.bed ] && python /home/daray/gitrepositories/bioinfo_tools/RM2bed.py -d . -sp class -p ${GENOME} -o lower_divergence ${GENOME}.fa.align.gz' + '\n' )
    f.write('module load gcc/10.1.0 bedtools2/2.29.2' + '\n')
    f.write('RMPATH=/lustre/work/daray/software/RepeatMasker-4.1.0/util' + '\n')
    f.write('gunzip -c ${GENOME}.fa.out.gz > {GENOME}.fa.out' + '\n')
    f.write('perl $RMPATH/rmOutToGFF3.pl {GENOME}.fa.out > {GENOME}.gff' + '\n')
    f.write('bedtools maskfasta -soft -fi {GENOME}.fa - BED {GENOME}.gff -fo {GENOME}.masked.fa' + '\n')
    f.write( 'echo "Masking genome done!"'+ '\n')
   
subprocess.check_call (['sbatch', 'Repeatmasker_slurm.sh']) 
```


The final output of this genome will generate a masked genome: ***prefix.masked.fa***, that will be used for TOGA annotation.

:::info
:pencil: <span style="color:red">**Note**</span>: In this script we are only using the repeatmodeler library as evidence for RepeatMasker, to match our TOGA annotation with Ariadna's. For future annotations, is better to mask the genomes with the customize libraries generated by David Ray lab. 
:::


###### tags: `scripts` `project_tyd`


Sending doblastz jobs to login1 node
localhost issue 
cd ~/.ssh 
create a public/private key
add the public key to your account on unity as @login1 and @login2
login to unity and try ssh login1 and ssh login2

