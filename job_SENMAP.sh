#!/bin/bash

#SBATCH -o SENMAP.out
#SBATCH -e SENMAP.err
#SBATCH --mail-type END
#SBATCH --mail-user simon.orozco.arias@gmail.com
#SBATCH -J SENMAP
#SBATCH --time 3-00:00:00
#SBATCH --partition gpu
#SBATCH -n 2
#SBATCH -N 1
#SBATCH --gres=gpu:1
#SBATCH --mem 200GB

source ~/.bashrc
source activate tf23

time python3 SENMAP.py -l GCA_002050065.1_CA_8.1_PBcR_genomic.fna.mod.EDTA.TElib.fa.ltr -o Dmel_EDTA_raw

time python3 SENMAP.py -l Dmelanogaster_refTEs.fa_short.ltr -o Dmel_REPET_raw
