#!/bin/bash
#SBATCH -J raxml
#SBATCH -n 10
#SBATCH -N 5
#SBATCH -p normal
#SBATCH -o raxml.o%j
#SBATCH -e raxml.e%j
#SBATCH -t 12:00:00
#SBATCH -A Recombining-sex-chro
#SBATCH --mail-type=ALL
#SBATCH --mail-user=grovesdixon@gmail.com

module load launcher





export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=raxml

$LAUNCHER_DIR/paramrun
