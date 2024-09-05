#PBS -A UPSU0068
#PBS -N vr_1986-2014
#PBS -l walltime=15:00:00
#PBS -M ajb8224@psu.edu
#PBS -l select=1:ncpus=10:mem=20GB
#PBS -l job_priority=economy
#PBS -m ae
#PBS -q casper

# Make sure to match numcores to ncpus above
syr=1986
eyr=2014
numcores=10
cmd_file="commands_${syr}_${eyr}.txt"

source ~/.bashrc
conda activate anapy

# Create the command file
touch ${cmd_file}

# Loop through years and append commands to the command file
for yr in `seq $syr $eyr`
do
    echo "python3 python_wrapper.py $yr" >> ${cmd_file}
done

# Execute commands in parallel using GNU Parallel
parallel --jobs ${numcores} --workdir $PWD < ${cmd_file}

# Clean up command file
rm -v ${cmd_file}
