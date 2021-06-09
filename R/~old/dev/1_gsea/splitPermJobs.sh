#Preparing 20 plink files with 5 permutations per file
for i in `seq 5 5 100`;
        do cat plinkJobs.txt | head -${i} | tail -5 > plinkJobs${i}.txt;

        jobFile=plinkJobs${i}.txt.sh
        jobDir=/scratch/g/gbader/cmross/out_160411/gwas

#Writing each plinkJob file       
        cat /dev/null > $jobFile;
echo "#!/bin/bash" >> $jobFile
echo "#PBS -l nodes=1:ppn=8,walltime=01:00:00" >> $jobFile;
echo "module load gnu-parallel/20140622" >> $jobFile;
echo "parallel -j 2 < ${jobDir}/plinkJobs${i}.txt" >> $jobFile;
echo ""                         >> $jobFile;

#Submitting the files to the debug queue
chmod u+x $jobFile
qsub $jobFile -q debug

done
