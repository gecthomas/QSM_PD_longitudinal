flag=$1			# group (PD)
design=$2		# e.g. fup_CCog_age_sex_fupinterval
V=$3			# V1 / V3 / V1 + V3
Nperm=$4		# 10000

d0=$PWD

echo +++ Permutation statistics

# brainmask for analysis
msk=brainmask_${flag}_${V}
# 4D smoothed data file
data=all_${flag}_cs3_abs_brain_wQSM_${V}
# design matrix names
designMat=${flag}_${design}
# output directory / input directories
inputDir=${d0}/../data/
designDir=${inputDir}/design_matrices/
outputDir=${d0}/../results/${flag}_${V}_${design}/

if [ ! -d $outputDir ]; then
	mkdir $outputDir
fi

if [ -f ${designDir}/${designMat}.grp ] 
then
	randomise -i ${inputDir}/${data}.nii.gz \
		  -o ${outputDir}/stats_${data}_${design} \
		  -m ${inputDir}${msk}.nii.gz \
		  -d ${designDir}/${designMat}.mat \
		  -t ${designDir}/${designMat}.con \
		  -e ${designDir}/${designMat}.grp \
		  -n $Nperm \
		  -T -D --uncorrp
else
	randomise -i ${inputDir}/${data}.nii.gz \
		  -o ${outputDir}/stats_${data}_${design} \
		  -m ${inputDir}${msk}.nii.gz \
		  -d ${designDir}/${designMat}.mat \
		  -t ${designDir}/${designMat}.con \
		  -n $Nperm \
		  -T -D --uncorrp
fi

