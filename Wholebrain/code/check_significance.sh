designDir=$1	# which analysis to check 

d0=$PWD

# report FWE significance
echo $designDir
cd $d0/../results/$designDir
echo "+++ FWE corrected"
for ii in stats*corrp*; do echo $ii

	maxp=( $(fslstats $ii -R) )

	if (( $(echo "${maxp[1]} > 0.95" |bc -l) )); then
		echo "FWE p<0.05 "
	else
		echo "ns"
	fi

done

cd $d0
