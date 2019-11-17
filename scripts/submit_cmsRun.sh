#!/bin/bash
mkdir -p log
mkdir -p submit

echo `pwd`

cd ../
RazorAnalyzerDir=`pwd`
cd -

job_script=${CMSSW_BASE}/src/cms_lpc_llp/llp_ntupler/scripts/runRazorJob_llp_vH.sh
echo $job_script

#year=RunIISummer16
year=Fall18
#samples=WminusH_HToSSTobbbb_ms55_pl10000_ev150000_${year}
samples=WplusH_HToSSTobbbb_ms55_pl10000_ev150000_RunII${year}
filesPerJob=1
for sample in ${samples}
do
	echo "Sample " ${sample}
	version=/v9/
	echo "${version}"
	#output=/store/group/phys_exotica/privateProduction/ntuple/RunIISummer16/WminusH_HToSSTobbbb_ms55_pl10000/${version}/${sample}
	output=/store/group/phys_exotica/privateProduction/ntuple/RunII${year}/${sample}/${version}/${sample}

	echo "output ${output}"
	inputfilelist=/src/cms_lpc_llp/llp_ntupler/lists/${sample}.txt
	nfiles=`cat ${CMSSW_BASE}$inputfilelist | wc | awk '{print $1}' `
        maxjob=`python -c "print int($nfiles.0/$filesPerJob)"`
	mod=`python -c "print int($nfiles.0%$filesPerJob)"`
        if [ ${mod} -eq 0 ]
        then
                maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
        fi
	#config=displacedJetMuon_step2_rechit_studies_cfg_condor
	config=displacedJetMuon_step2_ntupler_MC_${year}_cfg_condor
	echo "${config}"
	rm -f submit/${config}_{sample}_Job*.jdl
	rm -f log/${config}_${sample}_Job*
	for jobnumber in `seq 0 1 ${maxjob}`
	do
		echo "job " ${jobnumber} " out of " ${maxjob}
		jdl_file=submit/${config}_${sample}_Job${jobnumber}_Of_${maxjob}.jdl
		echo "Universe = vanilla" > ${jdl_file}
		echo "Executable = ${job_script}" >> ${jdl_file}
                echo "Arguments = ${config}.py ${inputfilelist} ${filesPerJob} ${jobnumber} ${sample}_Job${jobnumber}_Of_${maxjob}.root ${output} ${CMSSW_BASE} ${HOME}/ ${year}" >> ${jdl_file}
		
		# option should always be 1, when running condor
		echo "Log = log/${config}_${sample}_Job${jobnumber}_Of_${maxjob}_\$(Cluster).\$(Process).log" >> ${jdl_file}
		echo "Output = log/${config}_${sample}_Job${jobnumber}_Of_${maxjob}_\$(Cluster).\$(Process).out" >> ${jdl_file}
	        echo "Error = log/${config}_${sample}_Job${jobnumber}_Of_${maxjob}_\$(Cluster).\$(Process).err" >> ${jdl_file}

		#echo "Requirements=(TARGET.OpSysAndVer==\"CentOS7\" && regexp(\"blade.*\", TARGET.Machine))" >> ${jdl_file}
		echo "RequestMemory = 2000" >> ${jdl_file}
		echo "RequestCpus = 1" >> ${jdl_file}
		echo "RequestDisk = 4" >> ${jdl_file}
		echo "+RunAsOwner = True" >> ${jdl_file}
		echo "+InteractiveUser = true" >> ${jdl_file}
		if [ ${year} == "Summer16" ]
		then
			echo "+SingularityImage = \"/cvmfs/singularity.opensciencegrid.org/bbockelm/cms:rhel7\"" >> ${jdl_file}
		else
			echo "+SingularityImage = \"/cvmfs/singularity.opensciencegrid.org/bbockelm/cms:rhel6\"" >> ${jdl_file}
		fi
		echo '+SingularityBindCVMFS = True' >> ${jdl_file}
#		echo "transfer_input_files = tarball/${sample}.tar" >> ${jdl_file}
		echo "run_as_owner = True" >> ${jdl_file}
		echo "x509userproxy = ${HOME}/x509_proxy" >> ${jdl_file}
		echo "should_transfer_files = YES" >> ${jdl_file}
		echo "when_to_transfer_output = ON_EXIT" >> ${jdl_file}
		echo "Queue 1" >> ${jdl_file}
		echo "condor_submit ${jdl_file}"
		condor_submit ${jdl_file}
	done
done
