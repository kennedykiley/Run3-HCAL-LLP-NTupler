#!/bin/sh
mkdir -p log
mkdir -p submit

echo `pwd`

cd ../
RazorAnalyzerDir=`pwd`
cd -

job_script=${CMSSW_BASE}/src/cms_lpc_llp/llp_ntupler/scripts/runRazorJob_llp_vH.sh
echo $job_script
samples=WplusH_HToSSTobbbb_ms55_pl10000_ev150000
filesPerJob=1
for sample in ${samples}
do
	echo "Sample " ${sample}
	version=/v6/
	echo "${version}"
	output=/store/group/phys_exotica/privateProduction/ntuple/RunIIFall17/WplusH_HToSSTobbbb_ms55_pl10000/${version}/${sample}
	inputfilelist=/src/cms_lpc_llp/llp_ntupler/lists/${sample}.txt
	nfiles=`cat ${CMSSW_BASE}$inputfilelist | wc | awk '{print $1}' `
        maxjob=`python -c "print int($nfiles.0/$filesPerJob)"`
	mod=`python -c "print int($nfiles.0%$filesPerJob)"`
        if [ ${mod} -eq 0 ]
        then
                maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
        fi
	#config=displacedJetMuon_step2_rechit_studies_cfg_condor
	config=displacedJetMuon_step2_ntupler_cfg_condor
	echo "${config}"
	rm -f submit/${config}_{sample}_Job*.jdl
	rm -f log/${config}_${sample}_Job*
	for jobnumber in `seq 0 1 ${maxjob}`
	do
            jdl_file=submit/${config}_${sample}_Job${jobnumber}_Of_${maxjob}.jdl
	    outRoot="/mnt/hadoop/${output}/${sample}_Job${jobnumber}_Of_${maxjob}.root"

            minimumsize=1000
            actualsize=0
            if [ -f ${outRoot} ]
            then
                    actualsize=$(wc -c <"${outRoot}")
            fi
            if [ $actualsize -ge $minimumsize ]
            then
                    finished=yes
            else
                    echo "job ${config}_${sample}_Job${jobnumber}_Of_${maxjob} failed, now being resubmitted"
                    rm -f log/${config}_${sample}_Job${jobnumber}_Of_${maxjob}*
            	    condor_submit ${jdl_file}
	    fi




	done
done
