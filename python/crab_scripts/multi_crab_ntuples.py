if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    from WMCore.Configuration import Configuration
    config = Configuration()

    config.section_("General")
    config.General.workArea = 'crab'
    config.General.transferOutputs = True
    config.General.transferLogs = True

    config.section_("JobType")
    config.JobType.pluginName = 'Analysis'
    config.JobType.psetName = '/afs/cern.ch/work/c/christiw/public/LLP/CMSSW_9_4_4/src/cms_lpc_llp/llp_ntupler/python/llp_ntupler_MC_Summer16.py'
    config.JobType.numCores = 1
    config.section_("Data")
    config.Data.inputDBS = 'phys03'
    config.Data.splitting = 'FileBased'
    config.Data.unitsPerJob = 10 #when splitting is 'Automatic', this represents jobs target runtime(minimum 180)
    config.Data.publication = True
    config.Data.ignoreLocality = True

    config.section_("Site")
    config.Site.storageSite = 'T2_US_Caltech'
    config.Site.whitelist = ['T2_US_Caltech']
    config.Site.ignoreGlobalBlacklist = True
    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################
    ev = 100000
    mh_list = [125,300, 500, 1000, 2000]
    pl_list = [500, 1000, 10000]
    mode_list = ["ppTohToSS1SS2_SS1Tobb_SS2Toveve_ggh_withISR","ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR","ppTohToSS1SS2_SS1Tobb_SS2Toveve_vbfh_withISR","ppTohToSS1SS2_SS1Tobb_SS2Tobb_ggh_withISR","ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR","ppTohToSS1SS2_SS1Tobb_SS2Tobb_vbfh_withISR"]
    mode_list = ["ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR","ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR"]
    for i in range(len(mode_list)):
        mode = mode_list[i]
        index = mode.find("with")
        name = mode[:index]+mode[index+4:]
        print(mode,name)
        for mh in mh_list:
            mx = (mh-50)/2
            if mh == 125:
                mx = 50
            for pl in pl_list:
                spec = mode+"_mh{}_mx{}_pl{}_ev{}".format(mh,mx,pl,ev)
                config.General.requestName = 'CMSSW_9_4_4_'+name+"_mh{}_mx{}_pl{}_ev{}".format(mh,mx,pl,ev)+'_CaltechT2_v1'
                config.Data.inputDataset = '/'+spec +'/christiw-crab_CMSSW_8_0_21_'+name+'_mh{}_mx{}_pl{}_ev{}'.format(mh,mx,pl,ev)+'_AOD_CaltechT2_v1-b1a4edca9adfa7a2e4059536bf605cd7/USER'
                config.Data.outLFNDirBase = '/store/group/phys_exotica/delayedjets/llpntuple/V1p0/MC_Summer16/v1/christiw/'+spec+'/'
                #print(config.Data.outLFNDirBase)
                submit(config)

                                                                        
