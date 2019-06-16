if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    #from CRABClient.UserUtilities import config
    #config = config()
    from WMCore.Configuration import Configuration
    def resubmit(directory):
        try:
            crabCommand('resubmit', directory)
        except HTTPException as hte:
            print "Failed resubmitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed resubmitting task: %s" % (cle)


    def kill(directory):
        try:
            crabCommand('kill', directory)
        except HTTPException as hte:
            print "Failed killing task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed killing task: %s" % (cle)



    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################
    ev = 100000
    mh_list = [125,300, 500, 1000, 2000]
    pl_list = [500, 1000, 10000]
    mode_list = ["ppTohToSS1SS2_SS1Tobb_SS2Toveve_ggh_withISR","ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR","ppTohToSS1SS2_SS1Tobb_SS2Toveve_vbfh_withISR","ppTohToSS1SS2_SS1Tobb_SS2Tobb_ggh_withISR","ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR","ppTohToSS1SS2_SS1Tobb_SS2Tobb_vbfh_withISR"]
    mode_list = ["ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR", "ppTohToSS1SS2_SS1Tobb_SS2Tobb_vbfh_withISR"]

    for i in range(len(mode_list)):
	mode = mode_list[i]
        index = mode.find("with")
        name = mode[:index]+mode[index+4:]
	for mh in mh_list:
	    mx = (mh-50)/2
	    if mh == 125:
		mx = 50
	    for pl in pl_list:
    		directory = 'crab/crab_CMSSW_9_4_4_'+name+"_mh{}_mx{}_pl{}_ev{}".format(mh,mx,pl,ev)+'_CaltechT2'
		resubmit(directory)
