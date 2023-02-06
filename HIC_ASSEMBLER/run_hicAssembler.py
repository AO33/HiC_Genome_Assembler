###############
### Imports ###
import argparse
import time
import sys
##########


def readConfigFileToVariables(configFile):
	'''This function reads in the configFile and fills in the values of a variable dictionary
	   that is then used as input to the differnt parts of the pipeline. All variables must have a non '' value (empty value)
	   Input = full file path to configuration file
	   Output = { variableName:value }'''
	varDict = { "resolution":'',
				"saveFilesDirectory":'',
				"savePlotsDirectory":'',
				"hicProBedFile":'',
				"hicProBiasFile":'',
				"hicProMatrixFile":'',
				"hicProScaffSizeFile":'',
				"dendrogramOrderFile":'',
				"avgClusterPlot":'',
				"avgClusterPlot_outlined":'',
				"binGroupFile":'',
				"assessmentFile":'',
				"minSize":20,
				"modularity":.05,
				"convergenceRounds":8,
				"louvainRounds":20,
				"chromosomeGroupFile":'',
				"chromosomeOrderFile":'',
				"chromosomePlotSuffix":'',
				"fullGenomePlot":'',
				"fullGenomePlotTitle":'',
				"plotOrderFile":'',
				"nScaffolds":6,
				"scanScaffolds":5,
				"lengthCutoff":500000,
				"restrictionSiteFile":'',
				"validPairFile":'',
				"finalOrderingsFile":'',
				"originalFastaFile":'',
				"assembledFastaFile":''}
	#############################   
	inFile = open(configFile,'r')
	for line in inFile:
		line = line.strip('\r').strip('\n')
		##############
		if line == '':
			continue
		##################
		if line[0] == "#":
			continue
		#############################################################
		arg,val = str(line.split(' = ')[0]),str(line.split(' = ')[1])
		##################################################################
		### Variables that are used in 2 or more parts of the pipeline ###
		if arg == "resolution" and val:
			try:
				val = int(val)
				varDict["resolution"] = val
			except:
				print("ERROR... resolution must be a an integer value equal to the resolution of the contact map used. Exiting...")
				sys.exit()
		if arg == "saveFilesDirectory" and val:
			varDict["saveFilesDirectory"] = val
		if arg == "savePlotsDirectory" and val:
			varDict["savePlotsDirectory"] = val
		if arg == "hicProBedFile" and val:
			varDict["hicProBedFile"] = val
		if arg == "hicProBiasFile" and val:
			varDict["hicProBiasFile"] = val
		if arg == "hicProMatrixFile" and val:
			varDict["hicProMatrixFile"] = val
		if arg == "hicProScaffSizeFile" and val:
			varDict["hicProScaffSizeFile"] = val
		if arg == "chromosomeGroupFile" and val:
			varDict["chromosomeGroupFile"] = varDict["saveFilesDirectory"]+'/'+val
		if arg == "chromosomeOrderFile" and val:
			varDict["chromosomeOrderFile"] = varDict["saveFilesDirectory"]+'/'+val
		if arg == "finalOrderingsFile" and val:
			varDict["finalOrderingsFile"] = varDict["saveFilesDirectory"]+'/'+val
		####################################################
		### Variables only used in part1 of the pipeline ###
		if arg == "dendrogramOrderFile" and val:
			varDict["dendrogramOrderFile"] = varDict["saveFilesDirectory"]+'/'+val
		if arg == "avgClusterPlot" and val:
			varDict["avgClusterPlot"] = varDict["savePlotsDirectory"]+'/'+val
		if arg == "avgClusterPlot_outlined" and val:
			varDict["avgClusterPlot_outlined"] = varDict["savePlotsDirectory"]+'/'+val
		if arg == "binGroupFile" and val:
			varDict["binGroupFile"] = varDict["saveFilesDirectory"]+'/'+val
		if arg == "assessmentFile" and val:
			varDict["assessmentFile"] = varDict["saveFilesDirectory"]+'/'+val
		if arg == "minSize" and val:
			try:
				val = int(val)
				varDict["minSize"] = val
			except:
				print("WARNING... minSize must be an integer value... setting minSize=20 which is the default value")
		if arg == "modularity" and val:
			try:
				val = float(val)
				if val > 1.:
					print("WARNING... modularity must be a value between 0.0 and 1.0... setting modularity=.05 which is the default value")
					val = .05
				varDict["modularity"] = val
			except:
				print("WARNING... modularity must be a floating point value... setting modularity=.05 which is the default value")
		if arg == "convergenceRounds" and val:
			try:
				val = int(val)
				varDict["convergenceRounds"] = val
			except:
				print("WARNING... convergenceRounds must be an integer value... setting convergenceRounds=8 which is the default value")
		if arg == "louvainRounds" and val:
			try:
				val = int(val)
				varDict["louvainRounds"] = val
			except:
				print("WARNING... louvainRounds must be an integer value... setting louvainRounds=20 which is the default value")
		####################################################
		### Variables only used in part2 of the pipeline ###
		if arg == "chromosomePlotSuffix" and val:
			varDict["chromosomePlotSuffix"] = val
		if arg == "fullGenomePlot" and val:
			varDict["fullGenomePlot"] = varDict["savePlotsDirectory"]+'/'+val
		if arg == "fullGenomePlotTitle" and val:
			varDict["fullGenomePlotTitle"] = val
		if arg == "plotOrderFile" and val:
			varDict["plotOrderFile"] = varDict["saveFilesDirectory"]+'/'+val
		if arg == "nScaffolds" and val:
			try:
				val = int(val)
				varDict["nScaffolds"] = val
			except:
				print("WARNING... nScaffolds must be an integer value... setting nScaffolds=6 which is the default value")
		if arg == "scanScaffolds" and val:
			try:
				val = int(val)
				varDict["scanScaffolds"] = val
			except:
				print("WARNING... scanScaffolds must be an integer value... setting scanScaffolds=5 which is the default value")
		####################################################
		### Variables only used in part3 of the pipeline ###
		if arg == "lengthCutoff" and val:
			try:
				val = int(val)
				varDict["lengthCutoff"] = val
			except:
				print("WARNING... lengthCutoff must be an integer value... setting lengthCutoff=500000 which is the default value")
		if arg == "restrictionSiteFile" and val:
			varDict["restrictionSiteFile"] = val
		if arg == "validPairFile" and val:
			varDict["validPairFile"] = val
		####################################################
		### Variables only used in part4 of the pipeline ###
		if arg == "originalFastaFile" and val:
			varDict["originalFastaFile"] = val
		if arg == "assembledFastaFile" and val:
			varDict["assembledFastaFile"] = varDict["saveFilesDirectory"]+'/'+val
	##############
	inFile.close()
	##############
	return varDict

def ensureAllVariablesAreSet(varDict):
	exitFail = 0
	unsetVariables = []
	for key,val in varDict.items():
		if val == '':
			exitFail = 1
			unsetVariables.append(key)
	#################
	if exitFail == 1:
		print("The following variable(s) do not have any value assossicated with them. Please set this variables to continue...")
		for v in unsetVariables:
			print(v)
		print("Exiting...")
		return True
	else:
		return False

if __name__ == "__main__":
	###############
	## ARG PARSE ##
	###############
	def mm():
		parser = argparse.ArgumentParser(description='Runs the various parts of the HiC assembly pipeline. \
													  Each Part requires the previous Part(s) to be run before hand. \
													  Each Part can be run independantly or sequentially and any combination of Part(s)1-4 is allowed. \
													  Example, -part1 -part2 -part3 -part4 will run the entire pipeline. \
													  Example, -part2 -part3 will just run parts2 and 3 (assuming part1 was run beforehand)')
		##################################################################################
		parser.add_argument("-part1",help="Run part1 of the pipeline",action='store_true')
		parser.add_argument("-part2",help="Run part2 of the pipeline",action='store_true')
		parser.add_argument("-part3",help="Run part3 of the pipeline",action='store_true')
		parser.add_argument("-part4",help="Run part4 of the pipeline",action='store_true')
		parser.add_argument("-config",help="Full file path to the config file. All arguments must have a value in the config file or the program will exit", required=True, type=str)
		return parser.parse_args()
	###########
	args = mm()
	startTime = time.time()
	################################################
	varDict = readConfigFileToVariables(args.config)
	exitStatus = ensureAllVariablesAreSet(varDict)
	if exitStatus == True:
		sys.exit()
	##############
	if args.part1:
		import scaffoldToChromosomes as part1
		part1.runPipeline(varDict["hicProBedFile"],varDict["hicProBiasFile"],varDict["hicProMatrixFile"],varDict["hicProScaffSizeFile"],
						  varDict["dendrogramOrderFile"],varDict["avgClusterPlot"],varDict["avgClusterPlot_outlined"],
						  varDict["binGroupFile"],varDict["assessmentFile"],varDict["chromosomeGroupFile"],
						  varDict["minSize"],varDict["modularity"],varDict["convergenceRounds"],varDict["louvainRounds"],varDict["resolution"])
	##############
	if args.part2:
		import orderGenome as part2
		part2.runPipeline(varDict["hicProBedFile"],varDict["hicProBiasFile"],varDict["hicProMatrixFile"],varDict["chromosomeGroupFile"],varDict["chromosomeOrderFile"],
						  varDict["savePlotsDirectory"],varDict["chromosomePlotSuffix"],varDict["fullGenomePlot"],varDict["fullGenomePlotTitle"],varDict["plotOrderFile"],
						  varDict["nScaffolds"],varDict["scanScaffolds"],varDict["resolution"])
	##############
	if args.part3:
		import orientSmallScaffolds as part3
		part3.runPipeline(varDict["chromosomeOrderFile"],varDict["hicProScaffSizeFile"],
						  varDict["restrictionSiteFile"],varDict["validPairFile"],varDict["finalOrderingsFile"],
						  varDict["lengthCutoff"],varDict["resolution"])
	##############
	if args.part4:
		import writeAssembledFasta as part4
		part4.runPipeline(varDict["originalFastaFile"],
						  varDict["finalOrderingsFile"],
						  varDict["assembledFastaFile"])
	################################################################
	print("Total run-time = "+str(time.time()-startTime)+" seconds")
