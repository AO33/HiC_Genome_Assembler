###############
### Imports ###
import time
import math
###########

class RestrictionScaffold:
	'''This class of object just stores information about each scaffold we have in our genome assembly.
	   Maily, the restriction enzyme cutsite coordinates for each scaffold so we can normalize read pair counts accordingly.'''
	def __init__(self,name,orientation,size,resCoords,binCount):
		self.name = name
		self.orientation = orientation
		self.size = size
		self.resCoords = resCoords
	def getBinCount(self,resolution):
		self.binCount = math.ceil( float(self.size) / float(resolution) )
	def getResCounts(self,lengthCutoff):
		left,right = 0,0
		for c in self.resCoords:
			if c <= lengthCutoff:
				left += 1
			elif c > (self.size - lengthCutoff):
				right += 1
		#############
		if left == 0:
			left = 1
		if right == 0:
			right = 1
		###################
		self.resLeft = left
		self.resRight = right

def readPreliminaryOrientationFromFile(orientationFile):
	'''This function reads into memory the order and orientation of each scaffold (preliminary orientation for scaffolds containing only a single bin)
	   Input = full file path to the chromosome ordering file produced in the second step of the pipeline (one of the two ouput files)
	   Output = [ [ RestrictionScaffoldObject , ... ] , ... ] list of RestrictionScaffold Objects for each chromosome
				{ scaffoldName:RestrictionScaffoldObject } This dictionary points to the objects contained within the lists for each chromosome found'''
	chrGroups,chrGroup = [],[]
	scaffoldDict = {}
	orientFile = open(orientationFile,'r')
	orientFile.readline()
	for line in orientFile:
		line = line.strip('\r').strip('\n')
		if line[0] == "#":
			chrGroups.append(chrGroup)
			chrGroup = []
		else:
			cols = line.split('\t')
			scaffoldMan = RestrictionScaffold(cols[0],cols[1],0.,[],0)
			scaffoldDict.update({ cols[0]:scaffoldMan })
			chrGroup.append(scaffoldMan)
	##################
	orientFile.close()
	chrGroups.append(chrGroup)
	return chrGroups,scaffoldDict

def readScaffSizeFile(scaffSizeFile,scaffDict,resolution):
	'''Simple function that reads in the HiCpro required file that has one column as the scaffold name,
	   and the other column the size of the scaffold. It then adds the size as an attribute to each RestrictionScaffoldObject.
	   It also adds the binCount attribute to each object by calling the getBinCount method of each object ( scaffoldSize / resolution )
	   Input = full file path to scaffold size file, { scaffoldName:RestrictionScaffoldObject }, resolution of the contact map (size of each genomic loci)
	   Output = updated version of the { scaffoldName:RestrictionScaffoldObject } where size has been added as an attribute'''
	sizeFile = open(scaffSizeFile,'r')
	for line in sizeFile:
		cols = line.strip('\r').strip('\n').split('\t')
		scaffName,scaffSize = cols[0],float(cols[1])
		if scaffName in scaffDict:
			scaffDict[scaffName].size = scaffSize
			scaffDict[scaffName].getBinCount(resolution)
	################
	sizeFile.close()
	return scaffDict

def readRestrictionsFile(restrictionFile,scaffDict):
	'''This function takes in the restriction enzyme file cut site from HiCpro and adds the coordinates to each RestrictionScaffoldObject.
	   Input = full file path the restriction cut site file, { scaffoldName:RestrictionScaffoldObject }
	   Output = { scaffoldName:RestrictionScaffoldObject } with restriction cut site coordinates added to each object (as a list)'''
	resFile = open(restrictionFile,'r')
	for line in resFile:
		cols = line.strip('\r').strip('\n').split('\t')
		scaffName,resCoord = cols[0],int(cols[2])
		if scaffName in scaffDict:
			scaffDict[scaffName].resCoords.append(resCoord)
	###############
	resFile.close()
	###################################
	for scaffObj in scaffDict.values():
		scaffObj.resCoords = sorted(scaffObj.resCoords,key=lambda x: x)
	################
	return scaffDict

def initiateScaffoldObjects(orientationFile,scaffSizeFile,restrictionFile,resolution):
	'''This function initiates our data structures that we will then use to orient the scaffolds <= resolution size.
	   Input = preliminary scaffold order file from previous step,
			   scaffold size file used by HiCPro,
			   restriction enzyme cutsite file used by HiCPro,
			   resolution size (contact map resolution)
	   Output = [ [RestrictionScaffoldObject , ... ] , ... ] list of RestrictionScaffold objects for each chromosome,
				{ scaffoldName:RestrictionScaffoldObject } this dictionary just points to the RestrictionScaffold objects we put in the list'''
	scaffGroups,scaffDict = readPreliminaryOrientationFromFile(orientationFile)
	scaffDict = readScaffSizeFile(scaffSizeFile,scaffDict,resolution)
	scaffDict = readRestrictionsFile(restrictionFile,scaffDict)
	############################
	return scaffGroups,scaffDict

def pullTriplets(scaffoldList):
	'''This function takes in a list of ordered / oriented scaffolds (for a given chromosome),
	   and pulls out scaffolds whose binCount == 1 as well as the adjacent scaffolds.
	   These are what I refer to as triplets, and each triplet is a list. So the return is a list of triplets.
	   Note that edge cases are taken into account (ie if the scaffold is located on the very edge of the chromosome)
	   Input = [ RestrictionScaffoldObject , ... ]
	   Output = [ triplet, triplet , ... ] where each triplet is a list of RestrictionScaffoldObjects '''
	triplets = []
	for i,scaffObj in enumerate(scaffoldList):
		if scaffObj.binCount == 1:
			s1 = scaffoldList[i]
			##########
			if i != 0:
				s0 = scaffoldList[i-1]
			else:
				s0 = "NA"
			################################
			if i <= (len(scaffoldList) - 2):
				s2 = scaffoldList[i+1]
			else:
				s2 = "NA"
			#############################
			if s0 != "NA" and s2 != "NA":
				triplets.append([s0,s1,s2])
			###############################
			elif s0 == "NA" and s2 != "NA":
				triplets.append([s1,s2])
			###############################
			elif s0 != "NA" and s2 == "NA":
				triplets.append([s0,s1])
	###############
	return triplets

def produceReadPairKeys(allChromosomeTriplets):
	'''This function takes in our list of triplets (adjacent ResctrictionScaffold Objects) for all chromosomes,
	   and generates a dictionary where the key is (scaffoldName,scaffoldName) for all the adjacent scaffolds found
	   in our triplets. This is for use in quantifying how many read pairs there are between each triplet.
	   Input = [ [ triplets, ... ], ... ] where each triplet is a list of RestrictionScaffoldObjects
	   Output = { (scaffoldName,scaffoldName):[] } where the key is scaffoldNames of adjacent scaffolds'''
	keyDict = {}
	for chromosomeTriplets in allChromosomeTriplets:
		for triplet in chromosomeTriplets:
			if len(triplet) == 3:
				keyDict.update({ (triplet[0].name,triplet[1].name) : [] })
				keyDict.update({ (triplet[1].name,triplet[0].name) : [] })
				keyDict.update({ (triplet[1].name,triplet[2].name) : [] })
				keyDict.update({ (triplet[2].name,triplet[1].name) : [] })
			else:
				keyDict.update({ (triplet[0].name,triplet[1].name) : [] })
				keyDict.update({ (triplet[1].name,triplet[0].name) : [] })
	##############
	return keyDict

def readValidPairFile(pairFile,pairDict):
	'''This function takes in the all valid pairs file output by HiCPro and queries our dictionary of adjacent scaffold names.
	   It then adds the readPair coordinates of the R1 and R2 reads to any key in the dictionary.
	   Input = full file path the all valid pairs file, { (scaffoldName,scaffoldName):[] } generated by the produceReadPairKeys function
	   Output = { (scaffoldName,scaffoldName):[ [R1Coordinate,R2Coordinate] , ... ] } where we add in the readPair coordinates for our adjacent scaffolds'''
	pFile = open(pairFile,'r')
	pairsExamined = 0
	for line in pFile:
		cols = line.strip('\r').strip('\n').split('\t')
		s1,s2 = cols[1],cols[4]
		if (s1,s2) in pairDict:
			pairDict[(s1,s2)].append([s1,s2,int(cols[2]),int(cols[5])])
		##################
		pairsExamined += 1
		if pairsExamined % 10000000 == 0:
			print("Read pairs looked at "+str(pairsExamined)+"...")
	#############
	pFile.close()
	return pairDict

def orientTrueTriplet(triplet,pairDict,lengthCutoff):
	'''This function takes in a list of three scaffolds (triplet) and will orient the middle scaffold with respect to the other two.
	   The middle scaffold should only consist of a single loci from hiCPro ( size <= resolution of contact map ).
	   We calculate the number of restriction sites of the other two scaffolds with respect to their orientation such that the distance
	   of a restriction site with respect to the middle scaffolds edge must be <= lengthCutoff.
	   Input = [ RestrictionScaffoldObjec0,RestrictionScaffoldObjec1,RestrictionScaffoldObjec2 ] must be 3 Restriction scaffold Objects
			   { (scaffoldName,scaffoldName): [ coordinates of readPairs , ... ] This allows us to query the number of links between any of the 3 scaffolds
			   lengthCutoff = distance we move into an adjacent scaffold before we no longer count read pairs or resctriction cut sites
	   Output = scaffoldName,+ or -, the name and orientation of the scaffold we just oriented'''
	for s in triplet:
		s.getResCounts(lengthCutoff)
	###########################################
	s0,s1,s2 = triplet[0],triplet[1],triplet[2]
	p,m = 0,0
	###########################################
	### Deal with s1 with respect to s2 (+) ###
	key = "NA"
	if len(pairDict[(s1.name,s2.name)]) != 0:
		key = (s1.name,s2.name)
		coordDict = { s1.name:2, s2.name:3 }
	elif len(pairDict[(s2.name,s1.name)]) != 0:
		key = (s2.name,s1.name)
		coordDict = { s2.name:2, s1.name:3 }
	if key != "NA":
		for rPair in pairDict[key]:
			s1Coord,s2Coord = rPair[coordDict[s1.name]],rPair[coordDict[s2.name]]
			if s2.orientation == "+":
				if s2Coord <= lengthCutoff:
					p += 1
			elif (s2.size)-s2Coord <= lengthCutoff:
				p += 1
		#########################
		if s2.orientation == "+":
			p = ( float(p) / float(s1.resRight+s2.resLeft) )
		else:
			p = ( float(p) / float(s1.resRight+s2.resRight) )
	###########################################
	### Deal with s1 with respect to s0 (-) ###
	key = "NA"
	if len(pairDict[(s1.name,s0.name)]) != 0:
		key = (s1.name,s0.name)
		coordDict = { s1.name:2, s0.name:3 }
	elif len(pairDict[(s0.name,s1.name)]) != 0:
		key = (s0.name,s1.name)
		coordDict = { s0.name:2, s1.name:3 }
	if key != "NA":
		for rPair in pairDict[key]:
			s1Coord,s0Coord = rPair[coordDict[s1.name]],rPair[coordDict[s0.name]]
			if s0.orientation == "-":
				if s0Coord <= lengthCutoff:
					m += 1
			elif (s0.size)-s0Coord <= lengthCutoff:
				m += 1
		#########################
		if s0.orientation == "-":
			m = ( float(m) / float(s1.resRight+s0.resLeft) )
		else:
			m = ( float(m) / float(s1.resRight+s0.resRight) )
	##########
	if p >= m:
		return s1.name,"+"
	else:
		return s1.name,"-"

def orientLeftEdgeCase(scaffLeft,scaffRight,pairDict,lengthCutoff):
	'''This function orients a scaffold positioned on the very left edge (with respect to the positive strand) of a chromosome, where the scaffold
	   is smaller than the contact map resolution. We cut the left most scaffold in half and calculate the number of restriction cut sites on each half,
	   as well as the number of links from the respective half there are to the immediate scaffold to the right within the give lengthCutoff distance.
	   Input = leftScaffoldObject,rightScaffoldObject,
			   { (scaffoldName,scaffoldName): [ coordinates of readPairs , ... ] This allows us to query the number of links between the two scaffolds,
			   lengthCutoff = distance we move into an adjacent scaffold before we no longer count read pairs or resctriction cut sites
	   Output = scaffoldName,+ or -, the name and orientation of the scaffold we just oriented'''
	scaffLeft.getResCounts(float(scaffLeft.size/2.))
	scaffRight.getResCounts(lengthCutoff)
	#########
	p,m = 0,0
	################################################################
	### Deal with scaffLeft with respect to scaffRight (+ and -) ###
	key = "NA"
	if len(pairDict[(scaffLeft.name,scaffRight.name)]) != 0:
		key = (scaffLeft.name,scaffRight.name)
		coordDict = { scaffLeft.name:2, scaffRight.name:3 }
	elif len(pairDict[(scaffRight.name,scaffLeft.name)]) != 0:
		key = (scaffRight.name,scaffLeft.name)
		coordDict = { scaffRight.name:2, scaffLeft.name:3 }
	if key != "NA":
		#################################
		if scaffRight.orientation == "+":
			minRight,maxRight = 0,lengthCutoff
		else:
			minRight,maxRight = scaffRight.size-lengthCutoff,scaffRight.size
		###########################    
		for rPair in pairDict[key]:
			scaffLeftCoord,scaffRightCoord = rPair[coordDict[scaffLeft.name]],rPair[coordDict[scaffRight.name]]
			if (scaffLeftCoord >= float(scaffLeft.size/2.)) and (minRight <= scaffRightCoord <= maxRight):
				p += 1
			elif (minRight <= scaffRightCoord <= maxRight):
				m += 1
	#################################
	if scaffRight.orientation == "+":
		p =  float(p) / ( float(scaffLeft.resRight + scaffRight.resLeft) )
		m = float(m) / ( float(scaffLeft.resLeft + scaffRight.resLeft) )
	else:
		p =  float(p) / ( float(scaffLeft.resRight + scaffRight.resRight) )
		m = float(m) / ( float(scaffLeft.resLeft + scaffRight.resRight) )
	##########
	if p >= m:
		return scaffLeft.name,"+"
	else:
		return scaffLeft.name,"-"

def orientRightEdgeCase(scaffLeft,scaffRight,pairDict,lengthCutoff):
	'''This function orients a scaffold positioned on the very right edge (with respect to the positive strand) of a chromosome, where the scaffold
	   is smaller than the contact map resolution. We cut the right most scaffold in half and calculate the number of restriction cut sites on each half,
	   as well as the number of links from the respective half there are to the immediate scaffold to the left within the give lengthCutoff distance.
	   Input = leftScaffoldObject,rightScaffoldObject,
			   { (scaffoldName,scaffoldName): [ coordinates of readPairs , ... ] This allows us to query the number of links between the two scaffolds,
			   lengthCutoff = distance we move into an adjacent scaffold before we no longer count read pairs or resctriction cut sites
	   Output = scaffoldName,+ or -, the name and orientation of the scaffold we just oriented'''
	scaffLeft.getResCounts(lengthCutoff)
	scaffRight.getResCounts(float(scaffRight.size/2.))
	#########
	p,m = 0,0
	################################################################
	### Deal with scaffRight with respect to scaffLeft (+ and -) ###
	key = "NA"
	if len(pairDict[(scaffLeft.name,scaffRight.name)]) != 0:
		key = (scaffLeft.name,scaffRight.name)
		coordDict = { scaffLeft.name:2, scaffRight.name:3 }
	elif len(pairDict[(scaffRight.name,scaffLeft.name)]) != 0:
		key = (scaffRight.name,scaffLeft.name)
		coordDict = { scaffRight.name:2, scaffLeft.name:3 }
	if key != "NA":
		#################################
		if scaffLeft.orientation == "+":
			minLeft,maxLeft = scaffLeft.size-lengthCutoff,scaffLeft.size
		else:
			minLeft,maxLeft = 0,lengthCutoff
		###########################    
		for rPair in pairDict[key]:
			scaffLeftCoord,scaffRightCoord = rPair[coordDict[scaffLeft.name]],rPair[coordDict[scaffRight.name]]
			if (scaffRightCoord < float(scaffRight.size/2.)) and (minLeft <= scaffLeftCoord <= maxLeft):
				p += 1
			elif (minLeft <= scaffLeftCoord <= maxLeft):
				m += 1       
	#################################
	if scaffLeft.orientation == "+":
		p =  float(p) / ( float(scaffLeft.resRight + scaffRight.resLeft) )
		m = float(m) / ( float(scaffLeft.resRight + scaffRight.resRight) )
	else:
		p =  float(p) / ( float(scaffLeft.resLeft + scaffRight.resLeft) )
		m = float(m) / ( float(scaffLeft.resLeft + scaffRight.resRight) )
	##########
	if p >= m:
		return scaffRight.name,"+"
	else:
		return scaffRight.name,"-"

def orientTriplet(triplet,scaffList,pairDict,lengthCutoff):
	'''In this function we orient the scaffold in question with respect to the adjacent scaffolds.
	   There are 3 scenarios that can occur so there are 3 different functions (1 for each).
	   The first is that the scaffold to orient has scaffolds that flank it to the right and left (with respect to the + strand).
	   The second is that the scaffold to orient is on the very left edge of the chromosome.
	   The third is that the scaffold to orient is on the very right edge of the chromosome.
	   In the first case we simple normalize the links to the left and right scaffolds by the number of restriction sites
	   contained on the middle scaffold + the left restriction sites or + the right restriction sites. The number of links is calculated by,
	   the number of links between the middle scaffold and left right respectively that are within a distance of lengthCutoff from the middle
	   scaffolds respective edge.
	   In the other two scenarios, we calculate the number of cut sites on the left and right half of the scaffold to orient and we calculate
	   the number of links with respect to each half, and again normalize by the number of cutSites on the other scaffold + the number of cutSites
	   on the left or right half of the scaffold to orient.
	   Input = triplet which is a list of 2 or 3 adjacent scaffolds. Only one will be oriented.
	   Output = scaffoldName,+ or -, the name and orientation of the scaffold we just oriented'''
	### Deal with cases where the scaffold of interest to be oriented is surrounded by one or more scaffold on either side ###
	if len(triplet) == 3:
		scaffName,scaffOrientation = orientTrueTriplet(triplet,pairDict,lengthCutoff)
	else:
		s0,s1 = triplet[0],triplet[1]
		########################################################################################################################
		### Deal with the edge case where the scaffold of interest to be oriented is at the very beggining of the chromosome ###
		if s0.name == scaffList[0].name:
			scaffName,scaffOrientation = orientLeftEdgeCase(s0,s1,pairDict,lengthCutoff)
		##################################################################################################################
		### Deal with the edge case where the scaffold of interest to be oriented is at the very end of the chromosome ###
		else:
			scaffName,scaffOrientation = orientRightEdgeCase(s0,s1,pairDict,lengthCutoff)
	#################################
	return scaffName,scaffOrientation

def giveFinalChromOrdering(trips,scaffGroups,scaffDict,validPairs,resolution,lengthCutoff=500000):
	'''This function orients scaffolds that are smaller than the contact resolution size.
	   Input = list of RestrictionScaffolds that are adjacent to each other [ [adjacentScaffolds , ... ] ],
			   list of RestrictionScaffold objects belonging to a given chromosome genome wide,
			   { scaffoldName:RestrictionScaffoldObject } a dictionary that points to the scaffGroups objects,
			   { (scaffoldName,scaffoldName): [[coordinates of read pairs] , ... ] } a dictionary of adjacent scaffold pairs and their read pair coordinates,
			   resolution = the resolution of the contact map
			   lengthCutoff = the distance into an adjacent scaffold that we consider read pairs and restriction cut sites to determine orientation
	   Output = [ [ [scaffoldName,orientation] , ... ] ] lists of scaffolds within each chromosome group describing the order and orientations'''
	if lengthCutoff < resolution:
		print("lengthCutoff variable is set too low... Setting equal to resolution variable")
		lengthCutoff = resolution
	##########################
	chromosomeGroupOrders = []
	for chromosomeTriplets,chromosomeScaffs in zip(trips,scaffGroups):
		if len(chromosomeTriplets) != 0:
			chromGroup = []
			for trip in chromosomeTriplets:
				scaffName,scaffOrientation = orientTriplet(trip,chromosomeScaffs,validPairs,lengthCutoff=lengthCutoff)
				scaffDict[scaffName].orientation = scaffOrientation
		##################################################################################
		chromosomeGroupOrders.append([ [s.name,s.orientation] for s in chromosomeScaffs ])
	############################
	return chromosomeGroupOrders

def writeScaffoldOrderingsToFile(sOrderings,outFile):
	'''Simple function that writes out our final ordering/orientations of all scaffolds to the final output file!
	   Input = output from the giveFinalChromOrdering function, full file path to the outfile that we want to write'''
	outFile = open(outFile,'w')
	chromCount,scaffsWritten = 0,0
	for scaffGroup in sOrderings:
		chromCount += 1
		outFile.write("### Chromosome grouping "+str(chromCount)+" ###"+'\n')
		for s in scaffGroup:
			if (type(s) == type(RestrictionScaffold('','','','',''))):
				outFile.write(s.name+'\t'+s.orientation+'\n')

			elif type(s) == type([]):
				outFile.write(s[0]+'\t'+s[1]+'\n')
			else:
				print("- WARNING invalid output type {}... Expecting list or Scaffold like class...".format(type(s)))
			scaffsWritten += 1
	###############
	outFile.close()
	print("Chromosome groups written to file "+str(chromCount))
	print("Scaffolds written to file "+str(scaffsWritten))


########################################################
### Function that runs the entire pipeline for Part3 ###
def runPipeline(chromosomeOrderFile,scaffSizeFile,restrictionSiteFile,validPairFile,finalOrderingFile,lengthCutoff,resolution):
	#################################################
	print("########################################")
	print("### Working on Part3 of the pipeline ###")
	startTime = time.time()
	scaffGroups,scaffDict = initiateScaffoldObjects(chromosomeOrderFile,scaffSizeFile,restrictionSiteFile,resolution)
	trips = [ pullTriplets(sGroup) for sGroup in scaffGroups ]
	validPairs = produceReadPairKeys(trips)
	validPairs = readValidPairFile(validPairFile,validPairs)
	finalChromGroups = giveFinalChromOrdering(trips,scaffGroups,scaffDict,validPairs,resolution=resolution,lengthCutoff=lengthCutoff)
	writeScaffoldOrderingsToFile(finalChromGroups,finalOrderingFile)
	print("Total run-time  for Part3 = "+str(time.time()-startTime))
	print("- Part 3 (optional orientation of scaffolds smaller than resulution size) completed successfully")
