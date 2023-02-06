###############
### Imports ###
import math
import numpy
import scipy
import time
import numpy
from numpy import trace as numpyTrace
import math
import copy
import numba
from numba import jit
import plotContactMaps as plotModule
####################################


################################################################################################
### Contact map initiation, and matrix manipulation code (the exact same code used in Part1) ### 
class Bin:
	'''This object class represents a single genomic loci within the genome.
	   It contains the base information from the output of HiC-Pro'''
	def __init__(self,ID,chrom,start,stop,bias,rowSum):
		self.ID = ID ### ID number give by HiC-Pro
		self.chrom = chrom ### Chromosome and or Scaffold the loci belongs to
		self.start = start ### Start coodinate (5'-3' with respect to positive strand)
		self.stop = stop ### Stop coordinate (5'-3' with respect to positive strand)
		self.bias = bias ### Bias value assigned by HiC-Pro
		self.rowSum = rowSum ### The total sum of interaction values (will be the same for all valid bin objects and zero for non valid objects due to the ICE Normalization step performed by HiCPro)

def initiateLoci(bedFile,biasFile,binID_dict=False):
	''' This initiates our inititial objects (Bin) of genomic loci.
		Input = bedFile is the full file path to the bed file output by HiCPro
				biasFile is the full file path to the biasFile output by HiCPro
		Output = list of Bin objects'''
	binList = []
	bedFile = open(bedFile,'r')
	biasFile = open(biasFile,'r')
	while True:
		bedLine = bedFile.readline()
		biasLine = biasFile.readline()
		if not bedLine:
			break
		#####################################################
		bedCols = bedLine.strip('\r').strip('\n').split('\t')
		chrom,start,stop,bID = bedCols[0],int(bedCols[1]),int(bedCols[2]),int(bedCols[3])
		bias = biasLine.strip('\r').strip('\n')
		#######################
		if binID_dict != False:
			if bID not in binID_dict:
				continue
		#################
		if bias != "nan":
			try:
				biasValue = float(bias)
			except:
				biasValue = 0.
			binMan = Bin(bID,chrom,start,stop,biasValue,0.)
			binList.append(binMan)
	###############
	bedFile.close()
	biasFile.close()
	print("Genomic loci found"+'\t'+str(len(binList)))
	return binList

def buildAdjacencyMatrix(matrixFile,binList,binID_dict=False):
	'''This builds our initial adjacency matrix from the interaction values provided by HiCPro
	   Input = matrixFile is the ICE normalized matrix file produced by HiCPro
			   binList is the list of Bin objects provided by the initiateLoci function
	   Output = adjacency matrix of interaction values'''
	matrixFile = open(matrixFile,'r')
	adjacencyMatrix = [ [ 0. for i in range(0,len(binList)) ] for i in range(0,len(binList)) ]
	binIndexDict = { binMan.ID:i for i,binMan in enumerate(binList) }
	edgeCount = 0
	print("Adjacency matrix initiated")
	print("Adding edges...")
	for line in matrixFile:
		cols = line.strip('\r').strip('\n').split('\t')
		bID1,bID2,val = int(cols[0]),int(cols[1]),float(cols[2])
		if bID1 not in binIndexDict:
			continue
		if bID2 not in binIndexDict:
			continue
		adjacencyMatrix[ binIndexDict[bID1] ][ binIndexDict[bID2] ] = val
		adjacencyMatrix[ binIndexDict[bID2] ][ binIndexDict[bID1] ] = val
		edgeCount += 1
		if edgeCount % 1000000 == 0:
			print("Edges processed "+str(edgeCount))
	##################
	matrixFile.close()
	############################################################
	print("Edges added to adjacency matrix"+'\t'+str(edgeCount))
	print("Rows in adjacency matrix "+str(len(adjacencyMatrix)))
	return adjacencyMatrix

def removeRows(matrix,binList,zeroRows=True,biasVals=False):
	''' Function to remove rows from the matrix 
		and remove the corresponding binObject from the binList. 
		Input is the matrix, and corresponding list of bin objects.
		If zeroRows==True then we remove rows that sum to zero
		If biasVals!=False then we only keep rows that satisify
		   biaseVals[0] < rowBiasValue < biasVals[1]
		Input = adjacencyMatrix produced by the buildAdjacencyMatrix function
				list of Bin objects produced by the initiateLoci function
		Output = adjacencyMatrix and binList, both with removed rows/elements that sum to zero or have invalid biasValue'''
	indicesToRemove = []
	###################################################################
	for i,rSum in enumerate(numpy.sum(numpy.asmatrix(matrix),axis=-1)):
		if zeroRows == True:
			if numpy.asarray(rSum) == 0:
				indicesToRemove.append(i)
				continue
		#######################
		if (biasVals != False):
			if (binList[i].bias > biasVals[1]) or (binList[i].bias < biasVals[0]):
				indicesToRemove.append(i)
	##########################################################
	print("Rows/columns to remove "+str(len(indicesToRemove)))
	count = 0
	for i in indicesToRemove:
		del matrix[i-count]
		del binList[i-count]
		count += 1
	#############################
	for i,r in enumerate(matrix):
		count = 0
		for ii in indicesToRemove:
			del r[ii-count]
			count += 1
		binList[i].rowSum = sum(numpy.asarray(r))
	#####################
	return matrix,binList

def convertMatrix(adjacencyMatrix,binList,distance=True,similarity=False):
	'''This function is a general function to convert back and forth from a distance to a similarity matrix
	   Input = adjacencyMatrix and corresponding Bin object list (we use the preComputed rowSum attribute of the binObjects to speed up the process )
	   Output = transformed adjacencyMatrix'''
	startTime = time.time()
	adjacencyMatrix = scipy.asmatrix(adjacencyMatrix)
	########################################
	for i,row in enumerate(adjacencyMatrix):
		if distance == True:
			adjacencyMatrix[i] = (1. - (row/row.sum())) + 1.
		elif similarity == True:
			adjacencyMatrix[i] = binList[i].rowSum*(1.-(row-1.))
	#################################################################################
	#condensedMatrix = scipy.spatial.distance.squareform(adjacencyMatrix,checks=False)
	######################
	stopTime = time.time()
	print("Time to convert matrix "+str(stopTime-startTime))
	return adjacencyMatrix

def reorderMatrix(matrix,binList,newOrder):
	'''Fancy little function to reorder the rows and columns of our matrix according to new indices (also the list of Bin objects)
	   Input = adjacencyMatrix, list of Bin objects, indices of newOrder
	   Output = reordered matrix and list of Bin objects'''
	matrix = matrix[:, newOrder][newOrder]
	binList = [ binList[i] for i in newOrder ]
	return matrix,binList

def logTransformMatrix(matrix,logBase=10,reverse=False):
	'''Simple function to convert a matrix to and from a logTransform
	   Input = matrix and logBase or reverse=True/False
	   Output = transformed matrix'''
	matrix = scipy.asmatrix(matrix)
	newMat = [ [ 0.0 for ii in range(0,len(matrix)) ] for i in range(0,len(matrix)) ]
	if reverse == False:
		for i,row in enumerate(matrix):
			for ii,v in enumerate(numpy.asarray(row)[0]):
				if v != 0.:
					newMat[i][ii] = math.log(float(v),logBase)
	#####
	else:
		for i,row in enumerate(matrix):
			for ii,v in enumerate(numpy.asarray(row)[0]):
				if v != 0.:
					newMat[i][ii] = (float(logBase)**v)
	#############################
	return numpy.asmatrix(newMat)


#############################################################################
### This is the cost/score function that we need to precompile with numba ###
### We use numba to speed it up significantly ###############################
@jit(nopython=True)
def costFunction_numba(matrix,total):
	cumulitiveSum,cost = 0.,0.
	for i in range(1,len(matrix)):
		cumulitiveSum += numpyTrace(matrix,offset=i)
		cost += (cumulitiveSum/total/float(i))
	###########
	return cost

costFunctionPreCompile = costFunction_numba(numpy.asmatrix([ [ float(i) for i in range(0,10) ] for ii in range(0,10) ]),10.)
### Matrix entries MUST be floating points
##########################################


##############################################################
### Code to order and orient scaffolds within a chromosome ###
def readGroupingsToValidBins(chromosomeGroupFile):
	'''This function simply reads in the file that contains all the valid bins from the chromosome assignment step.
	   It reads them into a dictionary so that when we initiate our adjacency matrix and list of binObjects, we are only including the valid ones.
	   Input = full file path to the chromosome group file
	   Output = { binID:'' }'''
	binID_dict = {}
	inFile = open(chromosomeGroupFile,'r')
	for line in inFile:
		line = line.strip('\r').strip('\n')
		if line[0] != "#":
			bID = int(line.split('\t')[0])
			binID_dict.update({ bID:'' })
	##############
	inFile.close()
	return binID_dict

def readChromsFromFile(inFile):
	'''This function takes in the final output from the scaffold to chromosome assignment part of the algorithm.
	   It reads the input file into a list data structure, where each element in the list is another list whose elements are composed of the lines corresponding to the chromosome.
	   Input = full file path to the chromosome groups file (this is the file that should have the properly named chromosomes along with wich binIDs correspond to the chromosome)
	   Output = [ [ [binID,scaffoldName] ,... ] ,... ] where each [ [binID,scaffoldName] ,... ] belongs to any given chromosome'''
	cList,chrom = [],[]
	inFile = open(inFile,'r')
	inFile.readline()
	for line in inFile:
		line = line.strip('\r').strip('\n')
		if line[0] != "#":
			binID = int(line.split('\t')[0])
			scaff = line.split('\t')[1]
			chrom.append([ binID,scaff ])
		else:
			cList.append(chrom)
			chrom = []
	###################
	cList.append(chrom)
	print("Chromosomes found "+str(len(cList)))
	print("Nodes found "+str(sum([ len(c) for c in cList ])))
	return cList

class Scaffold:
	'''This is one of the core data structures we use in the ordering and orienting of scaffolds within a chromosome.
	   It contains all the useful information about the given scaffold. By default, the orientation is "+"'''
	def __init__(self,name,binList,orientation):
		self.name = name
		self.binList = binList
		self.orientation = orientation
	def flipOrientation(self):
		'''Flips the orientation from "+" to "-" and vice versa.
		   It also switches the order of the bins contained within the scaffold accordingly'''
		if self.orientation == "+":
			self.orientation = "-"
		else:
			self.orientation = "+"
		#################################
		self.binList = self.binList[::-1]

def initiateBinsAndScaffolds(nodeList):
	'''This function takes in a list of the nodes/bins from any give chromosome group as input.'''
	'''Each entry in the nodeList is [ binID, scaff ]. It then creates a { scaffoldName:ScaffoldObj } for each scaffold found within the nodeList.
	   It also creates a size sorted list of scaffoldObjects. These objects point to the scaffold objects in the scaffoldDictionary.
	   Each scaffold object contains a list of the nodes/bins that is sorted from 5'-->3' end
	   Input = [ [binID,scaffoldName],...] The bins corresponding to a chromosome group
	   Output = [ ScaffoldObject,... ] , { scaffoldName:ScaffoldObject } The objects point to each other'''
	scaffDict = { n[1]:'' for n in nodeList }
	print("Scaffolds to order for this chromosome "+str(len(scaffDict)))
	scaffDict = { s:Scaffold(s,[],"+") for s in scaffDict }
	##################
	for n in nodeList:
		scaffDict[n[1]].binList.append(n[0])
	#####################################################################################
	### The node/binIDs are ordered from smallest to biggest --> 5'--->3' end of scaffold
	for s in scaffDict:
		scaffDict[s].binList = sorted(scaffDict[s].binList,key=lambda binID: binID)
	##############
	scaffList = []
	for sObj in scaffDict.values():
		sObj.nodeCount = len(sObj.binList)
		scaffList.append(sObj)
	#############################################################################
	scaffList = sorted(scaffList,key=lambda sObj: len(sObj.binList),reverse=True)
	return scaffList,scaffDict

def pullScaffolds(puller,pullee,scaffsToPull):
	'''Pulls scaffold objects from one list (pullee) and places them in another (puller)
	   Input = number of scaffolds to pull. Will pull scaffolds from the very beginning of the list first 
			   list of scaffold objects to pull from
			   list of scaffoldNames to pull
	   Output = the two input lists with the corresponding scaffolds moved to and frow'''
	for i in range(0,scaffsToPull):
		if len(pullee) == 0:
			break
		sObj = pullee.pop(0)
		puller.append(sObj)
	####################
	return puller,pullee

def giveNewAdjMat(matrix,scaffList,binList):
	'''This function takes in full genome adjacency matrix and generates a subgraph/sub-adjacencyMatrix from the selected scaffolds.
	   It keeps track of the indices that each bin has within the current order of the chromosome the matrix when necessary.
	   Input = full genome adjacency matrix, list of scaffold objects, list of Bin objects
	   Output = scaffold specific adjacency matrix, { binID:index in chromosome specific list of bins }'''
	nList = [ n for sObj in scaffList for n in sObj.binList ]   
	orderDict = { gID:i for i,gID in enumerate(nList) }
	binIndices = { b.ID:i for i,b in enumerate(binList) }
	indices = [ binIndices[n] for n in nList ]
	####################################################################
	adjacencyMatrix = numpy.asmatrix(matrix[numpy.ix_(indices,indices)])
	################################
	return adjacencyMatrix,orderDict

def reorderScaffList(orderList,orientationList,scaffDict):
	'''Reorders the list of scaffolds to a new order / orientation
	   Input = new order of the scaffolds, new order of the orientation of each scaffold, { scaffoldDict:scaffoldObj }
	   Output = newelyOrdered list of scaffolds, newelyOrdered list of bins'''
	newScaffList,newNodeList = [],[]
	for scaffoldName,orientation in zip(orderList,orientationList):
		if scaffDict[scaffoldName].orientation != orientation:
			scaffDict[scaffoldName].flipOrientation()
		newScaffList.append(scaffDict[scaffoldName])
		newNodeList += scaffDict[scaffoldName].binList
	###############################
	return newScaffList,newNodeList

def costFunction(matrix,total):
	''' The none numba version (slower) but I am leaving it in for now'''
	cumulitiveSum,cost = 0.,0.
	for i in range(1,len(matrix)):
		cumulitiveSum += numpyTrace(matrix,offset=i)
		cost += (cumulitiveSum/total/float(i))
	###########
	return cost

def checkAllScores(adjMat,orderDict,orderedScaffs,scaffToCheck):
	'''This function takes in an input of an adjacency matrix where one scaffold is unordered
	   and tests all possible scores of where that scaffold could go. For N ordered scaffolds there are 2(N+1) possible scores to check.
	   Input = adjacencyMatrix where one scaffolds ordering is unkown, { binID:currentIndex in the scaffold }
			   [ list of ordered scaffolds ], scaffoldName to check all possible scores for.
	   Output = [ list of newely ordered scaffolds ], best score that we found after checking all possible combos '''
	bestCost = 0.
	bestOrd = "NA"
	bestInsertionPoint = 0
	bestOrientation = "+"
	length = len(orderedScaffs)+1
	total = sum([ numpy.trace(adjMat,offset=i) for i in range(1,len(adjMat)) ]) ### Calculate the total weight of the upper triangle of the adjacencyMatrix
	for i in range(0,length):
		### Insert and check original orientation of scaff to check ###
		orderedScaffs.insert(i,scaffToCheck)
		newNodeOrder = [ orderDict[n] for sObj in orderedScaffs for n in sObj.binList ]
		cost = costFunction_numba(adjMat[numpy.ix_(newNodeOrder,newNodeOrder)],total)
		if cost > bestCost:
			bestCost = cost
			bestOrd = [ [sObj.name,sObj.orientation] for sObj in orderedScaffs ]
			bestInsertionPoint = i
			bestOrientation = orderedScaffs[i].orientation
		####################################################
		### Now we want to test the opposite orientation ###
		orderedScaffs[i].flipOrientation()
		newNodeOrder = [ orderDict[n] for sObj in orderedScaffs for n in sObj.binList ]
		cost = costFunction_numba(adjMat[numpy.ix_(newNodeOrder,newNodeOrder)],total)
		if cost > bestCost:
			bestCost = cost
			bestOrd = [ [sObj.name,sObj.orientation] for sObj in orderedScaffs ]
			bestInsertionPoint = i
			bestOrientation = orderedScaffs[i].orientation
		###################################
		scaffToCheck = orderedScaffs.pop(i)
	###############################################
	if scaffToCheck.orientation != bestOrientation:
		scaffToCheck.flipOrientation()
	#####################################################
	orderedScaffs.insert(bestInsertionPoint,scaffToCheck)
	#############################
	return orderedScaffs,bestCost

def calcPossiblePerms(N):
	'''Calculates the number of permutations of order and orientation of N number of scaffolds.
	   Input = number of scaffolds
	   Output = number of possible permutations'''
	p = math.factorial(N) * (2**N) / 2
	return p

def permutations(elementList,paths,k=0):
	'''Fancy little algorithm I stole from a python algorithms book at some point, and modified slightly for my purposes.
	   It returns all permutations of a list. Will contain reverse duplicates, which need to be removed for our purposes.
	   Input = list of elements to permute, permutations that we have found so far, a number to progress the recursion
	   Ouput = at the final step of recursion we will have all permutaitons of the element list'''
	if k == len(elementList):
		paths.append([ copy.copy(e) for e in elementList ])
	else:
		for i in range(k, len(elementList)):
			elementList[k],elementList[i] = elementList[i],elementList[k]
			permutations(elementList,paths,k+1)
			elementList[k],elementList[i] = elementList[i],elementList[k]
		############
		return paths

def removeReverseDuplicates(permList):
	'''Removes reverse duplicates from the permutation list. Cuts our compute time in half!
	   Input = permutations containing all
	   Output = permutations containing only the foreward ones'''
	forewardDict = {}
	counter = 0
	for i in range(0,len(permList)):
		revMan = tuple(permList[i-counter][::-1])
		if revMan not in forewardDict:
			forewardDict.update({tuple(permList[i-counter]):''})
		else:
			del forewardDict[revMan]
			permList.pop(i-counter)
			counter += 1
	###############
	return permList

def plusMinusPerms(elementList):
	'''Finds all possible permutations of the "+" and "-" orientation for the scaffolds / elements in the list of things to permute.
	   Input = number of elements to permute
	   Output = list of all possible orientations for all the elements '''
	plusMinusPerms = [ ["+"] * len(elementList) ]
	for i in range(0,len(elementList)):
		permList = ["+"] * i
		permList += ["-"] * (len(elementList) - i)
		perms = permutations(permList,[],k=0)
		for p in perms:
			plusMinusPerms.append(p)
	#############
	seenDict = {}
	for p in plusMinusPerms:
		if tuple(p) not in seenDict:
			seenDict.update({ tuple(p):'' })
	########################################
	return [ list(key) for key in seenDict ]

def bruteForceBestScore(sObjList,scaffDict,matrix,orderDict):
	'''This function generates best ordering/orientation of the scaffold objects passed in as well as the adjacency matrix passed in.
	   We first find all possible permutations of scaffold orderings/orientations and calculate the scores for all of them.
	   Then we need to calculate the total weight of the upper triangle of the adjMatrix so we don't have to do it over and over again...
	   Input = list of scaffold objects to find the best ordering/orientation for
			   { scaffoldName:ScaffoldObject } this points to the objects in the list that is passed in
			   adjacencyMatrix of the scaffolds in the list
			   { binID:currentIndex in the scaffold }
	   Output = [ list of scaffoldNames in the best order found ], [ list of orientations of each scaffold in the best order found ]
				best score from cost function'''
	elementsToPermute = [ sObj.name for sObj in sObjList ] ### Generate the list of scaffold names to find permutations for
	differentOrderings = permutations(elementsToPermute,[],k=0) ### Find all permutations of scaffold orderings 
	differentOrderings = removeReverseDuplicates(differentOrderings) ### Remove reverse duplicates of the scaffold orderings (so we don't double our work)
	orientations = plusMinusPerms(elementsToPermute) ### Add in all the possible orientations of the scaffolds to mix
	############################################################
	count,bestCost,bestOrder,bestOrientation = 0, 0., "NA", "NA"
	total = sum([ numpy.trace(matrix,offset=i) for i in range(1,len(matrix)) ]) ### Calculate the total weight of the upper triangle of the adjacencyMatrix
	print("Initial permutations to test "+str(len(differentOrderings)*len(orientations))+"...")
	startTime = time.time()
	for dOrder in differentOrderings:
		for dOrient in orientations:
			nScaffs,nNodes = reorderScaffList(dOrder,dOrient,scaffDict)
			nOrder = [ orderDict[n] for n in nNodes ]
			#cost = costFunction(matrix[:, nOrder][nOrder],total)
			#cost = costFunction_numba(matrix[:, nOrder][nOrder],total)
			cost = costFunction_numba(matrix[numpy.ix_(nOrder,nOrder)],total) ### This speeds up computation dramatically
			if cost > bestCost:
				bestOrder,bestOrientation = dOrder,dOrient
				bestCost = cost
			count += 1
			if count % 10000 == 0:
				stopTime = time.time()
				print("Time for "+str(count)+" calculations = "+str(stopTime-startTime))
				startTime = time.time()
	#########################################
	return bestOrder,bestOrientation,bestCost

def orderRemainderScaffolds(orderedScaffolds,scaffoldList,orderDict,matrix,binList):
	'''This functions takes in the scaffolds that have already been ordered/oriented and iteratively adds in new scaffolds.
	   For each unordered scaffold we test the score at each possible position and orientation.
	   For N already ordered scaffolds, there are 2(N+1) possible scores to check. 
	   Input = [ list of scaffold names that have already been properly ordered ], [ list of scaffoldNames left to order ],
			   { binID:currentIndex in the scaffold }, list of Bin objects
	   Output = [ list of ordered scaffolds (all of them) ], best cost found for full adjacency matrix'''
	scaffsAdded = 0
	while True:
		orderedScaffolds,scaffoldList = pullScaffolds(orderedScaffolds,scaffoldList,1) ### Pull one scaffold from the unordered list and add it to the end of the ordered list
		adjMat,orderDict = giveNewAdjMat(matrix,orderedScaffolds,binList) ### Generate new adjacency matrix with new scaffold values included
		newScaff = orderedScaffolds.pop(-1) ### Remove the newely added scaffold from the ordered list
		orderedScaffolds,bestCost = checkAllScores(adjMat,orderDict,orderedScaffolds,newScaff) ### Check all the possible scores from the input ajacency matrix and add scaffold
		scaffsAdded += 1
		##########################
		if len(scaffoldList) == 0:
			break
	################################
	return orderedScaffolds,bestCost

def scanOrdering(orderedScaffolds,scaffoldDict,orderDict,matrix,binList,bestCost,scanScaffolds=5):
	'''This function takes the input matrix of the entire chromosome in its best found order, and slides along the scaffold list in a sliding window.
	   For each step in the sliding window, all permutations of scaffolds within the window are calculated with the entire chromosome being taken
	   into account for each permutation. The best score / ordering is taken at each step of the sliding window, and the window slides with the newest order found.
	   All scores / permutations are calcutated again until the window reaches the end of a chromosome. If the best order found is not the same as when the algorithm first started,
	   the algorithm repeats until it converges on a single best ordering. In theory, at the beggining of the algorithm, a scaffold found within the first window could migrate
	   all the way to the end of the chromosome if at each step a new better ordering is found (unlikely, but possible).
	   Input = [ list of ordered scaffold objects ], { scaffold:scaffoldObject }, { binID:current index of matrix / binList }, 
			   full adjacency matrix, [ list of Bin objects], best score from starting best order, size of the sliding window (in number of scaffolds)
	   Output = [ final ordered scaffoldObjects ], final best score from the algorithm'''
	adjMat,orderDict = giveNewAdjMat(matrix,orderedScaffolds,binList) ### Generate chromosome specific adjacency matrix
	total = sum([ numpy.trace(adjMat,offset=i) for i in range(1,len(adjMat)) ]) ### Calculate the total weight of the upper triangle of the adjacencyMatrix
	bestOrder = [ sObj.name for sObj in orderedScaffolds ]
	bestOrientation = [ sObj.orientation for sObj in orderedScaffolds ]
	roundNumber = 0
	while True:
		stop = 0
		print("Working on round "+str(roundNumber+1)+" of final step...")
		for i in range(0,len(orderedScaffolds)-scanScaffolds+1):
			elementsToPermute = [ sObj.name for sObj in orderedScaffolds[i:i+scanScaffolds] ]
			differentOrderings = permutations(elementsToPermute,[],k=0)
			differentOrderings = removeReverseDuplicates(differentOrderings)
			orientations = plusMinusPerms(elementsToPermute)
			count = 0
			for dOrder in differentOrderings:
				for dOrient in orientations:
					### Getting the new order of scaffolds ###
					begginingOrder = [ sObj.name for sObj in orderedScaffolds[0:i] ] ### Beggining
					windowedScaffoldsOrder,windowedNodeList = reorderScaffList(dOrder,dOrient,scaffoldDict) ### Middle
					endOrder = [ sObj.name for sObj in orderedScaffolds[i+scanScaffolds:] ] ### End
					newOrder = begginingOrder+[ sObj.name for sObj in windowedScaffoldsOrder ]+endOrder
					################################################
					### Getting the new orientation of scaffolds ###
					newNodes,newOrient = [],[]
					for scaffName in newOrder:
						newNodes += scaffoldDict[scaffName].binList
						newOrient.append(scaffoldDict[scaffName].orientation)
					nOrder = [ orderDict[n] for n in newNodes ]
					#################################################################
					cost = costFunction_numba(adjMat[numpy.ix_(nOrder,nOrder)],total)
					if cost > bestCost:
						bestOrder,bestOrientation = newOrder,newOrient
						bestCost = cost
						stop = 1
			##################################################################################
			orderedScaffolds,nNodes = reorderScaffList(bestOrder,bestOrientation,scaffoldDict)
			adjMat,orderDict = giveNewAdjMat(matrix,orderedScaffolds,binList)
		################
		roundNumber += 1
		if stop == 0:
			break
	####################################################################
	print("Sliding window conversion after "+str(roundNumber)+" rounds")
	print("Best cost at the end of the final step = "+str(bestCost))
	return orderedScaffolds,bestCost

def orderChromosome(chromGroup,matrix,binList,nScaffolds=6,scanScaffolds=5):
	'''This function ordereds/orients the scaffolds contained within a given chromosome. In essence, we define a scoring/cost function and order/orient
	   the top largest nScaffolds by testing all permutations. Then we iteratively add in the smaller scaffolds one by one by testing all positions that they could go.
	   For N already ordered scaffolds, there are 2(N+1) number of possible scores to check for each subsequent scaffold we add in. Then in a final step after all
	   scaffolds have been ordered, we use a sliding window to test all permutations of scanScaffolds number of scaffolds within the window (but using the full adjacency matrix for the scaffold).
	   This in theory helps us overcome local optima, but using the whole matrix, but only permuting a select portion of it iteratively until we reach the end, taking the best score / permutation
	   at each step. Until we converge on a single score (the global optima approximation for our case (or the best we can do given the heuristics))
	   Input = [ [binID,scaffoldName],... ] For all bins contained within the chromosome,
			   the full adjacency matrix for the data (full genome), [ list of binObjects ],
			   number of scaffolds to initially use for bruteForce computations, number of scaffolds to use in the sliding window of the final stage of the algorithm.
	   Output = [ list of ordered scaffoldObjects ]'''
	if nScaffolds >= 9:
		print("Number of initial scaffolds to order by brute force method is set too high...")
		print(str(calcPossiblePerms(nScaffolds))+" Different permutations would need to be calculated with current setting")
		print("Setting number of initial scaffolds to 8")
		nScaffolds = 8
		print(str(calcPossiblePerms(nScaffolds)) +" Different permutations will be calculated initially. This could take a while still...")
	if scanScaffolds > nScaffolds:
		scanScaffolds = nScaffolds
	############################################################################################
	### Steps of algorithm #####################################################################
	scaffoldList,scaffoldDict = initiateBinsAndScaffolds(chromGroup) ### Initiate dataStructures
	orderedScaffolds,scaffoldList = pullScaffolds([],scaffoldList,nScaffolds) ### Pull n largest scaffolds into new list to be ordered/oriented
	adjMat,orderDict = giveNewAdjMat(matrix,orderedScaffolds,binList) ### Generate adjacency matrix containing scaffolds that need to be ordered
	bfBestOrder,bfBestOrientation,bfBestScore = bruteForceBestScore(orderedScaffolds,scaffoldDict,adjMat,orderDict) ### Test all permutations of scaffold orders/orientations
	### Special ammendment for chr1 of gularis ###

	chr1Dict = { s:'' for s in bfBestOrder }
	if "ScHqM4t_47_G" in chr1Dict:
		bfBestOrder = ['ScHqM4t_47_G','ScHqM4t_65_G', 'ScHqM4t_29_G', 'ScHqM4t_29.1_G', 'ScHqM4t_25.1_G', 'ScHqM4t_28_G']
		bfBestOrientation = ['-','+', '-', '-', '-', '-']
		#################################################

	orderedScaffolds,orderedNodes = reorderScaffList(bfBestOrder,bfBestOrientation,scaffoldDict) ### Reorder/orient scaffold list according to best score from the bruteforce method
	orderedScaffolds,bestCost = orderRemainderScaffolds(orderedScaffolds,scaffoldList,orderDict,matrix,binList) ### Pull remaining unordered scaffolds into matrix and find their best placement one at a time until all are used
	print("BestCost at the end of first two steps "+str(bestCost))
	if len(orderedScaffolds) > nScaffolds:
		orderedScaffolds,bestCost = scanOrdering(orderedScaffolds,scaffoldDict,orderDict,matrix,binList,bestCost,scanScaffolds=scanScaffolds)
	##########################
	print("Final ordering...")
	for s in orderedScaffolds: ### Print out results
		print(s.name,s.orientation)
	#######################
	return orderedScaffolds


###############################################################
### Code to wrap everything together across all chromosomes ###
def orderGenome(matrix,chromList,binList,resolution,nScaffolds=6,scanScaffolds=5,plotChrom=True,showPlot=True,savePlotDir=False,plotTitleSuffix=False):
	'''This function will order and orient all scaffolds with respect to their chromosome for all chromosomes according to heuristics layed out by the orderChromosome function.
	   This function will also plot the contact map for each chromosome according to the best order and orientation found.
	   Input = full genome adjacency matrix, [ [ [binID,scaffold], ... ] , ... ] list of binIDs/Scaffold for each scaffold within a chromosome,
			   list of Bin objects corresponding to all genomic loci in the contact map,
			   resolution of contact map (should be an integer value),
			   nScaffolds = number of largest scaffolds within a chromosome to brute force all permutations of their ordering
			   scanScaffolds = number of scaffolds within sliding window at the end of the algorithm to find all perumations for,
			   plotChrom = True or False whether or not we want to plot the contact map,
			   showPlot = True or False whether or not we want to visualize the plot in realtime,
			   savePlotDir = False or filePath of directory we would like to save our chromosome plots,
			   plotTitleSuffix = False or suffix we want to add to the end of each plot title
		Output = [ [scaffoldObj , ... ] , ... ] list of ordered/oriented scaffold objects for each chromosome (ordered from largest to smallest presumably)'''
	startTime = time.time()
	fullGenomeOrder = []
	chromosomeGroupOrder = []
	matrix = numpy.asmatrix(matrix)
	for i,chromGroup in enumerate(chromList):
		print("#####################"+'\n'+"#####################")
		print("Working on Chr_"+str(i+1)+"...")
		chromOrder = orderChromosome(chromGroup,matrix,binList,nScaffolds=nScaffolds,scanScaffolds=scanScaffolds)
		fullGenomeOrder.append(chromOrder)
		chromosomeGroupOrder.append([ [scaff.name,scaff.orientation] for scaff in chromOrder ])
		#####################
		if plotChrom != True:
			continue
		chrName = "Chr_"+str(i+1)
		adjMat,orderDict = giveNewAdjMat(matrix,chromOrder,binList)
		plotModule.plotContactMap(adjMat,resolution=resolution,tickCount=11,highlightChroms=False,
								  wInches=24,hInches=24,lP=1,hP=98,reverseColorMap='',
								  showPlot=showPlot,savePlot=savePlotDir+"/"+chrName+".png",
								  title=chrName,titleSuffix=plotTitleSuffix)

	######################
	stopTime = time.time()
	print("RunTime for total genome with plotting and saving .pngs = "+str(stopTime-startTime))
	######################
	return fullGenomeOrder

def writeScaffoldOrderingsToFile(sOrderings,outFile):
	'''This function writes out the order / orientation of each scaffold for each chromosome.
	   Input = [ [scaffoldObject , ...] , ... ] list of scaffold objects for each chromosome, full file path to output file'''
	outFile = open(outFile,'w')
	chromCount,scaffsWritten = 0,0
	for scaffGroup in sOrderings:
		chromCount += 1
		outFile.write("### Chromosome grouping "+str(chromCount)+" ###"+'\n')
		for sObj in scaffGroup:
			outFile.write(sObj.name+'\t'+sObj.orientation+'\n')
			scaffsWritten += 1
	###############
	outFile.close()
	print("Chromosome groups written to file "+str(chromCount))
	print("Scaffolds written to file "+str(scaffsWritten))

def writeBinIDsOrderingToFile(scaffoldList,outFile):
	'''This function writes out the order of binIDs contained within a list of scaffoldObjects.
	   Mainly just to produce the ordering of a contact map without having to do all the previous steps.
	   Input = [ scaffoldObject , ... ] , full file path to output file '''
	outFile = open(outFile,'w')
	outFile.write("#ScaffoldID"+'\t'+'HiCPro-BinID')
	binsWritten = 0
	for sObj in scaffoldList:
		for bId in sObj.binList:
			outFile.write('\n')
			outFile.write(sObj.name+'\t'+str(bId))
			binsWritten += 1
	###############
	outFile.close()
	print("BinIDs written to file "+str(binsWritten))

def getChromosomeOutlineCoords(orderedChromosomes):
	'''This function takes in the ordered chromosome groups output from the orderGenome function.
	   It figures out what the outline coordinates of each chromosome is for a contact map to be plotted.
	   Input = [ [scaffoldObj , ...] , ... ]
	   Output = [ chromosome1Coordinate, chromosom2Coordinate, etc... ]'''
	coords,index = [],0
	for scaffGroup in orderedChromosomes:
		for scaffoldObj in scaffGroup:
			index += len(scaffoldObj.binList)
		####################
		coords.append(index)
	#############
	return coords


########################################################
### Function that runs the entire pipeline for Part2 ###
def runPipeline(hicProBedFile,hicProBiasFile,hicProMatrixFile,chromosomeGroupFile,chromosomeOrderFile,
				savePlotsDirectory,chromosomePlotSuffix,fullGenomePlot,fullGenomePlotTitle,plotOrderFile,
				nScaffolds,scanScaffolds,resolution):
	#################################################
	print("########################################")
	print("### Working on Part2 of the pipeline ###")
	startTime = time.time()
	### Part1 ##########################
	### Initiate our data structures ###
	binDict = readGroupingsToValidBins(chromosomeGroupFile) ### Generate a dictionary of { binIDs:'' } to use in the analysis based on assignment from the first step
	binList = initiateLoci(hicProBedFile,hicProBiasFile,binID_dict=binDict) ### Initiate our list of Bin objects (only including those that are in binDict)
	adjMat = buildAdjacencyMatrix(hicProMatrixFile,binList) ### Build our adjacency matrix
	chromosomeList = readChromsFromFile(chromosomeGroupFile) ### Initiate our list of chromosomes and the bins contained within each chromosome based on assignment from the first step
	#################################################################################
	### Part2 #######################################################################
	### We want to order and orient scaffolds within a chromosome for all chromosomes
	orderedChromosomes = orderGenome(adjMat,chromosomeList,binList,resolution,
									 nScaffolds=nScaffolds,scanScaffolds=scanScaffolds,plotChrom=True,
									 showPlot=False,savePlotDir=savePlotsDirectory,
									 plotTitleSuffix=chromosomePlotSuffix)
	################################################################################
	### Part3 ######################################################################
	### We want to plot the final results genome wide and write out results to files
	outLineCoords = getChromosomeOutlineCoords(orderedChromosomes)
	adjMat,orderedDict = giveNewAdjMat(numpy.asmatrix(adjMat),[ scaffObj for scaffList in orderedChromosomes for scaffObj in scaffList ],binList)
	plotModule.plotContactMap(adjMat,resolution=resolution,tickCount=11,highlightChroms=outLineCoords,
							  wInches=32,hInches=32,lP=2,hP=98,reverseColorMap='',
							  showPlot=False,savePlot=fullGenomePlot,
							  title=fullGenomePlotTitle,titleSuffix=False)
	writeScaffoldOrderingsToFile(orderedChromosomes,chromosomeOrderFile)
	writeBinIDsOrderingToFile([ sObj for sObjs in orderedChromosomes for sObj in sObjs ],plotOrderFile)
	################################################################
	print("Total run-time  for Part2 = "+str(time.time()-startTime))
