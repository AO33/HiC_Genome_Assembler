###############
### Imports ###
import math # Needed in the Chromosome finding section
import numpy # Needed in the Chromosome finding section
import scipy # Needed in the Chromosome finding section
import scipy.spatial # Needed in the chromosome finding section (clustering/dendrogram)
import scipy.cluster # Needed in the chromosome finding section (clustering/dendrogram)
import sys
sys.setrecursionlimit(100000) ### The higher the value the more recursion calls can be made when clustering
from hmmlearn import hmm # Needed in the Chromosome finding section
import time # Just generally useful and is used in most sections
import copy # For ordering and orienting scaffolds
import networkx as nx # Needed in the Chromosome finding section
import community # Needed in the Chromosome finding section
import gzip # Needed in the Scaffold assignment section
import plotContactMaps as plotModule
####################################


############################################################
### Contact map initiation, and matrix manipulation code ### 
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
					newMat[i][ii] = math.log(float(v)+1.,logBase)
	#####
	else:
		for i,row in enumerate(matrix):
			for ii,v in enumerate(numpy.asarray(row)[0]):
				if v != 0.:
					newMat[i][ii] = (float(logBase)**v)-1.
	#############################
	return numpy.asmatrix(newMat)


#############################
### UPGMA clustering code ### 
def averageClusterNodes(adjacencyMatrix,nodeLabels,noPlot=True):
	'''Performs UPGMA (average) heirarchical clustering for a distance transformed adjacency matrix
	   The actual matrix passed in is not altered in any way. It must be roered in another function
	   Input = distance transformed adjacencyMatrix
			   labels of the nodes of the dendrogram that is produced
	   Output = dendrogram object produced by scipy.cluster.dendrogram'''
	startTime = time.time()
	condensedMatrix = scipy.spatial.distance.squareform(adjacencyMatrix,checks=False)
	stopTime = time.time()
	print("Time to condensed matrix "+str(stopTime-startTime))
	clustMan = scipy.cluster.hierarchy.average(condensedMatrix)
	if noPlot == False:
		fig = plt.figure()
		fig.set_size_inches(32,8)
		dendro = scipy.cluster.hierarchy.dendrogram(clustMan,labels=nodeLabels,leaf_rotation=90,no_plot=noPlot,get_leaves=True,count_sort ='descendent') ### This is what you use for clustering bins to chroms
		plt.show()
	else:
		dendro = scipy.cluster.hierarchy.dendrogram(clustMan,labels=nodeLabels,leaf_rotation=90,no_plot=noPlot,get_leaves=True,count_sort ='descendent') ### This is what you use for clustering bins to chroms
	######################
	stopTime = time.time()
	print("Time to cluster "+str(stopTime-startTime))
	return dendro

def dendrogramLeafOrder_toFile(dendrogramObj,outFile):
	'''This is to simply write out the order of the nodes of the dendrogram produced by the averageClusterNodes function
	   Input = dendrogram object from scipy, full file path to output file'''
	outFile = open(outFile,'w')
	for i,l in enumerate(dendrogramObj['ivl']):
		if i != (len(dendrogramObj['ivl'])-1):
			outFile.write(l+'\t'+str(dendrogramObj['leaves'][i])+'\n')
		else:
			outFile.write(l+'\t'+str(dendrogramObj['leaves'][i]))
	###############
	outFile.close()

def readDengrogramLeavesFromFile(dendrogramFile):
	'''This simply reads in the output file from the dendrogramLeafOrder_toFile
	   Input = full file path to dendrogramFile from dendrogramLeafOrder_toFile function
	   Output = a dictionary with the same notation as the scipy dendrogram object, but only "ivl" and "leaves" are used" '''
	inFile = open(dendrogramFile,'r')
	dendoLeaves = { 'ivl':[],'leaves':[] }
	for line in inFile:
		cols = line.strip('\r').strip('\n').split('\t')
		dendoLeaves['ivl'].append(cols[0])
		dendoLeaves['leaves'].append(int(cols[-1]))
	##############
	inFile.close()
	return dendoLeaves


#################################################################
### Chromosome assignment of bins via HMM and modularity code ###
def identifyBoundry(hiddenStates,cutIndices,switchCount=10):
	'''Simple heuristic function that identifies when a TRUE switch in states occured from the hidden markov model
	   While each state has technically already predicted at this point, there is still some very small level of noise.
	   The idea here is that when we see a dramatic switch from the start state to the next state we return the index where that happens.
	   Input = decoded states of observations, matrix cutIndices
			   switchCount defines how many consecutive non-startState observations we must see to then define the boundary
	   Output = newely defined boundary'''
	### We define our start state by the majority observation state of the first switchCount number of observations ###
	countDict = { 0:0, 1:0 }
	for s in hiddenStates[0:switchCount]:
		countDict[s] += 1
	startState = sorted([ [s,c] for s,c in countDict.items() ], key=lambda x: x[1], reverse=True)[0][0]
	################################################################################################
	states = [ hiddenStates[ind:ind+switchCount] for ind in range(0,len(hiddenStates)-switchCount) ]
	cutInd = 0
	###############################
	for ind,s in enumerate(states):
		switchState = sum([ 1 for hS in s if hS != startState ])
		if switchState == switchCount:
			cutInd = ind+cutIndices[-1]
			break
	#############
	return cutInd

def hmmChromosomes(adjacencyMatrix,cutIndices,binList,minSize=15,convergenceRounds=8,lookAhead=False):
	'''This function will be called to define a chromosome boundary in the UPGMA clustered matrix (but the logTransformed version NOT the distance version)
	   In essence, we want to use a 2-state GaussianHMM where the first state (starting state) belongs to a single chromosome, and the 2nd state belongs to the remaining chromosome(s).
	   But, there is a caveat in that sometimes the HMM will define multiple chromosomes belonging to the first state. If this happens, we remove all data pertaining to the second state,
	   and then rerun the algorithm. We repeat this process until the algorithm/HMM converges on a single point in the matrix that separates the starting and second states.
	   Only matrix rows/column indices of the matrix will be used corresponding the last entry of cutIndices and greater.
	   Input = logTransformed adjacencyMatrix, list of chromosomeBoundaries (cutIndices) and the list of Bin objects
			   convergenceRounds defines how many times we are allowed to run the algorithm (so we don't get stuck in an infinite loop, although this shouldn't happen)
			   minSize defines the minimum chromosome size (in the number of genomic loci from HiCPro)
			   lookAhead is the fraction of the matrix that we want to initially use for our HMM. The full matrix is recommended, but smaller values can improve runtime.
	   Output = list of chromosomeBoundaries (cutIndices) with a newely appeneded boundary)'''
	adjacencyMatrix = numpy.asarray(adjacencyMatrix)
	######################
	if lookAhead != False:
		lookAhead = int(float(len(adjacencyMatrix)) * lookAhead) + cutIndices[-1]
	else:
		lookAhead = len(adjacencyMatrix)
	###################################
	prevCutInd,roundCount = lookAhead,1
	while roundCount <= convergenceRounds:
		### This is for the termination of the algorithm if moudularity is set to zero percent of the data ###
		if (len(adjacencyMatrix)-cutIndices[-1])/2 < minSize:
			cutInd = prevCutInd
			cutIndices.append("NA")
			break
		##############################################################################
		X = [ r[cutIndices[-1]:prevCutInd] for r in adjacencyMatrix[cutIndices[-1]:] ]
		startMat = numpy.asarray([  1. , .0 ] )
		transMat = numpy.asmatrix([ [ .9 , .1 ], [ 0.0001 , .9999 ] ])
		##############################################################
		print("Input matrix size = "+str(len(X))+" x "+str(len(X[0])))
		print("HMM round = "+str(roundCount))
		model = hmm.GaussianHMM(n_components=2,covariance_type="diag",n_iter=1000,init_params="cm",params="cmt") # tested tried and true
		model.startmat_ = startMat
		model.transmat_ = transMat
		model.fit(X)
		hiddenStates = model.predict(X)
		cutInd = identifyBoundry(hiddenStates,cutIndices,switchCount=minSize)
		########################
		if cutInd != prevCutInd:
			prevCutInd = cutInd
			roundCount += 1
			continue
		else:
			print("HMM convergence rounds = "+str(roundCount))
			cutIndices.append(int(cutInd))
			break
	##################################
	if roundCount > convergenceRounds:
		cutIndices.append(int(cutInd))
		print("WARNING... HMM failed to converge after "+str(roundCount)+" rounds...")
		print("Proceeding with last found cutIndex of "+str(cutInd)+"...")
	#################
	return cutIndices

def modularityChromosomes(adjacencyMatrix,binList,cutIndices,louvainRounds=20):
	'''This function will identify communities (ie chromosomes) by using the Louvain method to optimize the modularity score of the input matrix
	   This is used to define the chromosome boundaries of the smaller chromosomes usually.
	   We convert our matrix to a network/graph and then identify wich matrix rows belong to which community/chromosome.
	   Input = fullSize adjacency matrix, list of Bin objects, chromosomeBoundaries/cutIndices
			   The last entry in cutIndices will be the starting index of the adjacencyMatrix generated from the full adjacencyMatrix
	   Output = list of chromosomeBoundaries with the newely found ones via modularity maximization'''
	startTime = time.time()
	startIndex = cutIndices[-1]
	adjacencyMatrix = numpy.asarray(adjacencyMatrix)
	subGraph = nx.Graph()
	indexDict = { i:binObj.ID for i,binObj in enumerate(binList[startIndex:]) }
	### Add nodes to graph ############
	for binObj in binList[startIndex:]:
		subGraph.add_node(binObj.ID)
	### Add edges to graph ############################
	for i,r in enumerate(adjacencyMatrix[startIndex:]):
		for ii,v in enumerate(r[startIndex:]):                 
			subGraph.add_edge(indexDict[i],indexDict[ii],weight=v)
	###########################################
	print("Maximizing so-called modularity...")
	print("Graph created with "+str(subGraph.number_of_nodes())+" nodes, and "+str(subGraph.number_of_edges())+" edges")
	print("Performing "+str(louvainRounds)+" rounds of the louvain method...")
	##################
	bestModScore = -2.
	for i in range(0,louvainRounds):
		chromosomePartition = community.best_partition(subGraph,randomize=True)
		modScore = community.modularity(chromosomePartition,subGraph)
		if modScore > bestModScore:
			bestModScore = modScore
			bestPartition = chromosomePartition
			print("Louvain round "+str(i+1)+" , "+"ModScore "+str(modScore))
	####################################################################################
	communityDict = { node:communityKey for node,communityKey in bestPartition.items() }
	chromCount = len({ communityKey:'' for communityKey in bestPartition.values() })
	prevComm = communityDict[indexDict[0]]
	for i in range(1,len(indexDict)):
		comm = communityDict[indexDict[i]]
		if comm != prevComm:
			cutIndices.append(i+startIndex)
			prevComm = comm
	#################################
	modTime = time.time() - startTime
	print("Chromosomes found via modularity maximization = "+str(chromCount))
	print("Modularity maximization total time = "+str(modTime))
	return cutIndices

def identifyChromosomeGroups(adjacencyMatrix,binList,minSize=15,modularity=.05,convergenceRounds=8,lookAhead=.5,louvainRounds=20):
	'''This function essentially wraps around the hmmChromosomes and modularityChromosomes functions
	   to iteratively define chromosome boundaries in the input adjacencyMatrix.
	   Input = logTransformed adjacencyMatrix, list of Bin objects
			   minSize defines the minimum chromosome size (in the number of genomic loci from HiCPro)
			   lookAhead is the fraction of the matrix that we want to initially use for our HMM. The full matrix is recommended, but smaller values can improve runtime
			   modularity is the fraction of the lower right hand side of the matrix to perform modularity calculations on. Can be set to zero, but this might give weird results...
			   louvainRounds is the number of initiations of the Louvain method to perform. 
							 We randimize the order in which we iterate through the graph nodes each round of the Louvain algorithm to avoid potential suboptimal Modularity scores
	   Output = list of chromosomeBoundaries (indices of the input adjacencyMatrix that define our chromosome boundaries)'''
	print("#########################"+'\n'+"#########################")
	print("Working on iterative 2 state HMMs to identify chromosome boundaries...")
	startTime = time.time()
	adjacencyMatrix = scipy.asmatrix(adjacencyMatrix)
	matrixLength = float(len(adjacencyMatrix))
	remainder = matrixLength - (modularity*matrixLength)
	modFallBack,cutIndices = 0,[0]
	##################################
	while cutIndices[-1] <= remainder:
		print("#########################"+'\n'+"#########################")
		cutIndices = hmmChromosomes(adjacencyMatrix,cutIndices,binList,
									minSize=minSize,
									convergenceRounds=convergenceRounds,
									lookAhead=lookAhead)
		########################################
		print("Cut indices =  "+str(cutIndices))
		if cutIndices[-1] == 0:
			print("Algorithm terminated. No obvious chromome boundry could be found... ")
			break
		##########################
		if cutIndices[-1] == "NA":
			modFallBack = 1
			cutIndices.pop(-1)
			break
	######################
	if cutIndices[0] == 0:
		cutIndices.pop(0)
	###################################################################
	print("#########################"+'\n'+"#########################")
	print("HMM rounds completed in "+str(time.time()-startTime)+" seconds")
	print("Chromosome groups found via HMMs "+str(len(cutIndices))+" / "+str(len(cutIndices)+1))
	####################
	if modFallBack == 0:
		print("#########################"+'\n'+"#########################")
		cutIndices = modularityChromosomes(adjacencyMatrix,binList,cutIndices,louvainRounds=louvainRounds)
		print("#########################"+'\n'+"#########################")
	#############################################################################################
	print("Total time to identify chromosome boundries = "+str(time.time()-startTime)+" seconds")
	print("Chromosomes found = "+str(len(cutIndices)+1))
	#################
	return cutIndices

def writeBinGroupingsToFile(coords,binList,outFile):
	'''This function writes the genomic loci (Bin objects) associated with each chromosome found to a file.
	   Just a simple intermediate file that is then used in the final step of scaffold to chromosome assignment. Chromosome names in this file are just place holders. NOT final
	   Input = list of chromosomeBoundaries/matrix cutIndices identifed by the identifyChromosomeGroups function
			   list of Bin objects that contain the actual informatin we will be writing out
			   full file path to the output file'''
	chromGroups = []
	prevIndex = 0
	for cutIndex in coords:
		chromGroups.append(binList[prevIndex:cutIndex])
		prevIndex = cutIndex
	chromGroups.append(binList[prevIndex:])
	###########################
	outFile = open(outFile,'w')
	for i,chrGroup in enumerate(chromGroups):
		outFile.write("### Chromosome group "+str(i+1)+" ###"+'\n')
		for binObj in chrGroup:
			outFile.write(str(binObj.ID)+'\t'+binObj.chrom+'\t'+str(binObj.start)+'\t'+str(binObj.stop)+'\t'+str(binObj.bias)+'\n')
	###############
	outFile.close()


########################################################################################
### Scaffold assignment to chromosome and renaming of chromosomes based on size code ###
def readSizeFileToDict(sizeFile):
	'''Takes the hicPro scaffold size file and reads it into memory"
	   Input = full file path to scaffold size file
	   Output = { scaffoldName:size }'''
	scaffDict = {}
	inFile = open(sizeFile,'r')
	for line in inFile:
		cols = line.strip('\r').strip('\n').split('\t')
		scaffDict.update({ cols[0]:int(cols[1]) })
	##############
	inFile.close()
	return scaffDict

def readBinGroupingsFromFile(binGroupingsFile):
	'''Reads in the bin/chromosome groups file from the initital UPGMA clustering-->HMM/Modularity chromosome group assignment.
	   Input = binGroupings full file path
	   Output = A list of chromosome groups, where each element in the list is another list, where each element is a line from the input file'''
	chromList,chromGroup = [],[]
	inFile = open(binGroupingsFile,'r')
	inFile.readline()
	for line in inFile:
		line = line.strip('\n').strip('\r')
		if line[0] != "#":
			chromGroup.append(line)
		else:
			chromList.append(chromGroup)
			chromGroup = []
	############################
	chromList.append(chromGroup)
	inFile.close()
	print(str(len(chromList))+" chromosomes read in from file")
	return chromList

def assessClusterList(cList,scaffDict,outFile,percentToAssign=51.):
	'''One of these function calls is used per chromosome group found, and assigns scaffolds to a particular chromosome group.
	   Input = cList is the list of lines read in from the binGroupings file (contains information about a particular bin, its ID and scaffold)
			   scaffDict = { scaffoldName:[ binIDs belonging to the scaffold ] }
			   outFile = variable name with which an output file has been opened for writing with. Note!!! The chromosome group names are not necessarily final.
						   They are only finalized in the writeChromosomeGroupingsToFile function
			   percentToAssign = percentage of bins belonging to a scaffold that are found in the cList input argument for confident assignment to this chromosome group
	   Output = List of binIDs assigned to chromosome
				List of false positive bin IDs that do not belong to the chromosome
				Number of scaffolds found assigned to chromosome'''
	scaffs = {}
	for c in cList:
		binID,scaff = int(c.split('\t')[0]),c.split('\t')[1]
		if scaff not in scaffs:
			scaffs.update({scaff:[binID]})
		else:
			scaffs[scaff].append(binID)
	####################
	finalClustering = []
	totalScaffoldsAssigned = 0
	falsePositives = 0
	outFile.write("#Scaffold"+'\t'+"NodesAssigend"+'\t'+"TotalNodes"+'\t'+"Assigned%"+'\n')
	for s,nodeList in scaffs.items():
		nodesAssigned,totalNodes = len(nodeList),len(scaffDict[s])
		percentAssigned = round(((float(nodesAssigned)/float(totalNodes)) * 100.),2)
		outFile.write(str(s)+'\t'+str(nodesAssigned)+'\t'+str(totalNodes)+'\t'+str(percentAssigned)+"%"+'\n')
		if percentAssigned >= percentToAssign:
			finalClustering += scaffDict[s]
			totalScaffoldsAssigned += 1
		else:
			falsePositives += nodesAssigned
	###############################################################################
	outFile.write("Total scaffolds clustered to chromosome "+str(len(scaffs))+'\n')
	outFile.write("Total scaffolds assigned to chromosome "+str(totalScaffoldsAssigned)+'\n')
	############################################################
	return finalClustering,falsePositives,totalScaffoldsAssigned

def assessChromosomeClustering(chromList,statsFile,percentToAssign=51.):
	'''Assignes scaffolds to chromosomes.
	   Input = output of the readBinGroupings file
			   percentToAssign = percentage of bins belonging to a scaffold that are found in the cList input argument for confident assignment to this chromosome group
			   statsFile = full filePath to stats regarding the clustering assessment. Note!!! The chromosome group names are not necessarily final.
						   They are only finalized in the writeChromosomeGroupingsToFile function
	   Output = List of lists of binIDs belonging to any given chromosome group'''
	### First figure out how many nodes are contained within each scaffold ###
	scaffolds = {}
	fullNodeList = [ cc for c in chromList for cc in c]
	for node in fullNodeList:
		binID,scaff = int(node.split('\t')[0]),node.split('\t')[1]
		if scaff not in scaffolds:
			scaffolds.update({scaff:[[binID,scaff]]})
		else:
			scaffolds[scaff].append([binID,scaff])
	for s in scaffolds:
		scaffolds[s] = sorted(scaffolds[s],key=lambda x: x[0])
	##############################################################################
	### For each chromosome found, assign scaffolds based on simple heuristics ###
	finalChromList = []
	falsePositives = 0
	totalScaffsAssigned = 0
	outFile = open(statsFile,'w')
	for i,c in enumerate(chromList):
		outFile.write("### Chromosome"+str(i+1)+" ###"+'\n')
		chromosomeNodes,falsePos,scaffsAssigned = assessClusterList(c,scaffolds,outFile,percentToAssign=percentToAssign)
		if len(chromosomeNodes) > 0:
			finalChromList.append(chromosomeNodes)
		falsePositives += falsePos
		totalScaffsAssigned += scaffsAssigned
		outFile.write("####################"+'\n')
	##############################
	totalNodes = len(fullNodeList)
	outFile.write("Total Nodes "+str(len(fullNodeList))+'\n')
	outFile.write("Properly clustered nodes "+str(totalNodes-falsePositives)+'\n')
	outFile.write("Falsely clustered nodes "+str(falsePositives)+'\n')
	outFile.write("Total scaffolds assigned to chromosomes "+str(totalScaffsAssigned)+'\n')
	outFile.write("Error rate ~"+str(round((float(falsePositives)/float(totalNodes)) * 100., 2))+"%"+'\n')
	return finalChromList

def writeChromosomeGroupingsToFile(chromList,scaffSizeDict,outFile):
	'''Writes out our chromosome groups to a new file that is sorted by the size of the chromosome (biggest chromosome = Chr_1)
	   Very similar to the binGroupings file, but ordered by size and doesnt contain the bias values.
	   Input = output from the assessChromosomeClustering function
			   { scaffoldName:sizeInBasePairs }
			   full file path of output file '''
	### Find the sizes of all the chromosome groups ###
	sizes = []
	for c in chromList:
		scaffs = { cc[1]:'' for cc in c }
		sizes.append(sum([ scaffSizeDict[cc] for cc in scaffs ]))
	###############################################################################################
	### Then sort our list of chromosome groups by size ###########################################
	chromList = [ c for c,s in sorted(zip(chromList,sizes),key=lambda pair: pair[1],reverse=True) ]
	###########################
	outFile = open(outFile,'w')
	for i,c in enumerate(chromList):
		outFile.write("### Chromosome group "+str(i+1)+" ###"+'\n')
		for cc in c:
			outFile.write(str(cc[0])+'\t'+str(cc[1])+'\n')
	###############
	outFile.close()


########################################################
### Function that runs the entire pipeline for Part1 ###
def runPipeline(hicProBedFile,hicProBiasFile,hicProMatrixFile,hicProScaffSizeFile,
				dendrogramOrderFile,avgClusterPlot,avgClusterPlot_outlined,
				binGroupFile,assessmentFile,chromosomeGroupFile,
				minSize,modularity,convergenceRounds,louvainRounds,resolution):
	#################################################
	print("########################################")
	print("### Working on Part1 of the pipeline ###")
	totalStartTime = time.time()
	### Part 1 ############
	startTime = time.time()
	binList = initiateLoci(hicProBedFile,hicProBiasFile)
	adjMat = buildAdjacencyMatrix(hicProMatrixFile,binList)
	adjMat,binList = removeRows(adjMat,binList,zeroRows=True,biasVals=False)
	adjMat = convertMatrix(adjMat,binList,distance=True,similarity=False)
	dendroLabels = [ binList[i].chrom+'_'+str(binList[i].ID) for i in range(0,len(binList)) ]
	dendrogram = averageClusterNodes(adjMat,dendroLabels,noPlot=True)
	dendrogramLeafOrder_toFile(dendrogram,dendrogramOrderFile)
	dendoLeaves = readDengrogramLeavesFromFile(dendrogramOrderFile)
	adjMat,binList = reorderMatrix(scipy.asmatrix(adjMat),binList,dendoLeaves['leaves'])
	plotModule.plotContactMap(adjMat,resolution=resolution,highlightChroms=False,showPlot=False,savePlot=avgClusterPlot)
	stopTime = time.time()
	print("Total run-time to cluster and plot = "+str(stopTime-startTime))
	#######################
	### Part 2 ############
	startTime = time.time()
	adjMat = convertMatrix(adjMat,binList,distance=False,similarity=True)
	adjMat = logTransformMatrix(adjMat,logBase=10)
	cutIndices = identifyChromosomeGroups(adjMat,binList,minSize=minSize,modularity=modularity,convergenceRounds=convergenceRounds,lookAhead=False,louvainRounds=louvainRounds)
	adjMat = convertMatrix(logTransformMatrix(adjMat,logBase=10,reverse=True),binList,distance=True)
###	cutIndices = [1387, 2596, 3724, 4794, 5824, 6625, 7611, 9683, 10423, 11281, 11837, 12468, 12908, 13132, 13356, 13626, 13843, 14090, 14275, 14436, 14562, 14661]
	plotModule.plotContactMap(adjMat,resolution=resolution,highlightChroms=cutIndices,showPlot=False,savePlot=avgClusterPlot_outlined)
	writeBinGroupingsToFile(cutIndices,binList,binGroupFile)
	stopTime = time.time()
	print("Total run-time to identify chromosome boundaries = "+str(stopTime-startTime))
	#######################
	#######################
	### Part 3 ############
	startTime = time.time()
	fastaSizeDict = readSizeFileToDict(hicProScaffSizeFile)
	binGroups = readBinGroupingsFromFile(binGroupFile)
	chrGroups = assessChromosomeClustering(binGroups,assessmentFile)
	writeChromosomeGroupingsToFile(chrGroups,fastaSizeDict,chromosomeGroupFile)
	stopTime = time.time()
	print("Total run-time to assign scaffolds to chromosomes = "+str(stopTime-startTime))
	############################
	totalStopTime = time.time()
	print("Total run-time of Part1 = "+str(totalStopTime-totalStartTime))
	print("CutIndices = "+str(cutIndices))
