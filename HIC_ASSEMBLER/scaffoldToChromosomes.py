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
import collections # Needed in tallying modularity partition data
import plotContactMaps as plotModule
####################################


############################################################
### Contact map initiation, and matrix manipulation code ### 
### Matrix Initialization and manipulation code ###
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
	adjacencyMatrix = numpy.asmatrix(adjacencyMatrix)
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
	matrix = numpy.asmatrix(matrix)
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
		dendro = scipy.cluster.hierarchy.dendrogram(clustMan,labels=nodeLabels,leaf_rotation=90,no_plot=noPlot,get_leaves=True,count_sort ='ascending')
		plt.show()
	else:
		dendro = scipy.cluster.hierarchy.dendrogram(clustMan,labels=nodeLabels,leaf_rotation=90,no_plot=noPlot,get_leaves=True,count_sort ='ascending')
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

####################################################################################
### Chromosome assignment of bins via Hyper geometrics, HMM, and modularity code ###
### Modularity assignment code ###
def modularity_rounds(graph, louvain_rounds=1, induce_comm_subgraph=False):
	"""
	Performs multiple rounds of the louvain method (with random start orderings) and returns the best partition (i.e. maximimum modularity score partion found)
	"""
	
	# Set value to impossible score for start condition
	best_mod_score = -2.
	
	# Expects a {community:nodeID} formatted partition 
	if induce_comm_subgraph != False:
		graph = community.induced_graph(induce_comm_subgraph, graph)
	
	# Random modularity rounds for best result
	for i in range(0, louvain_rounds):
		chromosome_partition = community.best_partition(graph, randomize=True)
		mod_score = community.modularity(chromosome_partition, graph)
		if mod_score > best_mod_score:
			borg = best_mod_score
			best_mod_score = mod_score
			best_partition = chromosome_partition
			print("Previous best modularity score {}, Current best found {}, Louvain round {}".format(borg, mod_score, i+1))
	
	return best_partition, best_mod_score
	
def modularity_remaining_data(adjacencyMatrix, binList, cutIndices, n_rounds=20):
	'''
	Attempts to take the remaining data in the adjacency matrix, presumably small scale groups and partions the subgraph
	via modularity maximization. The ordering of the passed in adjacency matrix is preserved EXCEPT for the remaininng portion,
	which is resolved via the modulity partition.
	
	Both input matrix and binList are reordered (the remaining portion), and updated cutIndices are returned.
	Note, the submatrix used for modularity will begin at the last cutindex found
	'''
	
	startTime = time.time()
	if len(cutIndices) == 0:
		print("- Attempting to resolve groupings by modularity alone... This could take a while if matrix size is large and n_rounds is set high as well...")
		cutIndices = [0]
	
	# Ensure proper cutindex formatting
	cutIndices = sorted(cutIndices)
	startIndex = cutIndices[-1]
	adjacencyMatrix = numpy.asarray(adjacencyMatrix)
	
	# Create local and global index lookup for bin objects
	indexDict = {i:binObj.ID for i, binObj in enumerate(binList[startIndex:])}
	bin_to_index = {binObj.ID:i for i, binObj in enumerate(binList)}
	index_to_bin = {i:binObj.ID for i, binObj in enumerate(binList)}
	
	# Create graph of remaining data
	subGraph = nx.Graph()
	
	# Add nodes to graph
	for binObj in binList[startIndex:]:
		subGraph.add_node(binObj.ID)
	
	# Add edges to graph
	for i,r in enumerate(adjacencyMatrix[startIndex:]):
		for ii,v in enumerate(r[startIndex:]):                 
			subGraph.add_edge(indexDict[i],indexDict[ii],weight=v)
			
	# Progress report
	print("- Maximizing so-called modularity...")
	print("- Graph created with "+str(subGraph.number_of_nodes())+" nodes, and "+str(subGraph.number_of_edges())+" edges")
	print("- Performing "+str(n_rounds)+" rounds of the louvain method...")
	
	# Perform louvain rounds for best partition
	partition, mod_score = modularity_rounds(subGraph, louvain_rounds=n_rounds, induce_comm_subgraph=False)
	
	# Process output partition to groups and cut indices for plotting
	node_to_group = partition
	group_sizes = collections.Counter(list(node_to_group.values()))
	group_count = len(group_sizes)
	
	# By using this method, if the ordering of the matrix doesn't match the ordering found by the modularity method,
	# then more cut indices than there are groups can be obtained. 
	
	# Arrange large to small (create and sort group sizes list, then take the transpose to get back just the groupID)
	remaining_groups = numpy.asarray(sorted([[k,v] for k, v in group_sizes.items()], key=lambda x: x[1], reverse=True)).T[0]
	
	# Make a new order for the matrix and binLists passed in (needs to be of original index)
	remaining_order = []
	for rg in remaining_groups:
		for remaining_bin in binList[startIndex:]:
			if node_to_group[remaining_bin.ID] == rg:
				remaining_order.append(bin_to_index[remaining_bin.ID])
		
		# Add a new cut index for plotting and group partitioning
		cutIndices.append(cutIndices[-1]+group_sizes[rg])
	
	# Now reorder the matrix and binList for plotting
	new_order = [bin_to_index[bb.ID] for bb in binList[0:startIndex]] + remaining_order
	adjacencyMatrix, binList = reorderMatrix(adjacencyMatrix, binList, new_order)
	
	# Remove first and last cut index if equal to bounds of data sizes for reporting and plotting purposes
	# Groups found should be equal to the number of cutIndices + 1
	if cutIndices[0] == 0:
		cutIndices.pop(0)
	
	if cutIndices[-1] == len(adjacencyMatrix):
		cutIndices.pop(-1)

	
	# Progress report
	modTime = time.time() - startTime
	total_groups = len(cutIndices)+1
	print("- Modularity maximization total time = "+str(modTime))
	print("- Chromosomes found via HMMs or Hyper geometrics = {}".format(total_groups - group_count))
	print("- Chromosomes found via modularity maximization = "+str(group_count))
	print("- Total chromosomes found {}".format(total_groups))
	return adjacencyMatrix, binList, cutIndices

### Hyper geometric "dendrogram" cutting strategy (note only the relative ordering of dengroram is used)
def hyper_geom(x, M, n, N):
	"""
	x = How many specific objects we want to sample (or more)
	M = How many objects are in the bag
	n = Number of objects matching specification in the bag
	N = How many times we draw from the bag
	global_index = Specific to this particular function for multiprocessing and retrieving relevant data later
	
	https://alexlenail.medium.com/understanding-and-implementing-the-hypergeometric-test-in-python-a7db688a7458
	"""
	
	# Frozen distribution
	rv = scipy.stats.hypergeom(M, n, N)
	
	# Inverse of the Cumulitive distribution function to get the probability (survival function)
	prob = scipy.stats.hypergeom.sf(x-1, M, n, N)
	return prob

def get_sliding_window_distance_metrics(input_array, window_size, global_index=0):
	"""
	The input array is scanned with a window size set to (min_size * 2). 
	The array within each window is then split in half, and the left half is subtracted from the right and the resulting array is summed,
	to yield a value bound between (+-window_size * maximum arrayvalue). For example, [1,1,1,1,10, 0,1,0,0,0] with a min_size of 3 would yield a maximum distance of (1,1,10) - (0,1,0) = 11.
	The distances for each window are recording (at a step size of 1) and the maximum (+) distance index is then computed.
	"""
	
	
	# Ensure proper formatting of input information
	if window_size >= len(input_array):
		return ["NA", "NA", "NA"]
	
	if type(input_array) != numpy.asarray([0]):
		input_array = numpy.asarray(input_array)
	
	# Distance is measured by taking the first half of data within the window and subtracting the second half, and then taking the sum of that
	# -+ values should be retained because we are moving from left to right and lower values are expected using this method to be on the right than the left.
	# This elevates the background removal signal and allows higher accuracy on smaller scales
	half_dist = int(window_size)
	break_sigs = []
	half_dist2x = half_dist+half_dist
	for i in range(0, len(input_array)-half_dist):
		a1 = input_array[i:i+half_dist]
		a2 = input_array[i+half_dist:i+half_dist2x]
		if a1.shape != a2.shape:
			break_sig = 0
		else:
			break_sig = sum(input_array[i:i+half_dist] - input_array[i+half_dist:i+half_dist2x])
			
		break_sigs.append(break_sig)
	
	# return the the distances, the window size / 2, and the index to break data on according to the maximum distance found
	# If more than one maximum is found (ie datapoints with same values) then the left most (minimum index) is returned
	### If more than one maximum is found, then a sorted list of the indices can be returned. By default we ruturn the first one here
	best_inds = numpy.argwhere(break_sigs == numpy.amax(break_sigs))
	best_ind = sorted(best_inds.flatten().tolist())[0]
	
	### Old / bugged way where first occurence of max is returned
	### half_dist + numpy.argmax(break_sigs)
	
	return break_sigs, half_dist, best_ind

def find_matrix_pvalue_breakpoints(argsorted_mat, start, min_size, world_size, psig=.05):
	"""
	The input to this algorithm is a square, rank ordered index matrix that retains the row ordering of some sort of clustering algorithm.
	For example, averge heirarchical clustering --> numpy argsort. The algorithm will identify rows along the square matrix from top to bottom (left to right) 
	where there is a jump in hyper geometric pvalues from 0. --> greater than .05 by default.
	
	The algorithm starts by scanning a single row, with only the first column included. The next iteration scans the second row but with the first two columns included,
	the third with the third row and first 3 columns etc...
	The number of indices within in a rows current window that are between our start and stop index are then calculated as parameter x (or K...)
	The population size of the remaining matrix yet to be processed after the current start index is used as parameter M.
	Both draw count and available success count are set to the current window size (n, N).
	A hyper geometric test is then calculated and all pvalues are recored.
	Note, that the row we start on is dynamically passed in so that we can iterate over this function.
	
	Next, the pvalues are converted to 1. if < .05, and 0. if >= .05 and are scanned with a window size set to (min_size * 2). 
	The array of 1s and 0s within each window is then split in half, and the left half is subtracted from the right and the resulting array is summed,
	to yield a maximum value of min_size. For example, [1,1,1,1,1, 0,1,0,0,0] with a min_size of 3 would yield a maximum distance of (1,1,1) - (0,1,0) = 2.
	Only places where the distance found is == min_size are returned
	"""
	
	
	M = world_size
	org_size = min_size
	ws = min_size
	break_sig = 0
	loop_count = 0
	
	# The double while True allows us to run the algorithm and check for breakpoints. 
	# If no breakpoints are found then we can repeat the algorithm with a smaller window size until we reach zero
	while True:

		# Find next cut site from most recent cut site found
		while True:

			# Loop through dataset to generate pvalues for windowing cuttoff
			dist_sigs = [0]
			for i in range(start, len(argsorted_mat)):

				if i == start:
					continue

				# Create params for hyper geom
				curr_size = i-start
				row = argsorted_mat[i]
				prow = row[:curr_size]
				xr = prow[(prow >= start) & (prow <= i)]
				x = len(xr)
				n, N = curr_size, curr_size

				# Perform hyper geom
				pval = hyper_geom(x, M, n, N)
				
				# This allows us to squash data so that our maximim distances fall between 0 and min_size
				if pval >= psig:
					dist_sigs.append(0)
				else:
					dist_sigs.append(1)

			# Figure out fraction of data that was found to be significant. At lower levels, we need to adjust world size to something smaller to get groups to pop out
			# If the world size is dynamically scaled with the remaining data size, then at lower resolutions (ie 100Kb vs 500Kb) the smaller M size can "choke" the signal.
			loop_count += 1
			if (sum(dist_sigs) / len(dist_sigs)) >= .9:
				morg = M
				M = int(M-start)
				print("- M value (world_size) changed to dynamic {} --> {}".format(morg, M))
			else:
				break_sig = 1

			# Break conditions
			if (break_sig == 1) or (loop_count >= 5):
				break

		# After looping through dataset, find indices where pvalues go from ~0 --> 1
		# These sites will either be noise or the true index we need to cut our matrix at to get optimal groupings
		wmetrics = get_sliding_window_distance_metrics(dist_sigs, window_size=ws)
		mdata = numpy.asarray([[wmetrics[0][ii], ii+min_size] for ii in range(0, len(wmetrics[0])) if wmetrics[0][ii] == min_size]).T
		
		if len(mdata) > 0:
			pre_cut_vals, pre_cut_inds = mdata[0], mdata[1]
		else:
			pre_cut_vals, pre_cut_inds = [],[]
		

		if len(pre_cut_inds) > 0:
			break

		# Repeat entire process with smaller window size
		else:
			mws = ws
			ws -= 1
			if ws == 0:
				print("- Warning - No cut index found after scanning through all window sizes between 1 and {}".format(ws, org_size))
				break

			else:
				print("- Warning - No cut index found with window size of {}, decreasing by one to {}".format(mws, ws))

	
	return pre_cut_vals, pre_cut_inds

def pre_process_all_matrix_breakpoints(argsorted_mat, min_size=5, min_frac=.05, psig=.05):
	"""
	This function is intented to scan through the argsorted matrix from left to right to find the left most breakpoint.
	The algorithm repeats the process starting at the previous left most breakpoint found until we run out of data
	to process.
	
	Breakpoints are found by calling the find_matrix_pvalue_breakpoints function which will calculate pvalues from the current start index
	to the remainder of the data.
	"""
	
	mat_size = len(argsorted_mat)
	stop_ind = int(mat_size - (mat_size * min_frac))
	ind = 0
	count = 0
	cinds = []

	# Means we want to use modularity instead (presumably, or a mistake)
	if min_frac == 1:
		return cinds
	
	# Find left most breakpoint and repeat. Note that we update our worldsize for the hyper geometric test at each new iteration of the algorithm
	while True:
		pre_cut_vals, pre_cut_inds = find_matrix_pvalue_breakpoints(argsorted_mat, ind, min_size, mat_size-ind, psig=.05)
		
		# No cut points found wich means we are at the end of the algorithm
		if len(pre_cut_inds) == 0:
			break
		
		# Otherwise, update our information and continue forward
		ind += pre_cut_inds[0]
		cinds.append(ind)
		print(ind, pre_cut_inds)

		# Also means we have reached the end of our algorithm (no more data to process essentially)
		if (ind >= stop_ind) or ((mat_size-ind) <= min_size):
			break
	
	print("- Breakpoints found {}".format(len(cinds)))
	return cinds

def filter_noisy_breakpoints(argsorted_mat, original_inds, psig=.05):
	"""
	This algorithm is designed to take an aggressive set of cutindices for a square matrix and smooth them out to the most probable set via hyper geometric tests,
	utilizing the rankordered matrix as a backboard against the current matrix ordering (presumably ordered via heirarchical clustering, but can be anything).
	In other words, we would like to test how good our input cut indices are by testing them against the rank ordered matrix via hyper geometric tests.
	
	We first calculate how many rows between the current cut index and previous start have a pvalue < psig via hyper geometic tests via the following setup...
	The number of indices within in the rankordered matrix rows current window that are between our start and stop index are calculated as parameter x (or K...)
	The population size of the remaining matrix yet to be processed after the current start index is used as parameter M.
	Both draw count and available success count are set to the current window size (n, N).
	A hyper geometric test is then calculated and all pvalues are recored.
	
	Next, we calculate how many of these rows we would expect to see significant pvalues for just by random chance between our previous and current / future cut indices. 
	We utilizie yet another hyper geometric test where the x parameter is calculated by simply summing the total number of significant pvalues we found from the first step between the desired indices.
	The population size of the remaining matrix yet to be processed after the current start index is used as parameter M.
	The draw count (n) is set to the difference between our start index, and the current cut index in question. 
	The number of success available in our population parameter (N) is set to the difference between the start index and the current cut index.
	Note, that unlike the find_matrix_pvalue_breakpoints function, N will always be bigger than n except for the first iteration.
	
	We take the right most significant pvalue cut index found, and remove all cut indices found to the left. The process then starts over using
	our right most significant index found as the start index.
	"""
	
	# Allows us to return data within the pipeline essentially without having to have a safegaurd outside of it
	
	if len(original_inds) == 0:
		return []

	# Predetermined parameters to prevent data from spinning to infinity if something goes completely awry
	# Note, that if data is massive (> 50,000 rows) or large number of cutsites is expected then these can be easily increased to accomodate
	MD = int(len(argsorted_mat) / 5)
	MAX_ROUNDS = 10 * len(original_inds)
	GLOBAL_MAX_ROUNDS = 10
	
	altered_inds = copy.copy(original_inds)
	prev_filtered_inds = {}
	
	# Allows us to "recurse" on this loop by simply calling it again until we get the same results back, or a maximum of 10 default rounds
	while True:
		
		# Set (or reset) variables for this iteration
		global_round_count = 0
		start = 0
		filtered_inds = {}
		
		# Take original indices and iteratively remove "noisy" indices until convergence
		round_count = 0
		while True:

			# Limit our while true loop in case something is wrong with data
			if round_count >= MAX_ROUNDS:
				print("- WARNING - Maximum number of rounds {} exceeded... Data appears to be extremely noisy or something went wrong".format(MAX_ROUNDS))
				break
			
			# Allows the algorithm to parse smaller groups by scaling the world size value back for geom tests
			M = len(argsorted_mat) - start
			pop_ind = "NA"
			noise_found = 0

			print("- Working on round {}... {} current filtered indices found".format(round_count, filtered_inds))
			for i, c in enumerate(altered_inds):

				# Set initial variables and data container for hyper geom pval calculations
				local_size = c - start
				n, N = local_size, local_size
				pvals = {}
				print("- Local_Size {}".format(local_size))

				# Loop through current and future cut indices to see if significant connections exist within truncated window of argsoted matrix
				for ii in range(start, len(argsorted_mat)):

					# Prevents pvalue switching at large scales past ~50% of our original dataset
					# In other words, this prevents us from looking too far forward where it doesn't make sense 
					if (ii - start) > MD:
						pvals.update({ii:0})
						continue

					# Otherwise we perform the hyper geometric test with smaller window for argsorted matrix for row in question
					x = sum([1 for v in argsorted_mat[ii][0:local_size] if (start <= v) and (v <= c)])
					pval = hyper_geom(x, M, n, N)
					if pval < psig:
						pvals.update({ii:1})
					else:
						pvals.update({ii:0})

				# Initialize variables for next section
				right_most = "NA"
				sigs = []
				fc_prev = start

				# Now figure out our right most coordinate that has signifcant connections to previous rows
				for ai_ind, ai in enumerate(altered_inds):
					ps = [pvals[ind] for ind in range(fc_prev, ai)]

					# Deal with edge case conditions for the algorithm
					# Means we've already determined this index is valid and we are at the start
					if ai == fc_prev:
						continue

					# Continue forward with new index and set previous
					else:
						fc_prev = ai

					# Means we have reached the very end and no new values can be obtained
					if len(ps) == 0:
						break

					# Calculate our pvalue here via another hyper geometeric based on how many we could expect by chance
					x = sum(ps)
					n = local_size
					N = len(ps)

					psp = round((sum(ps) / len(ps)) * 100, 2)
					noise_pval = hyper_geom(x, M, n, N)

					# Means our current index is not a good fit according to future indices that have significant high connections with older ones (scanning left to right)                
					if noise_pval < psig:
						right_most = ai
						sigs.append([ai, [x, M, n, N, noise_pval, psp]])
						right_most_ind = ai_ind

					print(ai, x, M, n, N, noise_pval, len(ps))

				# Means our current cut index row has significant ties with a cut index later on, so we want to jump forward to that one
				if len(sigs) > 0:
					start = right_most
					filtered_inds.update({right_most:''})
					noise_found = 1
					altered_select_inds = set([iii for iii in range(0, len(altered_inds)) if iii >= right_most_ind])
					print("- Right most sig pvalue coordinate found {}".format(right_most))
					print("- Sig coords found {}".format(sigs))
					break      

				# Means our current cut index row is the optimal one
				else:
					filtered_inds.update({c:''})
					altered_select_inds = set([iii for iii in range(0, len(altered_inds)) if iii >= i])

			# Update round count and check for convergence
			round_count += 1
			if noise_found == 0:
				print("- Exiting algorithm... No significant connections found between current inds")
				break

			# Otherwise, we need to update our indices and repeat the process (note we want to maintain our order from left --> right)
			else:
				altered_inds = [a for i, a in enumerate(altered_inds) if i in altered_select_inds]
				print("- New round with {}, {}".format(start, altered_inds[0]))
	
		
		# Repeat the process if convergence hasn't been found
		if prev_filtered_inds != filtered_inds:
			print("- Previous cutindices {} don't match current indices found... Repeating the algorithm for the {}th round".format(filtered_inds, global_round_count))
			
			# Reset variables
			altered_inds = sorted([v for v in filtered_inds])
			prev_filtered_inds = filtered_inds
			
			# Last exit condition to so we don't run to infinity
			if global_round_count >= GLOBAL_MAX_ROUNDS:
				print("- WARNING - Algorithm as not yet converged after {} rounds... Exiting prematurely due to limit of {} round".format(rc, GLOBAL_MAX_ROUNDS))
				break
			else:
				global_round_count += 1
		
		# Termination condition
		else:
			print("- Algorithm appears to have converged as previous cutindices match current cutindices. Exiting...")
			break
		
	# Return sorted (smallest first) list of cut indices
	return_inds = sorted(list(filtered_inds.keys()))
	print("- Original cut indices {}".format(original_inds))
	print("- Filtered cut indices {}".format(return_inds))
	return return_inds

### HMM --> modularity assigmnet based code ###
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

def hmmChromosomes(adjacencyMatrix,cutIndices,binList,minSize=20,convergenceRounds=8,lookAhead=False):
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
		#lookAhead = int(float(len(adjacencyMatrix)) * lookAhead) + cutIndices[-1]
		lookAhead = int((float(len(adjacencyMatrix) - cutIndices[-1]) * lookAhead) + cutIndices[-1])
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
		
		# Simply pass if too small of submatrix is created
		if len(X[0]) < minSize:
			cutInd = lookAhead
		
		# Otherwise, perform hmm and find new cut index
		else:
		
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
	
def identifyChromosomeGroupsHMM(adjacencyMatrix,binList,minSize=5,modularity=.05,convergenceRounds=5,lookAhead=.2,louvainRounds=20, prev_cutInds=False):
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
	adjacencyMatrix = numpy.asmatrix(adjacencyMatrix)
	matrixLength = float(len(adjacencyMatrix))
	remainder = matrixLength - (modularity*matrixLength)
	modFallBack,cutIndices = 0,[0]
	# Means we want to use modularity instead (presumably, or a mistake)
	if modularity == 1:
		return []
	
	# Deal with using pass in cut indices
	if prev_cutInds != False:
		cutIndices = prev_cutInds
	
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
	
	# Prevents an edge case where hmm can cut at the very end of the matrix, thus providing no data for further modularity calculation
	# or if the previous index is taken and the remaining set of data is massive, then we don't want to run modularity a whole bunch of times.
	# Instead, we can alter the convergenceRounds parameter by subtracting 1 and recursing on the algorithm. Not ideal behavior but as a fall back this works
	if cutIndices[-1] == len(adjacencyMatrix):
		print("- WARNING - Last cut index found to be length of current matrix removing index values of {}".format(cutIndices[-1]))
		cutIndices.pop(-1)
		
		# We don't want to run n rounds of louvain modularity calculation on a large data set. 
		# So this prevents that by setting the maximum default modularity calculation to 5 times the intended fraction
		if (len(adjacencyMatrix) - cutIndices[-1]) >= (5 * (len(adjacencyMatrix[0]) * modularity)):
			
			# Reduce convergence rounds and make sure it is still within bounds
			print("- convergenceRounds reduced from {} --> {}".format(convergenceRounds, convergenceRounds-1))
			if convergenceRounds - 1 == 0:
				print("- Failed to converge after reducing convergence rounds all the way to 1... Returning current indices")
				return cutIndices
			
			# Otherwise, recurse on this function, but reduce our convergenceRounds to minus one and start over at our current index
			else:
				print("- Recursing on identifyChromosomeGroupsHMM function, due to remaining fraction of data being greater than 5x than that of desired fraction")
				cutIndices = identifyChromosomeGroupsHMM(adjacencyMatrix,binList,minSize=5,modularity=.05,convergenceRounds=convergenceRounds-1,lookAhead=.5,louvainRounds=20, prev_cutInds=cutIndices)
	
	print("Total time to identify chromosome boundries = "+str(time.time()-startTime)+" seconds")
	#################
	return cutIndices

### Write out coordinate based results 
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
				hyperGeom,hmm,minSize,modularity,louvainRounds,
				psig,convergenceRounds,lookAhead,resolution):
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
	adjMat,binList = reorderMatrix(numpy.asmatrix(adjMat),binList,dendoLeaves['leaves'])
	plotModule.plotContactMap(adjMat,resolution=resolution,highlightChroms=False,showPlot=False,savePlot=avgClusterPlot)
	stopTime = time.time()
	print("Total run-time to cluster and plot = "+str(stopTime-startTime))
	############################## 
	### Part 2 Hyper geometric ###
	startTime = time.time()
	if hyperGeom == True:
		adjMat = convertMatrix(adjMat,binList,distance=False,similarity=True)
		argsorted_adjMat = numpy.asarray(numpy.argsort(adjMat,axis=1)[:, ::-1])
		initial_cut_inds = pre_process_all_matrix_breakpoints(argsorted_adjMat, min_size=minSize, min_frac=modularity, psig=psig)
		cutIndices = filter_noisy_breakpoints(argsorted_adjMat, initial_cut_inds, psig=psig)
		adjMat = logTransformMatrix(adjMat, logBase=10)
	###################
	### Part 2  HMM ###
	elif hmm == True:
		adjMat = convertMatrix(adjMat,binList,distance=False,similarity=True)
		adjMat = logTransformMatrix(adjMat,logBase=10)
		cutIndices = identifyChromosomeGroupsHMM(adjMat,binList,minSize=minSize,modularity=modularity,convergenceRounds=convergenceRounds,lookAhead=lookAhead,louvainRounds=louvainRounds)

	
	### Note, at this stage in the algorithm we have a log transformed similarity matrix...###
	### This is the DEFAULT matrix that is input into the modularity calculations... ###
	### Scale matters here where if the fraction of the data is very large to compute modularity for, then we might be better off using the non-log transformed similarity matrix...
	### For visual purposes, I like the untransformed distance matrix so this is default but should be easy enough to change
	########################
	###	Part2 modulartiy ###
	if modularity != False:
		if modularity > 0.0:
			adjMat, binList, cutIndices = modularity_remaining_data(adjMat,binList,cutIndices,n_rounds=louvainRounds)
	####################################################################################################################
	### Results of part 2 plotting and writing to file #############################################
	adjMat = convertMatrix(logTransformMatrix(adjMat,logBase=10,reverse=True),binList,distance=True)
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
	print("- Part 1 (grouping bins to groups) completed successfully")
