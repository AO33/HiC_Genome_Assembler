###############
### Imports ###
import time
import gzip ### In case the original fasta file is gzipped
###########


##################################################################################################################################
### Code to read in the original fastaFile, the final ordering/orientations for each scaffold and write out the final assembly ###
def readFastaIntoMem(fastaFile):
	'''This function simple reads a fasta into memory
	   Input = full file path to fasta file
	   Output = { fastaEntryName:sequence }'''
	if ".gz" in fastaFile:
		inFile = gzip.open(fastaFile,mode='rt')
	else:
		inFile = open(fastaFile,'r')
	############
	chrDict = {}
	for line in inFile:
		line = line.strip('\r').strip('\n')
		if line[0] == '>':
			chrName = line[1:]
			chrDict.update({chrName:[]})
		else:
			chrDict[chrName].append(line)
	##############
	inFile.close()
	for c in chrDict:
		chrDict[c] = ''.join(chrDict[c])
	##############
	return chrDict

def readChromosomeOrderingFile(chrOrderFile):
	'''This function takes in the chromosome ordering file and reads it into memory
	   Input = full file path to chromosomeOrdering file (from part2 or part3 of the pipeline)
	   Output = [ [ [scaffold,orientation] , ... ] , ... ] list of scaffolds and orientations per chromosome'''
	orderFile = open(chrOrderFile,'r')
	orderFile.readline()
	groupList,currentGroup = [],[]
	for line in orderFile:
		line = line.strip('\r').strip('\n')
		if line[0] != "#":
			cols = line.split('\t')
			currentGroup.append([cols[0],cols[1]])
		else:
			groupList.append(currentGroup)
			currentGroup = []
	##############################
	groupList.append(currentGroup)
	#################
	orderFile.close()
	return groupList

def reverseTranscribeSeq(seq):
	'''Simple function that takes in a sequence and reverse transcribes it
	   Input = a string
	   Output = the reverse transcript of the input string'''
	oppDict = { "A":"T", "T":"A", "a":"t", "t":"a",
				"G":"C", "C":"G", "g":"c", "c":"g", 
				"N":"N", "n":"n" }
	###################################################
	return ''.join([ oppDict[s] for s in seq[::-1] ])

def writeSeqToFile(fileMan,seq,charsPerLine=50):
	'''Function that writes out a string to a file.
	   Input = an openedFile object, sequence (string) to write out, charsPerLine=how many characters there are per line for the input sequence
	   Output = the openedFile object with the sequence written to it'''
	prevNum = 0
	for num in range(charsPerLine,len(seq)+charsPerLine,charsPerLine):
		fileMan.write(seq[prevNum:num]+'\n')
		prevNum = num
	##############
	return fileMan

def writeNewFasta(chrGroups,oldFastaDict,outFile,charsPerLine=50,nGapLength=100):
	'''Function that writes out a fasta file by taking in an input fasta file and lists of entries to combine together for a new fasta entry.
	   This function will write out all entries in the fasta file that are not included in the chrGroups argument as well.
	   Gaps of Ns of size nGapLength will be introduced inbetween adjoined scaffolds.
	   Input = [ [ [scaffold,orientation] , ... ] , ... ] list of scaffolds and orientations per chromosome,
			   oldFastaDict = { scaffoldName:sequence },outFile = full file path to new fastaFile to be written,
			   charsPerLine = number of characters to write out per line in the fasta file
			   nGapLength = Gaps of Ns of size nGapLength will be introduced inbetween adjoined scaffolds
	   Output = standard out with some basic assembly statistics'''
	### Initiate basic assembly stats ###
	groupedLength,ungroupedLength = 0,0
	newNsWritten,gaps = 0,0
	scaffoldsGrouped,ungrouped = 0,0
	################
	writtenDict = {}
	outFile = open(outFile,'w')
	for i,chromosomeGroup in enumerate(chrGroups,1):
		outFile.write(">Chr_"+str(i)+'\n')
		seqToWrite = []
		for scaffNum,scaffoldInfo in enumerate(chromosomeGroup):
			scaffoldsGrouped += 1
			scaffName,orientation = scaffoldInfo[0],scaffoldInfo[1]
			#if scaffName == "Scpiz6a_72.1": ### Specific to a particular genome assembly I was working on. Not useful except for record keeping for me...
			#    continue
			writtenDict.update({scaffName:''})
			if orientation == "+":
				seqToWrite += [oldFastaDict[scaffName]]
			else:
				seqToWrite += [reverseTranscribeSeq(oldFastaDict[scaffName])]
			if scaffNum != len(chromosomeGroup)-1:
				newNsWritten += nGapLength
				gaps += 1
				seqToWrite += ["N" for n in range(0,nGapLength) ]
		################################
		seqToWrite = ''.join(seqToWrite)
		groupedLength += len(seqToWrite)
		outFile = writeSeqToFile(outFile,seqToWrite,charsPerLine=charsPerLine)
	#########################################################################
	### Write the remainder of scaffolds that were excluded from analysis ###
	for c in oldFastaDict:
		if c not in writtenDict:
			outFile.write(">"+c+'\n')
			seqToWrite = oldFastaDict[c]
			ungroupedLength += len(seqToWrite)
			ungrouped += 1
			outFile = writeSeqToFile(outFile,seqToWrite,charsPerLine=charsPerLine)
	###############
	outFile.close()
	print("Total scaffolds grouped into chromosomes"+'\t'+str(scaffoldsGrouped))
	print("Total genome length grouped into chromosomes"+'\t'+str(groupedLength-newNsWritten))
	print("Total new gaps introduced"+'\t'+str(gaps))
	print("Total ungrouped scaffolds"+'\t'+str(ungrouped))
	print("Total genome length ungrouped "+'\t'+str(ungroupedLength))


########################################################
### Function that runs the entire pipeline for Part4 ###
def runPipeline(originalFastaFile,finalOrderingFile,assembledFastaFile):
	#################################################
	print("########################################")
	print("### Working on Part4 of the pipeline ###")
	startTime = time.time()
	fastaDict = readFastaIntoMem(originalFastaFile)
	chrGroups = readChromosomeOrderingFile(finalOrderingFile)
	writeNewFasta(chrGroups,fastaDict,assembledFastaFile,charsPerLine=50,nGapLength=100)
	print("Total run-time  for Part4 = "+str(time.time()-startTime))
	print("- Part 4 (writing of new super-scaffolded genome .fasta) completed successfully")
