###############################################
#
# file utility functions
# 
# some useful file functions for fastq/a
# processing
#
###############################################

def grouped(iterator, size):
	""" 
	return tuple of a certain size
	
	"""
	yield tuple(next(iterator) for _ in range(size))
	
def grabLine(file, num):
	""" 
	return a list containing num lines from
	a file. 
	
	file = file to grab lines from
	num = number of lines to grab
	"""
	fileList = []
	append = fileList.append
	for i in num:
		line = str(file[i]).strip('\n')
		append(line.strip('\n'))
	return(fileList)
	
def addLine(finalList, filename):
	"""
	Add lines to a file
	
	finalList = list of lines to add
	filename = file to add to
	"""
	with open(filename, "a") as f:
  	 	f.writelines(finalList)
  	 	
def headerRep(fastqFile, outFile):
	""" 
	Replace header of a fastq file with a simpler
	name for compatibility with some tools
	
	@seq1 .... @seqN
	
	fastqFile = input fastq file to clean
	outFile = output file name to write
	"""
	
	# get number of lines first
	with open(fastqFile) as f:
		size=sum(1 for _ in f)
	f.close()
	
	# replace header
	with open(fastqFile) as f:
		fastq = f.readlines()
		lineNum = (0,1,2,3)
		countSeq=0
		while lineNum[0] < size:
			seqHeader, seq, extraName, seqQuality = grabLine(fastq, lineNum)
			countSeq += 1			
			finalList = ['@seq'+str(countSeq)+'\n', seq+'\n', extraName+'\n', seqQuality+'\n']
			addLine(finalList, outFile)	
			lineNum = (lineNum[0] + 4, lineNum[1] + 4, lineNum[2] + 4, lineNum[3] + 4)		 	
	f.close()	
