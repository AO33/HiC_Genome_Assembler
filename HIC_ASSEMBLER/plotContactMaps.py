###############
### Imports ###
import time
import numpy
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg') ### This is to suprress the Xwindows backend. Won't work on the IMB compute server otherwise
import xarray
#############


####################################################################################################################
### Plot contact map (basically a glorified heatMap where we can outline chromosomes given a set of coordinates) ###
### Default coordinates are reported in terms of MegaBases (Mb) on the axis but this is simple enough to change ####
def plotContactMap(adjMat,resolution=100000,tickCount=11,highlightChroms=False,wInches=32,hInches=32,lP=1,hP=98,reverseColorMap='_r',showPlot=False,savePlot=False,title=False,titleSuffix=False):
	print("- Attempting to plot matrix of shape {}".format(adjMat.shape))
	startTime = time.time()
	matrixLength = len(adjMat)
	fig,ax = plt.subplots()
	fig.set_size_inches(wInches,hInches)
	xarray.plot.pcolormesh(xarray.DataArray(adjMat[::-1]),
						 ax=ax,
						 add_colorbar=False,
						 cmap="plasma"+reverseColorMap, #plasma_r for the reverse which is best for distance matricies. Normal or '' is best for similarity matrices
						 robust=True,
						 vmin=numpy.percentile(adjMat,lP),
						 vmax=numpy.percentile(adjMat,hP))
	##########################
	if highlightChroms!=False:
		prevIndex = 0
		for index in highlightChroms:
			###################################################################################################
			### Top of chromosome horizontalLine ##############################################################
			ax.plot([prevIndex,index],[matrixLength-prevIndex,matrixLength-prevIndex],color='white',figure=fig)
			ax.plot([prevIndex,index],[matrixLength-prevIndex,matrixLength-prevIndex],color='white',figure=fig)
			###########################################################################################
			### Bottom of chromosome horizontalLine ###################################################
			ax.plot([prevIndex,index],[matrixLength-index,matrixLength-index],color='white',figure=fig)
			ax.plot([prevIndex,index],[matrixLength-index,matrixLength-index],color='white',figure=fig)
			###################################################################################################
			### Left of chromosome verticallLine ##############################################################
			ax.plot([prevIndex,prevIndex],[matrixLength-index,matrixLength-prevIndex],color='white',figure=fig)
			ax.plot([prevIndex,prevIndex],[matrixLength-index,matrixLength-prevIndex],color='white',figure=fig)
			###########################################################################################
			### Right of chromosome verticalLine ######################################################
			ax.plot([index,index],[matrixLength-index,matrixLength-prevIndex],color='white',figure=fig)
			ax.plot([index,index],[matrixLength-index,matrixLength-prevIndex],color='white',figure=fig)
			#################
			prevIndex = index
		##########################################################################################################
		### For the very last chromosome highlight, Top of chromosome horizontalLine #############################
		ax.plot([prevIndex,matrixLength],[matrixLength-prevIndex,matrixLength-prevIndex],color='white',figure=fig)
		ax.plot([prevIndex,matrixLength],[matrixLength-prevIndex,matrixLength-prevIndex],color='white',figure=fig)
		##################################################################################
		### Left of chromosome verticallLine #############################################
		ax.plot([prevIndex,prevIndex],[0,matrixLength-prevIndex],color='white',figure=fig)
		ax.plot([prevIndex,prevIndex],[0,matrixLength-prevIndex],color='white',figure=fig)
	##################################
	tickDist = len(adjMat) / tickCount
	xTicks,tickMan = [0],0
	for i in range(0,tickCount-1):
		tickMan += tickDist
		xTicks.append(tickMan)
	##########################
	xTicks.append(len(adjMat))
	ax.set_xticks(xTicks)
	ax.set_xticklabels([ str(int((t*resolution)/1000000))+" Mb" for t in xTicks ],size=18)
	ax.set_xlabel('')
	xTicks.pop(0)
	ax.set_yticks(xTicks)
	ax.set_yticklabels([ str(int((t*resolution)/1000000))+" Mb" for t in xTicks ],size=18)
	ax.set_ylabel('')
	##################
	if title != False:
		if titleSuffix != False:
			title = title+titleSuffix
		ax.set_title(title,size=25)
	#####################
	if savePlot != False:
		plt.savefig(savePlot)
	#####################
	if showPlot != False:
		plt.show()
	plt.close(fig)
	#####################
	stopTime = time.time()
	print("Time to rearrange matrix and plot "+str(stopTime-startTime))

