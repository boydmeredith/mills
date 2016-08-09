*Full overview of the pipeline*
———————————————
1. Register each frame of the motion-corrected movie to reference z-stack by finding the x, y, z position and rotation angle (r) that maximize correlations between (36) partially overlapping blocks from the movie and the slices of the stack.
2a. Divide movie into sub-movies—series of frames with correlated intensities and x, y, z, r estimates
2b. Run cell finding on these sub-movies.
3a. Map found cells back into reference space.
3b. Determine whether rois found in each sub-movie came from the same cell or not based on position in reference.


*What does mills do ?*
Most of the mills package exists to perform #1 and then to visualize the results. The main function is blockwiseMovieStackCorr.m which was extensive documentation for its many optional arguments. There are also many functions which operate on the output localization summary file that generate useful visualizations and quality controls, including to help generate figures to guide manual definition of the of sub-movie divisions for #2a (although it would be nice to write a gui that gives the user feedback on this). Finally, a function called registerRois.m can accomplish #3a provided 2b has been accomplished by separate code. 

______________________________________________________________________


#1
=======
MAIN FUNCTION 
blockwiseMovieStackCorr - register movies to subject-specific reference stack 
				    (see bottom of page for more detailed description of algorithm)
IN: 		subject name
		movie
		optional arguments
OUT: 		xyzrcoPeak (best x,y,z,r values and their corresponding correlation values and oddball status)
         		params (parameters used to get this registration)
SAVED: 	summary.mat file with xyzrcoPeak and other info (see intaglio doc). 
		frameXXX/blockXXX.mat files containing sparse formatted correlation values and x,y,z,r indices 


VISUALIZATION TOOLS:
createPairsPlot - plot xy, xz, yz, xr, yr, zr correlation peaks for a given block
blocksByTimeHeatPlots - plot heat maps of x,y,z,r,c,o values
blocksByTimeLinePlots - plot line plots of x,y,z,r,c,o values. has useful optional args for mean centering, plotting mean lines, etc
blockRectsGifWrapper - create a gif of the block position in 2d space - helpful for seeing rotation
ballStickGifWrapper - create a gif of block position in 3d space 
blockComparisonGifs  - create gifs of block-stack montage, false color overlay, and 

unmaintained:
plotDriftAllSessions, inspectAllDatasets - these functions do something related to showing registration information related to all subjects (or all datasets for one subject) in one big plot
plotBlockShifts - at some point this allowed us to plot a vector field of block position in movie relative to stack; might still work?


#2a - useful for breaking up clusters
========
blocksByTimeHeatPlots, blocksByTimeLinePlots - see above visualization tools section

#3b
========
MAIN FUNCTION
registerRois - maps found cells from sub-movies into subject-specific reference stack space by using the peaks found in #1 to guess where the block should line up for these time points. It localizes the mean image for that block/sub-movie in the stack space (similar to above) and then uses that localization to put the rois in the stack space.
IN: 		subject name
		movie date
		location
		optional arguments		
OUT:		rois 	-  struct array of size nClusters, nBlocks with fields w (roi weights), x, y, z (position in stack). rotation is 				applied to w, and w is padded with nan
		xyzrcoClusterPeaks - array of size 6 x nClusters x nBlocks containing localization of each mean image from 				blocks in time clusters  
		
SAVED:	nothing yet!

VISUALIZATION TOOLS:
testRegRoiInVol - a script that will plot a stack slice and clusters from different dates. not a stand alone function yet

unmaintained:
The following functions have not been updated to work with current roi format, but might be helpful:
plotRegRoiInVol, plotRoisInVol, plotRegisteredRoisInSection
______________________________________________________________________

RUNNING THIS STUFF ON THE CLUSTER
=================================
findAndReportPeaks.sh  <subject> <movieDate> will run the localization and produce many of the desired graphics including:
block comparison gifs (montage, overlay, difference), the ball and stick gif, and the rotating rects gif. These will all be saved in 
the reference localization folder.

______________________________________________________________________



A more detailed description of the algorithm used in blockwiseMovieStackCorr

For frame 1:

search widely to narrow down search range
	find best XYZ for R = 0 
	optionally: use xy values from another date and use search range specific leash (xMargin, yMargin; default 50 pix)
	figure showing maximum correlation vs z in frame001/zFit.pdf 
	find XY outliers using robust planar fit, use inliers to fit XYZ
		outliers are defined as points with residuals greater than 8*robustSTD
	find best R in (-10:.5:10) for Z = best fit Z, XY = best fit XY ± 10 pix
	figure showing maximum correlation vs angle in frame001/rFit.pdf 
	fit XY with robust planar fit
	use inliers from XY fit to fit R with a plane, fit Z with local piecewise function (loess)
	figure showing fits and outliers in searchRangeFig.pdf
	save matrix of fits for future use in xyzrSearchRange.mat
	
use XYZR fit values to define center point of search range
	R -2 : 0.25 : 2 around fit
	Z -6 : 6 around fit
	XY search within 20 pixels
		XY location must have 80% overlap with reference
within search range, compute all correlation values, find best
	keep at most 400 positive correlation values 
	save kept values in frame001/bockXXX.mat
	difference image showing block - reference match in blockDiffs.gif 
	side by side comparison of the block and reference match in blockMontage.gif


for frame N+1:

fit XYZR values from frame N
	same fitting as above
	outliers are noted
use XYZR fit values to define center point of search range for frame N+1 as above
compute all correlation values within search range and find best as above 
computes in ~8 minutes (about 9.5 hrs for 70 frames)



to run everything at once on the cluster - use findAndReportPeaks.sh or findAndReportPeaksWrapper.sh

after computing:
	run xyzrcBallStickFig to get 3d scatter of block registration for a given frame, or run warpGif to get a movie
	run createPairPlots to see 6 plots of correlation peaks 
	run fitXYZRSearchRange with a frame’s xyzr peak as input to get the search range for the following frame

