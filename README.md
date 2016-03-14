We would like to put together a suite of functions that helps us
understand where the imaged slice from any given day of calcium
imaging falls as a manifold within the 3 dimensional stack acquired at
the beginning of the experiment. We expect that this surface will
be almost perflectly aligned to one of the slices toward the middle of
the volume at the beginning of the experiment, but that during the
collection of each dataset and over the days of the experiment the
shape of the brain will change. In order to determine which neurons we
are looking at, we'll need to estimate their most likely x,y, and z 
coordinates within the stack.

The goal is to produce a matrix C for each day of the experiment that
describes the normalized 2d cross correlation between an image from
that day for a range of z, x, y and rotation values.
