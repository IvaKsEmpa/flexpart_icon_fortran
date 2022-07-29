require(ncdf4)
require(myRplots)


#Converts geographical to cartesian coordinates.
#
#@par Revision History
#Developed  by Luis Kornblueh  (2004).
#
gc2cc = function(p_pos){
  
	z_sln = sin(p_pos$lon)
	z_cln = cos(p_pos$lon)
	z_slt = sin(p_pos$lat)
	z_clt = cos(p_pos$lat)

	p = data.frame( 
		x = z_cln*z_clt,
		y = z_sln*z_clt,
		z = z_slt)

	return(p)
}

#	p, v1, v2, v3: points cartesian (fields x,y,z)
inside.triangle = function(p, v1, v2, v3){
  
      c1  = ccw.spherical(v1, v2, p)
      c2  = ccw.spherical(v2, v3, p)
      c3  = ccw.spherical(v3, v1, p)

      inside.triangle = ((  c1) && (  c2) && (  c3))  ||
                        (( !c1) && ( !c2) && ( !c3))

	return(inside.triangle)
}

ccw.spherical.q128 = function(v1,v2,v3){
	require(Rmpfr)
	
	vh1 = mpfr(c(v1$x, v1$y, v1$z),128)
	vh2 = mpfr(c(v2$x, v2$y, v2$z),128)
	vh3 = mpfr(c(v3$x, v3$y, v3$z),128)

	ccw = vh3[1]*(vh1[2]*vh2[3] - vh2[2]*vh1[3]) -
          vh3[2]*(vh1[1]*vh2[3] - vh2[1]*vh1[3]) +
          vh3[3]*(vh1[1]*vh2[2] - vh2[1]*vh1[2])

	return(ccw <= 0.)
}

ccw.spherical = function(v1,v2,v3) {
  
   # where a is the angle between v3 and the normal to the plane
   # defined by v1 and v2.
  
	ccw =           v3$x*(v1$y*v2$z - v2$y*v1$z) -
                    v3$y*(v1$x*v2$z - v2$x*v1$z) +
                    v3$z*(v1$x*v2$y - v2$x*v1$y)
  
      #! we apply a static error of
      #!   | e - e'| <= 3*2^-48
      #! to decide if a floating-point evaluation e' of an expression e
      #! has the correct sign, see Section 2.2 of
      #!
      #! Burnikel, C.; Funke, S. & Seel, M. 
      #! "Exact geometric computation using Cascading"
      #! International Journal of Computational Geometry & Applications, 
      #! World Scientific, 2001, 11, 245-266
	if (abs(ccw) <= 1.1e-14){
		ccw.spherical = ccw.spherical.q128(v1,v2,v3)
	} else {
		ccw.spherical = ( ccw <= 0.)
	}
	return(ccw.spherical)
}


#	barycentric weights / coordinate transformation 
#	All inputs in cartesian coordinates
bic.weights = function(p, v1, v2, v3){
	a2 = v1$x * (v2$y - v3$y) + v2$x * (v3$y - v1$y) + v3$x * (v1$y - v2$y)	

	w = vector("numeric", 3)
	w[1] = 1/a2 * (
		(v2$x*v3$y - v3$x*v2$y) + 
		(v2$y - v3$y) * p$x + 
		(v3$x - v2$x) * p$y )
	w[2] = 1/a2 * (
		(v3$x*v1$y - v1$x*v3$y) + 
		(v3$y - v1$y) * p$x + 
		(v1$x - v3$x) * p$y )
	w[3] = 1/a2 * (
		(v1$x*v2$y - v2$x*v1$y) + 
		(v1$y - v2$y) * p$x + 
		(v2$x - v1$x) * p$y )

	return(w)	
}


####################################################################
#	SETTINGS
####################################################################

#	icon grid (any icon grid should work, but routine does not zoom in, but
#		plots whole grid, which is increasingly slow for large grids)
grid.fn = "coarse_DOM01.nc"
grid.fn = "fine_DOM01.nc"
grid.fn = "/input/ICON/ivme/grid/ICON-1E_DOM01_tri.nc"

plot.grid = FALSE
use.random = TRUE
p.longlat = data.frame(lon=8, lat=47)

####################################################################
#	MAIN 
####################################################################

#	read grid definition data
#################################
cn = nc_open(grid.fn)

#	locations of grid centers
lon = ncvar_get(cn, "lon_cell_centre") 
lat = ncvar_get(cn, "lat_cell_centre")
cell = data.frame(lon=lon, lat=lat)

#	locations of grid verticies
lon = ncvar_get(cn, "vlon") 
lat = ncvar_get(cn, "vlat") 
vertex = data.frame(lon=lon, lat=lat)

#	for each cell give the 3 vertex indices
vertex.of.cell = ncvar_get(cn, "vertex_of_cell")

#	for each triangle (between grid centers) give the cell indices
cell.of.tri = ncvar_get(cn, "cc_delaunay")
nn.tri = dim(cell.of.tri)[1]

nc_close(cn)


#	vertex and cell locations to cartesian
####################################################
vertex.cc = gc2cc(vertex)
cell.cc = gc2cc(cell)


if (plot.grid){
#	start a two panel plot (long/lat on the left; cartesian on the right)
##################################3
par(mfcol=c(1,2))

#	plot icon cells/centers + triangulation (all in long/lat system)
####################################################
plot.new()
plot.window(range(vertex$lon)*180/pi, range(vertex$lat)*180/pi, asp=1)
axis(1)
axis(2)
box()
title(xlab="Longitude (°E)", ylab="Latitude (°N)")
nn.cells = dim(vertex.of.cell)[1]
blues = brewer.ramp(nn.cells, "Blues", col0=NA)
for (ii in 1:nn.cells){
	lons = vertex$lon[vertex.of.cell[ii,]]*180/pi
	lats = vertex$lat[vertex.of.cell[ii,]]*180/pi 
#	polygon(lons, lats, col=blues[ii])
	lines(lons[c(1,2,3,1)], lats[c(1,2,3,1)], col=1, lwd=2)
}

for (ii in 1:dim(cell.of.tri)[1]){
	lons = cell$lon[cell.of.tri[ii,]]*180/pi
	lats = cell$lat[cell.of.tri[ii,]]*180/pi
	lines(lons[c(1:3,1)], lats[c(1:3,1)], col=2, lty=1)
}
points(cell$lon*180/pi, cell$lat*180/pi, pch=19, cex=1.4)
map(add=TRUE)

#	plot icon cells/centers + triangulation (all in cartesian system)
####################################################
plot.new()
plot.window(range(vertex.cc$x), range(vertex.cc$y), asp=1)
axis(1)
axis(2)
box()
title(xlab="x", ylab="y")

nn.cells = dim(vertex.of.cell)[1]
blues = brewer.ramp(nn.cells, "Blues", col0=NA)
for (ii in 1:nn.cells){
	lons = vertex.cc$x[vertex.of.cell[ii,]] 
	lats = vertex.cc$y[vertex.of.cell[ii,]] 
#	polygon(lons, lats, col=blues[ii])
	lines(lons[c(1,2,3,1)], lats[c(1,2,3,1)], col=1, lwd=2)
}

for (ii in 1:dim(cell.of.tri)[1]){
	lons = cell.cc$x[cell.of.tri[ii,]] 
	lats = cell.cc$y[cell.of.tri[ii,]] 
	lines(lons[c(1:3,1)], lats[c(1:3,1)], col=2, lty=1)
}
points(cell.cc$x, cell.cc$y, pch=19, cex=1)

}


#	select a random point in plot
####################################################
if (use.random){
	p = data.frame(
		lon = runif(1, min=min(cell$lon), max=max(cell$lon)),
		lat = runif(1, min=min(cell$lat), max=max(cell$lat)))
} else {
	p = p.longlat/180*pi
}
#	convert to cartesian
p.cc = gc2cc(p)

if (plot.grid){
points(p.cc$x, p.cc$y, col=1, pch=4, cex=3, lwd=2*par("lwd"))
}

#	find ICON grid containing point (simple loop)
####################################################
for (ii in 1:nn.cells){
	is.in = inside.triangle(p.cc, 
		vertex.cc[vertex.of.cell[ii,1],], 
		vertex.cc[vertex.of.cell[ii,2],], 
		vertex.cc[vertex.of.cell[ii,3],] )
	if (is.in) break
}

if (plot.grid){
#	highlight ICON cell
if (is.in){
	lons = vertex.cc$x[vertex.of.cell[ii,c(1:3,1)]] 
	lats = vertex.cc$y[vertex.of.cell[ii,c(1:3,1)]] 
	lines(lons, lats, col="gray", lwd=3)
	points(p.cc$x, p.cc$y, col="blue", pch=4, cex=3, lwd=2*par("lwd"))
}
}

#	find triangulation triangle containing point (simple loop)
####################################################
for (ii in 1:nn.tri){
	is.in = inside.triangle(p.cc, 
		cell.cc[cell.of.tri[ii,1],], 
		cell.cc[cell.of.tri[ii,2],], 
		cell.cc[cell.of.tri[ii,3],] )
	if (is.in) break
}

#	highlight triangle
if (plot.grid){
if (is.in){
	lons = cell.cc$x[cell.of.tri[ii,c(1:3,1)]] 
	lats = cell.cc$y[cell.of.tri[ii,c(1:3,1)]] 
	lines(lons, lats, col="red", lwd=3)
	points(p.cc$x, p.cc$y, col="blue", pch=4, cex=3, lwd=2*par("lwd"))
}
}

#	calculate weights for barycentric interpolation
############################################################
if (is.in){
	w = bic.weights(p.cc, v1=cell.cc[cell.of.tri[ii,1],], 
		v2=cell.cc[cell.of.tri[ii,2],], 
		v3=cell.cc[cell.of.tri[ii,3],])

	if (plot.grid){
		text(lons[1:3], lats[1:3], signif(w,3), col="blue", font=2, cex=1.2)
	}
	#	check if sum of weights is equal to 1
	print(sum(w)-1)
	cat("w", w)
}

