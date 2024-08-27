####IMPORTS####

import numpy as np
import astropy.io.fits
from astropy import units as u
from astropy.coordinates import SkyCoord

####DEFINITIONS####

#from a skycoord object gives the coordinates in the correct frame
def correctFrame(coords,frame):
	if frame=='galactic_cartesian':
		out=coords.galactic.cartesian
		x=out.x.value
		y=out.y.value
		z=out.z.value
		out=np.array((x,y,z))
		return out
	if frame=='equatorial_cartesian':
		out=coords.equatorial.cartesian
		x=out.x.value
		y=out.y.value
		z=out.z.value
		out=np.array((x,y,z))
		return out

#for given coordinates arrays and given axis of a map, returns a function that retuns the indexes of these coordinates.
def coordsToIndexFastArray(xAxis,yAxis,zAxis):
	xLen=len(xAxis)-1
	yLen=len(yAxis)-1
	zLen=len(zAxis)-1

	xMin,xMax=xAxis[0],xAxis[xLen]
	yMin,yMax=yAxis[0],yAxis[yLen]
	zMin,zMax=zAxis[0],zAxis[zLen]
	
	xSpan=xMax-xMin
	ySpan=yMax-yMin
	zSpan=zMax-zMin
	
	xDelta=(xSpan)/(xLen)
	yDelta=(ySpan)/(yLen)
	zDelta=(zSpan)/(zLen)	
	
	def out(array): #NOTE: the imput for this function is [xs,ys,zs] where xs, ys and zs are all numpy arrays
		xs,ys,zs=array
		
		xIndex=(xs-xMin)/xDelta
		yIndex=(ys-yMin)/yDelta
		zIndex=(zs-zMin)/zDelta

		xIndex=np.rint(xIndex).astype(int)
		yIndex=np.rint(yIndex).astype(int)
		zIndex=np.rint(zIndex).astype(int)
		
		return tuple(np.array((xIndex,yIndex,zIndex)))
	
	return out	

####the map3D class
class map3D:
	def __init__(self,xAxis,yAxis,zAxis,mapName,frame):

		extinctionMap=astropy.io.fits.open(mapName)
		extinctionMap=extinctionMap[0].data.T
		self.extinctionMap=extinctionMap #the actual 3D numpy array
	
		self.indexFunction=coordsToIndexFastArray(xAxis,yAxis,zAxis) #the function that gives the index of where to look 
		
		self.xAxis=xAxis
		self.yAxis=yAxis
		self.zAxis=zAxis
		
		self.mapName=mapName
		
		self.frame=frame #the frame in which the map is situated
	def __repr__(self):
		out='map3D object'
		out+='\nMap used: '+self.mapName
		out+='\n\nDimensions:'
		out+=f'\nx axis: from {min(self.xAxis)} to {max(self.xAxis)} with {len(self.xAxis)} points'
		out+=f'\ny axis: from {min(self.yAxis)} to {max(self.yAxis)} with {len(self.yAxis)} points'
		out+=f'\nz axis: from {min(self.zAxis)} to {max(self.zAxis)} with {len(self.zAxis)} points'	
		out+='\n\nFrame:'+str(self.frame)
		
		return out

#gives the extinction for a given star
def extinction( coords, map3D, output='full', steps=1000, observer=np.array((0,0,0)) ):
	frame=map3D.frame
	indexFunction=map3D.indexFunction
	extinctionMap=map3D.extinctionMap
	
	coordsInCorrectFrame=correctFrame(coords,frame) #getting the coordinates
	path=np.linspace(coordsInCorrectFrame,observer,steps).T #a path from the end to the start

	distance=np.sqrt(sum((coordsInCorrectFrame-observer)**2))	

	try:
		indexes=indexFunction(path) 
		values=extinctionMap[indexes]
		totalV=sum(values)
	except:
		raise IndexError('The object (or the observer) is outside the extinction map you used')
	
	totalV*=distance
	totalV/=steps
	
	totalG=totalV*0.835
	totalBp=totalV*1.139
	totalRp=totalV*0.650		
	
	out={'A_V':totalV,'A_G':totalG,'A_Bp':totalBp,'A_Rp':totalRp,'E(Bp_Rp)':totalBp-totalRp}
	if output=='full':
		return out
	else:
		return out[output]

#from the extinction in the V band, returns the extinction for given wavelength
def extinctionSpectroscopic( A_V , lambd, unit='meters' , fluxOutput = False):
	out=np.array(lambd)
	if unit in ['Å','A','Ångström','Angstrom','angstrom']:
		out = out * 1e-10
	elif unit in ['nm','nanometers','nanometer']:
		out = out * 1e-9	
	
	waveNumberBase=1/551e-9
	
	out = 1 / out
	out =  out / waveNumberBase
	
	if fluxOutput:
		return 10 ** (- out / 2.5 )
	return out * A_V
	
	
####LOADING MAPS AND CUSTOM FUNCTIONS####

#loading the small map
try:
	smallMap=map3D(np.linspace(-1500,1500,601),np.linspace(-1500,1500,601),np.linspace(-400,400,161),'explore_cube_density_values_010pc_v2.fits',frame='galactic_cartesian')	
	small=True
except:
	try:
		smallMap=map3D(np.linspace(-1500,1500,601),np.linspace(-1500,1500,601),np.linspace(-400,400,161),'smallMap.fits',frame='galactic_cartesian')	
		small=True
	except:	
		small=False
if small:
	def extinctionSmallMap( coords, output='full', steps=1000 ):
		return extinction(coords,smallMap,output,steps)

#loading the medium map
try:
	mediumMap=map3D(np.linspace(-3000,3000,601),np.linspace(-3000,3000,601),np.linspace(-400,400,81),'explore_cube_density_values_025pc_v2.fits',frame='galactic_cartesian')	
	medium=True
except:
	try:
		mediumMap=map3D(np.linspace(-3000,3000,601),np.linspace(-3000,3000,601),np.linspace(-400,400,81),'mediumMap.fits',frame='galactic_cartesian')	
		medium=True
	except:
		medium=False
if medium:
	def extinctionMediumMap( coords, output='full', steps=1000 ):
		return extinction(coords,mediumMap,output,steps)

#loading the large map
try:
	largeMap=map3D(np.linspace(-5000,5000,501),np.linspace(-5000,5000,501),np.linspace(-400,400,41),'explore_cube_density_values_050pc_v2.fits',frame='galactic_cartesian')	
	large=True
except:
	try:
		largeMap=map3D(np.linspace(-5000,5000,501),np.linspace(-5000,5000,501),np.linspace(-400,400,41),'largeMap.fits',frame='galactic_cartesian')	
		large=True	
	except:
		large=False
if large:
	def extinctionLargeMap( coords, output='full', steps=1000 ):
		return extinction(coords,largeMap,output,steps)


###HELP FUNCTIONS####

def info():
	out='Welcome to the 3D extinction software!\nCurrent status:\n'
	if small:
		out+='\nsmallMap loaded and ready. Feel free to print it. Use the function extinctionSmallMap to use it\n'
	else:
		out+='\nSmall map not found\n'
	if medium:	
		out+='\nmediumMap loaded and ready. Feel free to print it. Use the function extinctionMediumMap to use it\n'
	else:
		out+='\nMedium map not found\n'
	if large:
		out+='\nlargeMap loaded and ready. Feel free to print it. Use the function extinctionLargeMap to use it\n'
	else:
		out+='\nLarge map not found\n'
	
	out+='\nCredit for these maps: .L. Vergely, Rosine Lallement, and N.J.L. Cox. Three-dimensional extinction maps: Inverting inter-calibrated extinction catalogues. Astronomy Astrophysics, 664, 05 2022.'
	out+='\nYou can also upload and use your own 3D map by defining a map3D class object. Note that you need a .fits file containing the map for this to happen.\n'
	print(out)

def status():
	info()
def help():
	info()
