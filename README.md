# 3D-extinciton-map-software-for-astronomy
3-Dimensional maps of space dust in the galaxy can be used to get the extinction and reddening of stars and other objects with unprecedented accuracy. This code caclulates it.

The software is initially meant for using the public 3-Dimensional dust maps available at 
https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/664/A174#/article

They were published in the following paper:
J.L. Vergely, Rosine Lallement, and N.J.L. Cox. Three-dimensional extinction maps: Inverting inter-calibrated extinction
catalogues. Astronomy Astrophysics, 664, 05 2022.

The maps are in the unit of magnitudes/parsec. To get the extinciton of a star, it is necessary to integrate over the line of sight from the star to the observer. Further updates may allow for customisation with your own maps.

# USAGE:
## Dowload the extinction map
Go to 
https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=J/A%2bA/664/A174/list
Submit the query.
Dowload your desired map(s). If you wish, you can rename them as following to make it easier to understand what they are.

**smallMap.fits**: the file called explore_cube_density_values_010pc_v2.fits

**mediumMap.fits**: the file called explore_cube_density_values_025pc_v2.fits

**largeMap.fits**: the file called explore_cube_density_values_050pc_v2.fits

Note that the small map is also the one of with the _largest file size_. This is due to the fact that while the small map describes a smaller region of space, it has a higher resolutions, leading to a larger total nunber of pixels.

## Setup
Start by having the extinction3D.py file and an extinction map in the same directory as your code. If you use one of the pre-determined maps, the maps will will be loeaded automatically. You will also need a SkyCoord object containing the 3-Dimensional position of your object. Note that three values are required to place the star in 3-Dimensional space. Example of what would work here is right ascension, declination and distance.

## Execution
From the astropy coordinates object coords, the formatting to get the extinciton is as follows:
## extinction( coords, 3Dmap, output='full', steps=100, observer=(0,0,0) )
coords is a SkyCoord class object.

3Dmap is the map to be used when calculating extinciton.

output='full' denotes for what observational bands the output will be given. So far, extinction in four bands is avaliable: visual V 'A_V', gaia G band 'A_G', gaia Bp band 'A_Bp', and gaia Rp band 'A_Rp'. Also given is the reddening 'E(Bp_Rp)' = A_Bp - A_Rp . If the standard string 'full' is passed, all of these are returned as a python dictionary.

steps=1000 is the number of integration steps

observer=(0,0,0) is the coordinates of the observer from which photometry of a star was taken. For many currently available extinction maps and astrophysics research, the maps are centered at Earth, the same place where the telescopes are situated.

This function uses numerical integration over the map to get the extinction, and is as such computationally expensive. If desired, the following functions can be used instead to simplify the process:

## extinctionSmallMap( coords, output=full, steps=1000 )
This function uses a 3kpc x 3kpc x 0.8 kpc extinction map and is centered at Earth, where it is assumed that the observer also is. This function is equal to **extinction(coords,smallMap)**.
## extinctionMediumMap( coords, output=full, steps=1000 )
This function uses a 6kpc x 6kpc x 0.8 kpc extinction map and is centered at Earth, where it is assumed that the observer also is. This function is equal to **extinction(coords,mediumMap)**.
## extinctionLargeMap( coords, output=full, steps=1000 )
This function uses a 10kpc x 10kpc x 0.8 kpc extinction map and is centered at Earth, where it is assumed that the observer also is. This function is equal to **extinction(coords,largeMap)**.

The advantages of the larger maps is that more stars are contained within them. The disadvantage is that precision is reduced and error is increased.

To get values of extinction for spectroscopic observation, the following function can be used:
## extinctionSpectroscopic( A_V, lambd, unit='meters', fluxOutput = False ):
A_V is the extinction in the V band, obtained i.e. by using the **extinction function.

lambd is the wavelenght for which extinction is desired. Can be a single value or can be a list (or array) of wavelengths.

unit = 'meters' is the unit in which the wavelength is given. Currently supported units are Ångström and nanometers. If no unit is given, it is assumed that the wavelength is given in meters.

fluxOutput = False denotes wheter the output should be given in magnitudes of extinction, or, if set as True, in the fraction of the light that reaches the observer. For example, a flux output of 0.3 means that 30% of the light from a star is not scattered by space dust and reaches the observer.

## info()
Returns the current status of the script, showing which maps are loaded.
