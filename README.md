# 3D-extinciton-map-software-for-astronomy
3-Dimensional map of space dust in the galaxy can be used to get the total darkening this dust causes with unprecedented accuracy. This code caclulates the total extinction and magnitude.

14/05/2024.
The software is initially meant for using the public 3-Dimensional dust maps available at 
https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/664/A174#/article

They were published in the following paper:
J.L. Vergely, Rosine Lallement, and N.J.L. Cox. Three-dimensional extinction maps: Inverting inter-calibrated extinction
catalogues. Astronomy Astrophysics, 664, 05 2022.

The maps are in the unit of magnitudes/parsec. To get the extinciton of a star, it is necessary to integrate over the line of sight from the star to the observer. Further updates may allow for customisation with your own maps.

# USAGE:
To get the extinction of a star or other object, the program requires, at first, the 3-Dimensional positon of the star. This is to be given by you in the form of astropy coordinates. Note that three values are required to place the star in 3-Dimensional space. Example of what would work here is right ascension, declination and parallax.

From the astropy coordinates object coords, the formatting to get the extinciton is as follows:
## out = extinction( coords, 3Dmap, output='full', steps=100, observer=(0,0,0) )
3Dmap=map10kpcX10kpcX800pc is the map to be used when calculating extinciton. The standard is to use a 10kpc x 10kpc x 0.8kpc map.

output='full' denotes for what observational bands the output will be given. So far, extinction in four bands is avaliable: visual V ('A_V'), gaia G band 'A_G', gaia Bp band 'A_Bp', and gaia Rp band 'A_Rp'. Also given is the reddening 'E(Bp_Rp)' = A_Bp - A_Rp . If the standard string 'full' is passed, all of these are returned as a python dictionary.

steps=100 is the number of integration steps

observer=(0,0,0) is the coordinates of the observer from which photometry of a star was taken. For many currently available extinction maps and astrophysics research, the maps are centered at Earth, the same place where the telescopes are situated.

To help simplify the project, the following streamlined functions have been defined:
## extinctionSmallMap( coords, output=full, steps=100 )
This function uses a 3kpc x 3kpc x 0.8 kpc extinction map and is centered at Earth, where it is assumed that the observer also is. This function is equal to **extinction(coords,smallMap)**.
## extinctionMediumMap( coords, output=full, steps=100 )
This function uses a 6kpc x 6kpc x 0.8 kpc extinction map and is centered at Earth, where it is assumed that the observer also is. This function is equal to **extinction(coords,mediumMap)**.
## extinctionLargeMap( coords, output=full, steps=100 )
This function uses a 10kpc x 10kpc x 0.8 kpc extinction map and is centered at Earth, where it is assumed that the observer also is. This function is equal to **extinction(coords,largeMap)**.

The advantages of the larger maps is that more stars are contained within them. The disadvantage is that precision is reduced and error is increased.


# POSSIBLE ERRORS:

Apart from errors given by in-built functions of the package, custom errors that can be received are:
## OUT_OF_BOUNDS_ERROR:
The position of the star (or of the observer) is outside of the extinction map used. The only solution to this problem is finding a new map that contains both.

