This directory contains histograms of Av, as well as color-magnitude relations, for redMaGiC galaxies 
divided into 2 redshift bins. The string of numbers in each filename represent the mean redshift of galaxies in that bin.

Files beginning with "av_redmagic" show the histograms of Av, computed in that bin. The "variance" on those figures is actually
the weight of the Av estimate in that tranche, as defined by the extinction MLE.

Files beginning with "redmagic_colormag" are color-magnitude diagrams for galaxies in that redshift tranche. 
The "variance" on those figures is simply np.std(mag1-mag2).
