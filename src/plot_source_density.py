import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

hidens = Table.read('redmagic_hidens_y3_GOLD_JOINED_catalog.fits', memmap=True)

ra = hidens['ra_redmagic_hidens']

dec = hidens['dec_redmagic_hidens']


# Let's assume you have arrays 'ra' and 'dec' containing the RA and Dec of your sources in radians.
# First, convert them to degrees for easier interpretation.
ra_deg = np.rad2deg(ra)
dec_deg = np.rad2deg(dec)

# Define the size of your bins in arcminutes.
bin_size_arcmin = 100.0  # for example, you can adjust this to your preferred bin size
bin_size_deg = bin_size_arcmin / 60.0  # convert arcminute to degree

# Calculate the number of bins along each axis.
num_bins_ra = int(np.ptp(ra_deg) / bin_size_deg)
num_bins_dec = int(np.ptp(dec_deg) / bin_size_deg)

# Create a 2D histogram of the source counts.
hist, ra_edges, dec_edges = np.histogram2d(ra_deg, dec_deg, bins=[num_bins_ra, num_bins_dec])

# Calculate the area of each bin in square degrees.
area_per_bin = bin_size_deg ** 2

# Convert the histogram to sources per square arcminute.
density_map = hist / area_per_bin / 3600.0  # because there are 3600 square arcminutes in a square degree

# Now you can plot the density map.
fig, ax = plt.subplots(1, 1, figsize=(8, 8))

im = ax.imshow(
    density_map.T, origin='lower', aspect='auto'
    extent=[ra_edges[0], ra_edges[-1], dec_edges[0], dec_edges[-1]]
)
cax = divider.append_axes("right", size="5%", pad=0.07)
fig.colorbar(im, cax=cax, label='Number of sources per square arcminute')
ax.set_xlabel('RA (degrees)')
ax.set_ylabel('Dec (degrees)')
ax.set_title('Source Density Map')
fig.savefig('test_src_density.png')
