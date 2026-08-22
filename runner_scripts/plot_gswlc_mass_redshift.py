import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = False
plt.rcParams['font.family'] = 'sans-serif'
import numpy as np
from astropy.table import Table

print('Loading GSWLC catalog...')
cat = Table.read('/projects/mccleary_group/dusty_halos/catalogs/GSWLC-X2_in_SDSS_z_lt_0.18.fits')
good = (cat['log_M_star'] > -99) & (cat['log_SFRsed'] > -99)
t = cat[good]

ssfr = t['log_SFRsed'] - t['log_M_star']

idx = np.random.choice(len(t), size=80000, replace=False)
z = t['z'][idx]
mass = t['log_M_star'][idx]
color = ssfr[idx]

fig, ax = plt.subplots(figsize=(11, 8), tight_layout=True)

sc = ax.scatter(z, mass, c=color, s=1, alpha=0.4, cmap='RdYlBu', vmin=-13, vmax=-9)

cbar = plt.colorbar(sc, ax=ax)
cbar.set_label(r'log sSFR (yr$^{-1}$)', fontsize=13)

# Mass split lines
ax.axhline(10.7, color='black', linestyle='--', linewidth=1.2, label='log M = 10.7')
ax.axhline(11.3, color='purple', linestyle='--', linewidth=1.2, label='log M = 11.3 (LRG)')
ax.axhline(9.5,  color='gray',   linestyle='--', linewidth=1.0, label='log M = 9.5 (dwarfs)')

ax.set_xlabel('Redshift', fontsize=14)
ax.set_ylabel(r'log $M_*$ ($M_\odot$)', fontsize=14)
ax.set_title('GSWLC — Stellar Mass vs Redshift\ncolored by log sSFR', fontsize=14)
ax.set_xlim(0, 0.18)
ax.set_ylim(7, 13)
ax.legend(fontsize=11, loc='upper left')

fig.savefig('/projects/mccleary_group/hsia.i/dusthalos_output/gswlc_sdss/gswlc_mass_redshift.png', dpi=150)
print('Done!')
