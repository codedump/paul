#!/usr/bin/python

from matplotlib.collections import LineCollection
from matplotlib.transforms import Affine2D
import matplotlib.colors

import itertools
import numpy as np

import paul.toolbox.arplot as ap

# BZ projection along 001, as a (X,Y) line; repeat by kx*ky.
uru2si2_001 = [
 ( .9004750999, -.6241998107),
 ( .9004750999,  .6241998107),
 ( .6241998107,  .9004750999),
 (-.6241998107,  .9004750999),
 (-.9004750999,  .6241998107),
 (-.9004750999, -.6241998107),
 (-.6241998107, -.9004750999),
 ( .6241998107, -.9004750999),
 ( .9004750999, -.6241998107),
]

# kx, ky and ky vectors in 3D; when using,
# ignore dimensions as required
uru2si2_kvec = [
  (2*1.5246749107, 0),
  (1.5246749107, 1.5246749107)
]

# distance: 1.5246749107

# Paul: plot parameters begin
plot_params = {'rot': 45, # 57,
               'repeat': 1,
               'xlim': (-0.3, 1.1),
               'ylim': (-0.5, 0.5),
               'cmap_lim': 0.1,
               }
# Paul: plot parameters end

def populate (*args, **kwargs):

    p = plot_params

    kwargs['can'].reset()
    wav = kwargs['wav']
    ax = kwargs['can'].axes

    # custom color map:
    # we want blue for the lower 10% of intensity,
    # and red/green hues for the rest
    lim = p['cmap_lim']
    cmap_nodes = {'red':   [(0.0, 0.0, 0.0),
                            (lim, 0.0, 0.0),
                            (1.0, 1.0, 0.0),
                            ],
                  'green': [(0.0, 0.0, 0.0),
                            (lim, 1.0, 0.7),
                            (1.0, 0.2, 0.0),
                            ],
                  'blue':  [(0.0, 0.0, 0.0),
                            (lim, 1.0, 0.8),
                            (1.0, 0.1, 0.0),
                            ],
                  }
    cmap = matplotlib.colors.LinearSegmentedColormap ('URS-BZ', cmap_nodes)
    ax.imshow (wav[0], extent=wav[0].imlim, interpolation='none', cmap=cmap)

    decorate (*args, axes=ax, wav=wav)


def decorate (*args, **kwargs):
    p = plot_params

    ax = kwargs['axes']
    ax.set_aspect (1)
    ax.axhline (0, ls=':', color='black')
    ax.axvline (0, ls=':', color='black')

    #ax.set_xlim (p['xlim'])
    #ax.set_ylim (p['ylim'])

    bz = ap.plot_bz2d (ax, uru2si2_001, repeat=p['repeat'], rotation=p['rot'],
                       kvec=[(v[0], v[1]) for v in uru2si2_kvec[0:2]])
