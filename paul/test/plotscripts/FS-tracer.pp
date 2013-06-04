#!/usr/bin/python

#
# Simple script to assist the creation of Fermi surface traces.
#
# The idea is roughly this:
#  . load the selected file (an FS map) in a 2D image on the left
#  . slice the 2D image radially at equidistant angles
#  . display radial slices in a lines diagram on the right
#  . if graph is clicked in the lines diagram, plot
#    the corresponding position on top of the 2D image
#  . save list of points in a file to permit later processing
#

import logging, math
log = logging.getLogger (__name__)

from paul.toolbox.pscr import *
from paul.toolbox.mpltrix import *

import paul.base.wave as wave
import numpy as np
import scipy.ndimage.interpolation as spni
import matplotlib.cm as cm

plot_params = {
    'slice': {
        'no': 24,
        'center': (0.0, 0.0),
        'sheet file': "FS-sheet.ibw",
    },

    'lines': {
        'offs': (0, 0.3),
    },
}

slice_data = { }


def init(*av, **kav):
    log.debug ("registering mouse events...")
    global slice_data
    slice_data['pick_id'] = kav['can'].mpl_connect ('pick_event', on_pick)
    slice_data['click_id'] = kav['can'].mpl_connect ('button_press_event', on_click)

    
def exit(*av, **kav):
    log.debug ("done.")
    global slice_data
    kav['can'].mpl_disconnect (slice_data['pick_id'])
    kav['can'].mpl_disconnect (slice_data['click_id'])

    
def populate (*av, **kav):
    fsmap      = kav['wav'][0]
    kav['can'].clear()
    p  = plot_params

    global slice_data
    
    sd = slice_data

    ax_fs = kav['fig'].add_subplot (gridplot ((1,3), (0,0), colspan=2))
    ax_ln = kav['fig'].add_subplot (gridplot ((1,3), (0,2), colspan=1))

    ax_fs.imshow (fsmap, extent=fsmap.imlim, interpolation='none', cmap=cm.hot)

    sd['sheets'] = []  # list of FS sheets
    sd['sheets'].append ({
        'pts-pol': np.zeros([2, p['slice']['no']]),  # polar point coordinates
        'pts-car': np.zeros([2, p['slice']['no']]),  # cartesian point coordinates
        'wfile': None, 
        })

    # slicing...
    tmp = slice2d_fsmap_radial (fsmap, center=p['slice']['center'], num=p['slice']['no'])
    kav['slices'] = tmp['slices']
    kav['coords'] = tmp['coords']

    ## Obsolete... LineCollection doesn't work nice with pick events,
    ## for the purpose of selecting single lines :-(
    #kav['slices'] = plotwater (ax_ln, tmp['slices'], offs=p['lines']['offs'])
    kav['coords'] = plotwater (ax_fs,
                               wlist=[c[0] for c in tmp['coords']],
                               xlist=[c[1] for c in tmp['coords']])

    # Slices have to be drawn manually to be pickable :-(
    max_data = max([np.nanmax(l) for l in kav['slices']]) 
    offset = 0
    for sl, i in zip (kav['slices'], range(len(kav['slices']))):
        ax_ln.plot (sl.dim[0].range, sl+offset, picker=10)
        offset += max_data*p['lines']['offs'][1]

        
    #ax_ln.set_ylim ((0, max_data * len(kav['coords'])+1))
    
    decorate_fs (*av, axes=ax_fs, **kav)
    decorate_ln (*av, axes=ax_ln, **kav)

    kav['can'].draw()

    
def on_pick (event):
    '''
    Called when user selects object on the canvas.
    '''
    ipos = event.ind[0] + round(0.5*(event.ind[-1]-event.ind[0]))
    xpos = event.artist.get_xdata()[ipos]
    print "Grr!"
    print "  Artist: ", event.artist, "xpos", xpos, "at", ipos, event.ind, event.artist.get_picker()

def on_click (event):
    '''
    Called when used clicks the canvas.
    '''
    #print "Boo!"
    pass

    
def decorate_fs (*args, **kwargs):
    '''
    Valid kwargs: can, wav, fig, axes. 
    '''
    fsmap = kwargs['wav'][0]
    coords = kwargs['coords']
    ax = kwargs['axes']
    
    ax.set_xlim (fsmap.dim[1].lim)
    ax.set_ylim (fsmap.dim[0].lim)
    

def decorate_ln (*args, **kwargs):
    '''
    Valid kwargs: can, wav, fig, axes. 
    '''
    slices = kwargs['slices']
    ax = kwargs['axes']

    
def slice2d_fsmap_radial (fsmap, center=(0.0, 0.0), num=12):
    '''
    Slices the specified FS map radially.
    Returns a dictionary with following fields:
      'angles': the slice angles
      'slices': the 1D slices
      '': 
    '''

    sd = {}
    
    # initialize slicing infrastructure
    sd['angles'] = []
    sd['coords'] = []
    sd['slices'] = []
    
    # center of the radial slicing
    center_pt = np.array(center)

    # furthest distance from the center to either of the borders
    max_radius = np.array([ abs(l-c) for l,c in zip([d.lim for d in fsmap.dim], center_pt) ]).max()
    
    for angle in np.array(np.linspace(start=0, stop=math.pi*2, num=num, endpoint=False)):
        radial_pt = np.array([math.cos(angle), math.sin(angle)]) * max_radius
        sl, ln = slice2d_line (fsmap, xstart=center_pt, xstop=radial_pt)
        sd['angles'].append (angle)
        sd['slices'].append (sl)
        sd['coords'].append (ln)

    return sd



def slice2d_line (wav2d, numpts='auto',
                  xstart=None, xstop=None, xstep='auto',
                  xline=None, iline=None,
                  xline_out=None, iline_out=None):
    '''
    Returns a slice of the 2D wave along a specified line.
    The line has a resolution specified by 'xstep', or
    of the smallest dim-delta, if xstep='auto'.

    The line may be specified in either of the following ways:

    - as a start / stop tuple pair in axis coordinates via `xstart` and `xstop`
    - as an explicitly specified line [(x0, x1, ... xn), (y0, y1, ... yn)] via `xline` in axis coordinates
    - as an explicitly specified line [(x0, x1, ... xn), (y0, y1, ... yn)] via `iline` in index coordinates

    Returns a (slice, coords) tuple of waves containing the slice
    starting with the first point, and the corresponding coordinates in
    the 2D map.
    The delta of the slice is re-calculated accordingly to 'xstep' or
    to the distance between the specified point coordinates.
    '''

    # line was not explicitly specified -- calculate from xstart / xstop
    if (iline is None) and (xline is None):
        if xstep == 'auto':
            xstep = min ([d.delta for d in wav2d.dim])

        if numpts == 'auto':
            to = np.array (xstop)
            fr = np.array (xstart)
            Diff = to-fr
            numpts = math.ceil(np.sqrt(np.sum(Diff*Diff)) / xstep)
        else:
            raise ValueError ("Specific 'numpts' is not yet implemented")

        # axis coordinates
        xline = np.array([ np.linspace (start=x1, stop=x0, num=numpts)
                  for x1, x0 in zip(xstart, xstop) ])
        
    elif xline is not None:
        # calculate an average 'dim-delta' parameter.
        # This is only accurate if xline is aequidistant.
        xstep = (xline[-1] - xline[0])

    
    # at this point, an explicit 'xline' or 'iline' exists
    if iline is None:
        iline = np.array([ d.x2i(l) for d, l in zip (wav2d.dim, xline) ])
        
    # the magic: this is where the slicing happens
    line = spni.map_coordinates (wav2d.view(np.ndarray), iline, mode='constant', cval=np.NaN).view(wave.Wave)

    # remove coordinate points that have no data
    # (e.g. because they point outside the 2D map)
    null_pts = np.isnan(line)
    xline[0,null_pts] = np.NaN
    xline[1,null_pts] = np.NaN
    line.dim[0].delta = xstep

    return line, xline
