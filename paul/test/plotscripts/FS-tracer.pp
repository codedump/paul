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
import paul.loader.igor as igor

import numpy as np
import scipy.ndimage.interpolation as spni
import matplotlib.cm as cm

plot_params = {
    'slice': {
        'no': 24,
        'center': (0.0, 0.0),
        'picker': 3,
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
    '''
    This function needs two waves to operate:
    - wav[0] is the FS map to operate
    - wav[1] is the FS trace to save the trace points.
      wav[1] is a 3 x num_slices wave, where the 1st dimension
      has the following meaning:
        i=0: x-coordinates of the FS trace points
        i=1: y-coordinates of the FS trace points
        i=2: radial coordinate of the FS trace points relative to info['FS trace']['center']
    '''

    global slice_data, plot_params
    p = plot_params

    # figure initialization...
    kav['can'].clear()
    fig = kav['fig']

    # data initialization: define fsmap and fstrace waves
    fsmap = kav['wav'][0]
    if len(kav['wav']) > 1 and kav['wav'][1].info.has_key('FS trace'):
        if kav['wav'][1].shape != (3, p['slice']['no']):
            log.debug ("Found slice array, but it needs resizing (shape: %s)" % str(kav['wav'][1].shape))
            fstrace = None
        else:
            fstrace = kav['wav'][1].copy()
    else:
        fstrace = None
            
    if fstrace is None:
        fstrace = np.zeros([3, p['slice']['no']]).view(wave.Wave)
    fstrace.info['FS trace'] = {'source': fsmap.infs('name'),
                                'slices': p['slice']['no'],
                                'center': p['slice']['center'], }
    slice_data['trace file'] = get_trace_file(fsmap, fstrace)

    fig.text (0.01, 0.01, "mouse left:  set trace points\nmouse right: save trace\ntrace file: %s" 
              % os.path.relpath(slice_data['trace file']))
    
    # slicing...
    tmp = slice2d_fsmap_radial (fsmap, center=p['slice']['center'], num=p['slice']['no'])
    slice_data['angles']  = tmp['angles']
    slice_data['slices']  = kav['slices']  = tmp['slices']
    slice_data['coords']  = kav['coords']  = tmp['coords']
    slice_data['fsmap']   = kav['fsmap']   = fsmap
    slice_data['fstrace'] = kav['fstrace'] = fstrace
    slice_data['trace plot'] = None
    slice_data['trace marks'] = None
    slice_data['center'] = p['slice']['center']

    # displaying...
    ax_fs = fig.add_subplot (gridplot ((1,3), (0,0), colspan=2))
    ax_ln = fig.add_subplot (gridplot ((1,3), (0,2), colspan=1))

    slice_data['plot axes'] = [ax_fs, ax_ln]

    # display FS map, radial lines, and FS trace
    ax_fs.imshow (fsmap.view(np.ndarray), extent=fsmap.imlim, interpolation='none', cmap=cm.hot)
    slice_data['trace coord'] = plotwater (ax_fs,
                                            wlist=[c[0] for c in tmp['coords']],
                                            xlist=[c[1] for c in tmp['coords']])

    # Slices have to be drawn manually to be pickable :-(
    offset     = max([np.nanmax(l) for l in kav['slices']]) * p['lines']['offs'][1]
    slice_data['slices offset'] = offset
    slice_data['slices plots'] = [ ax_ln.plot (kav['slices'][i].dim[0].range,
                                               kav['slices'][i]+offset*i,
                                               picker=p['slice']['picker'])
                                               for i in range(p['slice']['no']) ]

    update_trace (ax_fs, ax_ln, slice_data)
            
    decorate_fs (*av, axes=ax_fs, **kav)
    decorate_ln (*av, axes=ax_ln, **kav)

    slice_data['canvas'] = kav['can']
    kav['can'].draw()


def update_trace (ax_fs, ax_ln, slice_data):
    '''
    Updates trace marks on the FS map and the lines axes.
    '''
    num_slices = len(slice_data['slices'])
    slices = slice_data['slices']

    trace_x = slice_data['fstrace'][0]
    trace_y = slice_data['fstrace'][1]
    trace_rad = slice_data['fstrace'][2]
    
    offset  = slice_data['slices offset']

    if slice_data.setdefault('trace marks', None) is not None:
        print "deleting old trace marks"
        slice_data['trace marks'].remove()
    slice_data['trace marks'] = ax_ln.vlines (trace_rad,
                                              [slices[i](trace_rad[i]) + (i-0.3)*offset for i in range(num_slices)],
                                              [slices[i](trace_rad[i]) + (i+0.3)*offset for i in range(num_slices)],
                                              color='black', lw=1.0)

    if slice_data.setdefault('trace plot', None) is not None:
        #slice_data['trace plot'].remove()
        del slice_data['trace plot']
    
    trace_x = np.sin(slice_data['angles']) * trace_rad
    trace_y = np.cos(slice_data['angles']) * trace_rad
    slice_data['trace plot'] = ax_fs.plot (trace_x, trace_y)
    

def get_trace_file (fsmap, fstrace, stub='fstrace'):
    '''
    Returns a file path suitable to save the fstrace in
    by sequentially trying a number of options:
    
      1. if fstrace already is associated with a file on disk,
      (i.e. fstrace.info['tmp']['last path'] exists),
      that one is used

      2. if fsmap is associated with a file on disk
      (i.e. fsmap.info['tmp']['last path'] is defined),
      a new (non-existent) file name is constructed using this
      and 'stub' as a starting point

      3. else a new file name is constructed using only stub
      as a file name
    '''

    #print "fstrace:", fstrace.infs('tmp', 'last path')
    #print "fsmap:  ", fsmap.infs('tmp', 'last path')
    #print

    tmp = fstrace.infs('tmp', 'last path')
    if tmp is not None and len(tmp) > 0:
        return tmp

    tmp = fsmap.infs('tmp', 'last path')
    if (tmp is not None) and (len(tmp) > 0):
        i = tmp.rfind(".")
        if i < 0:
            i = len(tmp)
        prefix = tmp[:i]

    mkfn = lambda pre, stb, cnt: "%s_%s%s.ibw" % (pre, stb, ("%d" % cnt) if cnt >= 0 else "")
    cnt = -1
    while os.path.exists(mkfn(prefix, stub, cnt)):
        cnt += 1

    return mkfn(prefix, stub, cnt)

    
def on_pick (event):
    '''
    Called when user selects object on the canvas.
    '''

    global slice_data
    
    ipos = event.ind[0] + round(0.5*(event.ind[-1]-event.ind[0]))
    xpos = event.artist.get_xdata()[ipos]

    id_list = [id(s[0]) for s in slice_data['slices plots']]

    slice_index = id_list.index(id(event.artist))
    if (slice_index < 0):
        print "Oops, unknown slice", event.artist
    print "  Slice: ", slice_index, "xpos", xpos, "at", ipos
    slice_data['fstrace'][2][slice_index] = xpos
    
    update_trace(slice_data['plot axes'][0],
                 slice_data['plot axes'][1],
                 slice_data)
    
    slice_data['canvas'].draw()

    

def on_click (event):
    '''
    Called when used clicks the canvas.
    '''

    global slice_data
    
    if event.button == 1:
        pass
    elif event.button == 2:
        pass
    elif event.button == 3:
        log.debug ("Writing '%s'" % slice_data['trace file'])
        igor.wave_write (slice_data['fstrace'], slice_data['trace file'])
    else:
        print "Whoops!"

    
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
