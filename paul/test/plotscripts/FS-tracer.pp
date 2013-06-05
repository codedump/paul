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

    #
    # general strategy:
    #  - (we're working with a single trace per session, but in
    #     principle the code should work with multiple traces)
    #  - there is _definitely_ only one FS map per session.
    #    relevant data to the 2D FS map is saved in slice_data[] directly
    #  - for every trace (again, currently only one per session),
    #    relevant data to the trace is saved in slice_data['traces'][n][...]
    #

    global slice_data, plot_params
    p = plot_params
    

    # figure initialization...
    kav['can'].clear()
    fig = kav['fig']
    ax_fs = fig.add_subplot (gridplot ((1,3), (0,0), colspan=2))
    ax_ln = fig.add_subplot (gridplot ((1,3), (0,2), colspan=1))
    slice_data['plot axes'] = [ax_fs, ax_ln]
    slice_data['canvas'] = kav['can']
    slice_data['traces'] = []

    # data init: fsmap and general settings
    fsmap = kav['wav'][0]
    slice_data['fsmap']  = fsmap

    # data init: trace wave
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
                                'center': p['slice']['center'],
                                'type': 'radial',}

    trace_data = {}
    trace_data['slices'] = slice2d_fsmap_radial (fsmap, center=p['slice']['center'], num=p['slice']['no'])
    trace_data['wave'] = fstrace
    trace_data['center'] = p['slice']['center']
    trace_data['file'] = get_trace_file(fsmap, fstrace)
    trace_data['plot'] = None
    trace_data['marks'] = None
    trace_data['marks pos'] = fstrace[2] if fstrace.shape[0] > 2 \
      else np.zeros(p['slice']['no'])
    fig.text (0.01, 0.01, "mouse left:  set trace points\nmouse right: save trace\ntrace file: %s" 
              % os.path.relpath(trace_data['file']))

    # add 0-th (and only) trace to the traces list
    slice_data['traces'].append(trace_data)


    # display FS map, radial lines, and FS trace
    ax_fs.imshow (fsmap.view(np.ndarray), extent=fsmap.imlim, interpolation='none', cmap=cm.hot)
    trace_data['coord'] = plotwater (ax_fs,
                                     wlist=[c[0] for c in trace_data['slices']['coords']],
                                     xlist=[c[1] for c in trace_data['slices']['coords']])

    # Slices have to be drawn manually to be pickable :-(
    offset = max([np.nanmax(l) for l in trace_data['slices']['waves']]) * p['lines']['offs'][1]
    trace_data['slices']['offset'] = offset
    trace_data['slices']['plots']  = [ ax_ln.plot (trace_data['slices']['waves'][i].dim[0].range, 
                                                   trace_data['slices']['waves'][i]+offset*i,
                                                   picker=p['slice']['picker'])[0]
                                                   for i in range(p['slice']['no']) ]
                                               

    update_trace (ax_fs, ax_ln, trace_data)
    decorate_fs (*av, axes=ax_fs, **kav)
    decorate_ln (*av, axes=ax_ln, **kav)
    kav['can'].draw()


def load_trace (axes, wav=None):
    '''
    Loads and prepares the specified trace.
    Returns a 'trace_data' map.

    `wav` may be empty, in which case a new trace will
    be created from scratch.
    '''


def update_trace (ax_fs, ax_ln, trace_data):
    '''
    Updates trace marks on the FS map and the lines axes.
    '''
    num_slices = len(trace_data['slices']['waves'])
    slices = trace_data['slices']['waves']

    trace_x   = trace_data['wave'][0]
    trace_y   = trace_data['wave'][1]
    trace_rad = trace_data['marks pos']
    offset    = trace_data['slices']['offset']

    # display tracing dashes
    if trace_data.setdefault('marks', None) is not None:
        trace_data['marks'].remove()
    ymin = [slices[i](trace_rad[i]) + (i-0.3)*offset for i in range(num_slices)]
    ymax = [slices[i](trace_rad[i]) + (i+0.3)*offset for i in range(num_slices)]
    trace_data['marks'] = ax_ln.vlines (trace_rad, ymin, ymax, color='black', lw=1.0)

    # display FS trace ontop of the FS map
    if trace_data.setdefault('plot', None) is not None:
        for p in trace_data['plot']:
            p.remove()
    trace_data['plot'] = ax_fs.plot (trace_y, trace_x, 'o', color='white')
    

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

    trace_data = slice_data['traces'][0]
    slices = trace_data['slices']
    
    ipos = event.ind[0] + round(0.5*(event.ind[-1]-event.ind[0]))
    xpos = event.artist.get_xdata()[ipos]

    slice_index = slices['plots'].index(event.artist)
    if (slice_index < 0):
        print "Oops, unknown slice", event.artist
    print "Slice %d, position %d (%f)" % (slice_index, ipos, xpos)

    # save marker position
    trace_data['marks pos'][slice_index] = xpos

    # save corresponding FS map position in the trace map
    for i in range(slice_data['fsmap'].ndim):
        trace_data['wave'][i][slice_index] = slices['coords'][slice_index][i][ipos]

    # update optical clues (markers, FS map plot, etc...)
    ax_fsmap  = slice_data['plot axes'][0]
    ax_slices = slice_data['plot axes'][1]
    update_trace(ax_fsmap, ax_slices, trace_data)
    
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
        log.debug ("Writing '%s'" % slice_data['traces'][0]['file'])
        print "Writing", os.path.relpath(slice_data['traces'][0]['file'])
        igor.wave_write (slice_data['traces'][0]['wave'], slice_data['traces'][0]['file'])
    else:
        print "Whoops!"

    
def decorate_fs (*args, **kwargs):
    '''
    Valid kwargs: can, wav, fig, axes. 
    '''
    fsmap = kwargs['wav'][0]
    ax = kwargs['axes']
    
    ax.set_xlim (fsmap.dim[1].lim)
    ax.set_ylim (fsmap.dim[0].lim)
    

def decorate_ln (*args, **kwargs):
    '''
    Valid kwargs: can, wav, fig, axes. 
    '''
    #slices = kwargs['slices']
    ax = kwargs['axes']

    
def slice2d_fsmap_radial (fsmap, center=(0.0, 0.0), num=12):
    '''
    Slices the specified FS map radially.
    Returns a dictionary with following fields:
      'angles': the slice angles
      'slices': the 1D slices
      '': 
    '''

    # initialize slicing infrastructure
    sd = { 'angles': [],
           'coords': [],
           'waves':  [], }
    
    # center of the radial slicing
    center_pt = np.array(center)

    # furthest distance from the center to either of the borders
    max_radius = np.array([ abs(l-c) for l,c in zip([d.lim for d in fsmap.dim], center_pt) ]).max()
    
    for angle in np.array(np.linspace(start=0, stop=math.pi*2, num=num, endpoint=False)):
        radial_pt = np.array([math.cos(angle), math.sin(angle)]) * max_radius
        sl, ln = slice2d_line (fsmap, xstart=center_pt, xstop=radial_pt)
        sd['angles'].append (angle)
        sd['waves'].append (sl)
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
