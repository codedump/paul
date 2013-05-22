#!/usr/bin/python

'''
Creates the normalized and FDD-ed versions of the data
from the raw versions in specified subdirectories.
'''

import paul.base.wave as wave
import paul.loader.igor as igor
import paul.toolbox.arpes as arpes

import numpy as np
import sys, os, argparse

def cook_scienta_get_ef (data, axis=0):
    '''
    Returns the position of the Fermi level in 'data'
    by analysing Scienta metadata.
    '''
    # this is the total energy offset of the scale:
    # i.e. usually, the scale will have an energy relative to
    # E_F = 0, but the data was measured at an E_kin,
    # i.e. relative to the vacuum energy. We need to know
    # the absolute value of the kinetic energy for the
    # polar -> kspace conversion.
    return ((data.infv('SES', 'Low Energy')  - idata.dim[axis].lim[0]) +
            (data.infv('SES', 'High Energy') - idata.dim[axis].lim[1])) / 2
    
def cook_ax_shift (data, offsets=None):
    '''
    Shifts axis offset(s) by the amount specified by the tuple 'offsets'.
    '''

    if offsets == None:
        offsets = [0] * data.ndim

    odata = data.copy()

    print " * Axis offsets:",
    for dim, offs in zip (odata.dim, offsets):
        print offs,
        dim.offset += offs
    print ""
    
    return odata
    

def cook_norm_noise (data, axis=0, norm_range=None):
    '''
    Normalizes data along along 'not axis' by the intensity
    integrated along 'axis' in the range 'norm_range'.
    '''
    if norm_range is None:
        norm_range = (0, data.dim[axis].size)
        
    print " * Normalizing intensity along %d, range %d...%d" % (axis, norm_range[0], norm_range[1])
    return arpes.norm_by_noise (data, axis=axis, ipos=norm_range)
    

def cook_deg2ky (data, axis=0, tilt=0.0, eoffs=0.0):
    '''
    Conversion polar -> k coordinates, where 'axis' is the energy axis.
    'tilt' specifies the tilt angle at which the data was measured.
    '''
    print " * Converting to k-coordinates (energy offset relative to output wave: %f, tilt: %f)" % (eoffs, tilt)
    axtag = ('ed' if axis==0 else 'de')
    return arpes.deg2ky_single (data, axes=axtag, tilt=tilt, eoffs=eoffs)

    
def cook_norm_gnd(data, axis=0, norm_range=None):
    '''
    Calculates and substracts background along 'axis' in the range 'norm_range'.
    '''
    if norm_range is None:
        norm_range = (0, data.dim[axis].size)
        
    gnd = np.average(data[norm_range[0]:norm_range[1]]) if axis==0 \
      else np.average(data.swapaxes(0,1)[norm_range[0]:norm_range[1]])

    print " * Substracting background (%d...%d: %f)" % (norm_range[0], norm_range[1], gnd)
    return data-gnd

    
def cook_fdd(data, axis=0, kT=0, Ef=0):
    '''
    Normalizes data by the Fermi-Dirac distribution with specified
    parameters kT and E_F. 'axis' designates the energy axis.
    '''
    print " * Normalizing to Fermi-Dirac distribution (kT: %f eV, Ef = %f eV)" % (kT, Ef)
    return arpes.norm_by_fdd (data, axis=axis, kT=kT, Ef=Ef)
    

def cook_from_beamline (data, eax=0,
                        eoffs=None, doffs=None, tilt=None, kT=None, Ef=None,
                        norm_region=None, **kwargs):
    '''
    Performs a user-specified data preparation cycle on beam-line data.
    Which steps are performed depends on the input parameters,
    roughly speaking the function is capable of:
      - axis shifting (eoffs, doffs parameters)
      - intensity normalization along momentum axis (norm_noise, norm_region parameters)
      - background substraction (norm_gnd, norm_region parameters)
      - polar -> ky conversion (deg2ky, tilt parameters)
      - FDD normalization (norm_fdd, kT, Ef parameters)

    Returns a single output wave.
    '''

    
    # degree offset update (energy offset needs to be performed after the deg2ky transform)
    if doffs is not None:
        data = cook_ax_shift (data, offsets=((0, -doffs) if eax==0 else (-doffs, 0)))
    if eoffs is not None:
        data = cook_ax_shift (data, offsets=((-eoffs, 0) if eax==0 else (0, -eoffs)))
        if Ef is not None: # E_F / E-offsset interplay?
            # user specified E_F was relative to original energy axis -- needs adjusting!
            Ef -= eoffs
            print "   Fermi level adjusted by %f eV, now at Ef=%f." % (eoffs, Ef)
        else:
            # no E_F specified: we default to 0 eV on the _new_ axis!
            print "   Fermi level is af 0 eV on the new axis (default)"
            Ef = 0.0

    # intensity normalization across momentum axis
    # (i.e. by calculating noise along energy axis)
    if kwargs['norm_noise']:
        data = cook_norm_noise (data, axis=eax, norm_range=norm_region)

    # polar -> ky transformation
    if kwargs['deg2ky']:
        # Energy offset for deg2k transformation:
        # deg2ky() gets energy information from wave axis. This application
        # does not (yet) support an explicit pol-to-k energy offset,
        # so if the energy axis of the original wave was modified, you're
        # out of luck.
        # *However*, as a courtesy to tiresome scientist eyes (and to make
        # things more complicated :-) )the following exceptions apply:
        #
        #  - if the energy was modified just now using the 'eoffs'
        #    option, we are feeding this as an explicit offset parameter (eoffs).
        #
        #  - if energy axis contains a zero (i.e. limits have different sign),
        #    then we assume that zero is at the Fermi level. In this case,
        #    if 'Ef' is set, we use this as an offset for deg2ky.
        # 
        if eoffs is not None:
            d2k_offs = eoffs
        elif (sign(data.dim[eax].start) == -sign(data.dim[eax].end)) and (Ef is not None):
            d2k_offs = Ef
        else:
            d2k_offs = 0.0
        data = cook_deg2ky (data, axis=eax, tilt=tilt, eoffs=d2k_offs)

    # background substraction, if requested
    if kwargs['norm_gnd']:
        data = cook_norm_gnd (data, axis=eax, norm_range=norm_region)

    # fdd normalization
    if kwargs['norm_fdd']:
        data = cook_fdd (data, axis=eax, kT=kT, Ef=Ef)

    return data

    
def arpes_opt_parse (args):
    '''
    Returns the parsed command line options
    '''
    parser = argparse.ArgumentParser (description="Raw ARPES data post-processing.")

    # actions -- what to do (multiple specifications possible)
    parser.add_argument ('-E', '--eoffs', '--Eoffset', type=float,  help="Substract energy axis offset")
    parser.add_argument ('-D', '--doffs', '--Doffset', type=float,  help="Substract k-parallel (degree) axis offset")
    parser.add_argument ('-K', '--deg2ky', '--ky',        action='store_true', help="Perform polar -> k(parallel) coordinate transformation.")
    parser.add_argument ('-N', '--norm-noise', '--noise', action='store_true', help="Normalize intensity distribution along the momentum axis using thermally populated area")
    parser.add_argument ('-G', '--norm-gnd', '--gnd',     action='store_true', help="Substract background ('auto' calculates background from the thermally populated area)")
    parser.add_argument ('-F', '--norm-fdd', '--fdd',     action='store_true', help="Perform Fermi-Dirac normalization")
    

    # various input flags
    parser.add_argument ('-n', '--norm-region', '--norm-range', metavar='FROM:TO', type=str, help="Region along the energy axis suitable for normalisation calculations (i.e. a region which lies in the thermally populated area). Units are assumed as indices, if numbers are integers, or intrinsic axis units otherwise.")
    parser.add_argument ('-a', '--eax', '--axis', type=int, default=0, help="Which axis is the energy dimension.")
    parser.add_argument ('-k', '--kT', type=float,                help="k(Boltzmann) * T(effective) parameter for Fermi-Dirac normalization, in eV.")
    parser.add_argument ('-e', '--Ef', type=float,                help="The position of the Fermi level along the energy axis (index units if data is integer, intrinsic axis units otherwise).")
    parser.add_argument ('-t', '--tilt', type=float, default=0.0, help="The tilt angle at which the sample was measured.")

    # I/O flags
    parser.add_argument('-o', '--output', type=str,           help='Where to write binary output (-- for stdout, --- for the same as input). May overwrite input file!')
    parser.add_argument('--safe',   action='store_true',      help='If output would overwrite input file, save output in a temporaray file first, then move it over the input file to lower the risk of data corruption.', default=True)
    parser.add_argument('--unsafe', action='store_false',     help='Opposite of "--safe". May lead to data loss if anything goes wrong with saving the data.', dest='safe')
    parser.add_argument('input', type=argparse.FileType('r'), help='Input file (raw ARPES data, IBW format).')

    # parse arguments
    return parser.parse_args (args)
    

def arpes_save (wav, param):
    '''
    Saves (a processed) wave according to command line specifications from 'param'
    '''

    from tempfile import NamedTemporaryFile
    
    # output
    if (param.output == '--'):
        param.output = sys.stdout

    elif param.output == '---':
        param.output = param.input.name
        
    elif param.output is None:
        pindex = param.input.name.rfind('.')
        if pindex < 0:
            pindex = len(param.input.name)-1
        base = param.input.name[0:pindex]
        ext  = param.input.name[pindex:]
        param.output = base + ("_ky"*param.deg2ky) + ("_fdd"*param.norm_fdd) + ext
        
    print "Writing: %s" % param.output

    if (param.output == param.input.name) and param.safe:
        with NamedTemporaryFile (delete=False) as tmp:
            tmp.close()
            igor.wave_write (wav, tmp.name)
            os.rename (tmp.name, param.output)
    else:
        igor.wave_write (wav, param.output)

    

if __name__ == "__main__":

    p = arpes_opt_parse (sys.argv[1:])

    print "Loading:", p.input.name
    inw = igor.load (p.input)

    ## perform argument post-processing (i.e. parsing of string
    ## arguments, unit conversion where necessary).

    # eax is always defined (default: 0)
    print " - Energy axis:", p.eax

    # Ef is always eV
    if p.Ef is not None:
        p.Ef = eval(p.Ef)
    if isinstance(p.Ef, int):
        p.Ef = inw.dim[p.eax].i2x(eval(p.Ef))
    print " - Fermi level is", p.Ef, "eV"

    # norm-region is 'index' if specified
    if p.norm_region is not None:
        nr = [eval(i) for i in p.norm_region.split(':')]
        p.norm_region = [ i if isinstance(i, int) else int(inw.dim[p.eax].x2i_rnd(i)) for i in nr ]
    print " - Normalization region is", p.norm_region

    ## this is where the magic happens
    outw = cook_from_beamline (inw, **vars(p))
    
    arpes_save (outw, p)
