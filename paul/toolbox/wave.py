#!/usr/bin/python

'''
Command line meta data manipulation tool.

It is implemented as a series of sub-applications, called by argv[0]:

  - wave-note: manipulte the note section of the waves
  - wave-dim:  manipulate dimension information

'''

import sys, argparse, os, string
from pprint import pprint
import paul.loader.igor as ig

from subprocess import call
from tempfile import NamedTemporaryFile

##
## Wave-Note part
##

def wave_note_read (wfile):
    '''
    Returns the note section of the specified wave.
    '''
    text = bytearray()
    wav = ig.wave_read (wfile, note_text=text, note_parse=False)
    return str(text)


def wave_note_write (wfile, note, wav=None):
    '''
    Writes the note section of the specified wave.
    '''

    wfile_name = wfile.name if hasattr (wfile, 'name') else str(wfile)
    
    if wav is None:
        wav = ig.wave_read (wfile_name, note_parse=False)

    # write the wave (using an intermediary temporary file)
    with NamedTemporaryFile (delete=False) as tmp:
        tmp.close()
        ig.wave_write (wav, tmp.name, note=note)
        os.rename (tmp.name, wfile_name)


def wave_note_edit (wfile, output=None):
    '''
    Reads the specified wave, opens the Info text in a text editor,
    waits for the user to edit the text (i.e. close the file),
    and saves the new info text into the wave.
    '''

    note0 = wave_note_read(wfile)

    # create temp file with the note text
    with NamedTemporaryFile (delete=False) as tmp:
        tmp.write(str(note0))
        tmp_name = tmp.name

    # start editor
    editor = os.environ['EDITOR']
    if call([editor, tmp_name]) != 0:
        print "Editor died unexpectedly. Aborting."
        return None

    # retrieve edited notes
    with open(tmp_name) as tmp:
        note1 = ''.join(tmp.readlines())
        os.remove (tmp_name)

    wave_note_write (wfile.name, note1)


def wave_note_init (subparser=None, func=None, name='wave-note'):
    '''
    Initializes application for wave_note sub-application.
    (Mainly just concerns setting up the option parser.)
    '''

    parser = subparser.add_parser (name, help="IBW data note manipulation routines.")

    act = parser.add_mutually_exclusive_group()
    act.add_argument('-p', '--print',  action="store_true", help="Print the current wave note and exit.", dest='nprint', )
    act.add_argument('-P', '--parse',  action="store_true", help="Parse the current wave note, print the parsed tree and exit.")
    act.add_argument('-e', '--edit',   action="store_true", help="Edit the wave note in $EDITOR.")
    act.add_argument('-k', '--kill',   action="store_true", help="Kill (i.e. remove) wave notes.")
    act.add_argument('-a', '--append', type=str,            help="Append the specified string to wave notes.")

    parser.add_argument('input', nargs='+', type=argparse.FileType('r'),
                        help="IBW file(s) to process.")

    if func is not None:
        # call wave_note sub-application by default
        parser.set_defaults (func=func)
    


def wave_note_main (param):
    '''
    wave_note: main application (entry) routine.
    '''

    files = param.input

    # edit note
    if param.edit:
        wave_note_edit (param.input[1], param.output)

    # parse note (and print)
    elif param.parse:
        [ pprint(ig.wave_note_parse_simple(wave_note_read(f))) for f in files]

    # kill note
    elif param.kill:
        [ wave_note_write (f, "") for f in files ]
        
    # append to note
    elif param.append:
        [ wave_note_write (f, "%s\n%s" % (wave_note_read(f), param.append))
          for f in files ]

    # print note (param.nprint)
    else:
        for f in files:
            print wave_note_read (f)



##
## Wave-Dim part
##

def wave_dim_init (subparser=None, func=None, name='wave-dim'):
    '''
    Initializes application for wave_dim sub-application.
    (Intrinsic dimension scaling manipulation.)
    '''

    parser = subparser.add_parser (name, help="IBW intrinsic dimension information manipulation routines.")

    # for reading comma-seperated lists of numbers from strings
    d2ilist = lambda nr: [int(i)   for i in nr.split(',')]
    d2flist = lambda nr: [float(i) for i in nr.split(',')]

    act = parser.add_mutually_exclusive_group()
    parser.add_argument('-a', '--axis', '--axes', type=d2ilist, help="Comma-seperated list of dimensions to manipulate.", default=None)
    parser.add_argument('-p', '--print', action="store_true",   help="Print the current wave note and exit.", dest='nprint', )

    parser.add_argument('-s', '--set', nargs=2,                 help="Set specified axis attribute.")
    parser.add_argument('-i', '--inc', nargs=2,                 help="Increment specified axis attribute.")
    parser.add_argument('-r', '--read', type=str,               help="Read specified axis attribute.")

    parser.add_argument('input', nargs='+', type=argparse.FileType('r'),
                        help="IBW file(s) to process.")

    if func is not None:
        # call wave_note sub-application by default
        parser.set_defaults (func=func)
    


def wave_dim_main (param):
    '''
    wave_dim entry routine (intrinsic dimension scaling manipulation)
    '''

    for f in param.input:
        wav = ig.wave_read(f)

        if param.axis is None:
            param.axis = range(wav.ndim)
            #print "Auto-selecting axes:", param.axis

        # set attribute
        if param.set:
            pass
            #[ d.offset = param.offset for d in wav.dim[param.axis] ]

        # increase attribute
        elif param.inc:
            pass

        # read attribute
        elif param.read:
            pass
        
        # print note (param.nprint)
        else:        
            print string.join([ str(d) for d in wav.dim ], '\n')
                



##
## Main application shell
##    

if __name__ == "__main__":

    '''
    This is the main program shell. The idea is that we gather
    general data here (debug/verbose flags, input/output file names
    etc).

    Then we delegate the processing to corresponding sub-applications.
    For now, the sub-applications are all in this file, implemented in
    the *_main() functions.
    '''

    # initialize main parser (general options)
    main_p = argparse.ArgumentParser(description=r"IBW data note manipulation routines.")
    main_p.add_argument ('-C', '--copyright', action='store_true', help="Prints licensing information and exits")
    main_p.add_argument ('-o', '--output',    type=argparse.FileType('w'), help="Output wave. If not specified, input wave will be overwritten.", nargs='?', default=None)
    main_p.add_argument ('-d', '--debug',     action='store_true',         help="Set verbosity level to debug.")
    main_p.add_argument ('-q', '--quiet',     action='store_true',         help="Surpresses informational output on stdout.")
    sub_p = main_p.add_subparsers()
    
    # add sub-applications to subparser
    wave_note_init (name='wave-note', subparser=sub_p, func=wave_note_main)
    wave_dim_init  (name='wave-dim',  subparser=sub_p, func=wave_dim_main)


    # parse command line and call main application
    args = [os.path.basename(sys.argv[0])]+sys.argv[1:]
    params = main_p.parse_args (args)

    params.func(params)
