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


def wave_note_write (infile=None, note='', wav=None, outfile=None):
    '''
    Writes the note section of the specified wave.
    To this purpose, wave is read from 'wfile' if not an existing
    wave is already specified by the 'wav' parameter.
    The wave is written back to 'wfile', if not a different output
    file is specified by the 'outfile' parameter.
    '''

    if outfile is not None:
        outfile_name = outfile.name if hasattr(outfile, 'name') else outfile

    if infile is not None:
        infile_name = infile.name if hasattr (infile, 'name') else infile

    if wav is None:
        wav = ig.wave_read (infile_name, note_parse=False)

    # write the wave (using an intermediary temporary file)
    with NamedTemporaryFile (delete=False) as tmp:
        tmp.close()
        ig.wave_write (wav, tmp.name, note=note)
        os.rename (tmp.name, outfile_name)


def wave_note_edit (wfile_in, wfile_out):
    '''
    Reads the specified wave, opens the Info text in a text editor,
    waits for the user to edit the text (i.e. close the file),
    and saves the new info text into the wave.
    '''

    note0 = wave_note_read(wfile_in)

    # create temp file with the note text
    with NamedTemporaryFile (delete=False) as tmp:
        tmp.write(str(note0))
        tmp_name = tmp.name

    # start editor
    editor = os.environ['EDITOR'] if os.environ.has_key("EDITOR") else "nano"
    if call([editor, tmp_name]) != 0:
        print "Editor died unexpectedly. Aborting."
        return None

    # retrieve edited notes
    with open(tmp_name) as tmp:
        note1 = ''.join(tmp.readlines())
        os.remove (tmp_name)

    wave_note_write (infile=wfile_in, note=note1, outfile=wfile_out)


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

    parser.add_argument ('-o', '--output',  type=argparse.FileType('w'), help="Output wave. If not specified, input wave will be overwritten.", nargs='?', default=None)
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
        if param.output is None:
            param.output = param.input[0]
        wave_note_edit (param.input[0], param.output)

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
    act.add_argument('-p', '--print', action="store_true",   help="Print the current wave note and exit.", dest='nprint', )
    act.add_argument('-s', '--set', nargs=2,                 help="Set specified axis attribute.")
    act.add_argument('-i', '--inc', nargs=2,                 help="Increment specified axis attribute.")
    act.add_argument('-r', '--read', type=str,               help="Read specified axis attribute.")

    parser.add_argument ('-d', '--dim', '--dims',     type=d2ilist, help="Comma-seperated list of dimensions to manipulate.", default=None)
    parser.add_argument ('-o', '--output', nargs='*', type=argparse.FileType('w'), help="Output wave. If not specified, input wave will be overwritten.", default=None)
    parser.add_argument('input',           nargs='+', type=argparse.FileType('r'),
                        help="IBW file(s) to process.")

    if func is not None:
        # call wave_note sub-application by default
        parser.set_defaults (func=func)
    


def wave_dim_main (param):
    '''
    wave_dim entry routine (intrinsic dimension scaling manipulation)
    '''

    if param.output is None:
        param.output = param.input

    for fin, fout in zip(param.input, param.output):
        wav = ig.wave_read(fin)
        wav_save = False

        # prepare list of dimensions to be processed (defaults to: all)
        if param.dim is None:
            param.dim = range(wav.ndim)

        # set attribute / increase attribute
        if param.set is not None:
            wav_save = True
            [ setattr (wav.dim[i], param.set[0], eval(param.set[1]))
              for i in param.dim ]

        # increase attribute -- special case of "set"
        elif param.inc is not None:
            wav_save = True
            [ setattr (wav.dim[i], param.inc[0],
                       getattr(wav.dim[i], param.inc[0]) + eval(param.inc[1]))
              for i in param.dim ]

        # read attribute
        elif param.read is not None:
            [ pprint (getattr(wav.dim[i], param.read))
              for i in param.dim ]
        
        # print note (param.nprint)
        else:        
            print string.join([ str(d) for d in wav.dim ], '\n')


        # re-write the wave, if modifications were performed
        if wav_save:
            wav_name = fout.name if hasattr(fout, "close") else fout
            with NamedTemporaryFile (delete=False) as tmp:
                tmp.close()
                ig.wave_write (wav, tmp.name)
                os.rename (tmp.name, wav_name)


##
## Wave-Dump: dump wave data to ASCII format
##
def wave_dump_init (subparser=None, func=None, name='wave-dump'):
    '''
    Wave-Dump: initialization routine.

    Purpose of the sub-application is to dump wave data in an (x,y)
    manner suitable for further processing with other applications
    (i.e. gnuplot).
    '''

    parser = subparser.add_parser (name, help="IBW ASCII dumping.")

    # for reading comma-seperated lists of numbers from strings

    #act = parser.add_mutually_exclusive_group()
    #act.add_argument('-p', '--print', action="store_true",   help="Print the current wave note and exit.", dest='nprint', )

    parser.add_argument('input',           nargs='+', type=argparse.FileType('r'),
                        help="IBW file(s) to process.")

    if func is not None:
        # call wave_note sub-application by default
        parser.set_defaults (func=func)
    


def wave_dump_main(param):
    '''
    wave_dim entry routine (intrinsic dimension scaling manipulation)
    '''
    for f in param.input:
        wav = ig.wave_read(f)

        if wav.ndim == 1:
            for x,val in zip(wav, wav.dim[0].range):
                print x, "\t", val

        elif wav.ndim == 2:
            for block, x in zip(wav, wav.dim[0].range):
                for val, y in zip(block, block.dim[0].range):
                    print x, "\t", y, "\t", val
                print

        else:
            print "%s: dimension %d > 2, don't know what to do." % wav.ndim
            sys.exit(-1)


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
    #main_p.add_argument ('-C', '--copyright', action='store_true', help="Prints licensing information and exits")
    #main_p.add_argument ('-d', '--debug',     action='store_true',         help="Set verbosity level to debug.")
    #main_p.add_argument ('-q', '--quiet',     action='store_true',         help="Surpresses informational output on stdout.")
    sub_p = main_p.add_subparsers()
    
    # add sub-applications to subparser
    wave_note_init (name='wave-note', subparser=sub_p, func=wave_note_main)
    wave_dim_init  (name='wave-dim',  subparser=sub_p, func=wave_dim_main)
    wave_dump_init (name='wave-dump', subparser=sub_p, func=wave_dump_main)


    # parse command line and call main application
    args = [os.path.basename(sys.argv[0])]+sys.argv[1:]
    params = main_p.parse_args (args)

    params.func(params)
