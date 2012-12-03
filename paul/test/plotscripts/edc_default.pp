#!/usr/bin/python

def decorate (*args, **kwargs):
    ax = kwargs['can'].axes
    wav = kwargs['wav']
    
    ## try to auto-magically guess the correct y scaling
    ## to sane values (note that FDD data will diverge,
    ## so using min()/max() is not an option here.)
    #ax.set_ylim (0, wav[wav.dim[0].size/2]*2)
    ax.set_ylim (0, wav(0)*2)    
    ax.set_xlim (-0.06, 0.02)
    
    ax.axvline (0, ls=':')
