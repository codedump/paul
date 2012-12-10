#!/usr/bin/python



def decorate (*args, **kwargs):
    ax = kwargs['can'].axes
    ax.axhline (0, ls=':', color='black')
    ax.axvline (0, ls=':', color='black')
    
    
