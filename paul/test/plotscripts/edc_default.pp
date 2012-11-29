#!/usr/bin/python

def decorate (*args, **kwargs):
    ax = kwargs['can'].axes
    ax.set_ylim (0, 200000)
    ax.set_xlim (-0.04, 0.02)
    pass
