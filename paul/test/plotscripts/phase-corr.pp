#!/usr/bin/python

from scipy.fftpack import fftn, ifftn
import numpy as np
import paul.base.wave as wave
import paul.toolbox.arpes as arpes

def populate (*args, **kwargs):
    '''
    Valid kwargs: can, wav, fig, axes. 
    '''
    kwargs['can'].clear()
    fig = kwargs['fig']
    wav = kwargs['wav']

    axes = [ fig.add_subplot(131),
             fig.add_subplot(132),
             fig.add_subplot(133),]


    if len(wav) < 2:
        return

    efi = [ arpes.fermi_guess_efi(w) for w in wav ]
    ef  = [ arpes.fermi_guess_ef (w) for w in wav ]

    imgA = wave.regrid(wav[0][efi[0]-20:efi[0]+20,:], {'numpts': 400})
    imgB = wave.regrid(wav[1][efi[0]-20:efi[0]+20,:], {'numpts': 400})
    #imgA = wav[0]
    #imgB = wav[1]

    foo      = ifftn(fftn (imgA) * ifftn(imgB))
    foo_roll = np.roll(np.roll(foo, foo.shape[0]/2, 0), foo.shape[1]/2, 1)
    wfoo     = foo_roll.real.view(wave.Wave)


    print "argmax:", wfoo.argmax()
    print "foo shape:", wfoo.shape, foo.shape
    pos0 = ((wfoo.argmax() / foo.shape[1]),
            (wfoo.argmax() % foo.shape[1]))

    pos = (pos0[0] - foo.shape[0]/2,
           pos0[1] - foo.shape[1]/2)


    print "pos:", pos0
    print "shift:", pos[0]/10.0, pos[1]/10.0
    print "Fermi levels:", efi, "delta", efi[1]-efi[0]

    axes[0].imshow (wfoo, interpolation='none')
    axes[1].imshow (wav[0], extent=wav[0].imlim)
    axes[2].imshow (wav[0], extent=wav[1].imlim)
    #axes[1].imshow (imgB, extent=imgA.imlim)
    #axes[2].imshow (imgA, extent=imgB.imlim)

    for a, e in zip(axes[1:], ef):
        a.set_ylim (sorted(a.get_ylim()))
        a.axhline (y=e, ls=':', color='black')
