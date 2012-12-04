paul
====

Paul is a visual data navigation and analyis tool written on top of Python,
Qt4, IPython, Numpy and Matplotlib (just to name the most prominent
dependencies). It's basically a series of scripts for easier everyday
handling of "normal" scientific data. '"normal" scientific data'
in this context means what you would usually get out of spectral
analysers, tunneling microscope images etc. In my case, it means
data from angle-resolved photoelectron spectroscopy.

It allows point-and-view displaying of scientific spectra, integrates
a Matplotlib viewer (and soon editor :-) ) with an easy-to-use file
navigator and offers powerful live data introspection and processing
from an integrated IPython shell.

It features some sueful ideas borrowed from a commercial data analysis 
tool called IgorPro(tm) made by a great company called Wavemetrics,
a number of data reading/writing functions for various more or less
obscure scientific data formats, and some new ideas of my own.

Here's an incomplete feature list:


Supported data formats:
-----------------------

  . Reading and writing if IgroPro(tm) IBW files (binary waves)

  . Reading and unpacking of IgorPro(tm) PXP and PXT
    files (packed experiment files)

  . Reading of Elmitec's proprietary LEED/LEEM image format
    (used in low electron energy microscopy)

  . Generating of IgorPro(tm) UXP (Unpacked Experiment Files)
    for easier exchange of data between Paul and Igor.


Data processing:
----------------

  . Implementation a "Wave" class, i.e. a numpy.ndarray subclass
    that can handle intrinsic scaling information as typically
    used in IgorPro(tm)

  . Implementation of indexing/slicing algorithms in "Wave" based
    not the coordinate system of the intrinsic scale (i.e.
    data slicing can be done not only like foo[x,y,z], but also
    like foo(p,q,r), where p,q,r are index objects specifying
    position information in the coordinate system intrinsic to
    the data :-) )

  . Basic functions for common data analysis tasks (slicing
    of data, generating of waterfall diagrams etc) that take
    advantage of the "Wave" interface


Data visualization and interaction:
-----------------------------------

  . paul-viewer:  a scripting-based Matplotlib plotting application
                  similarly to the GUI displayed by Pylab (in fact,
                  it makes heavy use of the Matplotlib's FigureCanvas
		  object :-) ) with the ability to load and
                  automatically reload the plotting script when
                  modified on disk.

  . paul-browser: a click-and-display data browser application, allowing
                  to quickly skim through 1D or 2D data files on disk.

  . paul-shell:   a module to expose as much functionality of 
                  paul-browser and paul-viewer to the IPython command
                  line, making Paul integration with the Spyder
                  toolset easier.


The Paul package is released under GPL-v3.
                  
-- 
Florin Boariu <florin.p(at)rootshell.ro>
                 
    
