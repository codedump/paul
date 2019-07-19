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
a Matplotlib viewer with an easy-to-use file
navigator and offers powerful live data introspection and processing
from an integrated IPython shell.

It features some useful ideas borrowed from a commercial data analysis 
tool called *IgorPro(tm)* made by *Wavemetrics*,
a number of data reading/writing functions for various more or less
obscure scientific data formats, and some new ideas.

Here's an incomplete feature list:


Supported data formats:
-----------------------

  - Reading and writing if *IgroPro(tm)* `IBW` files (binary waves)

  - Reading and unpacking of *IgorPro(tm)* `PXP` and `PXT`
    files (packed experiment files)

  - Reading of *Elmitec*'s proprietary LEED/LEEM image format
    (used in low electron energy microscopy)

  - Generating of *IgorPro(tm)* `UXP` (Unpacked Experiment Files)
    for easier exchange of data between Paul and Igor.


Data processing:
----------------

  - Implementation of a `Wave` class, i.e. a `numpy.ndarray` subclass
    that can handle intrinsic scaling information as typically
    used in *IgorPro(tm)* -- arguably one if its most powerful features.

  - Implementation of indexing/slicing algorithms in `Wave` based
    not the coordinate system of the intrinsic scale (i.e.
    data slicing can be done not only like `foo[x,y,z]`, but also
    like `foo(p,q,r)`, where `p,q,r` are index objects specifying
    position information in the coordinate system intrinsic to
    the data.

  - Basic functions for common data analysis tasks (slicing
    of data, generating of waterfall diagrams etc) that take
    advantage of the `Wave` interface.
    
  - Ability to use plot specification scripts in order to define
    Matplotlib visualisation presets and data plots -- a good idea
    also borrowed from *IgorPro(tm)*.
    
  - The ability to modularly "blend into" a Unix environment, manipulate
    `IBW` data files (...or possibly others; interface is easily extensible)
    without the need to load / write large, opaque, binary data files (think
    *IgoroPro(tm)'s* `PXP` files) which cannot be manipulated except by
    a single, monolithic system.
    Its integrability is arguably one of the most powerful features of *Paul*; in fact,
    *Paul* was extensively used in scientific day-to-day business within
    Unix `Makefile` setups to automatically generate graphs from data
    files and lot specifitation scripts.


Data visualization and interaction:
-----------------------------------

  - **paul-viewer**:  a scripting-based Matplotlib plotting application
                  similarly to the GUI displayed by Pylab (in fact,
                  it makes heavy use of the Matplotlib's FigureCanvas
		  object :-) ) with the ability to load and
                  automatically reload the plotting script when
                  modified on disk.

  - **paul-browser**: a click-and-display data browser application, allowing
                  to quickly skim through 1D or 2D data files on disk.

  - **paul-shell**:   a module to expose as much functionality of 
                  paul-browser and paul-viewer to the IPython command
                  line, making Paul integration with the Spyder
                  toolset easier.


Project Status / Usability 
--------------------------

Consider the status a "usable prototype". Most of the code that appears
to be working does indeed work robustly and was used to analyze data
and generate graphs for various high-ranking publications.
Formal documentation is lacking,
but many portions of the code, in particular the science-heavy parts,
are very well explained & documented.
Writing few examples of how to rapidly start using it would probably
go a long way in helping to get a new user started.

Data exchange code, in particular reading / writing of `IBW` and `UXP`
files has received a large amount of testing and can be considered
fairly robust; writing of `UXP` could still be a little buggy. However,
data exchange is probably out of date with regards to newer versions
of the file formats. Code state is from around 2014, I assume
*IgorPro(tm)* will have moved along with its formats during that time.

Architecture of the system is in principle validated, but largerly
undocumented, and in an unfinished state.

On the plus side, the codebase is fairly small and still possible
to understand and modify/refactor by a single person or a small
group of developers.

The project is mostly orphaned, meaning that I don't actively work
on it anymore, although I would love to see it continued. If you
feel inclined, let me know, and I'll support to you the best of
my abilities. 

The Paul package is released under GPL-v3.
                  
-- 
Florin Boariu <florin.p(at)rootshell.ro>
