UPDATING PPM.PY WEBSITE ON GH_PAGES:

A quick warning about updating the documentation:
    
    This involves switching branches on your git repository!
    To switch branches all changes must be updated to the master,
    thus when updating the makefile will add and commit all changes
    in your local repository to the remote repository! 
    If there are things in your local repository you don't want pushed 
    you should remove them.
    
The docs are updated through the command:

    $ make gh-pages

This command must be run in the PyPPM/docs folder on the master branch

TROUBLESHOOTING:

    Common errors: github must be configured before running or else github 
    will ask you to input your username and email as globals.
    
    If you try updating twice in a row it may not recognize your command
    simply change out of the directory and back in and it will work.
    
ABOUT WRITING DOCSTRINGS:

    Sphinx is a python package that reads the docstrings and generates a 
    webpage out of them.
    
    The format of the docstrings is thus particular since it will be translated 
    into markdown by Sphinx.
    
    A copy of the style guidlines can be found here:
    
    https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt
    
ABOUT PLOTTING EXAMPLES IN DOCSTRINGS:

    Sphinx also has the ability to generate plots using a matplotlib extension
    the command for this is:
    
    .. plot::
        
        plot(...)
    
    or 
    
    .. ipython::
    
        @savefig your_fig_name width=4
        plot(...)
        
    Both which generate plots
    
    and the code used for the plot should be indented, and the whitespace before and
    after is important.
    
    For more info see:
    
    http://matplotlib.org/sampledoc/extensions.html
    http://matplotlib.org/sampledoc/ipython_directive.html
    https://ipython.org/ipython-doc/3/api/generated/IPython.sphinxext.ipython_directive.html