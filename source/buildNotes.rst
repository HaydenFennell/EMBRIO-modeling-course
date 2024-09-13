Notes on Current Build (0.1)
============================

**Updated: 03/19/23**

This site is currently a prototype, intended primarily to test functionality of the content delivery platform. We will begin populating module content soon.

**Known issues:**

* There are some issues with sidebar navigation on certain pages. This is due to some pages not being updated as we make rapid changes to content. A full source file rebuild will be performed (which will fix sidebar navigation) in the next build upload.
* Some of the example Jupyter notebook imports (see sections below modules) have broken code blocks. These are from interactive blocks from the notebook that are not currently supported by default in Sphinx. We are still looking for a possible workaround for embedding executable Jupyter cells into Sphinx documentation. We are also exploring methods for hosting pre-built worksheets to correspond to each module as needed.
* There are a number of equation and image rendering errors in some of the auto-imported Jupyter notebooks. These are due primarily to inflexibilities in Sphinx regarding blank lines and spaces that Jupyter is happy to overlook. These will need to be fixed by hand before importing notebooks to avoid rendering issues.