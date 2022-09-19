#!/usr/bin/env python
#
# ravenpy documentation build configuration file, created by
# sphinx-quickstart on Fri Jun  9 13:47:02 2017.
#
# This file is execfile()d with the current directory set to its
# containing dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

# If extensions (or modules to document with autodoc) are in another
# directory, add these directories to sys.path here. If the directory is
# relative to the documentation root, use os.path.abspath to make it
# absolute, like shown here.
#
import os
import sys

sys.path.insert(0, os.path.abspath(".."))

import ravenpy  # noqa

# -- General configuration ---------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "nbsphinx",
    "sphinx_click",
    "sphinx_codeautolink",
    "sphinx_copybutton",
]

autosectionlabel_prefix_document = True
autosectionlabel_maxdepth = 2

autosummary_generate = True
nbsphinx_execute = "auto"

# To ensure that underlined fields (e.g. `_field`) are shown in the docs.
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "private-members": False,
    "special-members": False,
}

# To avoid having to install these and burst memory limit on ReadTheDocs.
autodoc_mock_imports = [
    "affine",
    "cftime",
    "cf_xarray",
    "click",
    "climpred",
    "clisops",
    "fiona",
    "gdal",
    "netCDF4",
    "osgeo",
    "geopandas",
    "haversine",
    "holoviews",
    "hvplot",
    "lxml",
    "owslib",
    "pandas",
    "pyproj",
    "rasterio",
    "rioxarray",
    "scipy",
    "shapely",
    "statsmodels",
    "xarray",
    "xclim",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = [".rst"]

# The master toctree document.
master_doc = "index"

# General information about the project.
project = "RavenPy"
copyright = "2021, David Huard"
author = "David Huard"

# The version info for the project you're documenting, acts as replacement
# for |version| and |release|, also used in various other places throughout
# the built documents.
#
# The short X.Y version.
version = ravenpy.__version__
# The full version, including alpha/beta/rc tags.
release = ravenpy.__version__

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = "en"

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This patterns also effect to html_static_path and html_extra_path
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False

# Suppress "WARNING: unknown mimetype for ..." when building EPUB.
suppress_warnings = ["epub.unknown_project_files"]

# Avoid "configuration.rst:4:duplicate label configuration, other instance in configuration.rst"
autosectionlabel_prefix_document = False

# -- Options for HTML output -------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Theme options are theme-specific and customize the look and feel of a
# theme further.  For a list of options available for each theme, see the
# documentation.
#
html_theme_options = {"style_external_links": True}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = list()  # "_static"

# -- Options for HTMLHelp output ---------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = "ravenpydoc"

# -- Options for LaTeX output ------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',
    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',
    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',
    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, documentclass
# [howto, manual, or own class]).
latex_documents = [
    (master_doc, "ravenpy.tex", "RavenPy Documentation", "David Huard", "manual"),
]

# -- Options for manual page output ------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [(master_doc, "ravenpy", "RavenPy Documentation", [author], 1)]

# -- Options for Texinfo output ----------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (
        master_doc,
        "ravenpy",
        "RavenPy Documentation",
        author,
        "ravenpy",
        "One line description of project.",
        "Miscellaneous",
    ),
]
