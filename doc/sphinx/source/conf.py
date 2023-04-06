# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'gauss2dfit'
copyright = '2023, Dan Taranu'
author = 'Dan Taranu'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'breathe',
    # TODO: solve exhale issues below.
    # 'exhale',
    'sphinx_rtd_theme',
    'sphinx.ext.autodoc',
    'sphinx.ext.coverage',
]

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# -- Breathe configuration ---------------------------------------------------
breathe_default_project = "gauss2dfit"

# Setup the exhale extension
# TODO: exhale tries to parse non-existent cppreference tagged rst files
# e.g. doc/sphinx/source/doc/api/function_cpp/atomic/atomic_fetch_sub_1.rst
exhale_args = {
    # These arguments are required
    # TODO: exhale needs a path, but no env var from sphinx or otherwise points to build/
    # TODO: exhale insists on an in-source build (i.e. in docs/, not build/docs/)
    "containmentFolder":     "", #os.environ["EXHALE_API_DOC"],
    "rootFileName":          "library_root.rst",
    "doxygenStripFromPath":  "..",
    # Heavily encouraged optional argument (see docs)
    "rootFileTitle":         "Library API",
    # Suggested optional arguments
    "createTreeView":        True,
    "exhaleExecutesDoxygen": False,
}

# Tell sphinx what the primary language being documented is.
primary_domain = 'cpp'
