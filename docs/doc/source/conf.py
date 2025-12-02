# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'ModEM'
copyright = '2025, Gary Egbert, Anna Kelber, Naser Meqbel, Hao Dong'
author = 'Gary Egbert, Anna Kelber, Naser Meqbel, & Hao Dong'
release = 'develop'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'bizstyle'
html_static_path = ['_static']

html_sidebars = {
    "globaltoc_maxdepth" : 3
}


html_theme_options = {
    "globaltoc_maxdepth" : 3
}