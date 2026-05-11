import os,sys; sys.path.insert(0, os.path.abspath('..'))

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Library for Low Frequencies'  # Replace with your actual project name [web:2]
copyright = '2025, Your Name'   # Update with current year and author [web:2]
author = 'ULU Group - Francesco de Gasperin'            # Replace with your name [web:2]

# The full version, including alpha/beta/rc tags.
release = '1.4'                 # Set your project's version [web:2]

# The short X.Y version.
version = '1.0'                 # Short version number [web:2]

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['autodoc2']                 # List of Sphinx extensions to enable (add 'sphinx.ext.autodoc' if needed) [web:2]


templates_path = ['_templates'] # Path to custom templates relative to this file [web:2]
exclude_patterns = []           # Patterns to exclude from source files [web:2]

language = 'en'                 # Documentation language [web:2]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'        # Theme for HTML output (options: 'alabaster', 'sphinx_rtd_theme') [web:2]
html_static_path = ['_static']  # Path to static files [web:2]

# Theme options are theme-specific and customize the look and feel of a theme.
html_theme_options = {}         # Additional theme options [web:2]


autodoc2_packages = [
    {
        "path": "../LiLF",  # Relative from docs/ to root/LiLF/ – no abspath needed [web:62]
        "auto_mode": True,  # Scans all .py files like lib_util.py
    }
]
