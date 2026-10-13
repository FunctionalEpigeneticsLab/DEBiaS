from datetime import datetime
import os
import sys

sys.path.insert(0, os.path.abspath("../.."))

project = "DEBiaS"
author = "Xinran Song"
copyright = f"{datetime.now().year}, {author}"

release = "1.0.0"
version = release

extensions = [
    "myst_parser",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

master_doc = "index"
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
