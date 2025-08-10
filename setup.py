from pathlib import Path
from numpy.distutils.core import Extension, setup

README = Path(__file__).with_name("README.MD")
long_description = README.read_text(encoding="utf-8") if README.exists() else ""
library_path = "/opt/homebrew/Cellar/gcc/15.1.0/lib/gcc/15"  # Ruta donde Homebrew instaló libgfortran.dylib

ext_modules = [
    Extension(
        name="wmf_heavy.cu",
        sources=["wmf_heavy/cuencas_heavy.f90"],
        f2py_options=["--quiet"],
    ),
    Extension(
        name="wmf_heavy.models",
        sources=["wmf_heavy/Modelos_heavy.f90"],
        f2py_options=["--quiet"],
    ),
]

setup(
    name="wmf_heavy",
    version="0.1.0",
    description="Water Model Framework (WMF) hydrological modeling package migration python 3.10",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Julian Uran",
    author_email="uranzea@gmail.com",
    license="GPL-3.0-or-later",
    url="https://github.com/uranzea/wmf_heavy",
    packages=["wmf_heavy"],
    ext_modules=ext_modules,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Fortran",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
    ],
)
