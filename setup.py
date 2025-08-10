from pathlib import Path
from numpy.distutils.core import Extension, setup

# Read README for long description
README = Path(__file__).with_name("README.MD")
long_description = README.read_text(encoding="utf-8") if README.exists() else ""

extensions = [
    Extension(name="cu", sources=["wmf_heavy/cuencas_heavy.f90"]),
    Extension(name="models", sources=["wmf_heavy/Modelos_heavy.f90"]),
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
    ext_modules=extensions,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Fortran",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
    ],
)
