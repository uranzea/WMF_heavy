from pathlib import Path
from setuptools import setup


README = Path(__file__).with_name("README.md")
long_description = README.read_text(encoding="utf-8") if README.exists() else ""


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
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
    ],
)

