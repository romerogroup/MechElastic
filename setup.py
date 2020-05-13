import setuptools
from numpy.distutils.core import Extension, setup

setup(
    name="MechElastic",
    description="A Python library to calculate elastic properties of materials.",
    version="0.1",
    author="Sobhit Singh",
    author_email="smsingh@mix.wvu.edu",
    url="https://github.com/sobhitsinghphy/MechElastic",
    download_url="https://github.com/sobhitsinghphy/MechElastic/archive/0.1.tar.gz",
    license="LICENSE.txt",
    scripts=["scripts/MechElastic"],
    install_requires=["numpy", "spglib", "prettytable"],
)
