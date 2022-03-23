from setuptools import setup, find_packages
import glob

setup(
    name = "acdc",
    version = "1.0",
    description = "Perform an optimized COS/FUV dark correction",
    author = "Jo Taylor",
    author_email = "jotaylor@stsci.edu",
    keywords = ["astronomy", "hst", "cos", "ultraviolet", "uv", "dark", 
                "instrumentation"],
    classifiers = ['Programming Language :: Python',
                   'Programming Language :: Python :: 3',
                   'Development Status :: 1 - Planning',
                   'Intended Audience :: Science/Research',
                   'Topic :: Scientific/Engineering :: Astronomy',
                   'Topic :: Scientific/Engineering :: Physics',
                   'Topic :: Software Development :: Libraries :: Python Modules'],
    # Only uncomment the below line if you are creating a package 
    packages = ["acdc"],
    package_dir = {"acdc": "acdc"},
    package_data = {"acdc": ["data/*"]},
    install_requires = ["setuptools",
                        "numpy",
                        "astropy",
                        "pandas",
                        "sqlalchemy",
                        "pymysql",
                        "pyyaml",
                        "bokeh>=2.1.1",
                        "holoviews>=1.14.4",
                        "datashader>=0.13.0",
                        "panel>= 0.11.3",
                        "dask",
                        "matplotlib",
                        "asdf"]
    )


