from setuptools import setup, find_packages
import glob

setup(
    name = "cos_dark",
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
    py_modules = [x.split(".py")[0] for x in glob.glob("*py") if "setup.py" not in x],
    # Only uncomment the below line if you are creating a package 
    # (with a package directory scheme and __init__.py)
#	packages = find_packages(),
    install_requires = ["setuptools",
                        "numpy",
                        "astropy",
						"pandas",
                        "sqlalchemy",
                        "pymysql",
                        "pyyaml",
                        "bokeh>=2.1.1",
                        "holoviews>=1.13.4",
                        "datashader>=0.11.0",
                        "dask"]
    )


