import setuptools
from setuptools.command.sdist import sdist
from setuptools.command.install import install
import platform, subprocess, os

VERSION = '1.5.0'

instructions = """
(Courtesy of Morgan Reusch)
Install XCode from the Mac App Store
Install the XCode command line tools:
1. Open Terminal
2. Type this command and hit enter: xcode-select --install

Install a fortran compiler:
1. Go to www.hpc.sourceforge.net
2. Download the latest version of GCC. File should be: gcc-9.2-bin.tar.gz
3. Open terminal
4. Go to directory where file was just downloaded. (e.g.: cd /Downloads)
5. Unzip the file (may not be necessary). Type: gunzip gcc-8.1-bin.tar.gz
6. Install GCC. Type: sudo tar -xvf gcc-8.1-bin.tar -C /   (it will ask you for your computer’s password)
7. To check if it installed properly you can type: gfortran -v
    If GCC installed correctly then you shouldn’t get an error.
"""

# From https://stackoverflow.com/a/36902139
class PostInstallCommand(install):
    """Pre-installation for installation mode."""
    def run(self):
      # Notes to future me: When pip installed, the package is installed to a temporary directory, install.run() will export everything to the package directory
      #    so you must do any compiling and moving files before that.
      # You can see printed output in pip by using --verbose
      # Before install.run is called you have the original directory structure with setup.py in here and CEA_Wrap in the next directory down
      
      # Note to maintainers: I stopped including "cea2.f" and "cea.inc" in the git. You can get these from the NASA website.
      sys = platform.system()
      if sys in ["Linux", "Darwin"]: # If we need to compile CEA for this system
        print("posix system detected, compiling FCEA for use on this machine")
        if sys == "Linux":
          print("Ensuring gfortran is installed")
          subprocess.check_call("sudo apt-get install gfortran".split())
        else:
          print("Hoping that you have gfortran and xcode installed. Otherwise, follow these instructions:\n", instructions)
        print("Compiling FCEA2")
        os.chdir("CEA_Wrap/assets")
        subprocess.check_call("gfortran cea2.f -o FCEA2".split())
        print("Compiling thermo and trans libs")
        subprocess.run("./FCEA2", input="thermo_spg\n", text=True, stdout=subprocess.DEVNULL)
        subprocess.run("./FCEA2", input="trans\n", text=True, stdout=subprocess.DEVNULL)
        os.chdir("../..")
        print("Process complete!")
      
      install.run(self) # Call normal install function
        
with open("README.md", "r") as fh:
    long_description = fh.read()
    # For viewing on PyPI
    long_description = "**PyPI NOTE: This package is installed with 'pip install CEA_Wrap' with an underscore. PyPI doesn't show underscores in package names**\n\n" + long_description

setuptools.setup( name='CEA_Wrap',
                  version=VERSION,
                  description='A Python-Based wrapper for the NASA CEA Thermochemical Code',
                  long_description=long_description,
                  long_description_content_type="text/markdown",
                  url='https://github.com/civilwargeeky/CEA_Wrap',
                  author='Daniel Klinger',
                  author_email='klingerd@purdue.edu',
                  license_files=["LICENSE"],
                  install_requires=["appdirs"],
                  packages=setuptools.find_packages(),
                  python_requires=">=3.6",
                  include_package_data=True, # Add in pdfs and executables defined in the manifest
                  zip_safe=False,
                  cmdclass={
                    'install': PostInstallCommand,
                  },
                  )