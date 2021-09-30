import setuptools

setuptools.setup( name='CEA_Wrap',
                  version='1.0',
                  description='A Python-Based wrapper for the NASA CEA Thermochemical Code',
                  url='https://github.com/civilwargeeky/CEA_Wrap',
                  author='civilwargeeky',
                  author_email='klingerd@purdue.edu',
                  license_files=["LICENSE"],
                  install_requires=[],
                  packages=setuptools.find_packages(),
                  python_requires=">=3.6",
                  package_data={
                    "": ["FCEA2.exe", "thermo.lib", "trans.lib", "thermo_spg.inp", "thermo_spg.out"],
                  },
                  include_package_data=True, # Add in pdfs and executables
                  zip_safe=False)