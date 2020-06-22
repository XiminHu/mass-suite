# import os
# import sys
import setuptools

# Get appropriate version and release info
# ver_file = os.path.join(mss)

# Grabs description from the README file
# with open(ver_file) as f:
#     exec(f.read())

# Give setuptools a hit to complain if it's too old of a version
# 24.2.0 added the python_requires option
# Enables setuptools to install wheel on-the-fly
#SETUP_REQUIRES = ['setuptools >= 24.2.0']
#SETUP_REQUIRES += ['wheel'] if 'bdist_wheel' in sys.argv else []

setuptools.setup(name = 'mass-suite',
	  version = '1.0.1',
      description = 'A package for HRMS data analysis',
      long_description = open('README.md', 'r').read(),
      long_description_content_type = 'text/markdown',
      url = 'https://github.com/XiminHu/mass-suite',
      author = 'Ximin Hu, Derek Mar, Nozomi Suzuki, Bowei Zhang',
      author_email = 'xhu66@uw.edu',
      include_package_data = True,
      packages = setuptools.find_packages(),
      classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
        python_requires='>=3.6',
)
