# -*- coding: UTF-8 -*-
from distutils.core import setup, Extension

with open("README.md", "r") as fh:
    long_description = fh.read()

isingLib=Extension('isinglib',sources=['./mcsolver/isingLib.c'],language='c',extra_compile_args=['-std=c99','-fPIC','-O3'])
xyLib=Extension('xylib',sources=['./mcsolver/xyLib.c'],language='c',extra_compile_args=['-std=c99','-fPIC','-O3'])
heisenbergLib=Extension('heisenberglib',sources=['./mcsolver/heisenbergLib.c'],language='c',extra_compile_args=['-std=c99','-fPIC','-O3'])

setup(
    name="mcsolver",
    version="3.0.4",
    author="Liang Liu",
    author_email="liangliu@main.sdu.edu.cn",
    description="A user friendly program to do Monte Carlo sims for magnetic systems",
    long_description=long_description,
    #long_description_content_type="text/markdown",
    url="https://github.com/golddoushi/mcsolver",
    packages=["mcsolver"],
    #include_package_data=True,
    ext_package="mcsolver.lib",
    ext_modules=[isingLib,xyLib,heisenbergLib],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Development Status :: 4 - Beta",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
)

#from distutils.sysconfig import get_python_lib
#from subprocess import check_call
#mcs_path=get_python_lib()+'/mcsolver/'
#check_call(["rm -rf %sisinglib.so"%(mcs_path)],shell=True)
#check_call(["rm -rf %sxylib.so"%(mcs_path)],shell=True)
#check_call(["rm -rf %sheisenberglib.so"%(mcs_path)],shell=True)
#check_call(["ln -s %s/lib/isinglib.* %sisinglib.so"%(mcs_path,mcs_path)],shell=True)
#check_call(["ln -s %s/lib/xylib.* %sxylib.so"%(mcs_path,mcs_path)],shell=True)
#check_call(["ln -s %s/lib/heisenberglib.* %sheisenberglib.so"%(mcs_path,mcs_path)],shell=True)