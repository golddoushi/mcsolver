# -*- coding: UTF-8 -*-
import setuptools

with open("README.md", "r", encoding='UTF-8') as fh:
    long_description = fh.read()

isingLib=setuptools.Extension('isinglib.so',sources=['./isingLib.c'],language='c',extra_compile_args=['-std=c99','-fPIC'])
xyLib=setuptools.Extension('xylib.so',sources=['./xyLib.c'],language='c',extra_compile_args=['-std=c99','-fPIC'])
heisenbergLib=setuptools.Extension('heisenberglib.so',sources=['./heisenbergLib.c'],language='c',extra_compile_args=['-std=c99','-fPIC'])

setuptools.setup(
    name="mcsolver",
    version="1.0.2",
    author="Liang Liu",
    author_email="liangliu@main.sdu.edu.cn",
    description="A user friendly program to do Monte Carlo sims for magnetic systems",
    long_description=long_description,
    #long_description_content_type="text/markdown",
    url="https://github.com/golddoushi/mcsolver",
    packages=setuptools.find_packages(),
    include_package_data=True,
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
