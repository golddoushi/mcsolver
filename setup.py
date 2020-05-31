# -*- coding: UTF-8 -*-
import setuptools

with open("README.md", "r", encoding='UTF-8') as fh:
    long_description = fh.read()
 
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
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: GNU License v3",
        "Operating System :: OS Independent",
    ],
)