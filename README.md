# mcsolver
A user friendly and efficient tool implementing Monte Carlo simulations to estimate Curie/Neel temperature

Support multiple ocassions, e.g. standard ferromangetic/anti-ferromagnetic systems, DMI, Kiteav non-diagonal exchange interactions, dipole-dipole long-range couplings, with external fields.

Original version contributor: Dr. Liang Liu* 1.Shenzheng University 2.Shandong University
Email: liangliu@mail.sdu.edu.cn

You can download the packed .exe (only tested in Windows 10 platform) from the following link. Wish it can find something helpful for you. And if it was used for publication, please cite:
[1] Magnetic switches via electric field in BN nanoribbons. Applied Surface Science 480(2019)

Link for exe: https://pan.baidu.com/s/1EaDqOOdB7AP9WXrwEIEaxQ
passwd: 52ze


Brief tutorial:

A. using mcsolver via .exe, e.g., in Windows platform

  NOTE: the mcsolver.exe maybe reported to be virus and removed by some anti-virus software. I still have no ideal about this and maybe you need add it into white list. Otherwise you can use mcsolver as a python package (see Section B below).

  Download and extract file from upper link (or download the .zip from https://github.com/golddoushi/mcsolver/raw/master/mcsolver.20.10.10update.zip), then open .exe (maybe wait 10+ sec.), fill out all parameters, click startMC Btn, then wait for the results.

  How to define parameters?

  Style I Define parameters via GUI:

    1. Define the three lattice vectors of primitive cell.

    2. Define all basic spins in primitive cell, note that the fractional coordinates are supposed. Ani represents the single-ion anisotropies in xyz directions (It is useless in Ising model, and only former two are used in XY model). As well as, note that the units of anisotropies are in Kelvin. 

    3. Define all exchange interactions (bonds). There are nine matrix elements for one J including Jxx, Jyy, Jzz, Jxy, Jxz, Jyz, Jyx, Jzx, Jzy, respectively. Each element discribes the coupling between two components of spins. For example, a basic bond term can be expressed as 
    $$S1\dot J\dot S2 = S_{1x}J_{xx}S_{2x} + S_{1y}J_{yy}S_{2y} + S_{1z}J_{zz}S_{2z} + S_{1x}J_{xy}S_{2y} + ...$$
    For Ising model, since only one component is available for each spin, only the first element Jxx is used. As well as for XY model, only Jxx, Jyy, Jxy, Jyz are used. 
    In this step, you can click one of the bonds to review the actual linking in lattice on view pannel. Activated bond is depicted with bold and yellow line while others are green. You may drag left/right mouse Btn to rotate and expand/shrink the model shown in view pannel. 
    By the way, these parameters can be obtained via fitting experimental spin-wave spectra or DFT calculations. One of the general method for calculating anisotropic bonds and orbitals was introduced in dio: 10.1021/acs.jpclett.0c01911.

    4. Define other parameters, including 
    the start and end temperatures, number of temperature interpolations (for the temperature scanning)
    the start, end and number of samplings for external field (for the magnetic field scanning)
    nthermal is the total steps to make system enter balanced states, nsweep is the total steps involved in mesearing, tau denotes the MC updates for each step
    xAxis denote the physical quantity put in x-axis of right-hand Result viewer, it can be either T(for illustration of M-T curv) or H (for illustration of hysteresis loop).
    model type, algorithm (only Metropolis and Wolff are supported now)
    nFrame is the num. of output spin configurations, using for illustrating spin configurations in equilibrium or non-equilibrium states.

    5. Set spin_i and spin_j and the lattice vector between them, for correlation mesearments. (if spin_i=spin_j and overLat=0 0 0, then you will get susceptibility for spin_i)

    6. Set the core resources for parallel calc.

    7. (Optional) Save current parameters into file.

    8. Click startMC Btn to start.

    9. Wait for the diagram update in right pannel. Afterwards, you can find a file result.txt in the root directory of mcsolver, there are many useful informations including the averaged spin (on spin_i and j defined in step 5), correlation between spin_i and j, internal energy, specific heat capacith, and Binder cumulant U4, etc. If you handle the sims with more than one cores then the results may not be ordered according to temperature, however, the correspondences in every line are ok.

  Style II Define parameters via loading file
  
    1. click load Btn to load settings, and here I prepared the setting for CrI3 with exchanges up to 2nd nearest neiboring. You can modify the sample file for your own purposes, with any txt editor. 

    2. You can define the Topological section to compute the (thermally averaged) topological charges. Every circuits are made by three orbitals enclosing the triangle anti-clockwisely. And all the circuits should cover the zone with exactly the same area of unit cell. 

    3. Click startMC Btn to start.

B. using mcsolver as a python package

  Style I Build form binary file in Windows platform

    Install the package via command:

    pip install mcsolver

    Note that python>=3, matplotlib, numpy, tkinter are prerequisite

    Afterwards, you can import mcsolver into your own python code and use function:

    mcsolver.loadMC("parameterfile")

    to start simulation. Preparation of parameterfile is the same as in section A.
    
    There are one sample file sim_XYmodel_under_Windows.py in sample folder. To use this, change your current path (in console) into sample folder, and type command: python sim_XYmodel_under_Windows.py

    NOTE: since mcsolver employ python-parallel you have to use freeze_support() before calling loadMC(...).

  Style II Build from source in Linux platform

    Download all codes presented here, use install command:

    python setup.py install

    Afterwards, you can import mcsolver into your own python code and use function:

    mcsolver.loadMC("parameterfile")

    to start simulation. Preparation of parameterfile is the same as in section A.

    There are one sample file sim_XYmodel_under_Linux.py in sample folder. To use this, cd to sample folder, and type command: python sim_XYmodel_under_Linux.py

    Note that the parallelization of mcsolver is not perfect. Now it cannot parallelize between multiple machines but amongst mutiple cores in a single machine (that is, only SMP mode is efficient). Therefore submit the job into one node if you are working with clusters.

C. using code mode (not recommend)
  
   Download all codes, compile all .c files inito .so files the dynamic libraries. And run cmd: "python win.py" to load GUI and go on.