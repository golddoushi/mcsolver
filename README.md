# mcsolver
A user friendly tools implementing Monte Carlo simulations to estimate Curie/Neel temperature

Original version contributor： Dr. Liang Liu* 1.Shenzheng University 2.Shandong University
Email: liangliu@mail.sdu.edu.cn

You can download the packed .exe (only tested in Windows 10 platform) from the following link. Wish it can find something helpful for you. And if it was used for publication, please cite:
[1] Magnetic switches via electric field in BN nanoribbons. Applied Surface Science 480(2019)

Link for exe：https://pan.baidu.com/s/1EaDqOOdB7AP9WXrwEIEaxQ
passwd: 52ze

Installation：
No need for installation

Brief tutorial:
Open exe (maybe wait 10 sec.), fill out all parameters, click startMC Btn, then wait for the results.

Style I Define parameters via GUI:
  1. Define the three lattice vectors of primitive cell.

  2. Define all basic spins in primitive cell, note that the fractional coordinates are supposed. Ani represents the single-ion anisotropies in xyz directions (It is useless in Ising model, and only former two are used in XY model). As well as, note that the units of anisotropies are in Kelvin. 

  3. Define all exchange interactions (bonds). Only Jz is used for Ising model, and Jz, Jx are used for XY model, and all three J for Heisenberg model. In this step, you can click one of the bonds to review the actual linking in lattice on view pannel. Activated bond is depicted with bold and yellow line while others are green. You may drag left/right mouse Btn to rotate and expand/shrink the model shown in view pannel. 

  4. Define other parameters, including the start and end temperatures, number of temperature interpolations, nthermal the total steps to make system enter balanced states, nsweep the total steps involved in mesearing, tau the MC updates for each step, and model type, algorithm (only Metropolis and Wolff are supported now)

  5. Set spin_i and spin_j and the lattice vector between them, for correlation mesearments.

  6. Set the core resources for parallel calc.

  7. (Optional) Save current parameters into file.

  8. Click startMC Btn to start.

  9. Wait for the diagram update in right pannel. Afterwards, you can find a file result.txt in the root directory of mcsolver, there are many useful informations including the averaged spin (on spin_i and j defined in step 5), correlation between spin_i and j, internal energy, specific heat capacith, and Binder cumulant U4, etc. If you handle the sims with more than one cores then the results may not be ordered according to temperature, however, the correspondences in every line are ok.

Style II Define parameters via load file
  1. click load Btn to load settings, and here I prepared the setting for CrI3 with exchanges up to 2nd nearest neiboring. You can modify the sample file for your own purposes, with any txt editor. 
  2. Click startMC Btn to start.