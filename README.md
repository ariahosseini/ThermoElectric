# Themoelectric.py — A Python Tool for Design of High ZT Nanoengineered Thermoelectrics
<div align="justify">
 
*Thermoelectric.py* is a computational framework that computes electron transport coefficients with unique features to design the nanoscale morphology of thermoelectrics (TEs) to obtain electron scattering that will enhance thermoelectric performance through electron energy filtering.
</div>

- [INSTALLATION](#INSTALATION)
- [THEORY](#THEORY)
  * [Model Electron Transport Coefficients In Bulk Thermoelectrics](#Model-Electron-Transport-Coefficients-In-Bulk-Thermoelectrics)
  * [Model Fermi Level](#Model-Fermi-Level)
  * [Model Electron Lifetime](#Model-Electron-Lifetime)
  * [Model Electron Transport Coefficients In Nanostructured Thermoelectrics](#Model-Electron-Transport-Coefficients-In-Nanostructured-Thermoelectrics)
- [CASE STUDY](#CASE STUDY)
  * [Mitigating the Effect of Nanoscale Porosity on Thermoelectric Power Factor of Si](#grid-world-generation-1)
    + [Complete Spec (key, lock, and pit)](#complete-spec--key--lock--and-pit-)
    + [Key Only](#key-only)
    + [Lock and Pit](#lock-and-pit)

- [CITATION](#Citation)


# INSTALLATION

Compatible with python 3.0 and upwards

```bash
git clone https://github.com/...
cd ...
pip install -e .
```


# THEORY

<div align="justify">
 
Thermoelectrics (TE) are a class of materials that convert heat directly into electricity. The performance of TE materials at a given temperature, 𝑇, is quantified by a dimensionless figure of merit ZT=(σS<sup>2</sup>)/κ T, where κ, σ and S are the material’s thermal conductivity, electrical conductivity and Seebeck coefficient, respectively. The power factor (σS<sup>2</sup>) in ZT depends on a combination of strongly interdependent electrical transport properties, that have a countervailing dependence of the charge carrier concentration. The tradeoff of these parameters is well understood, and it has become an accepted truth that optimal TE performance can only be obtained in semiconductors that are highly doped to a narrow window of optimized charge carrier concentration. If made sufficiently efficient and inexpensive, these materials could be used to recapturing low-grade waste heat from industrial process as useful electrical energy. The potential energy savings are vast. Recent studies have suggested that recuperating only 10% of heat lost into electricity can improve fuel energy efficiency by 20% while other studies has reported that more than 68% of U.S. energy consumption escaped as waste heat. 

</div>

## Model Electron Transport Coefficients In Bulk Thermoelectrics

<div align="justify">
 
<p> 
The electrical conductivity and thermopower of a population of independent charge carriers can be derived from the Boltzmann transport equation by integrating the contribution from all carrier states. In an isotropic system where the states can be enumerated by their energy, and using the single relaxation time approximation for the collision operator, these can be written as integrals over the carrier energy, E, so that σ, S, and κ<sub>e</sub> are given by
</p>

<p align="center">
<img src="https://latex.codecogs.com/png.latex?%5Cinline%20%5Cdpi%7B300%7D%20%5Csigma%20%3D%20%5Cdfrac%7B-1%7D%7B3%7De%5E2%5Cint%5Cchi%28E%2CT%29%5Ctau%28E%2CT%29dE%3D%5Cdfrac%7B-1%7D%7B3%7De%5E2%5CDelta_0" align="center" alt="drawing" width="500px">
</p>

<p align="center">
<img src="https://latex.codecogs.com/png.latex?%5Cinline%20%5Cdpi%7B300%7D%20S%20%3D%20%5Cdfrac%7B-1%7D%7BeT%7D%5Cdfrac%7B%5Cint%5Cgamma%28E%2CT%29%5Ctau%28E%2CT%29dE%7D%7B%5Cint%5Cchi%28E%2CT%29%5Ctau%28E%2CT%29dE%7D%3D%5Cdfrac%7B-1%7D%7BeT%7D%28%5CDelta_1-E_f%29" align="center" alt="drawing" width="500px">
</p>

<p align="center">
<img src="https://latex.codecogs.com/png.latex?%5Cinline%20%5Cdpi%7B300%7D%20%5Cbegin%7Balign*%7D%20%5Ckappa_e%20%26%3D%20%5Cdfrac%7B-1%7D%7B3T%7De%5E2%20%5Cleft%28%5Cint%5Czeta%28E%2CT%29%20%5Ctau%28E%2CT%29%20dE%20-%5Cdfrac%7B%28%5Cint%5Cgamma%28E%2CT%29%20%5Ctau%28E%2CT%29%20dE%29%5E2%7D%7B%5Cint%5Cchi%28E%2CT%29%5Ctau%28E%2CT%29%20dE%7D%5Cright%29%5C%5C%20%26%3D%20%5Cdfrac%7B-1%7D%7B3T%7D%5CDelta_0%28%5CDelta_2-%5CDelta_1%5E2%29%20%5Cend%7Balign*%7D" align="center" alt="drawing" width="700px">
</p>

<p>
Here the function χ(E,T)= ν(E)<sup>2</sup>∂f(E<sub>f</sub>,E,T)/∂E D(E), lumps together the materials density of carrier states, D(E), and group velocity, ν(E), with the energy derivative of the Fermi-Dirac occupancy, f(E,c,T), where <sub>f</sub> is the Fermi level. The functions γ(E,T)=(E-E<sub>f</sub>)χ(E,T) and ζ(E,T)=(E-E<sub>f</sub>)<sup>2</sup> χ(E,T). These equations also express the relationship between the transport properties and Δ<sub>n</sub>, the moments of the distribution of conductivity over carriers with different energy, defined as </p>

<p align="center">
<img src="https://latex.codecogs.com/png.latex?%5Cinline%20%5Cdpi%7B300%7D%20%5CDelta_n%20%3D%20%5Cleft%5C%7B%5Cbegin%7Bmatrix%7D%20%5Cint%20%5Cchi%5Ctau%20dE%20%26%20n%3D0%5C%5C%20%5Cfrac%7B1%7D%7B%5CDelta_0%7D%5Cint%20E%5En%20%5Cchi%20%5Ctau%20dE%20%26%20n%5Cneq%200%20%5Cend%7Bmatrix%7D%5Cright." align="center" alt="drawing" width="400px">
</p>
The Seebeck coefficient obtains its largest magnitude by maximizing the asymmetry of product Dτν<sup>2</sup> about the Fermi level to move its center of current, Δ<sub>1</sub>, away from the Fermi level.

## Model Fermi Level

<div align="justify">
  
<p>
The Fermi level depends strongly on the carrier concentration, which varies non-monotonically with temperature as the solubility of the dopant changes. For a given carrier concentration, a self-consistent approach is developed to compute E<sub>f</sub> by setting the conduction band edge as the reference frame and computing E<sub>f</sub> that gives the same carrier population in DFT computed band and the given carrier population. This circumvents the problem that DFT underestimates the band gap. In this method Joyce and Dixon approximation of E<sub>f</sub> for degenerate semiconductors ((E<sub>f</sub>-E<sub>c</sub>)/k<sub>B</sub> ≅ln⁡[(n/N<sub>c</sub>)+1/(√8)]n/N<sub>c</sub> -(3/16-√3/9) (n/N<sub>c</sub>)<sup>2</sup>) is used as the initial guess. The E<sub>f</sub> iterates to meet the relation between charge density and density of state, n=∫<sub>E<sub>c</sub></sub>D(E)f(E)dE.
</p>

<p align="center">
<img src="Figures/Figure_2.png" align="center" alt="drawing" width="500px"/> 
</p>
The experimental measurements are noisy and so for the transport model,  the  carrier  concentration  was  represented  with the continuous smoothing function fit through the experimental data. This panel shows the Bayesian interpolation fitted to the temperature dependence of the experimentally measured carrier concentration.  
 
<p align="center">
<img src="Figures/Figure_3.png" align="center" alt="drawing" width="500px"/> 
</p>
In this pane, the Fermi level is plotted for different carrier concentrations using self-consistent method described above.


## Model Electron Lifetime

<div align="justify">

<p>
Semiconductor TEs are generally doped to beyond their saturation level. In these materials, strongly screened Columbic force induced by ionized impurities is the main source of scattering. The transition rate between initial and final energy states has S(E<sub>i</sub>,E<sub>f</sub>)=(2πN<sub>i</sub> e<sup>4</sup> L<sub>D</sub><sup>4</sup>)/((4πϵϵ<sub>o</sub> )<sup>2</sup>ℏΩ)δ(E<sub>f</sub>-E<sub>i</sub>). In this case, the electron lifetime is defined as 
</p>
 
<p align="center">
<img src="https://latex.codecogs.com/png.latex?%5Cinline%20%5Cdpi%7B300%7D%20%5Ctau_%7Bim%7D%28E%29%3D%5Cfrac%7B%5Chbar%7D%7B%5Cpi%20N_i%20%5Cleft%28%5Cfrac%7Be%5E2%20L_D%5E2%7D%7B4%5Cpi%20%5Cepsilon%20%5Cepsilon_o%7D%5Cright%29%5E2%20D%28E%29%7D" align="center" alt="drawing" width="350px"> 
</p>

<p>
For the strongly screened Columbic potential L<sub>D</sub> is small so that 1/(L<sub>D</sub><sup>4</sup>) is dominant. In doped semiconductors the Debye length has generalized form of 
</p>

<p align="center">
<img src="https://latex.codecogs.com/png.latex?%5Cinline%20%5Cdpi%7B300%7D%20L_D%3D%5Cfrac%7Be%5E2%20N_c%7D%7B4%20%5Cpi%20%5Cepsilon%20%5Cepsilon_o%20k_B%20T%7D%20%5Cleft%20%5BF_%7B-%5Cfrac%7B1%7D%7B2%7D%7D%28%5Ceta%29&plus;%5Cfrac%7B15%20%5Calpha%20k_B%20T%7D%7B4%7D%20F_%7B%5Cfrac%7B1%7D%7B2%7D%7D%28%5Ceta%29%5Cright%5D" align="center" alt="drawing" width="450px"> 
</p>

<p>
where N<sub>c</sub>=2((m<sub>c</sub> k<sub>B</sub>T)/(2πℏ)<sup>2</sup>)<sup>(3/2)</sup>. While the electron lifetime in equation serves reasonably well for many semiconductors, one should note two shortcomings of the Born approximation failures for slow moving electrons in Coulomb potential and deficiency of simply computing scattering from a single impurity and then multiplying it by number of impurities in capturing interference effects occur as electron wave propagate through random distribution of impurities in deriving this equation. We model the conduction band effective mass variation with temperature using m<sub>c</sub>(T)=m<sub>c</sub><sup>*</sup>(1+5αk<sub>B</sub> T). This model assumes linear dependency on temperature and does not count for degeneracy in high carrier population. 
</p>

<p align="center">
<img src="Figures/Figure_15.png" align="center" alt="drawing" width="500px"/> 
</p>

This panel shows the Debye length for two different carrier concentrations. The solid lines show the model prediction using degenerate form while the dash lines are for the cases in which degeneracy is neglected.

<p>
The second important scattering mechanism specially at high temperature in nonpolar semiconductors like Si is the acoustic phonon deformation potential. For electron phonon interaction, Ravich defined the lifetime as 
 
 <p align="center">
<img src="https://latex.codecogs.com/png.latex?%5Cinline%20%5Cdpi%7B300%7D%20%5Ctau_p%28E%29%3D%5Cfrac%7B%5Crho%20%5Cnu%5E2%20%5Chbar%7D%7B%5Cpi%20D_A%5E2%20k_B%20T%20D%28E%29%7D%20%5Cleft%20%28%20%5Cleft%5B1-%5Cfrac%7B%5Calpha%20E%7D%7B1&plus;2%5Calpha%20E%7D%20%5Cleft%281-%5Cfrac%7BD_v%7D%7BD_A%7D%20%5Cright%29%5Cright%5D%5E2-%5Cfrac%7B8%7D%7B3%7D%20%5Cfrac%7B%5Calpha%20E%281&plus;%20%5Calpha%20E%29%7D%7B%281&plus;2%20%5Calpha%20E%29%5E2%7D%5Cfrac%7BD_v%7D%7BD_A%7D%20%5Cright%29%5E%7B-1%7D" align="center" alt="drawing" width="800px"/> 
</p>
 
<p> 
This equation accounts for both absorption and emission of phonons. Note that the electron lifetime is strongly dominated by ion scattering and has weak dependency on phonon scattering. The other scattering terms of electron-electron and electron intervalley scattering has negligible importance in determining the electron lifetime and are excluded in calculations without loss of accuracy.
</p>

<p>
The rate of electron scattering due to the disordered arrangement in alloys dielectrics model as
</p>

 <p align="center">
<img src="https://latex.codecogs.com/png.latex?%5Cinline%20%5Cdpi%7B300%7D%20%5Ctau_a%28E%29%3D%5Cfrac%7B8%5Csqrt%7B2%7D%20%5Cpi%5E2%20%5Chbar%5E4%7D%7B0.75%20%5Cmathrm%7Bx%7D%281-%5Cmathrm%7Bx%7D%293a%5E3%20%5Cpi%5E3%20U_A%5E2%20m%5E%7B*%5E%7B%5Cfrac%7B3%7D%7B2%7D%7D%7D%20%5Csqrt%7BE%7D%7D" align="center" alt="drawing" width="500px"/> 
</p>
<p>
where x is the atomic fraction of alloy, a is the lattice parameter and the term U<sub>A</sub> is the alloy scattering potential.
</p>

<p align="center">
<img src="Figures/Picture2.jpg" align="center" alt="drawing" width="500px"/> 
</p>
<p>
This panel show the magnitude of electrical conductivity and Seebeck coefficient in phosphorus-doped bulk silicon. The solid blue line shows the model prediction for electrical conductivity, and the red line shows the prediction for the Seebeck coefficient. The experimentally measured σ and S are marked with open circles.
</p>

## Model Electron Transport Coefficients In Nanostructured Thermoelectrics

<p>
In the nanostructured of interest in this study, there are two additional electron scattering processes that arise as a result of the morphology: electron scattering at grain boundaries, and scattering from pores. The rate of electron momentum relaxation due to elastic scattering from a uniform dispersion of pores can be modeled as
</p>

<p align="center">
<img src="https://latex.codecogs.com/png.latex?%5Cinline%20%5Cdpi%7B300%7D%20%5Ctau_%7Bnp%7D%5E%7B-1%7D%20%3D%20%5Cfrac%7BN%7D%7B8%5Cpi%5E3%7D%20%5Cint%20SR_%7Bkk%27%7D%20%281-%5Ccos%28%5Ctheta_%7Bkk%27%7D%29%29%5Cdelta%5Cleft%28E-E%27%5Cright%29dk%27" align="center" alt="drawing" width="500px"/> 
</p>

<p>
Here N is the number density of pores, and the term SR<sub>kk'</sub> is the rate of transition of an electron from an initial state with wave vector k and energy E to a state k' with energy E' due to a single pore. For a time-invariant potential, the transition rate SR is given by Fermi’s golden rule as
</p>

<p align="center">
<img src="https://latex.codecogs.com/png.latex?%5Cinline%20%5Cdpi%7B300%7D%20SR_%7Bkk%27%7D%3D%5Cfrac%7B2%5Cpi%7D%7B%5Chbar%7D%5Cleft%7CM_%7Bkk%27%7D%5Cright%7C%5E2%5Cdelta%28E-E%27%29" align="center" alt="drawing" width="500px"/> 
</p>

<p>
where the matrix element operator M describes the strength which the pore couples the initial and final states and the number of ways the transition can occur.  For Bloch waves, M is given by the integral of the overlap of the initial and final state with the pore potential U(r) so that
<p> 

<p align="center">
<img src="https://latex.codecogs.com/png.latex?%5Cinline%20%5Cdpi%7B300%7D%20M_%7Bkk%27%7D%3D%20%5Cint%20e%5E%7Bi%28k%27-k%29.r%7D%20U%28r%29dr" align="center" alt="drawing" width="500px"/> 
</p>

For energy conservative (elastic) scattering between eigenstates with the same energy equation scattering can be recast as a surface integral over the isoenergetic k-space contour that satisfies E(k')=E(k) 

<p align="center">
<img src="https://latex.codecogs.com/png.latex?%5Cinline%20%5Cdpi%7B300%7D%20%5Ctau_%7Bnp%7D%5E%7B-1%7D%28s%29%20%3D%20%5Cfrac%7BN%7D%7B%282%5Cpi%29%5E2%5Chbar%7D%5Coint_%7B%5CGamma%7D%5Cfrac%7B%5Cleft%7CM_%7Bkk%27%7D%5Cright%7C%5E2%7D%7B%5Cnabla%20E%28k%27%29%7D%281-%5Ccos%5Ctheta%29dS%28k%27%29" align="center" alt="drawing" width="500px"/> 
</p>

<p>
where dS is the incremental area of the isoenergetic k-space surface. In most indirect bandgap semiconductors such as Si, the contours of isoenergy states near to conduction band valley have ellipsoidal shape in momentum space that can be approximated as
</p>

<p align="center">
<img src="https://latex.codecogs.com/png.latex?%5Cinline%20%5Cdpi%7B300%7D%20E%28k%29%3D%5Chbar%5E2%20%5B%28%5Cfrac%7B%28k_l-k_%7Bol%7D%20%29%5E2%7D%7B2m_l%5E*%7D%20&plus;%5Cfrac%7B%28k_t-k_%7Bot%7D%20%29%5E2%7D%7Bm_t%5E*%7D%5D" align="center" alt="drawing" width="500px"/> 
</p>

<p>
The pore potential, U(r), is assumed to be 
</p>

<p align="center">
<img src="https://latex.codecogs.com/png.latex?%5Cinline%20%5Cdpi%7B300%7D%20U%28r%29%20%3D%5Cleft%5C%7B%5Cbegin%7Bmatrix%7D%20U_o%20%26%20%5Ctext%7Bfor%20r%20inside%20the%20pore%7D%20%5C%5C%200%26%20%5Ctext%7Botherwise%7D%20%5Cend%7Bmatrix%7D%5Cright." align="center" alt="drawing" width="500px"/> 
</p>

<p>
where U<sub>o</sub> is the electron affinity. For an infinitely long cylindrical pores with radius r<sub>o</sub>, and aligned with axis parallel to z, this gives the scattering matrix element operator
</p>

<p align="center">
<img src="https://latex.codecogs.com/png.latex?%5Cinline%20%5Cdpi%7B300%7D%20M_%7Bkk%27%7D%5E%7Bcylinder%7D%3D2%5Cpi%20r_o%20U_o%20l_z%20%5Cleft%28%20%5Cfrac%7BJ_1%20%28r_o%20q_r%29%7D%7Bq_r%7D%20%5Cright%29%5Cdelta_k%28q_z%29" align="center" alt="drawing" width="500px"/> 
</p>

or for the spherical pores/ nanoparticles we have

<p align="center">
<img src="https://latex.codecogs.com/png.latex?%5Cinline%20%5Cdpi%7B300%7D%20M_%7Bkk%27%7D%5E%7Bspherical%7D%3D%5Cfrac%7B4%5Cpi%20U_o%7D%7Bq%5E2%7D%5Cleft%28%20%5Cfrac%7B1%7D%7Bq%7D%5Csin%28r_oq%29-r_o%5Ccos%28r_oq%29%5Cright%29" align="center" alt="drawing" width="500px"/> 
</p>

<p>
We remark that for nanoparticels the band alignment should be used instead of electron affinity. A similar use of Fermi's Golden rule can be used to model the rate of electron scattering by grain boundaries. Minnich et al have suggested that grain boundaries provide a scattering potential of magnitude U<sub>GB</sub> that decays away from the grain boundary over distance z<sub>o</sub>. From this, they derived the scattering operator matrix element for a small disc of grain boundary with radius r<sub>o</sub> as
</p>

<p align="center">
<img src="https://latex.codecogs.com/png.latex?%5Cinline%20%5Cdpi%7B300%7D%20M_%7Bkk%27%7D%3D4%5Cpi%20U_g%20%5Cleft%5B%20%5Cfrac%7Bz_o%7D%7B1&plus;%28q_zz_o%29%5E2%7D%20%5Cright%5Dr_o%5E2%5Cleft%5B%20%5Cfrac%7BJ_1%28q_rr_o%29%7D%7Bq_rr_o%7D%20%5Cright%5D" align="center" alt="drawing" width="500px"/> 
</p>



</div>


## Case Study: Si based TE nanocomposite

### Silicon band structure

<img src="Figures/Figure_1.png" align="center" alt="drawing" width="600px"/> 
<div align="justify">
<p>The terms D(E), and ν(E) for Si were derived from the conduction band of Si computed with density functional theory (DFT) using the Vienna Ab initio Simulation Package (VASP) using generalized gradient approximation (GGA) with the Perdew-Burke-Erzerhof exchange correlation functional (PBE). Projector augmented wave (PAW) pseudopotentials is used represent the ion cores. The Kohm-Sham wave functions constructed using a planewave basis set with 700 eV energy cutoff. The Brillouin zone was sampled using 12×12×12 Monkhorst-Pack k-point grid. The forces on the atoms minimized to better than 10<sup>-6</sup> eV/Å to relax the Si primitive cell. The electronic band structure used to compute D(E) on a 45×45×45 k-point grid. The group velocity was obtained from the conduction band curvature, ν=1/ℏ|∇<sub>κ</sub> E| along the〈100〉directions on the Γ to X Brillouin zone path.</p>
<img src="Figures/Figure_0.png" align="center" alt="drawing" width="600px"/> 
</div>

### Model prediction for bulk Si
<div align="justify">
  
<p>We have validated the transport model in bulk materials against a set of phosphorous-doped Si. The experimentally measured values are marked with open circles and the solid lines show the prediction of the semiclassical BTE model informed by the experimentally measured carrier population. The model gives a good fit to the experimental data across the full range of temperatures. </p>

<img src="Figures/Picture2.jpg" align="center" alt="drawing" width="600px"/> 
  
</div> 

### Maximum enhancement in power factor 
<div align="justify">
  
<p> This panel shows the model prediction for the variation in larges achievable TE power factor with carrier concentration in Si based porous materials with optimal characteristic lengths at 300 K. The best power factor performance using ideal filtering model is plotted in green. The power factor in bulk Si is plotted in black. In the narrow carrier concentration window with the highest power factor, bulk Si shows slightly better performance. The large energy difference between the conduction band edge in the Si and the vacuum level in the pore (about 4.15 eV electron affinity of bulk Si) causes strong electron scattering with countervailing response of reduction in the electrical conductivity and enhancement of Seebeck coefficient that cancel out each other, leading to an overall unchanged PF value. The maximum PF in porous structures takes place at carrier concentrations higher than the optimal carrier concentration in bulk Si. This is a key insight for the design of thermoelectrics at room temperature: If one is planning to engineer porous thermoelectrics to reduce phonon conduction, then one should also plan to increase the carrier concentration above the optimal level for the bulk semiconductor. In the Si model the maximum power factor of porous materials takes place at 8×10<sup>19</sup> cm<sup>-3</sup> carrier concentration and is slightly less than the maximum power factor in bulk Si occurs at 6.3x10<sup>19</sup> cm<sup>-3</sup> (~25% increase in doing concentration is needed for the best performance in porous Si).  </p>
<img src="Figures/Figure_18.png" align="center" alt="drawing" width="600px"/> 
  
</div> 

# CITATION
[1] Mitigating the Effect of Nanoscale Porosity on Thermoelectric Power Factor of Si, Hosseini, S. Aria and Romano, Giuseppe and Greaney, P. Alex, ACS Applied Energy Materials,2021, https://doi.org/10.1021/acsaem.0c02640.



