## USER GUIDE

<div align="justify">
  
<p>
  
Following examples show how to use *Thermoelectric.py* to compute the ideal filtering in Si-based TEs and the electrical properties of Si composite containg nanoscale SiC particles. You need to include thermoelectricProperties.py and accum.py to your PYTHONPATH.
</p>

## DFT CALCULATIONS
  
<p>  
  
Finally, the Si band structure is computed The terms D(E), and ν(E) for Si were derived from the conduction band of Si computed with DFT using the Vienna Ab initio Simulation Package (VASP) using generalized gradient approximation (GGA) with the Perdew-Burke-Erzerhof exchange correlation functional (PBE). Projector augmented wave (PAW) pseudopotentials is used represent the ion cores. The Kohm-Sham wave functions constructed using a planewave basis set with 700 eV energy cutoff. The Brillouin zone was sampled using 12 by 12 by12 Monkhorst-Pack k-point grid. The forces on the atoms minimized to better than 10<sup>-6</sup> eV/Å to relax the Si primitive cell. The electronic band structure used to compute D(E) on a 45 by 45 by 45 k-point grid. The group velocity was obtained from the conduction band curvature, along the ⟨001⟩ directions on the Γ to X Brillouin zone path. 
</p>
