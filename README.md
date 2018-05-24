# thermoelectric.py

__init__(self, latticeParameter, dopantElectricCharge, electronEffectiveMass, energyMin, energyMax, kS, numKpoints, numBands, numQpoints, electronDispersian=None, kpoints=None, numEnergySampling=1e3)

energyRange(self): Energy range in conduction band, generally around 50 meV
electronBandStructure(self, path2eigenval, skipLines): Electron band structure from VASP Eigenval file
kpoints(self, path2kpoints): Electron Kpoints
temp(self, TempMin, TempMax, dT=10): Temperature range
bandGap(self, temp, bandGapEquation): Band gap as a function of temperature
carrierConcentration(self, effectiveDoSConductionBand, effectiveDoSValanceBand, path2extrinsicCarrierConcentration, bandGap, temp): Carrier concentration approximation from solid state formula
fermiLevel(self, effectiveDoSConductionBand, temp, carrierConcentration): Electron Fermi level
fermiDistribution(self, energyRange, temp, fermiLevel): Fermi-Dirac distribution
electronDoS(self, path2DoS, headerLines, numDoSpoints, valleyPoint, energyRange): Electron DoS from Vasp DOS file 
fermiLevelSelfConsistent(self, carrierConcentration, electronDoS, fermiDistribution, temp, energyRange, initialGuess=-0.1): Electron Fermi level from self consistant formula
electronGroupVelocity(self, electronBandStructure, conductionBandIndex, kpoints, KpathInitialpoint, KpathLastpoint, energyRange): Electron group velocity
matthiessen(self, energyRange, *args): Matthiessen rule to add up relaxation time
electricalProperties(self, energyRange, electronBandStructure, temp, electronDoS, electronGroupVelocity, matthiessen, bandGap, carrierConcentration, fermiLevel, fermiDistribution, factor, q): Seebeck and resistivity calculator from BTE formula
filteringEffect(self, U0, tau0, tauOff, energyRange, electronBandStructure, temp, electronDoS, electronGroupVelocity, bandGap, carrierConcentration, fermiLevel, fermiDistribution, factor, q, uIncrement=0.05, tauIncrement=1e-15, tempIndex=0): Relaxation time from electron selective scattering
qpoints(self): Phonon Qpoints
dynamicalMatrix(self, path2massWeightedHessian, path2atomsPositions, skipLines, numAtoms, baseLatticePoint, numAtomsInUnitCell, qpoints): Phonon Dynamical matrix from LAMMPS output file
phonopyQpointYamlInterface(self, path2QpointYaml): Extractin Phonopy outputs
gaussianDestribution(self, sigma, expectedValue, qpoints): Produce Gaussian distribution
singleWave(self, path2atomsPosition, numberOfAtomsInChain, skipLines): Phonon sigle wave
