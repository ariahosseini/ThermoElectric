import numpy as np
from numpy.linalg import norm
from os.path import expanduser
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib as mpl
from matplotlib import cm
from numpy.matlib import repmat
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import matplotlib.ticker
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import ScalarFormatter
from mpl_toolkits import mplot3d
from matplotlib.colors import LightSource
import seaborn as sns

sns.set()
sns.set_context("paper", font_scale=2, rc={"lines.linewidth": 4})
sns.set_style("ticks", {"xtick.major.size": 2, "ytick.major.size": 2})

class thermoelectricProperties:

    hBar = 6.582119e-16     # Reduced Planck constant in eV.s
    kB = 8.617330350e-5     # Boltzmann constant in eV/K
    e2C = 1.6021765e-19     # e to Coulomb unit change
    e0 = 8.854187817e-12    # Permittivity in vacuum F/m
    Ang2meter = 1e-10       # Unit conversion from Angestrom to meter

    def __init__(self, latticeParameter, dopantElectricCharge, electronEffectiveMass, dielectric, numKpoints, numBands=None, numQpoints=None, electronDispersian=None, kpoints=None, energyMin=0, energyMax=2, numEnergySampling=1000):

        self.latticeParameter = latticeParameter            # Lattice parameter in A
        self.dopantElectricCharge = dopantElectricCharge
        self.electronEffectiveMass = electronEffectiveMass
        self.energyMax = energyMax                          # Maximum energy in eV
        self.energyMin = energyMin                          # Minimum energy in eV
        self.dielectric = dielectric                                        # Relative permittivity
        self.numEnergySampling = numEnergySampling          # Number of energy space samples to generate in eV
        self.numKpoints = numKpoints
        self.numBands = numBands
        self.electronDispersian = electronDispersian
        self.numQpoints = numQpoints

    def energyRange(self):                                  # Create an array of energy space sampling
        energyRange = np.linspace(self.energyMin, self.energyMax, self.numEnergySampling)
        return np.expand_dims(energyRange, axis=0)

    def kpoints(self, path2kpoints, delimiter=None, skiprows=0):
        kpoints = np.loadtxt(expanduser(path2kpoints), delimiter=None, skiprows=0)
        return kpoints

    def temp(self, TempMin=300, TempMax=1301, dT=100):
        temperature = np.arange(TempMin, TempMax, dT)
        return np.expand_dims(temperature, axis=0)

    def bandGap(self, Eg_o, Ao, Bo, Temp=None):
        if Temp is None:
            T = self.temp()
        else:
            T = Temp
        Eg = Eg_o - Ao * np.divide(T**2, T + Bo)
        return Eg

    def analyticalDoS(self, energyRange, alpha):
        DoS_nonparabolic = 1/np.pi**2*np.sqrt(2*energyRange*(1+energyRange*np.transpose(alpha)))*np.sqrt(self.electronEffectiveMass/thermoelectricProperties.hBar**2)**3*(1+(2*energyRange*np.transpose(alpha)))/thermoelectricProperties.e2C**(3./2)
        DoS_parabolic = np.sqrt(energyRange)/np.pi**2*np.sqrt(2)/thermoelectricProperties.hBar**3*self.electronEffectiveMass**(3/2)/self.e2C**(3/2)
        DoS = [DoS_nonparabolic,DoS_parabolic]
        return DoS

    def carrierConcentration(self, path2extrinsicCarrierConcentration, bandGap, Ao=None, Bo=None, Nc=None, Nv=None, Temp=None):
        if Temp is None:
            T = self.temp()
        else:
            T = Temp
        if Ao is None and Nc is None:
            raise Exception("Either Ao or Nc should be defined")
        if Bo is None and Nv is None:
            raise Exception("Either Bo or Nv should be defined")
        if Nc is None:
            Nc = Ao * Temp**(3. / 2)
        if Nv is None:
            Nv = Bo * Temp**(3. / 2)
        exCarrierFile = np.loadtxt(expanduser(path2extrinsicCarrierConcentration), delimiter=None, skiprows=0)
        extrinsicCarrierConcentration_tmp = InterpolatedUnivariateSpline(exCarrierFile[0, :], exCarrierFile[1, :] * 1e6)
        extrinsicCarrierConcentration = extrinsicCarrierConcentration_tmp(T)
        intrinsicCarrierConcentration = np.multiply(np.sqrt(np.multiply(Nc, Nv)), np.exp(-(np.divide(bandGap, (2 * thermoelectricProperties.kB * T)))))
        totalCarrierConcentration = intrinsicCarrierConcentration + abs(extrinsicCarrierConcentration)
        return totalCarrierConcentration

    def fermiLevel(self, carrierConcentration, energyRange, DoS, Nc=None, Ao=None, Temp=None):
        if Temp is None:
            T = self.temp()
        else:
            T = Temp
        if Ao is None and Nc is None:
            raise Exception("Either Ao or Nc should be defined")
        if Nc is None:
            Nc = Ao * Temp**(3. / 2)
        JD_CC = np.log(np.divide(carrierConcentration, Nc)) + 1 / np.sqrt(8) * np.divide(carrierConcentration, Nc) - (3. / 16 - np.sqrt(3) / 9) * np.power(np.divide(carrierConcentration, Nc), 2)
        fermiLevelEnergy = thermoelectricProperties.kB * np.multiply(T, JD_CC)
        f, _ = self.fermiDistribution(energyRange=energyRange, fermiLevel=fermiLevelEnergy, Temp=T)
        n = np.trapz(np.multiply(DoS, f), energyRange, axis=1)
        return [fermiLevelEnergy,np.expand_dims(n,axis=0)]

    def fermiDistribution(self, energyRange, fermiLevel, Temp=None):
        if Temp is None:
            T = self.temp()
        else:
            T = Temp

        xi = np.exp((energyRange-fermiLevel.T)/T.T/thermoelectricProperties.kB)
        fermiDirac = 1/(xi+1)
        dfdE = -1*xi/(1+xi)**2/T.T/thermoelectricProperties.kB
        fermi = np.array([fermiDirac, dfdE])
        return fermi

    def electronBandStructure(self, path2eigenval, skipLines):
        with open(expanduser(path2eigenval)) as eigenvalFile:
            for _ in range(skipLines):
                next(eigenvalFile)
            block = [[float(_) for _ in line.split()] for line in eigenvalFile]
        eigenvalFile.close()
        electronDispersian = [range(1, self.numBands + 1)]  # First line is atoms id
        kpoints = np.asarray(block[1::self.numBands + 2])[:, 0:3]
        for _ in range(self.numKpoints):
            binary2Darray = []
            for __ in range(self.numBands):
                binary2Darray = np.append(binary2Darray, block[__ + 2 + (self.numBands + 2) * _][1])
            electronDispersian = np.vstack([electronDispersian, binary2Darray])
        dispersian = [kpoints, electronDispersian]
        return dispersian

    def electronDoS(self, path2DoS, headerLines, numDoSpoints, unitcell_volume, valleyPoint, energyRange):
        DoS = np.loadtxt(expanduser(path2DoS), delimiter=None, skiprows=headerLines, max_rows=numDoSpoints)
        valleyPointEnergy = DoS[valleyPoint, 0]
        DoSSpline = InterpolatedUnivariateSpline(DoS[valleyPoint:, 0] - valleyPointEnergy, DoS[valleyPoint:, 1] / unitcell_volume)
        DoSFunctionEnergy = DoSSpline(energyRange)  # Density of state
        return DoSFunctionEnergy

    def fermiLevelSelfConsistent(self, carrierConcentration, Temp, energyRange, DoS, fermilevel):
        fermi = np.linspace(fermilevel[0]-0.2, fermilevel[0]+0.2, 1000, endpoint=True).T
        result_array = np.empty((np.shape(Temp)[1], np.shape(fermi)[1]))
        idx_j = 0
        for j in Temp[0]:
            idx_i = 0
            for i in fermi[idx_j]:
                f, _ = self.fermiDistribution(energyRange=energyRange, fermiLevel=np.expand_dims(np.array([i]), axis=0), Temp=np.expand_dims(np.array([j]), axis=0))
                tmp = np.trapz(np.multiply(DoS, f), energyRange, axis=1)
                result_array[idx_j, idx_i] = tmp
                idx_i += 1
            idx_j += 1
        diff = np.tile(np.transpose(carrierConcentration), (1, np.shape(fermi)[1])) - abs(result_array)
        min_idx = np.argmin(np.abs(diff), axis=1)
        Ef = np.empty((1, np.shape(Temp)[1]))
        for Ef_idx in np.arange(len(min_idx)):
            Ef[0,Ef_idx] = fermi[Ef_idx,min_idx[Ef_idx]]
        elm = 0
        n = np.empty((1, np.shape(Temp)[1]))
        for idx in min_idx:
            n[0,elm] = result_array[elm, idx]
            elm += 1
        return [Ef,n]

    def electronGroupVelocity(self, kp, energy_kp, energyRange):
        dE = np.roll(energy_kp, -1, axis=0) - np.roll(energy_kp, 1, axis=0)
        dk = np.roll(kp, -1, axis=0) - np.roll(kp, 1, axis=0)
        dEdk = np.divide(dE, dk)
        dEdk[0] = (energy_kp[1] - energy_kp[0]) / (kp[1] - kp[0])
        dEdk[-1] = (energy_kp[-1] - energy_kp[-2]) / (kp[-1] - kp[-2])
        dEdkSpline = InterpolatedUnivariateSpline(energy_kp, np.array(dEdk))
        dEdkFunctionEnergy = dEdkSpline(energyRange)
        groupVel = dEdkFunctionEnergy / thermoelectricProperties.hBar
        return groupVel

    def matthiessen(self, *args):
        tau = 1. / sum([1. / arg for arg in args])
        tau[np.isinf(tau)] = 0
        return tau

    def tau_p(self, energyRange, alpha, Dv, DA, T, vs, D, rho):

        nonparabolic_term = (1-(alpha.T*energyRange)/(1+2*alpha.T*energyRange)*(1-Dv/DA))**2-8/3*(alpha.T*energyRange)/(1+2*alpha.T*energyRange)*(Dv/DA)
        tau = rho*vs**2*thermoelectricProperties.hBar/np.pi/thermoelectricProperties.kB/T.T/DA/DA*1e9/thermoelectricProperties.e2C/D
        tau_p = tau/nonparabolic_term
        return [tau,tau_p]

    def tau_ion(self,energyRange, LD, N):

        g = 8*self.electronEffectiveMass*LD.T**2*energyRange/thermoelectricProperties.hBar**2
        var_tmp = np.log(1+g)-g/(1+g)
        tau = 16*np.pi*np.sqrt(2*self.electronEffectiveMass)*(4*np.pi*self.dielectric*thermoelectricProperties.e0)**2/N.T*energyRange**(3/2)*thermoelectricProperties.e2C**(-5/2)
        print(tau)
        return tau

    def electricalProperties(self, E, DoS, vg, Ef, dfdE, Temp, tau):
        X = DoS * vg**2 * dfdE
        Y = (E - np.transpose(Ef)) * X
        Z = (E - np.transpose(Ef)) * Y
        Sigma = -1 * np.trapz(X * tau, E, axis=1) / 3 * thermoelectricProperties.e2C
        S = -1*np.trapz(Y * tau, E, axis=1)/np.trapz(X * tau, E, axis=1)/Temp
        PF = Sigma*S**2
        ke = -1*(np.trapz(Z * tau, E, axis=1) - np.trapz(Y * tau, E, axis=1)**2/np.trapz(X * tau, E, axis=1))/Temp/3 * thermoelectricProperties.e2C
        delta_0 = np.trapz(X * tau* E, E, axis=1)
        delta_1 = np.trapz(X * tau* E, E, axis=1)/ np.trapz(X * tau, E, axis=1)
        delta_2 = np.trapz(X * tau* E**2, E, axis=1)/ np.trapz(X * tau, E, axis=1)
        Lorenz = (delta_2-delta_1**2)/Temp/Temp
        coefficients = [Sigma, S[0], PF[0], ke[0], delta_1, delta_2, Lorenz[0]]
        return coefficients

    def filteringEffect(self, U0, tau0, tauOff, energyRange, electronBandStructure, temp, electronDoS, electronGroupVelocity, bandGap, carrierConcentration, fermiLevel, fermiDistribution, factor, q, uIncrement=0.05, tauIncrement=1e-15, tempIndex=0):
        n = 0
        m = 0
        sMatrix = np.array([])
        # rMatrix = np.array([])
        # pfMatrix = np.array([])

        for _ in np.arange(0.1, U0, uIncrement):
            m += 1
            for __ in np.arange(tauIncrement, tau0, tauIncrement):
                tauInc = np.ones(len(Erange))
                tauInc[np.where(Erange < _)] = __
                tau_on = self.matthiessen(energyRange, tauOff, tauInc)
                s, r, c, pf, X, Y = self.electricalProperties(energyRange, electronBandStructure, temp, electronDoS, electronGroupVelocity, tau_on, bandGap, carrierConcentration, fermiLevel, fermiDistribution, factor, q)
                sMatrix = np.append(sMatrix, s[tempIndex])
                n += 1
        sMatrix = np.reshape(sMatrix, (n, m))
        return sMatrix

    # def qpoints(self):
    #     qpoints = np.array([np.zeros(self.numQpoints), np.zeros(self.numQpoints), np.linspace(-math.pi / self.latticeParameter, math.pi / self.latticeParameter, num=self.numQpoints)])
    #     return qpoints

    def dynamicalMatrix(self, path2massWeightedHessian, path2atomsPositions, skipLines, numAtoms, baseLatticePoint, numAtomsInUnitCell, qpoints):
        with open(os.path.expanduser(path2massWeightedHessian)) as hessianFile:
            hessianMatrix = hessianFile.readlines()
        hessianMatrix = [line.split() for line in hessianMatrix]
        hessianMatrix = np.array([[float(_) for _ in __] for __ in hessianMatrix])
        hessianFile.close()
        hessianSymmetry = (np.triu(hessianMatrix) + np.tril(hessianMatrix).transpose()) / 2
        hessianMatrix = hessianSymmetry + np.triu(hessianSymmetry, 1).transpose()
        #
        with open(os.path.expanduser(path2atomsPositions)) as atomsPositionsFile:
            atomsPositions = atomsPositionsFile.readlines()
        atomsPositions = [line.split() for line in atomsPositions]
        [atomsPositions.pop(0) for _ in range(skipLines)]
        atomsPositions = np.array([[float(_) for _ in __] for __ in atomsPositions[0:numAtoms]])
        atomsPositions = atomsPositions[atomsPositions[:, 0].argsort()]
        # atomsPositions = np.sort(atomsPositions.view('i8,i8,f8,f8,f8,f8,f8,f8'), order=['f0'], axis=0).view(np.float)
        latticePoints = np.array([_[2:5] for _ in atomsPositions[::numAtomsInUnitCell]])
        latticePointsVectors = latticePoints - numpy.matlib.repmat(latticePoints[baseLatticePoint], len(latticePoints), 1)
        dynamicalMatrix = np.zeros((numAtomsInUnitCell * 3, numAtomsInUnitCell * 3))
        for _ in range(self.numQpoints):
            dynamMatPerQpoint = np.zeros((numAtomsInUnitCell * 3, numAtomsInUnitCell * 3))
            for __ in range(len(latticePointsVectors)):
                sumMatrix = hessianMatrix[__ * numAtomsInUnitCell * 3: (__ + 1) * numAtomsInUnitCell * 3, baseLatticePoint * numAtomsInUnitCell * 3: (baseLatticePoint + 1) * numAtomsInUnitCell * 3] * cmath.exp(-1j * np.dot(latticePointsVectors[__], qpoints[:, _]))
                dynamMatPerQpoint = dynamMatPerQpoint + sumMatrix
            dynamicalMatrix = np.append(dynamicalMatrix, dynamMatPerQpoint, axis=0)
        dynamicalMatrix = dynamicalMatrix[numAtomsInUnitCell * 3:]
        eigVal = np.array([])
        eigVec = np.zeros((numAtomsInUnitCell * 3, numAtomsInUnitCell * 3))
        for _ in range(self.numQpoints):
            dynmat = dynamicalMatrix[_ * numAtomsInUnitCell * 3:(_ + 1) * numAtomsInUnitCell * 3]
            eigvals, eigvecs, = np.linalg.eigh(dynmat)
            eigVal = np.append(eigVal, eigvals).reshape(-1, numAtomsInUnitCell * 3)
            eigVec = np.append(eigVec, eigvecs, axis=0)
        eigVec = eigVec[numAtomsInUnitCell * 3:]
        frequencies = np.sqrt(np.abs(eigVal.real)) * np.sign(eigVal.real)
        # conversion_factor_to_THz = 15.633302
        # frequencies = frequencies * conversion_factor_to_THz
        return eigVec

    def phonopyQpointYamlInterface(self, path2QpointYaml):
        qpointsData = yaml.load(open("qpoints.yaml"))
        nqpoint = qpointsData['nqpoint']
        natom = qpointsData['natom']
        qpoints = []
        qpoints = np.append(qpoints, [qpointsData['phonon'][_]['q-position'] for _ in range(nqpoint)]).reshape(-1, 3)
        frequency = []
        frequency = np.append(frequency, [[qpointsData['phonon'][_]['band'][__]['frequency'] for __ in range(3 * natom)] for _ in range(nqpoint)]). reshape(-1, 3 * natom)
        eigVal = np.array([])
        eigVec = np.zeros((natom * 3, natom * 3))
        for _ in range(nqpoint):
            dynmat = []
            dynmat_data = qpointsData['phonon'][_]['dynamical_matrix']
            for row in dynmat_data:
                vals = np.reshape(row, (-1, 2))
                dynmat.append(vals[:, 0] + vals[:, 1] * 1j)
            dynmat = np.array(dynmat)
            eigvals, eigvecs, = np.linalg.eigh(dynmat)
            eigVal = np.append(eigVal, eigvals).reshape(-1, natom * 3)
            eigVec = np.append(eigVec, eigvecs, axis=0)
        eigVec = eigVec[natom * 3:]
        frequencies = np.sqrt(np.abs(eigVal.real)) * np.sign(eigVal.real)
        conversion_factor_to_THz = 15.633302
        frequencies = frequencies * conversion_factor_to_THz
        return eigVec

    def gaussianDestribution(self, sigma, expectedValue, qpoints):
        gauss = (1.0 / np.sqrt(2 * pi) / sigma) * np.exp((-1.0 / 2) * np.power(((qpoints - expectedValue) / sigma), 2))
        return gauss

    # def singleWave(self, path2atomsPosition, numberOfAtomsInChain, skipLines)
    # def __str(self):

    # def __repr(self):


me = 9.109e-31
Si = thermoelectricProperties(latticeParameter=5.401803661945516e-10, dopantElectricCharge=1, electronEffectiveMass=1.08*me, energyMin=0.0, energyMax=2, dielectric=11.7, numKpoints=800, numBands=8, numQpoints=201, numEnergySampling=5000)


ml = 0.98*me
mt = 0.19*me
m_CB = 3/(1/ml+2/mt)


Lv = np.array([[1,1,0],[0,1,1],[1,0,1]])*Si.latticeParameter/2
a_rp = np.cross(Lv[1],Lv[2])/np.dot(Lv[0],np.cross(Lv[1],Lv[2]))
b_rp = np.cross(Lv[2],Lv[0])/np.dot(Lv[1],np.cross(Lv[2],Lv[0]))
a_rp = np.cross(Lv[0],Lv[1])/np.dot(Lv[2],np.cross(Lv[0],Lv[1]))
RLv = np.array([a_rp, b_rp, a_rp])


e = Si.energyRange()
g = Si.temp(TempMin=300, TempMax=1301, dT=100)
h = Si.bandGap(Eg_o=2, Ao=7.7e-4, Bo=600, Temp=g)
alpha = (1-m_CB/me)**2/h
dos_nonparabolic, dos_parabolic = Si.analyticalDoS(energyRange=e, alpha = alpha)
cc = Si.carrierConcentration(Nc=None, Nv=None, path2extrinsicCarrierConcentration='experimental-carrier-concentration-no-inc.txt', bandGap=h, Ao=5.3e21, Bo=3.5e21, Temp=g)
kp, band = Si.electronBandStructure(path2eigenval='EIGENVAL', skipLines=6)
kp_rl = 2*np.pi*np.matmul(kp,RLv)
kp_mag = norm(kp_rl, axis=1)
min_band = np.argmin(band[400:600, 4], axis=0)
max_band = np.argmax(band[400:600, 4], axis=0)
kp_vel = kp_mag[401 + max_band:401 + min_band]
energy_vel = band[401 + max_band:401 + min_band, 4] - band[401 + min_band, 4]
enrg_sorted_idx = np.argsort(energy_vel, axis=0)
gVel = Si.electronGroupVelocity(kp=kp_vel[enrg_sorted_idx], energy_kp=energy_vel[enrg_sorted_idx], energyRange=e)
DoS = 1/2*Si.electronDoS(path2DoS='DOSCAR', headerLines=6, unitcell_volume=19.70272e-30, numDoSpoints=2000, valleyPoint=1118, energyRange=e)
JD_f, JD_n = Si.fermiLevel(carrierConcentration=cc, energyRange=e, DoS= DoS, Nc=None, Ao=5.3e21, Temp=g)
fermi, cc_sc = Si.fermiLevelSelfConsistent(carrierConcentration=cc, Temp=g, energyRange=e, DoS=DoS, fermilevel=JD_f)
dis, dfdE = Si.fermiDistribution(energyRange=e, Temp=g, fermiLevel=fermi)

bulk_module = 98
rho = 2329
sp = np.sqrt(bulk_module/rho)

LD = np.sqrt(4*np.pi*Si.dielectric*thermoelectricProperties.e0*thermoelectricProperties.kB/thermoelectricProperties.e2C*g/cc)

tau_p = 220e-15/np.sqrt(np.transpose(g) * e)
tau_p_npb, tau_p_pb = Si.tau_p(energyRange=e, alpha=alpha, Dv=2.94, DA=-9.5, T=g, vs=sp, D=DoS, rho=rho)
tau_p_npb_type_2, tau_p_pb_type_2 = Si.tau_p(energyRange=e, alpha=alpha, Dv=2.94, DA=-9.5, T=g, vs=sp, D=dos_nonparabolic, rho=rho)
tau_p_npb_type_3, tau_p_pb_type_3 = Si.tau_p(energyRange=e, alpha=alpha, Dv=2.94, DA=-9.5, T=g, vs=sp, D=dos_parabolic, rho=rho)
# tau_ion = Si.tau_ion(energyRange=e, LD= LD, N= cc)

tau_p2 = 2500e-15/np.sqrt(e)/np.transpose(g)
tau = Si.matthiessen(e, tau_p)

Coeff = Si.electricalProperties(E=e, DoS=DoS, vg=gVel, Ef=fermi, dfdE=dfdE, Temp=g, tau=tau_p)

exp_r = np.loadtxt('f1_Resistivity')
exp_s = np.loadtxt('f1_Seebeck')

print("done")
fig_0 = plt.figure(figsize=(6.5,4.5))
ax_0 = fig_0.add_subplot(111)
ax_0.plot(g[0],h[0], 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_0.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax_0.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_0.tick_params(axis="x", labelsize=16)
ax_0.set_ylabel('Band gap (eV)', fontsize=16, labelpad=10)
ax_0.tick_params(axis="y", labelsize=16)
fig_0.tight_layout()

fig_1 = plt.figure(figsize=(6.5,4.5))
ax_1 = fig_1.add_subplot(111)
ax_1.plot(g[0],alpha[0], 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_1.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax_1.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_1.tick_params(axis="x", labelsize=16)
ax_1.set_ylabel('Nonparabolic term (eV$^{-1}$)', fontsize=16, labelpad=10)
ax_1.tick_params(axis="y", labelsize=16)
fig_1.tight_layout()

fig_2 = plt.figure(figsize=(6.5,4.5))
ax_2 = fig_2.add_subplot(111)
ax_2.plot(e[0],DoS[0], 'None', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_2.plot(e[0],dos_nonparabolic[0], 'None', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)
ax_2.plot(e[0],dos_parabolic[0], 'None', linestyle='-', color='olive',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='olive',
          markeredgewidth=1)
ax_2.yaxis.set_major_formatter(ScalarFormatter())
ax_2.set_xlabel('Energy (eV)', fontsize=16, labelpad=10)
ax_2.tick_params(axis="x", labelsize=16)
ax_2.set_ylabel('Density of state (#/eV/m$^3$)', fontsize=16, labelpad=10)
ax_2.tick_params(axis="y", labelsize=16)
ax_2.ticklabel_format(axis="y", style="sci", scilimits=None)
fig_2.tight_layout()

fig_3 = plt.figure(figsize=(6.5,4.5))
ax_3 = fig_3.add_subplot(111)
ax_3.plot(e[0],dos_nonparabolic[0], 'None', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_3.plot(e[0],dos_nonparabolic[-1], 'None', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)
ax_3.yaxis.set_major_formatter(ScalarFormatter())
ax_3.set_xlabel('Energy (eV)', fontsize=16, labelpad=10)
ax_3.tick_params(axis="x", labelsize=16)
ax_3.set_ylabel('Density of state (#/eV/m$^3$)', fontsize=16, labelpad=10)
ax_3.tick_params(axis="y", labelsize=16)
ax_3.ticklabel_format(axis="y", style="sci", scilimits=None)
fig_3.tight_layout()

fig_4 = plt.figure(figsize=(6.5,4.5))
ax_4 = fig_4.add_subplot(111)
ax_4.plot(e[0],-1*gVel[0]*1e-5, 'None', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_4.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax_4.set_xlabel('Energy (eV)', fontsize=16, labelpad=10)
ax_4.tick_params(axis="x", labelsize=16)
ax_4.set_ylabel('Group velocity (x10$^5$ m/s)', fontsize=16, labelpad=10)
fig_4.tight_layout()

fig_5 = plt.figure(figsize=(6.5,4.5))
ax_5 = fig_5.add_subplot(111)
ax_5.plot(band[1::,], 'None', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_5.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax_5.set_xlabel('BZ tour', fontsize=16, labelpad=10)
ax_5.tick_params(axis="x", labelsize=16)
ax_5.set_ylabel('Energy (eV)', fontsize=16)
ax_5.tick_params(axis="y", labelsize=16)
ax_5.set_xticks([0,199,399,599,799])
ax_5.set_xticklabels(["W", "L","$\Gamma$", "X", "W"])
fig_5.tight_layout()

fig_6 = plt.figure(figsize=(6.5,4.5))
ax_6 = fig_6.add_subplot(111)
ax_6.plot(e[0],dis[0], 'None', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_6.plot(e[0],dis[-1], 'None', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)
ax_6.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax_6.set_xlabel('Energy (eV)', fontsize=16, labelpad=10)
ax_6.tick_params(axis="x", labelsize=16)
ax_6.set_ylabel('Fermi distribution', fontsize=16, labelpad=10)
ax_6.tick_params(axis="y", labelsize=16)
fig_6.tight_layout()

fig_7 = plt.figure(figsize=(6.5,4.5))
ax_7 = fig_7.add_subplot(111)
ax_7.plot(e[0],dfdE[0], 'None', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_7.plot(e[0],dfdE[-1], 'None', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)
ax_7.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax_7.set_xlabel('Energy (eV)', fontsize=16, labelpad=10)
ax_7.tick_params(axis="x", labelsize=16)
ax_7.set_ylabel('Fermi window (eV$^{-1}$)', fontsize=16)
ax_7.tick_params(axis="y", labelsize=16)
fig_7.tight_layout()

fig_8 = plt.figure(figsize=(6.5,4.5))
ax_8 = fig_8.add_subplot(111)
ax_8.plot(g[0],cc[0], 'o', linestyle='None', color='black',
          markersize=12, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='black',
          markeredgewidth=1,zorder=0)
ax_8.plot(g[0],JD_n[0], 'o', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)

ax_8.plot(g[0],cc_sc[0], 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)

ax_8.yaxis.set_major_formatter(ScalarFormatter())
ax_8.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_8.tick_params(axis="x", labelsize=16)
ax_8.set_ylabel('Carrier concentration( #/m$^{3}$)', fontsize=16)
ax_8.tick_params(axis="y", labelsize=16)
ax_8.ticklabel_format(axis="y", style="sci", scilimits=None)
fig_8.tight_layout()

fig_9 = plt.figure(figsize=(6.5,4.5))
ax_9 = fig_9.add_subplot(111)
ax_9.plot(g[0],fermi[0], 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1,zorder=10)
ax_9.plot(g[0],JD_f[0], 'o', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)

ax_9.yaxis.set_major_formatter(ScalarFormatter())
ax_9.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_9.tick_params(axis="x", labelsize=16)
ax_9.set_ylabel('E$_f$ (eV)', fontsize=16)
ax_9.tick_params(axis="y", labelsize=16)
ax_9.ticklabel_format(axis="y", style="sci", scilimits=None)
fig_9.tight_layout()

fig_10 = plt.figure(figsize=(6.5,4.5))
ax_10 = fig_10.add_subplot(111)
ax_10.plot(g[0],Coeff[0]*1e-5, 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_10.plot(exp_r[0], np.divide(1, exp_r[1]) * 1e-5, '-o', linestyle='None', color='black',
        markersize=5, linewidth=4,
        markerfacecolor='white',
        markeredgecolor='black',
        markeredgewidth=1)
ax_10.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax_10.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_10.tick_params(axis="x", labelsize=16)
ax_10.set_ylabel('Conductivity(x10$^5$ S/m)', fontsize=16)
ax_10.tick_params(axis="y", labelsize=16)
fig_10.tight_layout()

fig_11 = plt.figure(figsize=(6.5,4.5))
ax_11 = fig_11.add_subplot(111)
ax_11.plot(g[0],Coeff[1]*1e6, 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_11.plot(exp_s[0], exp_s[1] * 1e6, '-o', linestyle='None', color='black',
          markersize=5, linewidth=4,
          markerfacecolor='white',
          markeredgecolor='black',
          markeredgewidth=1)

ax_11.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax_11.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_11.tick_params(axis="x", labelsize=16)
ax_11.set_ylabel('Seebeck($\mu$V/K)', fontsize=16)
ax_11.tick_params(axis="y", labelsize=16)
fig_11.tight_layout()

fig_12 = plt.figure(figsize=(6.5,4.5))
ax_12 = fig_12.add_subplot(111)
ax_12.plot(g[0],Coeff[2]*1e3, 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_12.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax_12.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_12.tick_params(axis="x", labelsize=16)
ax_12.set_ylabel('Power factor(mW/mK$^2$)', fontsize=16, labelpad=10)
ax_12.tick_params(axis="y", labelsize=16)
fig_12.tight_layout()


fig_13 = plt.figure(figsize=(6.5,4.5))
ax_13 = fig_13.add_subplot(111)
ax_13.plot(g[0],Coeff[3], 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)

ax_13.plot(g[0],Coeff[0]*g[0]*Coeff[6], 'o', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)

ax_13.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax_13.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_13.tick_params(axis="x", labelsize=16)
ax_13.set_ylabel('$\kappa_e$(W/mK)', fontsize=16, labelpad=10)
ax_13.tick_params(axis="y", labelsize=16)
fig_13.tight_layout()

fig_14 = plt.figure(figsize=(6.5,4.5))
ax_14 = fig_14.add_subplot(111)
ax_14.plot(g[0],Coeff[4], 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_14.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax_14.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_14.tick_params(axis="x", labelsize=16)
ax_14.set_ylabel('$\Delta_1$(eV)', fontsize=16, labelpad=10)
ax_14.tick_params(axis="y", labelsize=16)
fig_14.tight_layout()

fig_15 = plt.figure(figsize=(6.5,4.5))
ax_15 = fig_15.add_subplot(111)
ax_15.plot(g[0],Coeff[5], 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_15.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax_15.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_15.tick_params(axis="x", labelsize=16)
ax_15.set_ylabel('$\Delta_2$([eV]$^2$)', fontsize=16, labelpad=10)
ax_15.tick_params(axis="y", labelsize=16)
fig_15.tight_layout()


fig_16 = plt.figure(figsize=(6.5,4.5))
ax_16 = fig_16.add_subplot(111)
ax_16.plot(g[0],Coeff[6]*1e8, 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)

ax_16.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax_16.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_16.tick_params(axis="x", labelsize=16)
ax_16.set_ylabel('Lorenz number (x$10^{-8}$[V/K]$^2$)', fontsize=16, labelpad=10)
ax_16.tick_params(axis="y", labelsize=16)
fig_16.tight_layout()

plt.show()
exit()


# Qpoint = np.array([np.zeros(silicon.numQpoints), np.zeros(silicon.numQpoints), np.linspace(-math.pi / silicon.latticeParameter, math.pi / silicon.latticeParameter, num=silicon.numQpoints)])
# print len(Qpoint[1])
# dynamicalMatrix = thermoelectricProperties.dynamicalMatrix(silicon, '~/Desktop/Notes/Box_120a_Lambda_10a/Si-hessian-mass-weighted-hessian.d', '~/Desktop/Notes/Box_120a_Lambda_10a/data.Si-3x3x3', 15, 216, 14, 8, Qpoint)
# print dynamicalMatrix
# print(dynamicalMatrix)
# Eig = thermoelectricProperties.phonopyQpointYamlInterface(silicon, '~/Desktop/qpoints.yaml')
# np.savetxt('EigVec', Eig.real, fmt='%10.5f', delimiter=' ')
print('done')
