# Required Libs

"""
Thermoelectric.py is a computational framework that computes electron transport coefficients with unique features to design the nanoscale morphology of thermoelectrics (TEs) to obtain electron scattering that will enhance performance through electron energy filtering.

Author: S. Aria Hosseini
Email: shoss008@ucr.edu
"""

import numpy as np
from numpy.linalg import norm
from os.path import expanduser
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import PchipInterpolator
from scipy.special import jv
from accum import accum # This pakage is used to compute the cumulative average
import copy

# Post processing packages

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



class thermoelectricProperties:

    hBar = 6.582119e-16     # Reduced Planck constant in eV.s
    kB = 8.617330350e-5     # Boltzmann constant in eV/K
    e2C = 1.6021765e-19     # e to Coulomb unit change
    e0 = 8.854187817e-12    # Permittivity in vacuum F/m
    Ang2meter = 1e-10       # Unit conversion from Angestrom to meter
    me = 9.109e-31          # Electron rest mass in Kg

    def __init__(self, latticeParameter, dopantElectricCharge, electronEffectiveMass, dielectric, numKpoints, numBands=None, numQpoints=None, electronDispersian=None, kpoints=None, energyMin=0, energyMax=2, numEnergySampling=1000):

        self.latticeParameter = latticeParameter                # Lattice parameter in m
        self.dopantElectricCharge = dopantElectricCharge        # Impurity charge, generally equal to 1 for P dopand
        self.electronEffectiveMass = electronEffectiveMass      # Electron effective mass
        self.energyMax = energyMax                              # Maximum energy in eV, 1-2 eV sounds reasonable, defaul is 2
        self.energyMin = energyMin                              # Minimum energy in eV, 0 eV unless looking for high level contribution
        self.dielectric = dielectric                            # Relative permittivity
        self.numEnergySampling = numEnergySampling              # Number of energy space samples to generate, defaul is 1000
        self.numKpoints = numKpoints                            # Number of Kpints in DFT calculation
        self.numBands = numBands
        self.electronDispersian = electronDispersian
        self.numQpoints = numQpoints

    def energyRange(self):                                      # Create an 2D array of energy space sampling

        energyRange = np.linspace(self.energyMin, self.energyMax, self.numEnergySampling)
        return np.expand_dims(energyRange, axis=0)              # The array size is [1,numEnergySampling]

    def kpoints(self, path2kpoints, delimiter=None, skiprows=0):

        kpoints = np.loadtxt(expanduser(path2kpoints), delimiter=None, skiprows=0)
        return kpoints

    def temp(self, TempMin=300, TempMax=1301, dT=100):          # Create an 2D array of temperature sampling

        temperature = np.arange(TempMin, TempMax, dT)
        return np.expand_dims(temperature, axis=0)              # The array size is [1,(TempMax-TempMin)/dT]

    def bandGap(self, Eg_o, Ao, Bo, Temp=None):                 # general shap of bandgap: Eg(T)=Eg(T=0)-Ao*T**2/(T+Bo)
        if Temp is None:    # Use defaul temperature range (300K up to 1301K with step of 50K)
            T = self.temp()
        else:
            T = Temp
        Eg = Eg_o - Ao * np.divide(T**2, T + Bo)
        return Eg                                               # The array size is the same as the temp func.

    def analyticalDoS(self, energyRange, alpha):                # See the manual fot the eqution

        DoS_nonparabolic = 1/np.pi**2*np.sqrt(2*energyRange*(1+energyRange*np.transpose(alpha)))*np.sqrt(self.electronEffectiveMass/thermoelectricProperties.hBar**2)**3*(1+(2*energyRange*np.transpose(alpha)))/thermoelectricProperties.e2C**(3./2)
        DoS_parabolic = np.sqrt(energyRange)/np.pi**2*np.sqrt(2)/thermoelectricProperties.hBar**3*self.electronEffectiveMass**(3/2)/thermoelectricProperties.e2C**(3/2)
        DoS = [DoS_nonparabolic,DoS_parabolic]
        return DoS    # The array size is [2,numEnergySampling], first row is the with nonparabolic while the second row is parabolic DoS

    def carrierConcentration(self, path2extrinsicCarrierConcentration, bandGap, Ao=None, Bo=None, Nc=None, Nv=None, Temp=None):

        # This function computes the carrier concentration

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
        intrinsicCarrierConcentration = np.multiply(np.sqrt(np.multiply(Nc, Nv)), np.exp(-(np.divide(bandGap, (2 * thermoelectricProperties.kB * T)))))  # The formula for temperature dependent intrinsic carreir concentration is given in the manual
        totalCarrierConcentration = intrinsicCarrierConcentration + abs(extrinsicCarrierConcentration)
        return totalCarrierConcentration # Total carreir concentration is sum of intrinsic and extrinsic carrier concentration

    def fermiLevel(self, carrierConcentration, energyRange, DoS, Nc=None, Ao=None, Temp=None):

        # This function uses Joice Dixon approximation to predict Ef and thereby the carreir concentration at each temperature

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
        return [fermiLevelEnergy,np.expand_dims(n,axis=0)] # The list size is [2, size(temp)]

    def fermiDistribution(self, energyRange, fermiLevel, Temp=None):

        # This function compute the Fermi distribution

        if Temp is None:
            T = self.temp()
        else:
            T = Temp

        xi = np.exp((energyRange-fermiLevel.T)/T.T/thermoelectricProperties.kB)
        fermiDirac = 1/(xi+1)
        dfdE = -1*xi/(1+xi)**2/T.T/thermoelectricProperties.kB
        fermi = np.array([fermiDirac, dfdE])
        return fermi  # The array size is [2, size(temp)], The first row is the Fermi and the second row is the derivative(Fermi window)

    def electronBandStructure(self, path2eigenval, skipLines):

        # This function ead EIGENVAL file from VASP

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
        return dispersian # The array size is [(number of bands + 1) by (number of kpoints)]

    def electronDoS(self, path2DoS, headerLines, numDoSpoints, unitcell_volume, valleyPoint, energyRange):

        # This function ead DOSCAR file from VASP
        # The unitcell_volume is in unit of m

        DoS = np.loadtxt(expanduser(path2DoS), delimiter=None, skiprows=headerLines, max_rows=numDoSpoints)
        valleyPointEnergy = DoS[valleyPoint, 0]
        DoSSpline = InterpolatedUnivariateSpline(DoS[valleyPoint:, 0] - valleyPointEnergy, DoS[valleyPoint:, 1] / unitcell_volume)
        DoSFunctionEnergy = DoSSpline(energyRange)  # Density of state
        return DoSFunctionEnergy # The array size is [1, numEnergySampling]

    def fermiLevelSelfConsistent(self, carrierConcentration, Temp, energyRange, DoS, fermilevel):

        """
        A tool for self-consistent calculation of the Fermi level from a given carrier concentration to circumvent the problem that DFT underestimates the band gaps
        """

        fermi = np.linspace(fermilevel[0]-0.4, fermilevel[0]+0.2, 4000, endpoint=True).T
        # Joyce Dixon approx. is a good initial point for degenerate semiconductors. I considered 4000 points in energy range from Ef(JD)-0.4 eV up to Ef(JD)+0.2. This look reasonble in most cases.
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
        print("Fermi Level Self Consistent Index ",min_idx)
        # This print the index if it reaches the extreme index (0) or (4000), increase number of sampling points.

        Ef = np.empty((1, np.shape(Temp)[1]))
        for Ef_idx in np.arange(len(min_idx)):
            Ef[0,Ef_idx] = fermi[Ef_idx,min_idx[Ef_idx]]
        elm = 0
        n = np.empty((1, np.shape(Temp)[1]))
        for idx in min_idx:
            n[0,elm] = result_array[elm, idx]
            elm += 1
        return [Ef,n] # The array size is [2, size(temp)], The first row is the Fermi and the second row is the carrier concentration

    def electronGroupVelocity(self, kp, energy_kp, energyRange):
        dE = np.roll(energy_kp, -1, axis=0) - np.roll(energy_kp, 1, axis=0)
        dk = np.roll(kp, -1, axis=0) - np.roll(kp, 1, axis=0)
        dEdk = np.divide(dE, dk)
        dEdk[0] = (energy_kp[1] - energy_kp[0]) / (kp[1] - kp[0])
        dEdk[-1] = (energy_kp[-1] - energy_kp[-2]) / (kp[-1] - kp[-2])
        dEdkSpline = InterpolatedUnivariateSpline(energy_kp, np.array(dEdk))
        dEdkFunctionEnergy = dEdkSpline(energyRange)
        groupVel = dEdkFunctionEnergy / thermoelectricProperties.hBar
        return groupVel # This is the derivation of band structure from DFT. BTE needs single band data. See the example (si.py) or the manual for the details

    def analyticalGroupVelocity(self,energyRange, nk, m, valley, dk_len, alpha):

        """
        If no DFT calculation is availble this function approximate the group velocity near the conduction band edge. This works well up to few hundreds of mev.
        """

        meff = np.array(m)
        ko = 2 * np.pi / self.latticeParameter * np.array(valley)
        del_k = 2*np.pi/self.latticeParameter * dk_len * np.array([1, 1, 1])
        kx = np.linspace(ko[0], ko[0] + del_k[0], nk[0], endpoint=True)  # kpoints mesh
        ky = np.linspace(ko[1], ko[1] + del_k[1], nk[1], endpoint=True)  # kpoints mesh
        kz = np.linspace(ko[2], ko[2] + del_k[2], nk[2], endpoint=True)  # kpoints mesh
        [xk, yk, zk] = np.meshgrid(kx, ky, kz)
        xk_ = np.reshape(xk, -1)
        yk_ = np.reshape(yk, -1)
        zk_ = np.reshape(zk, -1)
        kpoint = np.array([xk_, yk_, zk_])
        mag_kpoint = norm(kpoint, axis=0)
        mc = 3/(1/meff[0]+1/meff[1]+1/meff[2])
        E = thermoelectricProperties.hBar**2 / 2 * ((kpoint[0] - ko[0])**2 / meff[0] + (kpoint[1] - ko[1])**2 / meff[1] + (kpoint[2] - ko[2]) ** 2 / meff[2]) * thermoelectricProperties.e2C
        vel = thermoelectricProperties.hBar*np.sqrt((kpoint[0]-ko[0])**2+(kpoint[1]-ko[1])**2+(kpoint[2]-ko[2])**2)/mc/(1+2*alpha*E)*thermoelectricProperties.e2C
        Ec, indices, return_indices = np.unique(E, return_index=True, return_inverse=True)
        vg = accum(return_indices, vel, func=np.mean, dtype=float)
        ESpline = PchipInterpolator(Ec, vg)
        velFunctionEnergy = ESpline(energyRange)
        return velFunctionEnergy  # The array size is [1, numEnergySampling]

    def matthiessen(self, *args):
        tau = 1. / sum([1. / arg for arg in args])
        tau[np.isinf(tau)] = 0
        return tau # Using  Matthiessen's rule to sum the scattering rates

    def tau_p(self, energyRange, alpha, Dv, DA, T, vs, D, rho):

        # Electron-phonon scattering rate using Ravich model

        nonparabolic_term = (1-((alpha.T*energyRange)/(1+2*alpha.T*energyRange)*(1-Dv/DA)))**2-8/3*(alpha.T*energyRange)*(1+alpha.T*energyRange)/(1+2*alpha.T*energyRange)**2*(Dv/DA)
        tau = rho*vs**2*thermoelectricProperties.hBar/np.pi/thermoelectricProperties.kB/T.T/DA/DA*1e9/thermoelectricProperties.e2C/D
        tau_p = tau/nonparabolic_term
        return [tau,tau_p] # The first row does not count for nonparabolicity, the second row does

    def tau_Screened_Coulomb(self,energyRange, m_c, LD, N):

        # Electron-ion scattering rate following Brook Heiring model

        g = 8*m_c.T*LD.T**2*energyRange/thermoelectricProperties.hBar**2/thermoelectricProperties.e2C
        var_tmp = np.log(1+g)-g/(1+g)
        tau = 16*np.pi*np.sqrt(2*m_c.T)*(4*np.pi*self.dielectric*thermoelectricProperties.e0)**2/N.T/var_tmp*energyRange**(3/2)/thermoelectricProperties.e2C**(5/2)
        where_are_NaNs = np.isnan(tau)
        tau[where_are_NaNs] = 0
        return tau  # The array size is [1, numEnergySampling]

    def tau_Unscreened_Coulomb(self,energyRange, m_c, N):

        # Electron-ion scattering rate for shallow dopants ~10^18 1/cm^3 (no screening effect is considered)

        g = 4*np.pi*(4*np.pi*self.dielectric*thermoelectricProperties.e0)*energyRange/N.T**(1/3)/thermoelectricProperties.e2C
        var_tmp = np.log(1+g**2)
        tau = 16*np.pi*np.sqrt(2*m_c.T)*(4*np.pi*self.dielectric*thermoelectricProperties.e0)**2/N.T/var_tmp*energyRange**(3/2)/thermoelectricProperties.e2C**(5/2)
        where_are_NaNs = np.isnan(tau)
        tau[where_are_NaNs] = 0
        return tau  # The array size is [1, numEnergySampling]

    def tau_Strongly_Screened_Coulomb(self, D, LD, N):
        tau = thermoelectricProperties.hBar/N.T/np.pi/D/(LD.T**2/(4*np.pi*self.dielectric*thermoelectricProperties.e0))**2*1/thermoelectricProperties.e2C**2
        return tau  # The array size is [1, numEnergySampling]

    def tau2D_cylinder(self,energyRange, nk, Uo, m, vfrac, valley, dk_len, ro, n=2000):

        """
        This is a fast algorithm that uses Fermi’s golden rule to compute the energy dependent electron scattering rate due cylindrical nanoparticles or pores extended perpendicular to the electrical current
        """

        meff = np.array(m) * thermoelectricProperties.me
        ko = 2 * np.pi / self.latticeParameter * np.array(valley)
        del_k = 2*np.pi/self.latticeParameter * dk_len * np.array([1, 1, 1])
        N = vfrac/np.pi/ro**2
        kx = np.linspace(ko[0], ko[0] + del_k[0], nk[0], endpoint=True)  # kpoints mesh
        ky = np.linspace(ko[1], ko[1] + del_k[1], nk[1], endpoint=True)  # kpoints mesh
        kz = np.linspace(ko[2], ko[2] + del_k[2], nk[2], endpoint=True)  # kpoints mesh
        [xk, yk, zk] = np.meshgrid(kx, ky, kz)
        xk_ = np.reshape(xk, -1)
        yk_ = np.reshape(yk, -1)
        zk_ = np.reshape(zk, -1)
        kpoint = np.array([xk_, yk_, zk_])
        mag_kpoint = norm(kpoint, axis=0)
        E = thermoelectricProperties.hBar**2 / 2 * ((kpoint[0, :] - ko[0])**2 / meff[0] + (kpoint[1, :] - ko[1])**2 / meff[1] + (kpoint[2, :] - ko[2]) ** 2 / meff[2]) * thermoelectricProperties.e2C
        t = np.linspace(0, 2*np.pi, n)
        a = np.expand_dims(np.sqrt(2 * meff[1] / thermoelectricProperties.hBar**2 * E / thermoelectricProperties.e2C), axis=0)
        b = np.expand_dims(np.sqrt(2 * meff[2] / thermoelectricProperties.hBar**2 * E / thermoelectricProperties.e2C), axis=0)
        ds = np.sqrt((a.T * np.sin(t))**2 + (b.T * np.cos(t))**2)
        cos_theta = ((a * kpoint[0]).T * np.cos(t) + (b * kpoint[1]).T * np.sin(t) + np.expand_dims(kpoint[2]**2, axis=1)) / np.sqrt(a.T**2 * np.cos(t)**2 + b.T**2 * np.sin(t)**2 + np.expand_dims(kpoint[2]**2, axis=1)) / np.expand_dims(mag_kpoint, axis=1)
        delE = thermoelectricProperties.hBar**2 * np.abs((a.T * np.cos(t) - ko[0]) / meff[0] + (b.T * np.sin(t) - ko[1]) / meff[1] + (np.expand_dims(kpoint[2]**2, axis=1) - ko[2] / meff[2]))
        qx = np.expand_dims(kpoint[0], axis=1) - a.T * np.cos(t)
        qy = np.expand_dims(kpoint[1], axis=1) - b.T * np.sin(t)
        qr = np.sqrt(qx**2 + qy**2)
        tau = np.empty((len(ro), len(E)))
        for r_idx in np.arange(len(ro)):
            J = jv(1, ro[r_idx] * qr)
            SR = 2 * np.pi / thermoelectricProperties.hBar * Uo**2 * (2 * np.pi)**3 * (ro[r_idx] * J / qr)**2
            f = SR * (1 - cos_theta) / delE * ds
            int_ = np.trapz(f, t, axis=1)
            tau[r_idx] = 1 / (N[r_idx] / (2 * np.pi)**3 * int_) * thermoelectricProperties.e2C
        Ec, indices, return_indices = np.unique(E, return_index=True, return_inverse=True)
        tau_c = np.empty((len(ro), len(indices)))
        tauFunctionEnergy = np.empty((len(ro), len(energyRange[0])))
        for r_idx in np.arange(len(ro)):
            tau_c[r_idx] = accum(return_indices, tau[r_idx], func=np.mean, dtype=float)
        for tau_idx in np.arange(len(tau_c)):
            ESpline = PchipInterpolator(Ec[30:], tau_c[tau_idx,30:])
            tauFunctionEnergy[tau_idx] = ESpline(energyRange)
        return tauFunctionEnergy

    def tau3D_spherical(self,energyRange, nk, Uo, m, vfrac, valley, dk_len, ro, n=32):

        """
        This is a fast algorithm that uses Fermi’s golden rule to compute the energy dependent electron scattering rate due spherical nanoparticles or pores.
        """

        meff = np.array(m) * thermoelectricProperties.me
        ko = 2 * np.pi / self.latticeParameter * np.array(valley)
        del_k = 2*np.pi/self.latticeParameter * dk_len * np.array([1, 1, 1])
        N = 3*vfrac/4/np.pi/ro**3
        kx = np.linspace(ko[0], ko[0] + del_k[0], nk[0], endpoint=True)  # kpoints mesh
        ky = np.linspace(ko[1], ko[1] + del_k[1], nk[1], endpoint=True)  # kpoints mesh
        kz = np.linspace(ko[2], ko[2] + del_k[2], nk[2], endpoint=True)  # kpoints mesh
        [xk, yk, zk] = np.meshgrid(kx, ky, kz)
        xk_ = np.reshape(xk, -1)
        yk_ = np.reshape(yk, -1)
        zk_ = np.reshape(zk, -1)
        kpoint = np.array([xk_, yk_, zk_])
        mag_kpoint = norm(kpoint, axis=0)
        E = thermoelectricProperties.hBar**2 / 2 * ((kpoint[0, :] - ko[0])**2 / meff[0] + (kpoint[1, :] - ko[1])**2 / meff[1] + (kpoint[2, :] - ko[2]) ** 2 / meff[2]) * thermoelectricProperties.e2C
        scattering_rate = np.zeros((len(ro), len(E)))
        nu = np.linspace(0, np.pi, n)
        z_ = -1 * np.cos(nu)
        r = np.sqrt(1.0 - z_**2)[:, None]
        theta = np.linspace(0, 2 * np.pi, n)[None, :]
        x_ = r * np.cos(theta)
        y_ = r * np.sin(theta)
        for u in np.arange(len(E)):

            Q = np.zeros((2 * (n-2) * (n - 1), 3))
            A = np.zeros((2 * (n-2) * (n - 1), 1))
            k = 0
            a_axis = np.sqrt(2 / (thermoelectricProperties.hBar**2 * thermoelectricProperties.e2C) * meff[0] * E[u])
            b_axis = np.sqrt(2 / (thermoelectricProperties.hBar**2 * thermoelectricProperties.e2C) * meff[1] * E[u])
            c_axis = np.sqrt(2 / (thermoelectricProperties.hBar**2 * thermoelectricProperties.e2C) * meff[2] * E[u])

            y = -1 * b_axis * y_ + ko[1]
            x = -1 * a_axis * x_ + ko[0]
            Z_ = c_axis * z_ + ko[2]
            z = np.tile(Z_[:, None], (1,n))
            for j in np.arange(1,n-1):
                for i in np.arange(2,n):
                    S = np.array(np.array([x[i,j],y[i,j],z[i,j]])+np.array([x[i-1,j],y[i-1,j],z[i-1,j]])+np.array([x[i-1,j-1],y[i-1,j-1],z[i-1,j-1]]))
                    Q[k] = S/3
                    a = norm(np.array([x[i,j],y[i,j],z[i,j]])-np.array([x[i-1,j],y[i-1,j],z[i-1,j]]))
                    b = norm(np.array([x[i-1,j],y[i-1,j],z[i-1,j]])-np.array([x[i-1,j-1],y[i-1,j-1],z[i-1,j-1]]))
                    c = norm(np.array([x[i-1,j-1],y[i-1,j-1],z[i-1,j-1]])-np.array([x[i,j],y[i,j],z[i,j]]))
                    s = a+b+c
                    s = s/2
                    A[k] = np.sqrt(s*(s-a)*(s-b)*(s-c))
                    k += 1
            for j in np.arange(1,n-1):
                for i in np.arange(1,n-1):
                    S = np.array([x[i,j-1],y[i,j-1],z[i,j-1]])+np.array([x[i,j],y[i,j],z[i,j]])+np.array([x[i-1,j-1],y[i-1,j-1],z[i-1,j-1]])
                    Q[k] = S/3
                    a = norm(np.array([x[i,j-1],y[i,j-1],z[i,j-1]])-np.array([x[i,j],y[i,j],z[i,j]]))
                    b = norm(np.array([x[i,j],y[i,j],z[i,j]])-np.array([x[i-1,j-1],y[i-1,j-1],z[i-1,j-1]]))
                    c = norm(np.array([x[i-1,j-1],y[i-1,j-1],z[i-1,j-1]])-np.array([x[i,j-1],y[i,j-1],z[i,j-1]]))
                    s = a+b+c
                    s = s/2
                    A[k] = np.sqrt(s*(s-a)*(s-b)*(s-c))
                    k += 1
            for i in np.arange(2,n):
                S = np.array([x[i,0],y[i,0],z[i,0]])+np.array([x[i-1,0],y[i-1,0],z[i-1,0]])+np.array([x[i-1,-2],y[i-1,-2],z[i-1,-2]])
                Q[k] = S/3
                a = norm(np.array([x[i,0],y[i,0],z[i,0]])-np.array([x[i-1,0],y[i-1,0],z[i-1,0]]))
                b = norm(np.array([x[i-1,0],y[i-1,0],z[i-1,0]])-np.array([x[i-1,-2],y[i-1,-2],z[i-1,-2]]))
                c = norm(np.array([x[i-1,-2],y[i-1,-2],z[i-1,-2]])-np.array([x[i,0],y[i,0],z[i,0]]))
                s = a+b+c
                s = s/2
                A[k] = np.sqrt(s*(s-a)*(s-b)*(s-c))
                k += 1
            for i in np.arange(1,n-1):
                S = np.array([x[i,-2],y[i,-2],z[i,-2]])+np.array([x[i,0],y[i,0],z[i,0]])+np.array([x[i-1,-2],y[i-1,-2],z[i-1,-2]])
                Q[k] = S/3
                a = norm(np.array([x[i,-2],y[i,-2],z[i,-2]])-np.array([x[i,0],y[i,0],z[i,0]]))
                b = norm(np.array([x[i,0],y[i,0],z[i,0]])-np.array([x[i-1,-2],y[i-1,-2],z[i-1,-2]]))
                c = norm(np.array([x[i-1,-2],y[i-1,-2],z[i-1,-2]])-np.array([x[i,-2],y[i,-2],z[i,-2]]))
                s = a+b+c
                s = s/2
                A[k] = np.sqrt(s*(s-a)*(s-b)*(s-c))
                k += 1
            qx = kpoint[0,u] - Q[:,0]
            qy = kpoint[1,u] - Q[:,1]
            qz = kpoint[2,u] - Q[:,2]
            q  = np.sqrt(qx**2+qy**2+qz**2)
            cosTheta = np.matmul(kpoint[:,u][None,:],Q.T)/norm(kpoint[:,u])/np.sqrt(np.sum(Q**2,axis=1))
            delE = np.abs(thermoelectricProperties.hBar**2*((Q[:,0]-ko[0])/meff[0]+(Q[:,1]-ko[1])/meff[1]+(Q[:,2]-ko[2])/meff[2]))
            for ro_idx in  np.arange(len(ro)):
              M = 4*np.pi*Uo*(1/q*np.sin(ro[ro_idx]*q)-ro[ro_idx]*np.cos(ro[ro_idx]*q))/(q**2)
              SR = 2*np.pi/thermoelectricProperties.hBar*M*np.conj(M)
              f = SR/delE*(1-cosTheta);
              scattering_rate[ro_idx,u] = N[ro_idx]/(2*np.pi)**3*np.sum(f*A.T)
        return scattering_rate


    def electricalProperties(self, E, DoS, vg, Ef, dfdE, Temp, tau):

        # This function returns a list of thermoelectric properties
        # See the manual for the detail of calculations

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
        return coefficients # The list is 7 by numEnergySampling

    def filteringEffect(self, U, E, DoS, vg, Ef, dfdE, Temp, tau_b):

        # This function returns list of electrical conductivity and Seebecl for the ideal filtering where all the electrons up to a cutoff energy level of U are completely filtered

        tauUo = np.ones(len(E[0]))
        _Conductivity = [np.empty([1, len(tau_b)])]
        _Seebeck = [np.empty([1, len(tau_b)])]
        for i in np.arange(len(U)):
            tau_idl = copy.copy(tauUo)
            tau_idl[E[0]<U[i]] = 0
            tau = self.matthiessen(E, tau_idl, tau_b)
            coefficients = self.electricalProperties(E=E, DoS=DoS, vg=vg, Ef=Ef, dfdE=dfdE, Temp=Temp, tau=tau)
            Sigma = np.expand_dims(coefficients[0],axis=0)
            S = np.expand_dims(coefficients[1],axis=0)
            _Conductivity = np.append(_Conductivity, [Sigma], axis=0)
            _Seebeck = np.append(_Seebeck, [S], axis=0)
            del tau_idl
        Conductivity = np.delete(_Conductivity, 0, axis = 0)
        Seebeck = np.delete(_Seebeck, 0, axis = 0)
        return [Conductivity, Seebeck] # The list is 2 by numEnergySampling

    def phenomenological(self, U, tauo, E, DoS, vg, Ef, dfdE, Temp, tau_b):

        # This function returns list of electrical conductivity and Seebecl for the phenomenological model where a frequency independent lifetime of tauo is imposed to all the electrons up to a cutoff energy level of U

        tauU = np.ones(len(E[0]))
        _Conductivity = [np.empty([1, 1])]
        _Seebeck = [np.empty([1, 1])]
        for _j in np.arange(len(tauo)):
            for _i in np.arange(len(U)):
                tau_ph = copy.copy(tauU)
                tau_ph[E[0]<U[_i]] = tauo[_j]
                tau = self.matthiessen(E, tau_ph, tau_b)
                coefficients = self.electricalProperties(E=E, DoS=DoS, vg=vg, Ef=Ef, dfdE=dfdE, Temp=Temp, tau=tau)
                Sigma = np.expand_dims(coefficients[0],axis=0)
                S = np.expand_dims(coefficients[1],axis=0)
                _Conductivity = np.append(_Conductivity, [Sigma], axis=0)
                _Seebeck = np.append(_Seebeck, [S], axis=0)
                del tau_ph
        __Conductivity = np.delete(_Conductivity, 0, axis = 0)
        __Seebeck = np.delete(_Seebeck, 0, axis = 0)
        Conductivity = np.reshape(__Conductivity,(len(tauo),len(U)))
        Seebeck = np.reshape(__Seebeck,(len(tauo),len(U)))
        return [Conductivity, Seebeck]  # The list is 2 by numEnergySampling
