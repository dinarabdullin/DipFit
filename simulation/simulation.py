'''
Simulation of RIDME spectra and corresponding time traces
'''

import sys
import math
import numpy as np
from constants import const
from math_func.randomPointsOnSphere import randomPointsOnSphere
from math_func.euler2rotationMatrix import euler2rotationMatrix
from math_func.spherical2cartesian import spherical2cartesian
from spin_func.gEffective import gEffective
from spin_func.quantisationAxis import quantisationAxis
from spin_func.flipProbabilities import flipProbabilities
from spin_func.gFactorHighSpinFe import gFactorHighSpinFe
from math_func.rmsd import rmsd

class Simulator:

    def __init__(self, calc):
        self.Ns = calc['Ns']
        self.r_distr = calc['r_distr']
        self.xi_distr = calc['xi_distr']
        self.phi_distr = calc['phi_distr']
        self.spc_max = calc['spc_max']
        self.gA = []
        self.gB = []
        self.qA = []
        self.qB = []
        self.weightB = []
        self.Nf = 0
        self.f = []
        self.fn = []
        self.spc = []
        self.spc_vs_theta = []
        self.spc_vs_xi = []
        self.spc_vs_phi = []
        self.spc_vs_E = []
        self.Nt = 0
        self.t = []
        self.sig = []
        self.modulation_depth = 0

    def set_faxis(self, exp, spinA = {}, spinB = {}, var = []):
        if not (exp['f'] == []):
            f = np.array(exp['f'])
        else:
            # Determine the min value of r
            rMin = 0.0
            if (var['r_width'] > 0):
                if (self.r_distr == 'n'):
                    r = var['r_mean'] + var['r_width'] * np.random.randn(self.Ns)
                elif (self.r_distr == 'u'):
                    r = var['r_mean'] + var['r_width'] * (0.5 - np.random.rand(self.Ns))
                rMin = np.amin(r) 
            else:
                rMin = var['r_mean']
            # Determine the max values of g factors of both spins
            gAMax = np.amax(spinA['g'])
            gBMax = np.amax(spinB['g'])
            # Estimate the max dipolar frequency
            fddMax = 2 * const['Fdd'] * gAMax * gBMax / rMin**3 
            # Set the frequency axis
            fMax = np.around(fddMax) + 5.0
            fMin = -fMax
            fStep = 0.1
            Nf = int(np.around((fMax - fMin) / fStep)) + 1
            f = np.linspace(fMin, fMax, Nf)
        return f

    def normalize_faxis(self, spinA = {}, spinB = {}, var = []):
        # Determine the reference frequency
        fddRef = const['Fdd'] * const['ge'] * const['ge'] / var['r_mean']**3 
        # Create an array with the normalized frequencies
        Nf = self.f.size
        fn = np.zeros(Nf)
        for i in range(Nf):
            fn[i] = self.f[i] / fddRef
        return fn 
        
    def set_taxis(self, exp, spinA = {}, spinB = {}, var = []):
        if not (exp['t'] == []):
            t = np.array(exp['t'])
        else:
            # Determine the max value of r
            if (var['r_width'] > 0):
                if (self.r_distr == 'n'):
                    r = var['r_mean'] + var['r_width'] * np.random.randn(self.Ns)
                elif (self.r_distr == 'u'):
                    r = var['r_mean'] + var['r_width'] * (0.5 - np.random.rand(self.Ns))
                rMax = np.amax(r) 
            else:
                rMax = var['r_mean']
            # Determine the min values g factors for both spins
            gAMin = np.amin(spinA['g'])
            gBMin = np.amin(spinB['g'])
            # Estimate the min dipolar frequency
            fddMin = const['Fdd'] * gAMin * gBMin / rMax**3
            # Estimate the max period of dipolar oscilation
            tddMax = 1 / fddMin
            # Set the time axis
            tMin = 0.0
            tMax = 3 * tddMax
            tStep = 0.008
            Nt = int(np.around((tMax - tMin) / tStep)) + 1
            t = np.linspace(tMin, tMax, Nt)
        return t       
        
    def set_modulation_depth(self, sigExp):
        if not (sigExp == []):
            depth = 1.0 - sigExp[-1]
        else:
            depth = 0
        return depth
        
    def precalculations(self, spinA, spinB):
        # Directions of the magnetic field in the lab frame
        eB0_L = randomPointsOnSphere(self.Ns)
        # Directions of the magnetic field in the frame of spin A
        if (spinA['gFrame'].all() == 0):
            eB0_A = eB0_L
        else:
            RM_L2A = euler2rotationMatrix(spinA['gFrame'])
            eB0_A = np.dot(RM_L2A.T, eB0_L.T).T
        # Directions of the magnetic field in the frame of spin B
        if (spinB['gFrame'].all() == 0):
            eB0_B = eB0_L
        else:
            RM_L2B = euler2rotationMatrix(spinB['gFrame'])
            eB0_B = np.dot(RM_L2B.T, eB0_L.T).T
        # Effective g factors of spin A 
        gA = gEffective(spinA['g'], eB0_A)
        # Effective g factors of spin B 
        gB = gEffective(spinB['g'], eB0_B)   
        # Quantisation axes of spin A
        qA = eB0_L
        # qA = quantisationAxis(spinA['g'], gA, eB0_L)
        # Quantisation axes of spin B
        qB = quantisationAxis(spinB['g'], gB, eB0_L)
        return [gA, gB, qA, qB]

    def dipolar_frequencies(self, var, gA, gB, qA, qB, calculateTheta=False):	
        # Set the spherical coordinates of the distance vector
        if (var['r_width'] > 0):
            if (self.r_distr == 'n'):
                r = var['r_mean'] + var['r_width'] * np.random.randn(self.Ns)
            elif (self.r_distr == 'u'):
                r = var['r_mean'] + var['r_width'] * (0.5 - np.random.rand(self.Ns))
        else:
            r = var['r_mean'] * np.ones(self.Ns)
        if (var['xi_width'] > 0):
            if (self.xi_distr == 'n'):
                xi = var['xi_mean'] + var['xi_width'] * np.random.randn(self.Ns)
            elif (self.xi_distr == 'u'):
                xi = var['xi_mean'] + var['xi_width'] * (0.5 - np.random.rand(self.Ns))
        else:
            xi = var['xi_mean'] * np.ones(self.Ns)	
        if (var['phi_width'] > 0):
            if (self.phi_distr == 'n'):
                phi = var['phi_mean'] + var['phi_width'] * np.random.randn(self.Ns)
            elif (self.phi_distr == 'u'):
                phi = var['phi_mean'] + var['phi_width'] * (0.5 - np.random.rand(self.Ns))
        else:
            phi = var['phi_mean'] * np.ones(self.Ns)	    
        # Directions of the distance vector
        ers = np.array([np.ones(self.Ns), xi, phi])
        ers = ers.T
        erc = spherical2cartesian(ers)
        # Calculate dipolar frequencies and theta
        fdd = np.zeros(self.Ns)
        theta = np.zeros(self.Ns)	
        for i in range(self.Ns):
            # Projection of the distance vector on the quantisation axis of spin A
            prA = np.dot(erc[i], qA[i].T)
            # Projection of the distance vector on the quantisation axis of spin B
            prB = np.dot(erc[i], qB[i].T)
            # Dipolar frequency
            fdd[i] = const['Fdd'] * gA[i] * gB[i] * (1.0 - 3.0 * prA * prB) / r[i]**3
            # Theta
            if (calculateTheta):
                theta[i] = np.arccos(prA) * const['rad2deg']
                if (theta[i] < 0.0):
                    theta[i] = -theta[i]
                if (theta[i] > 90.0):
                    theta[i] = 180.0 - theta[i]
        return [fdd, theta]

    def dipolar_spectrum(self, fdd, weights=[]):	
        # Set frequency bounds	
        fb = np.zeros(self.Nf + 1)
        fb[:-1] = self.f - 0.5 * (self.f[1] - self.f[0]) * np.ones(self.Nf)
        fb[self.Nf] = self.f[self.Nf-1] + 0.5 * (self.f[1] - self.f[0])
        # Set weights
        if (weights == []):
            weights = self.weightB
        # Calculate the spectrum
        spc, bin_edges = np.histogram(fdd, bins=fb, weights=weights)
        spc = spc.astype(float)
        spc = spc + spc[::-1]
        # Normalize the spectrum
        spc = spc * self.spc_max / np.amax(spc)
        return spc

    def dipolar_spectrum_vs_theta(self, fdd, theta, weights=[]):
        # Some parameters of the frequency axis
        Nf = self.Nf
        fMin = self.f[0]
        fStep = self.f[1] - self.f[0]
        # Some parameters of the theta axis
        Ntheta = const['theta'].size
        thetaMin = const['theta'][0]
        thetaStep = const['theta'][1] - const['theta'][0]
        if (weights == []):
            weights = self.weightB
        #  Calculate a dipolar spectrum vs theta
        spc_vs_theta = np.zeros((Ntheta,Nf))
        for i in range(self.Ns):
            idx_f1 = int(np.around((fdd[i] - fMin) / fStep))
            idx_f2 = int(np.around((-fdd[i] - fMin) / fStep))
            idx_theta = int(np.around((theta[i] - thetaMin) / thetaStep))
            spc_vs_theta[idx_theta][idx_f1] += weights[i]
            spc_vs_theta[idx_theta][idx_f2] += weights[i]
        # Normalize the spectrum
        spc_vs_theta = spc_vs_theta * self.spc_max / np.amax(spc_vs_theta)
        return spc_vs_theta
        
    def dipolar_spectrum_vs_xi(self, var):
        Nf = self.Nf
        Nxi = const['xi'].size
        spc_vs_xi = np.zeros((Nxi, Nf))
        for i in range(Nxi):
            sys.stdout.write('\r')
            status = int(float(i+1) / float(Nxi) * 100)
            sys.stdout.write("Calculating dipolar spectrum vs xi... %d%% " % (status))
            sys.stdout.flush()
            varNew = var
            varNew['xi_mean'] = const['xi'][i] * const['deg2rad']
            fdd, theta = self.dipolar_frequencies(varNew, self.gA, self.gB, self.qA, self.qB)
            spc = self.dipolar_spectrum(fdd)
            spc_vs_xi[i] = spc
        return spc_vs_xi

    def dipolar_spectrum_vs_phi(self, var):
        Nf = self.Nf
        Nphi = const['phi'].size
        spc_vs_phi = np.zeros((Nphi, Nf))
        for i in range(Nphi):
            sys.stdout.write('\r')
            status = int(float(i+1) / float(Nphi) * 100)
            sys.stdout.write("Calculating dipolar spectrum vs phi... %d%% " % (status))
            sys.stdout.flush()
            varNew = var
            varNew['phi_mean'] = const['phi'][i] * const['deg2rad']
            fdd, theta = self.dipolar_frequencies(varNew, self.gA, self.gB, self.qA, self.qB)
            spc = self.dipolar_spectrum(fdd)
            spc_vs_phi[i] = spc
        return spc_vs_phi	

    def dipolar_spectrum_vs_E(self, var, exp, spinA, spinB):
        Nf = self.Nf
        Ne = const['E'].size
        spc_vs_E = np.zeros((Ne, Nf))
        for i in range(Ne):
            sys.stdout.write('\r')
            status = int(float(i+1) / float(Ne) * 100)
            sys.stdout.write("Calculating dipolar spectrum vs E/D... %d%% " % (status))
            sys.stdout.flush()
            spinBNew = spinB
            spinBNew['D'][1] = spinBNew['D'][0] * const['E'][i]
            spinBNew['g'] = gFactorHighSpinFe(spinBNew, exp['mwFreq'])
            gA, gB, qA, qB = self.precalculations(spinA, spinBNew)
            weightB = flipProbabilities(spinBNew['g'], gB, exp['magnField'], var['temp'])
            fdd, theta = self.dipolar_frequencies(var, gA, gB, qA, qB)
            spc = self.dipolar_spectrum(fdd, weightB)
            spc_vs_E[i] = spc
        # # Normalize the integral
        # I = np.sum(spc_vs_E, axis=1)
        # Imin = np.amin(I)
        # for i in range(Ne):
            # spc_vs_E[i] = spc_vs_E[i] * Imin / I[i]
        return spc_vs_E

    def dipolar_timetrace(self, fdd):	
        # Calculate a time trace
        sig = np.zeros(self.Nt)
        for i in range(self.Ns):
            for j in range(self.Nt):
                pFlip = self.weightB[i]
                sig[j] += pFlip * (1.0 - np.cos(2.0 * np.pi * fdd[i] * self.t[j]))
        sigAmp = float(self.Ns)
        # Correction factor for a modulation depth
        if (self.modulation_depth):
            depth = sig[-1] / sigAmp
            depthCorr = self.modulation_depth / depth
        else:
            depthCorr = 1.0
        for j in range(self.Nt):
            sig[j] = (sigAmp - depthCorr * sig[j]) / sigAmp
        return sig

    def run_simulation(self, simPar, exp, spinA, spinB, calc):
        # Display status
        sys.stdout.write('Starting simulation...\n')
        sys.stdout.write('Running pre-calculations... ')
        # Set the frequency axis
        if (simPar['settings']['spc'] or simPar['settings']['spc_vs_theta'] or simPar['settings']['spc_vs_xi'] or simPar['settings']['spc_vs_phi'] or simPar['settings']['spc_vs_E']):
            self.f = self.set_faxis(exp, spinA, spinB, simPar['variables'])
            self.Nf = self.f.size
            if simPar['settings']['faxis_normalized']:
                self.fn = self.normalize_faxis(spinA, spinB, simPar['variables'])
        # Set the time axis
        if (simPar['settings']['timetrace']):
            self.t = self.set_taxis(exp, spinA, spinB, simPar['variables'])
            self.Nt = self.t.size
            self.modulation_depth = self.set_modulation_depth(exp['sig'])
        # Do some pre-calculations
        self.gA, self.gB, self.qA, self.qB = self.precalculations(spinA, spinB)
        # Calculate the Bolzmann population of the +1/2 energy level of spin B 
        # at a given temperature and a given magnetic field
        self.weightB = flipProbabilities(spinB['g'], self.gB, exp['magnField'], simPar['variables']['temp'])
        # Calculate the dipolar frequencies
        if (simPar['settings']['spc'] or simPar['settings']['spc_vs_theta'] or simPar['settings']['timetrace']):
            if simPar['settings']['spc_vs_theta']:
                fdd, theta = self.dipolar_frequencies(simPar['variables'], self.gA, self.gB, self.qA, self.qB, True)
            else:
                fdd, theta = self.dipolar_frequencies(simPar['variables'], self.gA, self.gB, self.qA, self.qB, False)
        sys.stdout.write('[DONE]\n')
        # Calculate the dipolar spectrum
        if (simPar['settings']['spc']):
            sys.stdout.write('Calculating dipolar spectrum... ')
            self.spc = self.dipolar_spectrum(fdd)
            sys.stdout.write('[DONE]\n') 
            if not (exp['f'] == []):
                fitness = rmsd(self.spc, exp['spc'], exp['f'], calc['f_min'], calc['f_max'])
                sys.stdout.write("RMSD = %f \n" % fitness)
        # Calculate the dipolar spectrum vs theta
        if (simPar['settings']['spc_vs_theta']):
            sys.stdout.write('Calculating dipolar spectrum vs theta... ')
            self.spc_vs_theta = self.dipolar_spectrum_vs_theta(fdd, theta)
            if not simPar['settings']['spc']:
                self.spc = self.dipolar_spectrum(fdd)
            sys.stdout.write('[DONE]\n')   
        # Calculate a time trace
        if simPar['settings']['timetrace']:
            sys.stdout.write('Calculating dipolar time trace... ')
            self.sig = self.dipolar_timetrace(fdd)
            sys.stdout.write('[DONE]\n')            
        # Calculate a dipolar spectrum vs xi
        if simPar['settings']['spc_vs_xi']:
            sys.stdout.write('Calculating dipolar spectrum vs xi... ')
            self.spc_vs_xi = self.dipolar_spectrum_vs_xi(simPar['variables'])
            sys.stdout.write('[DONE]\n')
        # Calculate a dipolar spectrum vs phi
        if simPar['settings']['spc_vs_phi']:
            sys.stdout.write('Calculating dipolar spectrum vs phi... ')
            self.spc_vs_phi = self.dipolar_spectrum_vs_phi(simPar['variables'])
            sys.stdout.write('[DONE]\n')
        # Calculate a dipolar spectrum vs E
        if (simPar['settings']['spc_vs_E'] and spinB['label'] == 'hs_fe'):
            sys.stdout.write('Calculating dipolar spectrum vs E/D... ')
            self.spc_vs_E = self.dipolar_spectrum_vs_E(simPar['variables'], exp, spinA, spinB)
            sys.stdout.write('[DONE]\n')
        sys.stdout.write('Simulation is finished\n\n')

    def init_fitting(self, fitPar, exp, spinA, spinB):
        # Set the frequency axis
        self.f = self.set_faxis(exp)
        self.Nf = self.f.size
        # Do some pre-calculations
        self.gA, self.gB, self.qA, self.qB = self.precalculations(spinA, spinB)
        # Calculate the Bolzmann population of the +1/2 energy level of spin B 
        # at a given temperature and a given magnetic field
        if (fitPar['variables']['indices']['temp'] == -1):
            self.weightB = flipProbabilities(spinB['g'], self.gB, exp['magnField'], fitPar['variables']['fixed']['temp'])