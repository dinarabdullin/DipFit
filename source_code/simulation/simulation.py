'''
Simulation of RIDME spectra and RIDME time traces
'''

import sys
import numpy as np
import time
import datetime
import copy
from mathematics.rand_points_on_sphere import rand_points_on_sphere
from mathematics.euler2RM import euler2RM
from mathematics.spherical2cartesian import spherical2cartesian
from mathematics.rmsd import rmsd
from spinphysics.effective_gfactor import effective_gfactor
from spinphysics.quantization_axis import quantization_axis
from spinphysics.flip_probabilities import flip_probabilities
from supplement.constants import const


class Simulator:

    def __init__(self, calcSettings):
        self.Ns = calcSettings['Ns']
        self.r_distr = calcSettings['r_distr']
        self.xi_distr = calcSettings['xi_distr']
        self.phi_distr = calcSettings['phi_distr']   
        self.field = []
        self.gA = []
        self.gB = []
        self.qA = []
        self.qB = []
        self.pB = []
        self.f = []
        self.fn = []
        self.Nf = 0
        self.t = []
        self.Nt = 0
        self.theta_bins = []
        self.xi_bins = []
        self.phi_bins = []
        self.temp_bins = []
        self.spc = []
        self.spc_max = calcSettings['spc_max']
        self.spc_vs_theta = []
        self.spc_vs_xi = []
        self.spc_vs_phi = []
        self.spc_vs_temp = []
        self.sig = [] 
        self.depth = 0
        self.g = []
        self.depth_vs_temp = []    

    def set_faxis(self, fExp, gA=[], gB=[], var={}):
        if not (fExp == []):
            f = np.array(fExp)
        else:
            # Determine the minimal value of r
            rMin = 0.0
            if (var['r_width'] > 0):
                if (self.r_distr == 'normal'):
                    r = var['r_mean'] + var['r_width'] * np.random.randn(self.Ns)
                elif (self.r_distr == 'uniform'):
                    r = var['r_mean'] + var['r_width'] * (0.5 - np.random.rand(self.Ns))
                rMin = np.amin(r) 
            else:
                rMin = var['r_mean']
            # Determine the maximal values of g-factors
            gAMax = np.amax(gA)
            gBMax = np.amax(gB)
            # Estimate the max dipolar frequency
            fddMax = 2 * const['Fdd'] * gAMax * gBMax / rMin**3 
            # Set the frequency axis
            fMax = np.around(fddMax) + 5.0
            fMin = -fMax
            fStep = 0.1
            Nf = int(np.around((fMax - fMin) / fStep)) + 1
            f = np.linspace(fMin, fMax, Nf)
        return f

    def normalize_faxis(self, var):
        # Determine the reference frequency
        fddRef = const['Fdd'] * const['ge'] * const['ge'] / var['r_mean']**3 
        # Create an array with the normalized frequencies
        Nf = self.f.size
        fn = np.zeros(Nf)
        for i in range(Nf):
            fn[i] = self.f[i] / fddRef
        return fn 
        
    def set_taxis(self, tExp, gA=[], gB=[], var={}):
        if not (tExp == []):
            t = np.array(tExp)
        else:
            # Determine the max value of r
            if (var['r_width'] > 0):
                if (self.r_distr == 'normal'):
                    r = var['r_mean'] + var['r_width'] * np.random.randn(self.Ns)
                elif (self.r_distr == 'uniform'):
                    r = var['r_mean'] + var['r_width'] * (0.5 - np.random.rand(self.Ns))
                rMax = np.amax(r) 
            else:
                rMax = var['r_mean']
            # Determine the mininimal values of g-factors
            gAMin = np.amin(gA)
            gBMin = np.amin(gB)
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

    def set_gaxis(self, g):
        gMin = np.amin(g)
        gMax = np.amax(g)
        gStep = 0.01
        Ng = int(np.around((gMax - gMin) / gStep) + 1)
        g = np.linspace(gMin, gMax, Ng)
        return g
          
    def set_modulation_depth(self, sigExp):
        if not (sigExp == []):
            depth = 1.0 - np.mean(sigExp[-10:])
        else:
            depth = 0
        return depth
        
    def spin_ensemble(self, spinA, spinB):
        # Directions of the magnetic field in the frame of spin B
        fieldB = rand_points_on_sphere(self.Ns)
        # Effective g-factors of spin B 
        gB = effective_gfactor(spinB['g'], fieldB)   
        # Quantization axes of spin B
        qB = quantization_axis(spinB['g'], gB, fieldB)
        # Directions of the magnetic field in the frame of spin A
        fieldA = []
        # Effective g-factors of spin A 
        gA = []
        # Quantization axes of spin A
        qA = []
        if spinA['type'] == "isotropic":
            fieldA = fieldB
            gA = effective_gfactor(spinA['g'], fieldA) 
            qA = fieldB
        return [fieldB, gA, gB, qA, qB]
    
    def geometric_parameters(self, var):
        r = []
        xi = []
        phi = []
        if (var['r_width'] > 0):
            if (self.r_distr == 'normal'):
                r = var['r_mean'] + var['r_width'] * np.random.randn(self.Ns)
            elif (self.r_distr == 'uniform'):
                r = var['r_mean'] + var['r_width'] * (0.5 - np.random.rand(self.Ns))
        else:
            r = var['r_mean'] * np.ones(self.Ns)
        if (var['xi_width'] > 0):
            if (self.xi_distr == 'normal'):
                xi = var['xi_mean'] + var['xi_width'] * np.random.randn(self.Ns)
            elif (self.xi_distr == 'uniform'):
                xi = var['xi_mean'] + var['xi_width'] * (0.5 - np.random.rand(self.Ns))
        else:
            xi = var['xi_mean'] * np.ones(self.Ns)	
        if (var['phi_width'] > 0):
            if (self.phi_distr == 'normal'):
                phi = var['phi_mean'] + var['phi_width'] * np.random.randn(self.Ns)
            elif (self.phi_distr == 'uniform'):
                phi = var['phi_mean'] + var['phi_width'] * (0.5 - np.random.rand(self.Ns))
        else:
            phi = var['phi_mean'] * np.ones(self.Ns)
        # Check that all distances are above 0
        for i in range(self.Ns):
            if (r[i] < 0) and (var['r_width'] > 0):
                while True:
                    if (self.r_distr == 'normal'):
                        r[i] = var['r_mean'] + var['r_width'] * np.random.randn()
                    elif (self.r_distr == 'uniform'):
                        r[i] = var['r_mean'] + var['r_width'] * (0.5 - np.random.rand())
                    if (r[i] > 1.0):
                        break
        return [r, xi, phi]    

    def dipolar_frequencies(self, var, spinA, spinB, gA=[], gB=[], qA=[], qB=[], calculateTheta=False):	
        # Parameters of the spin ensemble
        if gA == []:
            gA = self.gA
        if gB == []:
            gB = self.gB
        if qA == []:
            qA = self.qA
        if qB == []:
            qB = self.qB    
        # Geometric parameters
        r, xi, phi = self.geometric_parameters(var)	    
        # Distance vector
        rv_spherical = np.array([np.ones(self.Ns), xi, phi])
        rv_spherical = rv_spherical.T
        rv = spherical2cartesian(rv_spherical)
        # Dipolar frequencies and theta angles
        fdd = np.zeros(self.Ns)
        theta = np.zeros(self.Ns)	
        for i in range(self.Ns):
            if (spinA['type'] == "isotropic") and (spinB['type'] == "isotropic"):
                rv_field = np.dot(rv[i], self.field[i].T)
                rv_A = rv_field
                rv_B = rv_field
                fdd[i] = const['Fdd'] * gA[i] * gB[i] * (1.0 - 3.0 * rv_A * rv_B) / r[i]**3
            elif (spinA['type'] == "isotropic") and (spinB['type'] == "anisotropic"):
                rv_field = np.dot(rv[i], self.field[i].T)
                rv_A = rv_field
                rv_B = np.dot(rv[i], qB[i].T)
                fdd[i] = const['Fdd'] * gA[i] * gB[i] * (1.0 - 3.0 * rv_A * rv_B) / r[i]**3
            if (calculateTheta):
                theta[i] = np.arccos(rv_field) * const['rad2deg']
                if (theta[i] < 0.0):
                    theta[i] = -theta[i]
                if (theta[i] > 90.0):
                    theta[i] = 180.0 - theta[i]
        return [fdd, theta]

    def dipolar_spectrum(self, fdd, weights=[]):	
        # Set the frequency axis	
        fb = np.zeros(self.Nf + 1)
        fb[:-1] = self.f - 0.5 * (self.f[1] - self.f[0]) * np.ones(self.Nf)
        fb[self.Nf] = self.f[self.Nf-1] + 0.5 * (self.f[1] - self.f[0])
        # Set weights
        if (weights == []):
            weights = self.pB
        if (weights == []):
            weights = 1.0 * np.ones(self.Ns)
        # Calculate the spectrum
        spc, bin_edges = np.histogram(fdd, bins=fb, weights=weights) 
        spc = spc.astype(float)
        spc = spc + spc[::-1]
        # Normalize the spectrum
        spc = spc * self.spc_max / np.amax(spc)
        return spc

    def dipolar_timetrace(self, fdd, weights=[]):
        # Set weights
        if (weights == []):
            weights = self.pB
        if (weights == []):
            weights = 1.0 * np.ones(self.Ns)    
        # Set the modulation depth
        if self.depth:
            total_depth = self.depth
        else:
            total_depth = np.sum(weights) / float(self.Ns)
        individual_depths = total_depth * weights / np.sum(weights)      
        # Calculate the time trace
        sig = (1-total_depth) * np.ones(self.Nt)
        for i in range(self.Nt):
            for j in range(self.Ns):
                sig[i] += individual_depths[j] * np.cos(2.0 * np.pi * fdd[j] * self.t[i])
        return sig

    def dipolar_timetrace_fast(self, fdd, weights=[]):
        # Set weights
        if (weights == []):
            weights = self.pB
        if (weights == []):
            weights = 1.0 * np.ones(self.Ns) 
        # Create frequency bins
        fmax = max(np.amax(fdd), np.abs(np.amin(fdd)))
        df = 0.01
        Nf = int(fmax // df) + 1
        fb = np.linspace(-float(Nf)*df, float(Nf)*df, 2*Nf+1)
        fv = np.linspace(-(float(Nf)-0.5)*df, (float(Nf)-0.5)*df, 2*Nf)
        # Calculate the distribution of dipolar frequencies
        pf, f = np.histogram(fdd, bins=fb, weights=weights)
        # Set the modulation depth
        if self.depth:
            total_depth = self.depth
        else:
            total_depth = np.sum(weights) / float(self.Ns)
        individual_depths = total_depth * pf / np.sum(pf)  
        # Calculate the time trace
        sig = (1-total_depth) * np.ones(self.Nt)
        for i in range(self.Nt):
            for j in range(2*Nf):
                sig[i] += individual_depths[j] * np.cos(2.0 * np.pi * fv[j] * self.t[i])
        return sig

    def dipolar_spectrum_vs_theta(self, fdd, theta, weights=[]):
        # Set the frequency axis
        Nf = self.Nf
        fMin = self.f[0]
        fStep = self.f[1] - self.f[0]
        # Set parameters of the theta axis
        Ntheta = self.theta_bins.size
        thetaMin = self.theta_bins[0]
        thetaStep = self.theta_bins[1] - self.theta_bins[0]
        # Set weights
        if (weights == []):
            weights = self.pB
        if (weights == []):
            weights = 1.0 * np.ones(self.Ns)
        #  Calculate a dipolar spectrum vs theta
        spc_vs_theta = np.zeros((Ntheta,Nf))
        for i in range(self.Ns):
            idx_f1 = int(np.around((fdd[i] - fMin) / fStep))
            idx_f2 = int(np.around((-fdd[i] - fMin) / fStep))
            idx_theta = int(np.around((theta[i] - thetaMin) / thetaStep))
            spc_vs_theta[idx_theta][idx_f1] += weights[i]
            spc_vs_theta[idx_theta][idx_f2] += weights[i]
        spc_vs_theta = spc_vs_theta * self.spc_max / np.amax(spc_vs_theta)
        return spc_vs_theta
    
    def dipolar_spectrum_vs_xi(self, var, spinA, spinB):
        Nf = self.Nf
        Nxi = self.xi_bins.size
        spc_vs_xi = np.zeros((Nxi, Nf))
        for i in range(Nxi):
            sys.stdout.write('\r')
            status = int(float(i+1) / float(Nxi) * 100)
            sys.stdout.write("Calculating dipolar spectrum vs xi... %d%% " % (status))
            sys.stdout.flush()
            varNew = copy.deepcopy(var)
            varNew['xi_mean'] = self.xi_bins[i] * const['deg2rad']
            fdd, theta = self.dipolar_frequencies(varNew, spinA, spinB)
            spc = self.dipolar_spectrum(fdd)
            spc_vs_xi[i] = spc
        return spc_vs_xi

    def dipolar_spectrum_vs_phi(self, var, spinA, spinB):
        Nf = self.Nf
        Nphi = self.phi_bins.size
        spc_vs_phi = np.zeros((Nphi, Nf))
        for i in range(Nphi):
            sys.stdout.write('\r')
            status = int(float(i+1) / float(Nphi) * 100)
            sys.stdout.write("Calculating dipolar spectrum vs phi... %d%% " % (status))
            sys.stdout.flush()
            varNew = copy.deepcopy(var)
            varNew['phi_mean'] = self.phi_bins[i] * const['deg2rad']
            fdd, theta = self.dipolar_frequencies(varNew, spinA, spinB)
            spc = self.dipolar_spectrum(fdd)
            spc_vs_phi[i] = spc
        return spc_vs_phi	

    def dipolar_spectrum_vs_temp(self, fdd, spinA, spinB, calcSettings):
        Nf = self.Nf
        Ntemp = self.temp_bins.size
        spc_vs_temp = np.zeros((Ntemp, Nf))
        g = self.set_gaxis(spinB['g'])
        Ng = g.size
        depth_vs_temp = np.zeros((Ntemp, Ng))
        for i in range(Ntemp):
            sys.stdout.write('\r')
            status = int(float(i+1) / float(Ntemp) * 100)
            sys.stdout.write("Calculating dipolar spectrum vs temperature... %d%% " % (status))
            sys.stdout.flush()
            temp = self.temp_bins[i]
            pB = flip_probabilities(self.gB, spinB['g'], calcSettings['magnetic_field'], temp)
            spc = self.dipolar_spectrum(fdd, pB)
            spc_vs_temp[i] = spc
            pB = flip_probabilities(g, spinB['g'], calcSettings['magnetic_field'], temp)
            depth_vs_temp[i] = pB           
        return [spc_vs_temp, g, depth_vs_temp]

    def run_simulation(self, simSettings, expData, spinA, spinB, calcSettings):
        sys.stdout.write('Starting simulation...\n')
        sys.stdout.write('Running pre-calculations... ')
        # Set the frequency axis
        if simSettings['modes']['spc'] \
            or simSettings['modes']['spc_vs_theta'] or simSettings['modes']['spc_vs_xi'] or simSettings['modes']['spc_vs_phi'] or simSettings['modes']['spc_vs_temp']:
            self.f = self.set_faxis(expData['f'], spinA['g'], spinB['g'], simSettings['variables'])
            self.Nf = self.f.size
            # Normalize the frequency axis
            if simSettings['settings']['faxis_normalized']:
                self.fn = self.normalize_faxis(simSettings['variables'])
        # Set the time axis
        if (simSettings['modes']['timetrace']):
            self.t = self.set_taxis(expData['t'], spinA['g'], spinB['g'], simSettings['variables'])
            self.Nt = self.t.size
        # Set the modulation depth parameter
        if (simSettings['modes']['timetrace']):
            self.depth = self.set_modulation_depth(expData['sig'])
        # Generate the ensemble of spin pairs
        self.field, self.gA, self.gB, self.qA, self.qB = self.spin_ensemble(spinA, spinB)
        # Calculate the dipolar frequencies       
        if simSettings['modes']['spc'] or simSettings['modes']['timetrace'] or simSettings['modes']['spc_vs_theta'] or simSettings['modes']['spc_vs_temp']:
            fdd, theta = self.dipolar_frequencies(simSettings['variables'], spinA, spinB, self.gA, self.gB, self.qA, self.qB, simSettings['modes']['spc_vs_theta'])
        # Calculate weights for different g-values of spin B
        if calcSettings['g_selectivity'] and (spinB['type'] == "anisotropic"):
            self.pB = flip_probabilities(self.gB, spinB['g'], calcSettings['magnetic_field'], simSettings['variables']['temp'])
        sys.stdout.write('[DONE]\n')    
        # Calculate the spectrum
        if (simSettings['modes']['spc']):
            sys.stdout.write('Calculating dipolar spectrum... ')
            #time_start = time.time()
            self.spc = self.dipolar_spectrum(fdd)
            sys.stdout.write('[DONE]\n') 
            if not (expData['f'] == []):
                score = rmsd(self.spc, expData['spc'], expData['f'], calcSettings['f_min'], calcSettings['f_max'])
                sys.stdout.write("RMSD = %f \n" % score)
            #time_finish = time.time()
            #time_elapsed = str(datetime.timedelta(seconds = time_finish - time_start))
            #sys.stdout.write('Calculation time: %s\n' % (time_elapsed))
        # Calculate the time trace
        if simSettings['modes']['timetrace']:
            sys.stdout.write('Calculating dipolar time trace... ')
            #time_start = time.time()
            #self.sig = self.dipolar_timetrace(fdd)
            self.sig = self.dipolar_timetrace_fast(fdd) 
            sys.stdout.write('[DONE]\n')
            if not (expData['t'] == []):
                score = rmsd(self.sig, expData['sig'], expData['t'], calcSettings['t_min'], calcSettings['t_max'])
                sys.stdout.write("RMSD = %f \n" % score)
            #time_finish = time.time()
            #time_elapsed = str(datetime.timedelta(seconds = time_finish - time_start))
            #sys.stdout.write('Calculation time: %s\n' % (time_elapsed))   
        # Calculate the dipolar spectrum vs theta
        if (simSettings['modes']['spc_vs_theta']):
            sys.stdout.write('Calculating dipolar spectrum vs theta... ')
            self.theta_bins = simSettings['settings']['theta_bins']
            self.spc_vs_theta = self.dipolar_spectrum_vs_theta(fdd, theta)
            if not simSettings['modes']['spc']:
                self.spc = self.dipolar_spectrum(fdd)
            sys.stdout.write('[DONE]\n')      
        # Calculate a dipolar spectrum vs xi
        if simSettings['modes']['spc_vs_xi']:
            sys.stdout.write('Calculating dipolar spectrum vs xi... ')
            self.xi_bins = simSettings['settings']['xi_bins']
            self.spc_vs_xi = self.dipolar_spectrum_vs_xi(simSettings['variables'], spinA, spinB)
            sys.stdout.write('[DONE]\n')
        # Calculate a dipolar spectrum vs phi
        if simSettings['modes']['spc_vs_phi']:
            sys.stdout.write('Calculating dipolar spectrum vs phi... ')
            self.phi_bins = simSettings['settings']['phi_bins']
            self.spc_vs_phi = self.dipolar_spectrum_vs_phi(simSettings['variables'], spinA, spinB)
            sys.stdout.write('[DONE]\n')
        # Calculate a dipolar spectrum vs temperature
        #if simSettings['modes']['spc_vs_temp'] and calcSettings['g_selectivity'] and (spinB['type'] == "anisotropic"):
        if simSettings['modes']['spc_vs_temp']:
            sys.stdout.write('Calculating dipolar spectrum vs temperature... ')
            self.temp_bins = simSettings['settings']['temp_bins']
            self.spc_vs_temp, self.g, self.depth_vs_temp = self.dipolar_spectrum_vs_temp(fdd, spinA, spinB, calcSettings)
            sys.stdout.write('[DONE]\n')
        sys.stdout.write('Simulation is finished\n\n')

    def init_fitting(self, fitSettings, expData, spinA, spinB, calcSettings):
        if (fitSettings['settings']['fitted_data'] == 'spectrum'):
            # Set the frequency axis
            self.f = self.set_faxis(expData['f'])
            self.Nf = self.f.size     
        elif (fitSettings['settings']['fitted_data'] == 'timetrace'):
            # Set the time axis
            self.t = self.set_taxis(expData['t'])
            self.Nt = self.t.size
            # Set the modulation depth
            self.depth = self.set_modulation_depth(expData['sig'])
        # Generate the ensemble of spin pairs
        self.field, self.gA, self.gB, self.qA, self.qB = self.spin_ensemble(spinA, spinB)
        # Calculate weights for different g-values of spin B
        if (fitSettings['variables']['indices']['temp'] == -1) and calcSettings['g_selectivity'] and (spinB['type'] == "anisotropic"):
            self.pB = flip_probabilities(self.gB, spinB['g'], calcSettings['magnetic_field'], fitSettings['variables']['fixed']['temp'])
