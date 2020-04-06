'''
Simulation of RIDME spectra and RIDME time traces
'''

import sys
import numpy as np
import time
import datetime
import copy
from mathematics.random_points_on_sphere import random_points_on_sphere
from mathematics.spherical2cartesian import spherical2cartesian
from mathematics.rmsd import rmsd
from spinphysics.effective_gfactor import effective_gfactor
from spinphysics.quantization_axis import quantization_axis
from spinphysics.flip_probabilities import flip_probabilities
from supplement.constants import const


class Simulator:

    def __init__(self, calc_settings):
        self.Ns = calc_settings['Ns']
        self.r_distr = calc_settings['r_distr']
        self.xi_distr = calc_settings['xi_distr']
        self.phi_distr = calc_settings['phi_distr']  
        self.spc_max = calc_settings['spc_max']
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
        self.spc_vs_theta = []
        self.spc_vs_xi = []
        self.spc_vs_phi = []
        self.spc_vs_temp = []
        self.sig = [] 
        self.depth = 0
        self.g = []
        self.depth_vs_temp = []    

    def set_faxis(self, f_exp, gA=[], gB=[], parameters={}):
        if not (f_exp == []):
            f = np.array(f_exp)
        else:
            # Determine the minimal value of r
            r_min = 0.0
            if (parameters['r_width'] > 0):
                if (self.r_distr == 'normal'):
                    r = parameters['r_mean'] + parameters['r_width'] * np.random.randn(self.Ns)
                elif (self.r_distr == 'uniform'):
                    r = parameters['r_mean'] + parameters['r_width'] * (0.5 - np.random.rand(self.Ns))
                r_min = np.amin(r) 
            else:
                r_min = parameters['r_mean']
            # Determine the maximal values of g-factors
            gA_max = np.amax(gA)
            gB_max = np.amax(gB)
            # Estimate the max dipolar frequency
            fdd_max = 2 * const['Fdd'] * gA_max * gB_max / r_min**3 
            # Set the frequency axis
            f_max = np.around(fdd_max) + 5.0
            f_min = -f_max
            f_step = 0.1
            Nf = int(np.around((f_max - f_min) / f_step)) + 1
            f = np.linspace(f_min, f_max, Nf)
        return f

    def normalize_faxis(self, parameters):
        # Determine the reference frequency
        fdd_ref = const['Fdd'] * const['ge'] * const['ge'] / parameters['r_mean']**3 
        # Normalize the frequency axis by the reference frequency
        Nf = self.f.size
        fn = np.zeros(Nf)
        for i in range(Nf):
            fn[i] = self.f[i] / fdd_ref
        return fn 
        
    def set_taxis(self, t_exp, gA=[], gB=[], parameters={}):
        if not (t_exp == []):
            t = np.array(t_exp)
        else:
            # Determine the max value of r
            if (parameters['r_width'] > 0):
                if (self.r_distr == 'normal'):
                    r = parameters['r_mean'] + parameters['r_width'] * np.random.randn(self.Ns)
                elif (self.r_distr == 'uniform'):
                    r = parameters['r_mean'] + parameters['r_width'] * (0.5 - np.random.rand(self.Ns))
                r_max = np.amax(r) 
            else:
                r_max = parameters['r_mean']
            # Determine the mininimal values of g-factors
            gA_min = np.amin(gA)
            gB_min = np.amin(gB)
            # Estimate the min dipolar frequency
            fdd_min = const['Fdd'] * gA_min * gB_min / r_max**3
            # Estimate the max period of dipolar oscilation
            tdd_max = 1 / fdd_min
            # Set the time axis
            t_min = 0.0
            t_max = 3 * tdd_max
            t_step = 0.008
            Nt = int(np.around((t_max - t_min) / t_step)) + 1
            t = np.linspace(t_min, t_max, Nt)
        return t       

    def set_gaxis(self, g):
        gMin = np.amin(g)
        gMax = np.amax(g)
        gStep = 0.01
        Ng = int(np.around((gMax - gMin) / gStep) + 1)
        g = np.linspace(gMin, gMax, Ng)
        return g
          
    def set_modulation_depth(self, sig_exp, mod_depth = 0):
        if not (sig_exp == []):
            depth = 1.0 - np.mean(sig_exp[-10:])
        else:
            depth = mod_depth
        return depth
        
    def spin_ensemble(self, spinA, spinB):
        # Directions of the magnetic field in the frame of spin B
        fieldB = random_points_on_sphere(self.Ns)
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
    
    def geometric_parameters(self, parameters):
        r = []
        xi = []
        phi = []
        if (parameters['r_width'] > 0):
            if (self.r_distr == 'normal'):
                r = parameters['r_mean'] + parameters['r_width'] * np.random.randn(self.Ns)
            elif (self.r_distr == 'uniform'):
                r = parameters['r_mean'] + parameters['r_width'] * (0.5 - np.random.rand(self.Ns))
        else:
            r = parameters['r_mean'] * np.ones(self.Ns)
        if (parameters['xi_width'] > 0):
            if (self.xi_distr == 'normal'):
                xi = parameters['xi_mean'] + parameters['xi_width'] * np.random.randn(self.Ns)
            elif (self.xi_distr == 'uniform'):
                xi = parameters['xi_mean'] + parameters['xi_width'] * (0.5 - np.random.rand(self.Ns))
        else:
            xi = parameters['xi_mean'] * np.ones(self.Ns)	
        if (parameters['phi_width'] > 0):
            if (self.phi_distr == 'normal'):
                phi = parameters['phi_mean'] + parameters['phi_width'] * np.random.randn(self.Ns)
            elif (self.phi_distr == 'uniform'):
                phi = parameters['phi_mean'] + parameters['phi_width'] * (0.5 - np.random.rand(self.Ns))
        else:
            phi = parameters['phi_mean'] * np.ones(self.Ns)
        # Check that all distances are above 0
        for i in range(self.Ns):
            if (r[i] < 0) and (parameters['r_width'] > 0):
                while True:
                    if (self.r_distr == 'normal'):
                        r[i] = parameters['r_mean'] + parameters['r_width'] * np.random.randn()
                    elif (self.r_distr == 'uniform'):
                        r[i] = parameters['r_mean'] + parameters['r_width'] * (0.5 - np.random.rand())
                    if (r[i] > 1.0):
                        break
        return [r, xi, phi]    

    def dipolar_frequencies(self, parameters, spinA, spinB, gA=[], gB=[], qA=[], qB=[], calculate_theta=False):	
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
        r, xi, phi = self.geometric_parameters(parameters)	    
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
            if (calculate_theta):
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
        fdd_max = max(np.amax(fdd), np.abs(np.amin(fdd)))
        df = 0.01
        Nf = int(fdd_max // df) + 1
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
        fdd_min = self.f[0]
        fdd_step = self.f[1] - self.f[0]
        # Set parameters of the theta axis
        Ntheta = self.theta_bins.size
        theta_min = self.theta_bins[0]
        theta_step = self.theta_bins[1] - self.theta_bins[0]
        # Set weights
        if (weights == []):
            weights = self.pB
        if (weights == []):
            weights = 1.0 * np.ones(self.Ns)
        #  Calculate a dipolar spectrum vs theta
        spc_vs_theta = np.zeros((Ntheta,Nf))
        for i in range(self.Ns):
            idx_f1 = int(np.around((fdd[i] - fdd_min) / fdd_step))
            idx_f2 = int(np.around((-fdd[i] - fdd_min) / fdd_step))
            idx_theta = int(np.around((theta[i] - theta_min) / theta_step))
            spc_vs_theta[idx_theta][idx_f1] += weights[i]
            spc_vs_theta[idx_theta][idx_f2] += weights[i]
        spc_vs_theta = spc_vs_theta * self.spc_max / np.amax(spc_vs_theta)
        return spc_vs_theta
    
    def dipolar_spectrum_vs_xi(self, parameters, spinA, spinB):
        Nf = self.Nf
        Nxi = self.xi_bins.size
        spc_vs_xi = np.zeros((Nxi, Nf))
        for i in range(Nxi):
            sys.stdout.write('\r')
            status = int(float(i+1) / float(Nxi) * 100)
            sys.stdout.write("Calculating the dipolar spectrum vs xi... %d%% " % (status))
            sys.stdout.flush()
            new_parameters = copy.deepcopy(parameters)
            new_parameters['xi_mean'] = self.xi_bins[i] * const['deg2rad']
            fdd, theta = self.dipolar_frequencies(new_parameters, spinA, spinB)
            spc = self.dipolar_spectrum(fdd)
            spc_vs_xi[i] = spc
        return spc_vs_xi

    def dipolar_spectrum_vs_phi(self, parameters, spinA, spinB):
        Nf = self.Nf
        Nphi = self.phi_bins.size
        spc_vs_phi = np.zeros((Nphi, Nf))
        for i in range(Nphi):
            sys.stdout.write('\r')
            status = int(float(i+1) / float(Nphi) * 100)
            sys.stdout.write("Calculating the dipolar spectrum vs phi... %d%% " % (status))
            sys.stdout.flush()
            new_parameters = copy.deepcopy(parameters)
            new_parameters['phi_mean'] = self.phi_bins[i] * const['deg2rad']
            fdd, theta = self.dipolar_frequencies(new_parameters, spinA, spinB)
            spc = self.dipolar_spectrum(fdd)
            spc_vs_phi[i] = spc
        return spc_vs_phi	

    def dipolar_spectrum_vs_temp(self, fdd, spinA, spinB, calc_settings):
        Nf = self.Nf
        Ntemp = self.temp_bins.size
        spc_vs_temp = np.zeros((Ntemp, Nf))
        g = self.set_gaxis(spinB['g'])
        Ng = g.size
        depth_vs_temp = np.zeros((Ntemp, Ng))
        for i in range(Ntemp):
            sys.stdout.write('\r')
            status = int(float(i+1) / float(Ntemp) * 100)
            sys.stdout.write("Calculating the dipolar spectrum vs temperature... %d%% " % (status))
            sys.stdout.flush()
            temp = self.temp_bins[i]
            pB = flip_probabilities(self.gB, spinB['g'], calc_settings['magnetic_field'], temp)
            spc = self.dipolar_spectrum(fdd, pB)
            spc_vs_temp[i] = spc
            pB = flip_probabilities(g, spinB['g'], calc_settings['magnetic_field'], temp)
            depth_vs_temp[i] = pB           
        return [spc_vs_temp, g, depth_vs_temp]

    def run_simulation(self, sim_settings, exp_data, spinA, spinB, calc_settings):
        sys.stdout.write('Starting the simulation...\n')
        sys.stdout.write('Running the pre-calculations... ')
        # Set the frequency axis
        if (sim_settings['modes']['spc'] or sim_settings['modes']['spc_vs_theta'] or 
            sim_settings['modes']['spc_vs_xi'] or sim_settings['modes']['spc_vs_phi'] or 
            sim_settings['modes']['spc_vs_temp']):
            self.f = self.set_faxis(exp_data['f'], spinA['g'], spinB['g'], sim_settings['parameters'])
            self.Nf = self.f.size
            # Normalize the frequency axis
            if sim_settings['settings']['faxis_normalized']:
                self.fn = self.normalize_faxis(sim_settings['parameters'])
        # Set the time axis
        if sim_settings['modes']['timetrace']:
            self.t = self.set_taxis(exp_data['t'], spinA['g'], spinB['g'], sim_settings['parameters'])
            self.Nt = self.t.size
        # Set the modulation depth parameter
        if sim_settings['modes']['timetrace']:
            self.depth = self.set_modulation_depth(exp_data['sig'], sim_settings['settings']['mod_depth'])
        # Generate the ensemble of spin pairs
        self.field, self.gA, self.gB, self.qA, self.qB = self.spin_ensemble(spinA, spinB)
        # Calculate the dipolar frequencies       
        if (sim_settings['modes']['spc'] or sim_settings['modes']['timetrace'] or 
            sim_settings['modes']['spc_vs_theta'] or sim_settings['modes']['spc_vs_temp']):
            fdd, theta = self.dipolar_frequencies(sim_settings['parameters'], spinA, spinB, self.gA, self.gB, self.qA, self.qB, sim_settings['modes']['spc_vs_theta'])
        # Calculate weights for different g-values of spin B
        if (calc_settings['g_selectivity'] and spinB['type'] == "anisotropic"):
            self.pB = flip_probabilities(self.gB, spinB['g'], calc_settings['magnetic_field'], sim_settings['parameters']['temp'])
        sys.stdout.write('[DONE]\n')    
        # Calculate the spectrum
        if sim_settings['modes']['spc']:
            sys.stdout.write('Calculating the dipolar spectrum... ')
            self.spc = self.dipolar_spectrum(fdd)
            sys.stdout.write('[DONE]\n') 
            if not (exp_data['f'] == []):
                score = rmsd(self.spc, exp_data['spc'], exp_data['f'], calc_settings['f_min'], calc_settings['f_max'])
                sys.stdout.write("RMSD = %f \n" % score)
        # Calculate the time trace
        if sim_settings['modes']['timetrace']:
            sys.stdout.write('Calculating the dipolar time trace... ')
            #self.sig = self.dipolar_timetrace(fdd)
            self.sig = self.dipolar_timetrace_fast(fdd) 
            sys.stdout.write('[DONE]\n')
            if not (exp_data['t'] == []):
                score = rmsd(self.sig, exp_data['sig'], exp_data['t'], calc_settings['t_min'], calc_settings['t_max'])
                sys.stdout.write("RMSD = %f \n" % score) 
        # Calculate the dipolar spectrum vs theta
        if (sim_settings['modes']['spc_vs_theta']):
            sys.stdout.write('Calculating the dipolar spectrum vs theta... ')
            self.theta_bins = sim_settings['settings']['theta_bins']
            self.spc_vs_theta = self.dipolar_spectrum_vs_theta(fdd, theta)
            if not sim_settings['modes']['spc']:
                self.spc = self.dipolar_spectrum(fdd)
            sys.stdout.write('[DONE]\n')      
        # Calculate a dipolar spectrum vs xi
        if sim_settings['modes']['spc_vs_xi']:
            sys.stdout.write('Calculating the dipolar spectrum vs xi... ')
            self.xi_bins = sim_settings['settings']['xi_bins']
            self.spc_vs_xi = self.dipolar_spectrum_vs_xi(sim_settings['parameters'], spinA, spinB)
            sys.stdout.write('[DONE]\n')
        # Calculate a dipolar spectrum vs phi
        if sim_settings['modes']['spc_vs_phi']:
            sys.stdout.write('Calculating the dipolar spectrum vs phi... ')
            self.phi_bins = sim_settings['settings']['phi_bins']
            self.spc_vs_phi = self.dipolar_spectrum_vs_phi(sim_settings['parameters'], spinA, spinB)
            sys.stdout.write('[DONE]\n')
        # Calculate a dipolar spectrum vs temperature
        if (sim_settings['modes']['spc_vs_temp'] and calc_settings['g_selectivity']):
            sys.stdout.write('Calculating the dipolar spectrum vs temperature... ')
            self.temp_bins = sim_settings['settings']['temp_bins']
            self.spc_vs_temp, self.g, self.depth_vs_temp = self.dipolar_spectrum_vs_temp(fdd, spinA, spinB, calc_settings)
            sys.stdout.write('[DONE]\n')
        sys.stdout.write('The simulation is finished\n\n')

    def init_fitting(self, fit_settings, exp_data, spinA, spinB, calc_settings):
        if fit_settings['settings']['fitted_data'] == 'spectrum':
            # Set the frequency axis
            self.f = self.set_faxis(exp_data['f'])
            self.Nf = self.f.size     
        elif fit_settings['settings']['fitted_data'] == 'timetrace':
            # Set the time axis
            self.t = self.set_taxis(exp_data['t'])
            self.Nt = self.t.size
            # Set the modulation depth
            self.depth = self.set_modulation_depth(exp_data['sig'])
        # Generate the ensemble of spin pairs
        self.field, self.gA, self.gB, self.qA, self.qB = self.spin_ensemble(spinA, spinB)
        # Calculate weights for different g-values of spin B
        if (fit_settings['parameters']['indices']['temp'] == -1 and calc_settings['g_selectivity'] and spinB['type'] == "anisotropic"):
            self.pB = flip_probabilities(self.gB, spinB['g'], calc_settings['magnetic_field'], fit_settings['parameters']['fixed']['temp'])
