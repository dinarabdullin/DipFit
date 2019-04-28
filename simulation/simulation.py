'''
Simulation of RIDME spectra and corresponding time traces
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
from spinphysics.quantisation_axis import quantisation_axis
from spinphysics.flip_probabilities import flip_probabilities
from spinphysics.gfactor_hs_iron import gfactor_hs_iron
from supplement.constants import const

class Simulator:

    def __init__(self, calc):
        self.Ns = calc['Ns']
        self.spc_max = calc['spc_max'] 
        self.r_distr = calc['r_distr']
        self.xi_distr = calc['xi_distr']
        self.phi_distr = calc['phi_distr']   
        self.theta_bins = calc['theta_bins']
        self.xi_bins = calc['xi_bins']
        self.phi_bins = calc['phi_bins']
        self.E_bins = calc['E_bins'] 
        self.temp_bins = calc['temp_bins']
        self.gA = []
        self.gB = []
        self.qA = []
        self.qB = []
        self.weightB = []
        self.Nf = 0
        self.Nt = 0
        self.Ng = 0
        self.f = []
        self.fn = []
        self.t = []
        self.g = []
        self.sig = []
        self.pB = 0
        self.spc = []
        self.spc_vs_theta = []
        self.spc_vs_xi = []
        self.spc_vs_phi = []
        self.spc_vs_E = []
        self.spc_vs_temp = []
        self.pB_vs_temp = []

    def set_faxis(self, exp, spinA = {}, spinB = {}, var = []):
        if not (exp['f'] == []):
            f = np.array(exp['f'])
        else:
            # Determine the min value of r
            rMin = 0.0
            if (var['r_width'] > 0):
                if (self.r_distr == 'normal'):
                    r = var['r_mean'] + var['r_width'] * np.random.randn(self.Ns)
                elif (self.r_distr == 'uniform'):
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
                if (self.r_distr == 'normal'):
                    r = var['r_mean'] + var['r_width'] * np.random.randn(self.Ns)
                elif (self.r_distr == 'uniform'):
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

    def set_gaxis(self, gB):
        gBMin = np.amin(gB)
        gBMax = np.amax(gB)
        gBStep = 0.01
        Ng = int(np.around((gBMax - gBMin) / gBStep) + 1)
        g = np.linspace(gBMin, gBMax, Ng)
        return g
          
    def set_modulation_depth(self, sigExp):
        if not (sigExp == []):
            depth = 1.0 - np.mean(sigExp[-50:])
        else:
            depth = 0
        return depth
        
    def precalculations(self, spinA, spinB):
        # Directions of the magnetic field in the lab frame
        eB0_L = rand_points_on_sphere(self.Ns)
        # Directions of the magnetic field in the frame of spin A
        if (spinA['gFrame'].all() == 0):
            eB0_A = eB0_L
        else:
            RM_L2A = euler2RM(spinA['gFrame'])
            eB0_A = np.dot(RM_L2A.T, eB0_L.T).T
        # Directions of the magnetic field in the frame of spin B
        if (spinB['gFrame'].all() == 0):
            eB0_B = eB0_L
        else:
            RM_L2B = euler2RM(spinB['gFrame'])
            eB0_B = np.dot(RM_L2B.T, eB0_L.T).T
        # Effective g factors of spin A 
        gA = effective_gfactor(spinA['g'], eB0_A)
        # Effective g factors of spin B 
        gB = effective_gfactor(spinB['g'], eB0_B)   
        # Quantisation axes of spin A
        qA = eB0_L
        # qA = quantisation_axis(spinA['g'], gA, eB0_L)
        # Quantisation axes of spin B
        qB = quantisation_axis(spinB['g'], gB, eB0_L)
        return [gA, gB, qA, qB]

    def dipolar_frequencies(self, var, gA, gB, qA, qB, calculateTheta=False):	
        # Set the spherical coordinates of the distance vector
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
        # Set parameters of the frequency axis
        Nf = self.Nf
        fMin = self.f[0]
        fStep = self.f[1] - self.f[0]
        # Set parameters of the theta axis
        Ntheta = self.theta_bins.size
        thetaMin = self.theta_bins[0]
        thetaStep = self.theta_bins[1] - self.theta_bins[0]
        # Set weights
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
        Nxi = self.xi_bins.size
        spc_vs_xi = np.zeros((Nxi, Nf))
        for i in range(Nxi):
            sys.stdout.write('\r')
            status = int(float(i+1) / float(Nxi) * 100)
            sys.stdout.write("Calculating dipolar spectrum vs xi... %d%% " % (status))
            sys.stdout.flush()
            varNew = copy.deepcopy(var)
            varNew['xi_mean'] = self.xi_bins[i] * const['deg2rad']
            fdd, theta = self.dipolar_frequencies(varNew, self.gA, self.gB, self.qA, self.qB)
            spc = self.dipolar_spectrum(fdd)
            spc_vs_xi[i] = spc
        return spc_vs_xi

    def dipolar_spectrum_vs_phi(self, var):
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
            fdd, theta = self.dipolar_frequencies(varNew, self.gA, self.gB, self.qA, self.qB)
            spc = self.dipolar_spectrum(fdd)
            spc_vs_phi[i] = spc
        return spc_vs_phi	

    def dipolar_spectrum_vs_E(self, var, exp, spinA, spinB):
        Nf = self.Nf
        Ne = self.E_bins.size
        spc_vs_E = np.zeros((Ne, Nf))
        for i in range(Ne):
            sys.stdout.write('\r')
            status = int(float(i+1) / float(Ne) * 100)
            sys.stdout.write("Calculating dipolar spectrum vs E/D... %d%% " % (status))
            sys.stdout.flush()
            spinBNew = copy.deepcopy(spinB)
            spinBNew['D'][1] = spinBNew['D'][0] * self.E_bins[i]
            spinBNew['g'] = gfactor_hs_iron(spinBNew, exp['mwFreq'])
            gA, gB, qA, qB = self.precalculations(spinA, spinBNew)
            weightB = flip_probabilities(spinBNew['g'], gB, exp['magnField'], var['temp'])
            fdd, theta = self.dipolar_frequencies(var, gA, gB, qA, qB)
            spc = self.dipolar_spectrum(fdd, weightB)
            spc_vs_E[i] = spc
        # # Normalize the integral
        # I = np.sum(spc_vs_E, axis=1)
        # Imin = np.amin(I)
        # for i in range(Ne):
            # spc_vs_E[i] = spc_vs_E[i] * Imin / I[i]
        return spc_vs_E

    def dipolar_spectrum_vs_temp(self, fdd, exp, spinB):
        Nf = self.Nf
        Ntemp = self.temp_bins.size
        spc_vs_temp = np.zeros((Ntemp, Nf))
        g = self.set_gaxis(self.gB)
        Ng = g.size
        self.g = g
        self.Ng = Ng
        pb_vs_temp = np.zeros((Ntemp, Ng))
        gb = np.zeros(Ng + 1)
        gb[:-1] = g - 0.5*(g[1]-g[0])*np.ones(Ng)
        gb[Ng] = g[Ng-1] + 0.5*(g[1]-g[0])
        for i in range(Ntemp):
            sys.stdout.write('\r')
            status = int(float(i+1) / float(Ntemp) * 100)
            sys.stdout.write("Calculating dipolar spectrum vs temperature... %d%% " % (status))
            sys.stdout.flush()
            # Calculate the B-spin flip probability for each temperature
            temp = self.temp_bins[i]
            weightB = flip_probabilities(spinB['g'], self.gB, exp['magnField'], temp)
            # Calculate the spectrum
            spc = self.dipolar_spectrum(fdd, weightB)
            spc_vs_temp[i] = spc
            # Calculate the dependence of flip probability on the g-factor of the B-spin
            pb = flip_probabilities(spinB['g'], self.g, exp['magnField'], temp)
            pb_vs_temp[i] = pb           
        return [spc_vs_temp, pb_vs_temp]

    def dipolar_timetrace(self, fdd, weights=[]):
        # Set weights
        if (weights == []):
            weights = self.weightB
        # Set the modulation depth
        if (self.pB):
            depth = self.pB
            cdepth = (self.pB / np.sum(weights)) * weights
        else:
            depth = np.sum(weights) / float(self.Ns)
            cdepth = (1 / float(self.Ns)) * weights
        # Calculate the time trace
        sig = (1-depth) * np.ones(self.Nt)
        for i in range(self.Nt):
            for j in range(self.Ns):
                sig[i] += cdepth[j] * np.cos(2.0 * np.pi * fdd[j] * self.t[i])
        return sig

    def dipolar_timetrace_fast(self, fdd, weights=[]):
        # Set weights
        if (weights == []):
            weights = self.weightB
        # Set frequency bounds	
        fmax = max(np.amin(fdd), np.abs(np.amin(fdd)))
        df = 0.01
        Nf = int(fmax // df) + 1
        fb = np.linspace(-float(Nf)*df, float(Nf)*df, 2*Nf+1)
        fv = np.linspace(-(float(Nf)-0.5)*df, (float(Nf)-0.5)*df, 2*Nf)
        # Calculate the spectrum
        pf, f = np.histogram(fdd, bins=fb, weights=weights)
        # Set the modulation depth
        if (self.pB):
            depth = self.pB
            cdepth = (self.pB / np.sum(pf)) * pf
        else:
            depth = np.sum(pf) / float(self.Ns)
            cdepth = (1 / float(self.Ns)) * pf
        # Calculate the time trace
        sig = (1-depth) * np.ones(self.Nt)
        for i in range(self.Nt):
            for j in range(2*Nf):
                sig[i] += cdepth[j] * np.cos(2.0 * np.pi * fv[j] * self.t[i])
        return sig

    def run_simulation(self, simPar, exp, spinA, spinB, calc):
        # Display status
        sys.stdout.write('Starting simulation...\n')

        sys.stdout.write('Running pre-calculations... ')
        # Set the frequency axis
        if (simPar['settings']['spc'] or simPar['settings']['spc_vs_theta'] \
            or simPar['settings']['spc_vs_xi'] or simPar['settings']['spc_vs_phi'] \
            or simPar['settings']['spc_vs_E'] or simPar['settings']['spc_vs_temp']):
            self.f = self.set_faxis(exp, spinA, spinB, simPar['variables'])
            self.Nf = self.f.size
            # Normalize the frequency axis
            if simPar['settings']['faxis_normalized']:
                self.fn = self.normalize_faxis(spinA, spinB, simPar['variables'])
        # Set the time axis
        if (simPar['settings']['timetrace']):
            self.t = self.set_taxis(exp, spinA, spinB, simPar['variables'])
            self.Nt = self.t.size
            # Set the modulation depth parameter
            self.pB = self.set_modulation_depth(exp['sig'])
        # Do some pre-calculations
        self.gA, self.gB, self.qA, self.qB = self.precalculations(spinA, spinB)
        # Calculate the dipolar frequencies
        if (simPar['settings']['spc'] or simPar['settings']['spc_vs_theta'] \
            or simPar['settings']['timetrace'] or simPar['settings']['spc_vs_temp']):
            if simPar['settings']['spc_vs_theta']:
                fdd, theta = self.dipolar_frequencies(simPar['variables'], self.gA, self.gB, self.qA, self.qB, True)
            else:
                fdd, theta = self.dipolar_frequencies(simPar['variables'], self.gA, self.gB, self.qA, self.qB, False)
        # Calculate the weights that are assigned to the different components of the dipolar spectrum
        self.weightB = flip_probabilities(spinB['g'], self.gB, exp['magnField'], simPar['variables']['temp'])                
        sys.stdout.write('[DONE]\n')

        # Calculate the dipolar spectrum
        if (simPar['settings']['spc']):
            sys.stdout.write('Calculating dipolar spectrum... ')
            #time_start = time.time()
            self.spc = self.dipolar_spectrum(fdd)
            sys.stdout.write('[DONE]\n') 
            if not (exp['f'] == []):
                score = rmsd(self.spc, exp['spc'], exp['f'], calc['f_min'], calc['f_max'])
                sys.stdout.write("RMSD = %f \n" % score)
            #time_finish = time.time()
            #time_elapsed = str(datetime.timedelta(seconds = time_finish - time_start))
            #sys.stdout.write('Calculation time: %s\n' % (time_elapsed))

        # Calculate a time trace
        if simPar['settings']['timetrace']:
            sys.stdout.write('Calculating dipolar time trace... ')
            #time_start = time.time()
            #self.sig = self.dipolar_timetrace(fdd)
            self.sig = self.dipolar_timetrace_fast(fdd) 
            sys.stdout.write('[DONE]\n')
            if not (exp['t'] == []):
                score = rmsd(self.sig, exp['sig'], exp['t'], calc['t_min'], calc['t_max'])
                sys.stdout.write("RMSD = %f \n" % score)
            #time_finish = time.time()
            #time_elapsed = str(datetime.timedelta(seconds = time_finish - time_start))
            #sys.stdout.write('Calculation time: %s\n' % (time_elapsed))   

        # Calculate the dipolar spectrum vs theta
        if (simPar['settings']['spc_vs_theta']):
            sys.stdout.write('Calculating dipolar spectrum vs theta... ')
            self.spc_vs_theta = self.dipolar_spectrum_vs_theta(fdd, theta)
            if not simPar['settings']['spc']:
                self.spc = self.dipolar_spectrum(fdd)
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

        # Calculate a dipolar spectrum vs temperature
        if simPar['settings']['spc_vs_temp']:
            sys.stdout.write('Calculating dipolar spectrum vs temperature... ')
            self.spc_vs_temp, self.pb_vs_temp = self.dipolar_spectrum_vs_temp(fdd, exp, spinB)
            sys.stdout.write('[DONE]\n')

        sys.stdout.write('Simulation is finished\n\n')

    def init_fitting(self, fitPar, exp, spinA, spinB, calc):
        if calc['fitted_data'] == 'spectrum':
            # Set the frequency axis
            self.f = self.set_faxis(exp)
            self.Nf = self.f.size     
        elif calc['fitted_data'] == 'timetrace':
            # Set the time axis
            self.t = self.set_taxis(exp)
            self.Nt = self.t.size
            # Set the modulation depth parameter
            self.pB = self.set_modulation_depth(exp['sig'])
        # Do some pre-calculations
        self.gA, self.gB, self.qA, self.qB = self.precalculations(spinA, spinB)
        # Calculate the weights that are assigned to the different components of the dipolar spectrum
        if (fitPar['variables']['indices']['temp'] == -1):
            self.weightB = flip_probabilities(spinB['g'], self.gB, exp['magnField'], fitPar['variables']['fixed']['temp'])