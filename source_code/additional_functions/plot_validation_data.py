'''
Plots the validation data of DipFit
'''

import os
parent_directory = os.path.dirname(os.getcwd())
import sys
sys.path.append(parent_directory)
import io
import libconf
import wx
import numpy as np
from input.read_config import read_fitting_settings, read_validation_settings
from fitting.graphics.plot_score_vs_par import plot_score_vs_par
from supplement.keep_figures_live import keep_figures_live
from supplement.constants import const
from matplotlib import rcParams
rcParams['xtick.direction'] = 'out'
rcParams['ytick.direction'] = 'out'
rcParams['axes.facecolor']= 'white'
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = 'Arial'
rcParams['lines.linewidth'] = 2
rcParams['xtick.major.size'] = 8
rcParams['xtick.major.width'] = 1.5
rcParams['ytick.major.size'] = 8
rcParams['ytick.major.width'] = 1.5
rcParams['font.size'] = 24



validation_files = [
    {'filename': 'parameter_errors-r_mean-r_width.dat',     '2d': True, 'x1': const['variableNames'][0], 'x2': const['variableNames'][1]},
    {'filename': 'parameter_errors-xi_mean-xi_width.dat',   '2d': True, 'x1': const['variableNames'][2], 'x2': const['variableNames'][3]},
    {'filename': 'parameter_errors-phi_mean-phi_width.dat', '2d': True, 'x1': const['variableNames'][4], 'x2': const['variableNames'][5]},
    {'filename': 'parameter_errors-xi_mean-phi_mean.dat',   '2d': True, 'x1': const['variableNames'][2], 'x2': const['variableNames'][4]},
    {'filename': 'parameter_errors-xi_mean-r_width.dat',    '2d': True, 'x1': const['variableNames'][2], 'x2': const['variableNames'][1]},
    {'filename': 'parameter_errors-xi_width-r_width.dat',   '2d': True, 'x1': const['variableNames'][3], 'x2': const['variableNames'][1]},
    {'filename': 'parameter_errors-phi_mean-r_width.dat',   '2d': True, 'x1': const['variableNames'][4], 'x2': const['variableNames'][1]},
    {'filename': 'parameter_errors-phi_width-r_width.dat',  '2d': True, 'x1': const['variableNames'][5], 'x2': const['variableNames'][1]},
    {'filename': 'parameter_errors-temp.dat',               '2d': False,'x1': const['variableNames'][6], 'x2': ''                        },
]


def get_path(message):
    app = wx.App(None) 
    dialog = wx.FileDialog(None, message, wildcard='*.*', style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
    if dialog.ShowModal() == wx.ID_OK:
        path = dialog.GetPath()
    return path


def read_validation_data(dir):
    nfiles = len(validation_files)
    par = []
    score_vs_par = []
    for i in range(nfiles):
        filename = dir + validation_files[i]['filename']
        # check if the file exists
        if os.path.isfile(filename):
            # set parameters
            if validation_files[i]['2d'] == False:
                single_par = [validation_files[i]['x1']]
            else:
                single_par = [validation_files[i]['x1'], validation_files[i]['x2']]
            par.append(single_par)
            # read the data from the file
            data = np.genfromtxt(filename, skip_header=1)
            nx = data.shape[0]
            single_score_vs_par = {}
            if validation_files[i]['2d'] == False:
                name1 = validation_files[i]['x1']
                single_score_vs_par[name1] = [data[j][0] * const['variableScales'][name1] for j in range(nx)]
                single_score_vs_par['score'] = [data[j][1] for j in range(nx)]
            else:
                name1 = validation_files[i]['x1']
                name2 = validation_files[i]['x2']
                single_score_vs_par[name1] = [data[j][0] * const['variableScales'][name1]  for j in range(nx)]
                single_score_vs_par[name2] = [data[j][1] * const['variableScales'][name2]  for j in range(nx)]
                single_score_vs_par['score'] = [data[j][2] for j in range(nx)]
            score_vs_par.append(single_score_vs_par)
    return [par, score_vs_par]

	
if __name__ == '__main__':
	# Read the config file
    configPath = get_path("Open the config file...")
    mode = 1
    expData = {'t': [1], 'sig': [1], 'f': [1], 'spc': [1]}
    fitSettings = {}
    valSettings = {}
    with io.open(configPath) as file:
        config = libconf.load(file)
        fitSettings = read_fitting_settings(config, expData)
        valSettings = read_validation_settings(config, mode)
    # Read the validation files
    dataPath = get_path("Open the directory with validation data...")
    dataDir = os.path.dirname(dataPath) + '/'
    par, score_vs_par = read_validation_data(dataDir)
    # Determine the minimal score
    scores_min = []
    for i in range(len(par)):
        scores_min.append(np.amin(score_vs_par[i]['score']))
    scores_min = np.array(scores_min) 
    score_min = np.amin(scores_min)   
    score_threshold = valSettings['threshold'] * score_min
    # Plot the score in depedence of individual fitting parameters
    filename = dataDir + 'parameter_errors_new.png'
    plot_score_vs_par(score_vs_par, par, fitSettings, valSettings['display_threshold'], score_threshold, True, filename)	
    keep_figures_live()
    