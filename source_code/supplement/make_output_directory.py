'''
Generate an output directory
'''

import os
import errno
import datetime
import shutil


def make_output_directory(outputSettings, configPath):
    if outputSettings['save_data']:
        parentDir = outputSettings['directory']
        configDir, configName = os.path.split(os.path.abspath(configPath))
        if (parentDir):
            dir = parentDir
        else:
            dir = configDir
        now = datetime.datetime.now()
        folder = now.strftime("%Y-%m-%d_%H-%M")
        dir = dir + "/" + folder + "/"
        try:
            os.makedirs(dir)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        shutil.copy2(configPath, dir+configName)
        outputSettings['directory'] = dir