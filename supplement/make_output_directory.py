'''
Generates an output directory
'''

import os
import errno
import datetime
import shutil

def make_output_directory(output, configPath):
    if output['save_data']:
        # Set the parent directory
        parentDir = output['directory']
        configDir, configName = os.path.split(os.path.abspath(configPath))
        if (parentDir):
            dir = parentDir
        else:
            dir = configDir
        # Create a folder, which is named according to the current date and time, in the parent directory
        now = datetime.datetime.now()
        folder = now.strftime("%Y-%m-%d_%H-%M")
        dir = dir + folder + "/"
        try:
            os.makedirs(dir)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        # Copy the configuration file into the output directory
        shutil.copy2(configPath, dir+configName)
        output['directory'] = dir