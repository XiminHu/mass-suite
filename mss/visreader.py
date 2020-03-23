import pyteomics
from pyteomics import mzml, auxiliary
import matplotlib.pyplot as plt
import numpy as np
import math
import plotly.graph_objects as go
import re
from scipy.integrate import simps
import pandas as pd

import peakutils
from peakutils.plot import plot as pplot
from matplotlib import pyplot


#Base module
def mz_locator(input_list, mz, error):
    '''
    Find specific mzs from given mz and error range
    input list: mz list
    '''
    target_mz = []
    target_index = []
    
    lower_mz = mz - error
    higher_mz = mz + error

    for i, mzs in enumerate(input_list):
        if mzs < lower_mz:
            continue
        elif mzs >= lower_mz:
            if mzs <= higher_mz:
                target_mz.append(mzs)
                target_index.append(i)
        elif mzs > higher_mz:
                target_mz = 0
                target_index = 'NA'
                break
        
    return target_mz, target_index


#Visualization module 
def tic_plot(spectrum, interactive=True):
    '''
    Static tic plot function
    '''
    time=[]
    TIC=[]
    for i in range(len(spectrum)):
        time.append(spectrum[i]['scanList']['scan'][0]['scan start time'])
        TIC.append(spectrum[i]['total ion current'])
    
    if interactive == True:
        fig = go.Figure([go.Scatter(x=time, y=TIC,
                        hovertemplate = 'Int: %{y}' + '<br>RT: %{x}minute<br>')])

        fig.update_layout(
            template = 'simple_white',
            width = 1000,
            height = 600,
            xaxis = {'title':'Retention Time (min)'},
            yaxis = dict(
                showexponent = 'all',
                exponentformat = 'e',
                title = 'Intensity'))

        fig.show()
    
    elif interactive == False:
        plt.figure(figsize=(10,6))
        plt.plot(time,TIC)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.xlabel('RT (min)')
        plt.ylabel('TIC')
        plt.title('TIC spectrum')
        plt.show()
    
    return


def ms_plot(spectrum, time, interactive=False):
    '''
    Interactive spectrum plot with nearest retention time from the given time
    '''
    for i in range(len(spectrum)):
        if spectrum[i]['scanList']['scan'][0]['scan start time'] >= time:
            mz = f[i]['m/z array']
            ints = f[i]['intensity array']
            rt = spectrum[i]['scanList']['scan'][0]['scan start time']
            break
    
    if interactive == True:        
        fig = go.Figure([go.Bar(x=mz, y=ints, marker_color = 'red', width = 0.5,
                        hovertemplate =
                        'Int: %{y}'+
                        '<br>m/z: %{x}<br>')])
        fig.update_layout(
                title_text=str(round(rt, 3)) + ' min MS1 spectrum, input '+ str(time) + ' min',
                template = 'simple_white',
                width = 1000,
                height = 600,
                xaxis = {'title':'m/z ratio'},
                yaxis = dict(
                    showexponent = 'all',
                    exponentformat = 'e',
                    title = 'Intensity'))
        fig.show()
    
    elif interactive == False:
        plt.figure(figsize=(10,5))
        plt.bar(mz, ints, width = 1.0)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.xlabel('m/z')
        plt.ylabel('Intensity')
        plt.title('MS1 spectrum')
    
    return


def formula_mass(input_formula, mode = 'pos'):
    '''
    sudo code:
    convert input string into a list with element:number structure
    convert all the element into upper case
    match the string list into a given list of element weight
    add adduct/delete H according to mode -- also have neutral mode
    '''
    #Define a list
    elist = {'C': 12, 
            'H':1.00782,
            'N':14.0031,
            'O':15.9949,
            'S':31.9721,
            'P':30.973763,
            'e':0.0005485799}
    
    mol_weight = 0
    parsed_formula = re.findall(r'([A-Z][a-z]*)(\d*)', input_formula)
    for element_count in parsed_formula:
        element = element_count[0]
        count = element_count[1]
        if count == '':
            count = 1
        
        mol_weight += elist[element]*float(count)
    
    if mode == 'pos':
        mol_weight += elist['e'] + elist['H']
    elif mode == 'neg':
        mol_weight -= elist['e'] + elist['H']
    else:
        pass
    
    return mol_weight


def ms_chromatogram(ms_file, input_mz, error, smooth=False, mode='pos', interactive=True):
    '''
    Interactive chromatogram for selected m/z
    '''
    if type(input_mz) == float:
        pass
    elif type(input_mz) == int:
        pass
    elif type(input_mz) == str:
        input_mz = formula_mass(input_mz, mode)
    else:
        print('Cant recognize input type!')
    
    
    
    retention_time = []
    intensity = []
    for i in range(len(ms_file)):
        #print(i)
        retention_time.append(ms_file[i]['scanList']['scan'][0]['scan start time'])
        
        target_mz, target_index = mz_locator(ms_file[i]['m/z array'], input_mz, error)
        if target_index == 'NA':
            intensity.append(0)
        else:
            intensity.append(sum(ms_file[i]['intensity array'][target_index]))

    def peak_smooth(input_list, baseline=500):
        for i, int_ in enumerate(input_list):
            if i > 1 and i < len(input_list)-3:
                if int_ > baseline:
                    for index in np.arange(i+1,i+3):
                        if input_list[index] == 0:
                            input_list[index] = (input_list[index-1]+input_list[index+1])/2
                        else:
                            continue
            
    if smooth == True:
        peak_smooth(intensity)

    if interactive == True:
        fig = go.Figure([go.Scatter(x=retention_time, y=intensity,
                    hovertemplate = 'Int: %{y}' + '<br>RT: %{x}minute<br>')])

        fig.update_layout(
            title_text=str(round(input_mz, 2)) + ' chromatogram, error '+ str(error),
            template = 'simple_white',
            width = 1000,
            height = 600,
            xaxis = {'title':'Retention Time (min)'},
            yaxis = dict(
                showexponent = 'all',
                exponentformat = 'e',
                title = 'Intensity'))

        fig.show()
    elif interactive == False:
        plt.figure(figsize=(20,10))
        plt.plot(retention_time, intensity)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.xlabel('m/z')
        plt.ylabel('Intensity')
        plt.title('MS1 spectrum')
        plt.xlim(0,retention_time[-1])
        plt.ylim(0,)
        plt.show()
    return


#Intergration module
def ms_chromatogram_list(ms_file, input_mz, error, baseline = 5000):
    '''
    Generate a peak list for specific input_mz over whole rt period from the mzml file
    ***Most useful function!
    '''
    retention_time = []
    intensity = []
    for i in range(len(ms_file)):
        #print(i)
        retention_time.append(ms_file[i]['scanList']['scan'][0]['scan start time'])
        
        target_mz, target_index = mz_locator(ms_file[i]['m/z array'], input_mz, error)
        if target_index == 'NA':
            intensity.append(0)
        else:
            intensity.append(sum(ms_file[i]['intensity array'][target_index]))
    
    for i, ints in enumerate(intensity):
        if ints < baseline:
            intensity[i] = 0
            
    return retention_time, intensity


def mz_locator(input_list, mz, error):
    '''
    Find specific mzs from given mz and error range
    input list: mz list
    '''
    target_mz = []
    target_index = []
    
    lower_mz = mz - error
    higher_mz = mz + error

    for i, mzs in enumerate(input_list):
        if mzs < lower_mz:
            continue
        elif mzs >= lower_mz:
            if mzs <= higher_mz:
                target_mz.append(mzs)
                target_index.append(i)
        elif mzs > higher_mz:
                target_mz = 0
                target_index = 'NA'
                break
        
    return target_mz, target_index


def mz_locator(input_list, mz, error):
    '''
    Find specific mzs from given mz and error range
    input list: mz list
    '''
    target_mz = []
    target_index = []
    
    lower_mz = mz - error
    higher_mz = mz + error

    for i, mzs in enumerate(input_list):
        if mzs < lower_mz:
            continue
        elif mzs >= lower_mz:
            if mzs <= higher_mz:
                target_mz.append(mzs)
                target_index.append(i)
        elif mzs > higher_mz:
                target_mz = 0
                target_index = 'NA'
                break
        
    return target_mz, target_index