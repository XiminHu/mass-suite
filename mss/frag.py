import pandas as pd
import urllib
import urllib.request
import json
import pymzml
import numpy as np
import matplotlib.pyplot as plt
import plotly

def massbank(frag_list, minsimilarity=500):
    url = 'https://mona.fiehnlab.ucdavis.edu/rest/similarity/search'
    values = {"spectrum" : frag_list, "minSimilarity" : minsimilarity}

    data_str = json.dumps(values)
    data = data_str.encode()
    #print(data) 
    headers = {
      "content-type": "application/json",
    }
    req = urllib.request.Request(url, data, headers)
    response = urllib.request.urlopen(req)
    json_result = json.loads(response.read())
    
    result_list=[]
    for result in json_result:
        name = result['hit']['compound'][0]['names'][0]['name']
        cpd_meta = result['hit']['compound'][0]['metaData']
        try:
            formula = [i for i in cpd_meta if i['name']=='molecular formula' and i['category']=='none'][0]['value']
        except:
            formula = "NA"
        try:
            inchi = result['hit']['compound'][0]['inchiKey']
        except:
            inchi = "NA"
        ms_meta = result['hit']['metaData']
        try:
            instrument = [i for i in ms_meta if i['name']=='instrument type'][0]['value']
        except:
            instrument = "NA"
        mode = [i for i in ms_meta if i['name']=='ionization mode'][0]['value']
        try:
            energy = [i for i in ms_meta if i['name']=='collision energy'][0]['value']
        except:
            energy = "NA"
        try:
            p_type = [i for i in ms_meta if i['name']=='precursor type'][0]['value']
        except:
            p_type = "NA"
        score = result['score']
        link = 'https://mona.fiehnlab.ucdavis.edu/spectra/display/' + result['hit']['id']
        result_list.append([name, formula, inchi, instrument, energy, mode, p_type, score, link])
    
    d_result = pd.DataFrame(result_list,columns=['name', 'formula', 'inchiKey', 'instrument', 'CE', 'ionization mode', 'precursor type', 'score', 'url'])
    
    def make_clickable(val):
    # target _blank to open new window
        return '<a target="_blank" href="{}">{}</a>'.format(val, val)

    d_result = d_result.style.format({'url': make_clickable})
    
    return d_result

def frag_massbank(mzml_scans, precursor, error=20, noise_thr=50, scan_index=0):
    p_range_l = precursor * (1 - error * 1e-6)
    p_range_h = precursor * (1 + error * 1e-6)
    frag_scan = []
    for scan in mzml_scans:
        if scan.ms_level == 2:
            precursor = scan.selected_precursors[0]['mz']
            p_intensity = scan.selected_precursors[0]['i']
            if precursor < p_range_h and precursor > p_range_l:
                frag_scan.append([precursor, p_intensity, scan])
    frag_scan.sort(key=lambda x: x[1], reverse=True)
    if len(frag_scan) != 0:
        print('Now showing index', scan_index, 'of', str(len(frag_scan)), 'total found scans')
        plot_scan = frag_scan[scan_index][2]
        drop_index = np.argwhere(plot_scan.i <= noise_thr)
        plot_scan.i = np.delete(plot_scan.i, drop_index)
        plot_scan.mz = np.delete(plot_scan.mz, drop_index)
        mz = plot_scan.mz
        ints = plot_scan.i

    for i in range(len(mz)):
        if i == 0:
            list_string = str(round(mz[i],4)) + ':' + str(round(ints[i],1)) + ' '
        else:
            list_string += str(round(mz[i],4)) + ':' + str(round(ints[i],1)) + ' '

    d_result = massbank(list_string)

    return d_result

def frag_comp(mzml_scans, precursor, error=20, scan_index=0, mona_index=0, noise_thr = 50, interactive=False, search=False, source='MoNA'): 
    '''
    Interactive spectrum plot with nearest retention time from the given scan
    mzml_scans: mzfile
    time: selected time for the scan
    '''
    p_range_l = precursor * (1 - error * 1e-6)
    p_range_h = precursor * (1 + error * 1e-6)
    frag_scan = []
    for scan in mzml_scans:
        if scan.ms_level == 2:
            precursor = scan.selected_precursors[0]['mz']
            p_intensity = scan.selected_precursors[0]['i']
            if precursor < p_range_h and precursor > p_range_l:
                frag_scan.append([precursor, p_intensity, scan])
    frag_scan.sort(key=lambda x: x[1], reverse=True)
    if len(frag_scan) != 0:
        print('Now showing index', scan_index, 'of', str(len(frag_scan)), 'total found scans')
        plot_scan = frag_scan[scan_index][2]
        drop_index = np.argwhere(plot_scan.i <= noise_thr)
        plot_scan.i = np.delete(plot_scan.i, drop_index)
        plot_scan.mz = np.delete(plot_scan.mz, drop_index)
        mz = plot_scan.mz
        ints = plot_scan.i
        rt = plot_scan.scan_time[0]
        print('Precursor:', round(plot_scan.selected_precursors[0]['mz'],4), 'precursor intensity:', round(plot_scan.selected_precursors[0]['i'],1))
        print('Scan time:', round(plot_scan.scan_time[0], 2), 'minute')

        for i in range(len(mz)):
        	if i == 0:
        		list_string = str(round(mz[i],4)) + ':' + str(round(ints[i],1)) + ' '
        	else:
        		list_string += str(round(mz[i],4)) + ':' + str(round(ints[i],1)) + ' '

        url = 'https://mona.fiehnlab.ucdavis.edu/rest/similarity/search'
        values = {"spectrum" : list_string, "minSimilarity" : 500}

        data_str = json.dumps(values)
        data = data_str.encode()
        headers = {
          "content-type": "application/json",
        }
        req = urllib.request.Request(url, data, headers)
        response = urllib.request.urlopen(req)
        json_result = json.loads(response.read())

        mona_list = json_result[mona_index]['hit']['spectrum']
        spec1 = mona_list.split(' ')
        spec2 = [i.split(':') for i in spec1]
        mona_mz = [float(a) for a,b in spec2]
        mona_i = [float(b) for a,b in spec2]
        print(mona_mz)
        print(mona_i) #Integrate into the plot


        if interactive == True:
            plt.clf()
            fig = go.Figure([go.Bar(x=mz, y=ints, marker_color = 'red', width = 0.5,
                            hovertemplate =
                            'Int: %{y}'+
                            '<br>m/z: %{x}<br>')])
            fig.update_traces(marker_color='rgb(158,202,225)', marker_line_color='rgb(0,0,0)',
                      marker_line_width=0.5, opacity=1)
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
            plt.xlim(0,)

        if search==True:
            for i in range(len(mz)):
                if i == 0:
                    list_string = str(round(mz[i],4)) + ' ' + str(round(ints[i],1)) + '\r'
                else:
                    list_string += str(round(mz[i],4)) + ' ' + str(round(ints[i],1)) + '\r'
            pyperclip.copy(list_string)
            if source=='MoNA':
                webbrowser.open("https://mona.fiehnlab.ucdavis.edu/spectra/search")
            elif source=='metfrag':
                webbrowser.open("https://msbi.ipb-halle.de/MetFragBeta/")
        
    else:
        print('No MS2 spectrum found!')
    
    return