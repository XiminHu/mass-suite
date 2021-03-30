import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
import re
from scipy.integrate import simps
import peakutils
import webbrowser
import pyperclip
import pyisopach
import plotly.offline as py
from ipywidgets import interactive, HBox, VBox
import pandas as pd
from .mssmain import peak_pick, ms_chromatogram_list


# TIC plot
def tic_plot(mzml_scans, interactive=True, f_width=10, f_height=6):
    '''
    Static tic plot function
    '''
    time = []
    TIC = []
    for scan in mzml_scans:
        time.append(scan.scan_time[0])
        TIC.append(scan.TIC)

    if interactive is True:
        fig = go.Figure([go.Scatter(x=time, y=TIC,
                        hovertemplate='Int: %{y}' + '<br>RT: %{x}minute<br>')])

        fig.update_layout(
            template='simple_white',
            width=f_width * 100,
            height=f_height * 100,
            xaxis=dict(title='Retention Time (min)',
                       rangeslider=dict(visible=True)),
            yaxis=dict(
                showexponent='all',
                exponentformat='e',
                title='Intensity',
            ))

        fig.show()

    elif interactive is False:
        plt.figure(figsize=(f_width, f_height))
        plt.plot(time, TIC)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        plt.xlabel('RT (min)')
        plt.ylabel('TIC')
        plt.title('TIC spectrum')
        plt.show()

    return


def ms_plot(mzml_scans, time, interactive=False, search=False,
            f_width=10, f_height=6, source='MoNA'):
    '''
    Interactive spectrum plot with nearest retention time from the given scan
    mzml_scans: mzfile
    time: selected time for the scan
    '''
    for scan in mzml_scans:
        if scan.scan_time[0] >= time:
            mz = scan.mz
            ints = scan.i
            rt = scan.scan_time[0]
            break

    if interactive is True:
        plt.clf()
        fig = go.Figure([go.Bar(x=mz, y=ints, marker_color='red', width=0.5,
                         hovertemplate='Int: %{y}' + '<br>m/z: %{x}<br>')])
        fig.update_traces(marker_color='rgb(158,202,225)',
                          marker_line_color='rgb(0,0,0)',
                          marker_line_width=0.5, opacity=1)
        fig.update_layout(
                title_text=str(round(rt, 3)) +
                ' MS spectrum @ ' + str(time) + ' min',
                template='simple_white',
                width=f_width * 100,
                height=f_height * 100,
                xaxis={'title': 'm/z ratio'},
                yaxis=dict(
                    showexponent='all',
                    exponentformat='e',
                    title='Intensity'))
        fig.show()

    elif interactive is False:
        plt.figure(figsize=(f_width, f_height))
        plt.bar(mz, ints, width=1.0)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        plt.xlabel('m/z')
        plt.ylabel('Intensity')
        plt.title('MS spectrum @' + str(time) + 'min')

    if search is True:
        for i in range(len(mz)):
            if i == 0:
                list_string = str(round(mz[i], 4)) + ' ' +\
                    str(round(ints[i], 1)) + '\r'
            else:
                list_string += str(round(mz[i], 4)) + ' ' +\
                    str(round(ints[i], 1)) + '\r'
        pyperclip.copy(list_string)
        if source == 'MoNA':
            webbrowser.open("https://mona.fiehnlab.ucdavis.edu/spectra/search")
        elif source == 'metfrag':
            webbrowser.open("https://msbi.ipb-halle.de/MetFragBeta/")

    return


def frag_plot(mzml_scans, precursor, error=20, scan_index=0,
              noise_thr=50, interactive=False, search=False, source='MoNA'):
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
        print('Now showing index', scan_index, 'of',
              str(len(frag_scan)), 'total found scans')
        plot_scan = frag_scan[scan_index][2]
        drop_index = np.argwhere(plot_scan.i <= noise_thr)
        plot_scan.i = np.delete(plot_scan.i, drop_index)
        plot_scan.mz = np.delete(plot_scan.mz, drop_index)
        mz = plot_scan.mz
        ints = plot_scan.i
        rt = plot_scan.scan_time[0]
        print('Precursor:', round(plot_scan.selected_precursors[0]['mz'], 4),
              'precursor intensity:',
              round(plot_scan.selected_precursors[0]['i'], 1))
        print('Scan time:', round(plot_scan.scan_time[0], 2), 'minute')

        if interactive is True:
            plt.clf()
            fig = go.Figure([go.Bar(x=mz, y=ints,
                                    marker_color='red', width=0.5,
                                    hovertemplate='Int: %{y}' +
                                    '<br>m/z: %{x}<br>')])
            fig.update_traces(marker_color='rgb(158,202,225)',
                              marker_line_color='rgb(0,0,0)',
                              marker_line_width=0.5, opacity=1)
            fig.update_layout(
                    title_text=str(round(rt, 3)) +
                    ' min MS1 spectrum, input ' + str(rt) + ' min',
                    template='simple_white',
                    width=1000,
                    height=600,
                    xaxis={'title': 'm/z ratio'},
                    yaxis=dict(
                        showexponent='all',
                        exponentformat='e',
                        title='Intensity'))
            fig.show()

        elif interactive is False:
            plt.figure(figsize=(10, 5))
            plt.bar(mz, ints, width=1.0)
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            plt.xlabel('m/z')
            plt.ylabel('Intensity')
            plt.title('MS1 spectrum')
            plt.xlim(0,)

        if search is True:
            for i in range(len(mz)):
                if i == 0:
                    list_string = str(round(mz[i], 4)), ' ' +\
                        str(round(ints[i], 1)) + '\r'
                else:
                    list_string += str(round(mz[i], 4)) + ' ' +\
                        str(round(ints[i], 1)) + '\r'
            pyperclip.copy(list_string)
            if source == 'MoNA':
                (webbrowser.open
                 ("https://mona.fiehnlab.ucdavis.edu/spectra/search"))
            elif source == 'metfrag':
                webbrowser.open("https://msbi.ipb-halle.de/MetFragBeta/")

    else:
        print('No MS2 spectrum found!')

    return


def mz_locator(array, mz, error):
    '''
    Find specific mzs from given mz and error range out from a given mz array
    input list: mz list
    mz: input_mz that want to be found
    error: error range is now changed to ppm level
    all_than_close False only select closest one, True will append all
    '''
    # ppm conversion
    error = error * 1e-6

    lower_mz = mz - error * mz
    higher_mz = mz + error * mz

    index = (array >= lower_mz) & (array <= higher_mz)

    return array[index], np.where(index)[0]


def formula_mass(input_formula, mode='pos'):
    '''
    sudo code:
    convert input string into a list with element:number structure
    convert all the element into upper case
    match the string list into a given list of element weight
    add adduct/delete H according to mode -- also have neutral mode
    '''
    # Define a list -- expand?
    elist = {'C': 12,
             'H': 1.00782,
             'N': 14.0031,
             'O': 15.9949,
             'S': 31.9721,
             'P': 30.973763,
             'e': 0.0005485799}

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


def ms_chromatogram(mzml_scans, input_value, error,
                    fillgap=False, mode='pos', interactive=True,
                    search=False, source='pubchem',
                    f_width=20, f_height=10):
    '''
    Interactive chromatogram for selected m/z
    search is now only available on chemspider and pubchem
    '''
    if type(input_value) == float:
        input_mz = input_value
    elif type(input_value) == int:
        input_mz = input_value
    elif type(input_value) == str:
        input_mz = formula_mass(input_value, mode)
    else:
        print('Cant recognize input type!')

    retention_time = []
    intensity = []
    for scan in mzml_scans:
        # print(i)
        retention_time.append(scan.scan_time[0])

        _, target_index = mz_locator(scan.mz, input_mz, error)
        if target_index == 'NA':
            intensity.append(0)
        else:
            intensity.append(max(scan.i[target_index]))

    def fill_gap(input_list, baseline=500):
        for i, intens in enumerate(input_list):
            if i > 1 and i < len(input_list)-3:
                if intens > baseline:
                    for index in np.arange(i+1, i+3):
                        if input_list[index] == 0:
                            input_list[index] = (input_list[index-1] +
                                                 input_list[index+1])/2
                        else:
                            continue
        return

    if fillgap is True:
        fill_gap(intensity)

    if interactive is True:
        fig = go.Figure([go.Scatter(x=retention_time, y=intensity,
                        hovertemplate='Int: %{y}' + '<br>RT: %{x}minute<br>')])

        fig.update_layout(
            title_text=str(round(input_mz, 2)) +
            ' chromatogram, error ' + str(error),
            template='simple_white',
            width=f_width * 100,
            height=f_height * 100,
            xaxis={'title': 'Retention Time (min)'},
            yaxis=dict(
                showexponent='all',
                exponentformat='e',
                title='Intensity'))

        fig.show()
    elif interactive is False:
        plt.figure(figsize=(f_width, f_height))
        plt.plot(retention_time, intensity)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        plt.xlabel('Retention Time(min)')
        plt.ylabel('Intensity')
        plt.title('Chromatogram for m/z ' + str(round(input_mz, 3)))
        plt.xlim(0, retention_time[-1])
        plt.ylim(0, )
        plt.show()

    if search is False:
        pass
    if search is True:
        if type(input_value) == str:
            if source == 'chemspider':
                webbrowser.open("http://www.chemspider.com/Search.aspx?q="
                                + input_value)
            elif source == 'pubchem':
                webbrowser.open("https://pubchem.ncbi.nlm.nih.gov/#query="
                                + input_value)
        else:
            print('Please enter formula for search!')

    return


def integration_plot(mzml_scans, input_mz, error,
                     f_width=20, f_height=10):

    result_dict = peak_pick(mzml_scans, input_mz, error,
                            min_scan=5, peak_area_thres=0)
    
    rt = [i.scan_time[0] for i in mzml_scans]
    ints = ms_chromatogram_list(mzml_scans, input_mz, error)

    plt.figure(figsize=(f_width, f_height))
    plt.plot(rt, ints)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    plt.xlabel('Retention Time(min)')
    plt.ylabel('Intensity')
    plt.title('Integration result')
    plt.xlim(0, rt[-1])
    plt.ylim(0, )

    for index in result_dict:
        print(('Peak retention time: {:0.2f} minute, Peak area: {: 0.1f}'
               .format(rt[index], result_dict[index][2])))
        plt.fill_between(rt[result_dict[index][0]: result_dict[index][1]],
                         ints[result_dict[index][0]: result_dict[index][1]])

    return


def iso_plot(mzml_scan, input_mz, error, formula):
    '''
    Interactive spectrum plot with nearest retention time from the given scan
    mzml_scans: mzfile
    time: selected time for the scan
    '''
    def closest(lst, K):
        idx = np.abs(np.asarray(lst) - K).argmin()
        return idx

    select_intensity = ms_chromatogram_list(mzml_scan, input_mz, error)
    scan = mzml_scan[np.argmax(select_intensity)]

    mz = scan.mz
    ints = scan.i

    precursor_idx = closest(mz, input_mz)
    precursor_mz = mz[precursor_idx]
    precursor_ints = ints[precursor_idx]

    rel_abundance = [i / precursor_ints * 100 for i in ints]

    # Predicted isotope pattern
    mol = pyisopach.Molecule(formula)
    isotope_i = [-i for i in mol.isotopic_distribution()[1]]
    iso_mz = mol.isotopic_distribution()[0]

    wd = 0.05
    _, ax = plt.subplots(figsize=(12, 9))
    ax.bar(mz, rel_abundance, width=wd, label='scan spectrum')
    ax.bar(iso_mz, isotope_i, width=wd, label='predicted isotope pattern')
    ax.axhline(y=0, color='k')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    ticks = ax.get_yticks()
    ax.set_yticklabels([int(abs(tick)) for tick in ticks])
    plt.xlabel('m/z')
    plt.ylabel('Relative Intensity %')
    plt.title('Isotope pattern comparison')
    plt.legend()
    plt.xlim(precursor_mz - 5, precursor_mz + 10)

    return


def manual_integration(mzml_scans, input_mz, error, start, end):
    '''
    Area integration for selected mz and time
    '''
    rt_lst = []
    intensity = []
    for scan in mzml_scans:
        # print(i)
        rt_lst.append(scan.scan_time[0])

        _, target_index = mz_locator(scan.mz, input_mz, error)
        if target_index == 'NA':
            intensity.append(0)
        else:
            intensity.append(sum(scan.i[target_index]))

    def closest(lst, K):
        return lst[min(range(len(lst)), key=lambda i: abs(lst[i]-K))]

    s_index = rt_lst.index(closest(rt_lst, start))
    e_index = rt_lst.index(closest(rt_lst, end))

    integrated = simps(y=intensity[s_index:e_index], even='avg')

    return integrated


def overview_scatter(data):
    # Currently only use on MSS dataset
    # Original reference:
    # https://plotly.com/python/v3/selection-events/
    py.init_notebook_mode()

    df = data
    df['max area'] = df.iloc[:, 3:].max(1)

    f = go.FigureWidget([go.Scatter(x=df['Average Rt(min)'],
                                    y=df['Average Mz'],
                                    mode='markers')])
    f.layout.xaxis.title = 'Retention Time (min)'
    f.layout.yaxis.title = 'm/z Ratio'
    scatter = f.data[0]
    scatter.marker.opacity = 0.5

    data_col = ['Average Rt(min)', 'Average Mz', 'S/N average', 'max area']
    t = go.FigureWidget([go.Table(
        header=dict(values=data_col,
                    fill=dict(color='#C2D4FF'),
                    align=['left'] * 5),
        cells=dict(values=[df[col] for col in data_col],
                   fill=dict(color='#F5F8FF'),
                   align=['left'] * 5))])

    def selection_fn(trace, points, selector):
        t.data[0].cells.values =\
            [df.loc[points.point_inds][col] for col in data_col]

    scatter.on_selection(selection_fn)

    return VBox((HBox(), f, t))
