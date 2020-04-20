import pandas as pd
import matplotlib.pyplot as plt
from ast import literal_eval

path = input('Please input the file path:')
sample_size = input('Please input sample size:')

df = pd.read_excel(path, sheet_name='intensity')
d_rt = pd.read_excel(path, sheet_name='rt')
#Not necessary
df.columns = [1, 'mz', 'int']
df.drop(columns = [1],inplace = True)
d_rt.columns = ['drop', 'rt']
d_rt.drop(columns = ['drop'],inplace = True)


d_iter = df.sample(int(sample_size))
d_iter['label'] = 0
count = 1

for i, row in d_iter.iterrows():
    plt.figure(figsize = (12,6))
    plt.plot(d_rt.values, literal_eval(row['int']))
    plt.xlabel('RT')
    plt.ylabel('intensity')
    plt.show()
    label = input('1 for good, 2 for unsure, 3 for bad:')
    d_iter.at[i,'label'] = label
    print('finished', count, 'out of ', d_iter.shape[0])
    count += 1
    plt.clf()

outputpath = '../example_data/peakdata/labelled_output/'
outputname = input('output name, please use .csv:')
output = outputpath+outputname
print(output)
d_iter.to_csv(output)



