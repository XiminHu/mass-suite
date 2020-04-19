import pandas as pd
import matplotlib.pyplot as plt
from ast import literal_eval

path = input('Please input the file path:')
sample_size = input('Please input sample size:')

df = pd.read_csv(path)
#Not necessary
df.columns = [1, 'mz', 'int','rt']
df.drop(columns = [1],inplace = True)


d_iter = df.sample(int(sample_size))
d_iter['label'] = 0
count = 1

for i, row in d_iter.iterrows():
    plt.figure(figsize = (12,6))
    plt.plot(literal_eval(row['int']))
    plt.show()
    label = input('1 for bad, 2 for unsure, 3 for good:')
    d_iter.at[i,'label'] = label
    print('finished', count, 'out of ', d_iter.shape[0])
    count += 1
    plt.clf()

outputpath = '../example_data/peakdata/labelled_output/'
outputname = input('output name:')
output = outputpath+outputname
print(output)
d_iter.to_csv(output)



