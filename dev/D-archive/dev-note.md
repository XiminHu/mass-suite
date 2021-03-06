* Change out current mz_locator method, instead of using the sum of all i for all mzs, only use the i for the mz with least mz error

TAsk:
1. cloud option
2. dig into peak picking algorithm
3. peak boundary/validation algorithm
4. mz list binning algorithm for the looping peak list generation --> kde? check Craig et al 2020


Task for me 4-11~4-18:
refine peak boundary: not use loop until thres, calculate every 3 previous int R value, if > thres(indicate linear relationship) then proceed, if not then break and record as boundary -dev on peak picking function
update output value: considering include S/N(algorithm needed), half intensity rt point instead of full boundary for alignment -- dev on sample 1d array
check scipy.signal functions, could be helpful to get statistical value out from 1d array

small task:
1. peak width, 1/2rt, S/N calculation
2. add filter onto boundary detection
3. check scipy.signal and peakutils code

**alignment note: firstly read all features in seperate files, reorganize the features to order or rebin the file for different rt,mz bins to accelerate process time

check hyak
generate rt,int data with statistical numbers, dev ml model to do peak classification
write function for mz looping list generation - done

potentially switch all list out and instead use np.array to boost loop efficiency

need a looping script for the peak manual labeling, potentially logic in: one data-one plot- input label- next loop
the modeling classification will be used after peak detection to exclude bad peaks or before detection to only have good peak into detection -- need discussion and efficiency test

4.23 - coding note
1. hyak learning
2. pack pymzml notebook into py script and test through mzml files execution
3. start wrap up for modeling

5.1 coding note
1. finish up the modeling and close up functions
2. finish the test
3. formula assignment & isotopic check?

5.7 coding note
Get the complete script from mzml file to aligned worksheet -- done
1. Model refinement
2. Revise the workflow to boost up iteration speed

5.12 coding note
1. try different models to see which one works best, KNN/RF/etc.
2. hyak python management--transport to stf group folder instead of personal login node

5.15 coding note
1. decide whether use all data for clustering or only dilution series for clustering
2. plot all data as barplot at rows -- x axis = sample, y axis = peak area, then treat them the same way as spectrum
3. heatmap to do the clustering?

data approach:
1. non-parameter testing
2. unsupervised random forest

3. use optic clustering to fill up adduct algorithm

task for finals:
1. test
2. package test/release
3. travis CI
4. clean repo & documentation
5. showcase notebook & ppt for standup -- GUI if possible -- Django, Flask, Dash, tkinter

1. prediction of dilution rate -- collinearity -- firstly plot the correlation plot, find significant features using P value
2. identification of multilple sources

Multiple source note:
!!for now only the ID is possible, if further accurate prediction/calculation needed, the matrix effect must besolved first for the back calculations