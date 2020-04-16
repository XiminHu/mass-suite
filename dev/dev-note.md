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

**important!! mz_locator--add relative ppm error, now only used absolute value
**alignment note: firstly read all features in seperate files, reorganize the features to order or rebin the file for different rt,mz bins to accelerate process time

check hyak
generate rt,int data with statistical numbers, dev ml model to do peak classification
write function for mz looping list generation