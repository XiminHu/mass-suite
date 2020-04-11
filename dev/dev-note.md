* Change out current mz_locator method, instead of using the sum of all i for all mzs, only use the i for the mz with least mz error

TAsk:
1. cloud option
2. dig into peak picking algorithm
3. peak boundary/validation algorithm
4. mz list binning algorithm for the looping peak list generation --> kde? check Craig et al 2020


Task for me 4-11~4-18:
refine mz_locator: not use sum, use the closest mz
refine peak boundary: not use loop until thres, calculate every 3 previous int R value, if > thres(indicate linear relationship) then proceed, if not then break and record as boundary