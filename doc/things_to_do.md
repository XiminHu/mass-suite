1. Either finish own peak picking workflow or using the peakonly workflow(after totally understand it)
	** seperate task

	User case: input-raw mzml, output-csv/xlsx with m/z(centralized) rt(centralized) rtrange? peakarea
	May require super computer/aws for the intensive computing


2. Alignment based on either a) output spreadsheet(using points) or b) using chromatogram
	* main task

	for a) local work with selected algorithms
	b) may need computing assisstance and more complicated
	--Need discussion

3. Clustering from the alignment

	Unsupervised ML preferred

4. Prediction based on Clustering?

5. optimize the visualization part--give a handy GUI for users

assign works for pairs rather than single people, follow the interst, prioritize menu
one main--branches
set up process, QA/QC trail dataset, blinded elements, avoid garbage in/out
train/test samples


* Before next Wed:
- collect material to get people onboard:
-- MS background knowledge (youtube is better)
	https://www.youtube.com/watch?v=ZN7euA1fS4Y
	https://www.youtube.com/watch?v=CyJaxr6FXAE
-- peak picking review(one module)
	Arsenty et al 2019 Deep Learning for the Precise Peak Detection in High-Resolution MS data
-- alignment review(one module)
	rob et al 2013 LCMS alignment review
-- why clustering
-- clustering review(main module)
Update later


Make up a slide for overview, part of works, update as week post?
- overview
- part of works to assign(module assignment)
- main goal development(brain storm for the solution)
- train/test samples


Peak picking idea::
Peak range refinement: instead of loop for 0 point as boundary, loop for 1/2 intensity and expand the result by 2 to get the range?


Alignment Idea::
Two ways to go:
1. Use every first feature in first sample as ref and align subsequent ones according to setup, if no-align ref feature then append to the align feature list
2. Create a dataset with all feature presented then group them for the alignment -- better