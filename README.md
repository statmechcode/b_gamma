# Statistical Analysis of Magnitude Jumps in Earthquake Catalogs
Code to calculate the value of b and gamma of the p(dm) distribution

This script performs a statistical analysis of magnitude differences between nearby earthquakes in space, as a function of increasing distance thresholds (dth).

For each distance threshold dth, the script:
	1.	Reads an earthquake catalog (input file) with fields: time, latitude, longitude, depth, magnitude. \\
	2.	For each event, searches for close neighboring events in space (within dth km).\\
	3.	Computes magnitude differences (dm) between event pairs.\\
	4.	Selects positive dm greater than a threshold dmth = 0.5.\\
	5.	Fits a statistical model to the dm distribution: Finds optimal parameters b and gamma by maximum likelihood estimation.\\
	6.	Estimates uncertainties on b and gamma using bootstrap resampling.\\
	7.	Outputs results for each dth into output.dat

Dependencies:
	•	Python ≥ 3.6\\
	•	numpy\\
	•	math (standard library)\\

 For any queries or suggestions, please contact:
Giuseppe Petrillo - giuseppe.petrillo@ntu.edu.sg
