[![DOI](https://zenodo.org/badge/998813138.svg)](https://doi.org/10.5281/zenodo.15653633)

# Statistical Analysis of Magnitude Jumps in Earthquake Catalogs
Code to calculate the value of b and gamma of the p(dm) distribution with incompleteness

This script performs a statistical analysis of magnitude differences between nearby earthquakes in space.

## Overview

This repository contains a Python implementation of a maximum likelihood estimation method. It estimates parameters `b` and `gamma` of a two-branch model for magnitude differences (dm) between temporally-close earthquakes. The approach is based on:

- Selecting earthquake pairs separated by a small spatial distance (dth) and a magnitude threshold (dmth).
- Computing dm between each valid pair.
- Fitting the log-likelihood.
- Estimating parameter uncertainties via bootstrap.

## Files

- `like_bgamma.py`: main Python script (this file)
- `input_catalog.dat`: input catalog file with 5 columns:  
  `time`, `magnitude`, `latitude`, `longitude`, `depth`
- `output.dat`: output file containing the fitted parameters and uncertainties

## Usage

Run the main script with:
python main.py

Ensure your input_catalog.dat is in the correct format:

time    magnitude    latitude    longitude    depth

12345678  4.5          34.05       -118.25      12.0

...

Each line must have at least 5 floating point values.

## Usage

The output is written to output.dat and contains the following columns:

dth  b  sb  a  sa

Where:

	•	dth: distance threshold (fixed to 20.0 km)
 
	•	b, a: maximum likelihood estimates
 
	•	sb, sa: bootstrap standard deviations


## Notes

A small Gaussian noise is added to each magnitude to break ties:

q.append(qi + 0.1 * (np.random.rand() - 0.5))

## Requirements

Please install Python packages with:

- pip install numpy scipy

## Contact

For questions, suggestions, or contributions, please contact:

Giuseppe Petrillo

giuseppe51289@gmail.com


