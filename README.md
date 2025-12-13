# CMT Project: Algal Proliferation Prediction Model

## Project description 

This project aims to develop a predictive model for algal growth in Lake Zug using historical measurements of atmospheric CO₂, surface and bottom water temperatures, and
chlorophyll-a concentrations. By simulating algal biomass at different depths over time, this model provides insight into how these environmental changes may influence algal dynamics. While simplified, this approach allows us to simulate algal growth and visualize trends over decades, providing a useful tool for both researchers and policymakers.

### Input files

Surface and bottom (195 m depth) temperature measurements from Lake Zug, along with global atmospheric CO₂ records dating back to 1969.

### Output files

Graphical predictions showing:
- Algal growth over time (up to 2050) at depths of 0 m, 1 m, 2 m, 5 m, 10 m, and 30 m.
- Algal biomass over time at the same depths.

### Report

The full report is available in the `docs/` directory as **report.pdf**.

## Running the program

### Dependencies
This project runs on the SIE Linux VDI, running on Matlab 2021b, as well as on C (v1.29.3).
### Build
The C source file `projet.c` is compiled using `gcc` into a standalone executable, which is placed in the `bin/` directory and executed between two MATLAB stages: a preprocessing step and a final plotting step.
### Execute
git clone https://github.com/Lorenzo-Mezz/CMT-Project.git && cd CMT-Project && chmod +x shellscript.sh && ./shellscript.sh

## Contributors

Nolan Chappatte & Lorenzo Mezzanotte Amat

## Acknowledgments

### Data sources
Input data were obtained from Datalakes Eawag (surface and bottom water temperatures), Alplakes (chlorophyll-a concentrations), and the NOAA Global Monitoring Laboratory (atmospheric CO₂).
### Code
Copilot and ChatGPT were used occasionally to assist with parts of the C implementation of the differential equation and to help identify and fix some coding errors. All code was written, checked, modified and validated by the authors.
