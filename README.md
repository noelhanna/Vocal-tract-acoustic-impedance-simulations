# Vocal-tract-acoustic-impedance-simulations
1D plane wave propogation simulations of the impedance of concatenated cylinders

## Vocal_Tract_General_Model_2017_SVT
#### Input
This function takes three input structures
- VT is the supraglottal vocal tract
- glottis
- SG is the subglottal vocal tract
Each of these structures has radius, length, and alpha_multiplier which can be single values or arrays that define the airway geometry as the sum of cylinders
alpha_multiplier = 1 for standard visco-thermal losses. 5 is typical for human airways (Hanna et al., 2016)
#### Output
Structure result
- freq, frequency vector that covers the range of calculations
- tube, contains the input structures and all calculated impedances and transfer functions
- PARAMS, contains temperature, humidity, density of air, speed of sound etc.

## Vocal_Tract_Calculations 
This is a script that runs the above function and creates figures showing the airway profile as a function of radius, the impedance, and the transfer functions for rigid and non-rigid tracts

## Known Issues
- Vocal_Tract_Calculations contains hardcoded example values, these need to be commented out when loading a different airway geometry
- Vocal_Tract_Calculations does not clear the workspace
- Temperature, humidity and speed of sound are not transparent (may not be calculated properly)
- Flex calculations may not account for subglottal tract appropriately, currently just add inertance and compliance in parallel with Z_load acting on the glottis
- Airway radius profile is a bit ugly
