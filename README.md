# WSIGMA – evaluation of sound speed at given pressure, temperature and chemical composition

WSIGMA is a Python implementation for calculating sound speed in freshwater 
with dissolved ions based on experimentally retrieved coefficients.

The tool combines:

- Pure water sound speed (Belogolskii et al., 1999)
- Linear temperature dependence of ionic sound speed excess
- Additivity of ionic contributions (validated up to 0.05 mol/L total concentration)

The implementation follows the methodology described in:
[Journal reference here]


## List of available ions:
Na+, Cl-, K+, F-, Mg2+, Mn2+, Ca2+, Al3+, NH4+, CO3^-2, SO4^-2, Fe3+, HCO3-, OH-, H+, NO3-


## Validity Range

- Temperature: 1–30 °C
- Pressure: limnic range (Belogolskii formulation valid 0–40 °C, 0.1–60 MPa)
- Solute concentration: up to 0.05 mol/L (~5 mS/cm)


## Functions

### Wbel(T, P) calculates the speed of sound in pure water based on (Belogolskii et al., 1999).

Parameters:
- T: Temperature in [°C];
- P: Pressure in [MPa]
  
Returns:
- Sound speed in [m/s]

### calculate_sound_speed(molarities, T, P, save_file=None) calculates total sound speed

Parameters:
- Composition in [mol/L];
- T: temperature in [°C];
- P: pressure in [MPa];
- save_file (str, optional): path to save results as a text file

Example chemical composition dictionary:

sample =
    {
        'Na+': 0.01,
        'Cl-': 0.01,
        'Ca2+': 0.005
    }

Returns: 
a dictionary containing:
- temperature: temperature in [°C];
- pressure: pressure in [MPa];
- total_deltaW: Sound speed excess due to dissolved ions in [m/s];
- W_belogolskii: Sound speed of pure water (Belogolskii et al., 1999) in [m/s];
- W: Total sound speed in [m/s];
- sigma0: specific contribution of the dissolved electrolytes at 25°C in 	[(m/s) / (mS/cm)];
- sigma1:	temperature dependence coefficient in [(m/s) / (K·mS/cm)];

results = {
        'temperature': T,
        'pressure': P,
        'total_deltaW': total_deltaW,
        'W_belogolskii': W_bel,
        'W': W_total,
        'sigma0': sigma0,
        'sigma1': sigma1,
        'composition': composition
    }

## Example usage

sample = {'Na+': 0.01, 'Cl-': 0.01, 'Ca2+': 0.005}
result = calculate_sound_speed(sample, T=20, P=0.101325)
print(result)

or:

sample = {'Na+': 0.01, 'Cl-': 0.01, 'Ca2+': 0.005}
results = calculate_sound_speed(sample, T=20, P=0.101325)
print("Total sound speed:", results['W'], "m/s")

## Installation

pip install git+https://github.com/yourusername/wsigma.git
from wsigma import calculate_sound_speed
