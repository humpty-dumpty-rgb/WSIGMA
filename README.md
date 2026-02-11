## WSIGMA – evaluation of sound speed at given pressure, temperature and chemical composition

WSIGMA is a Python implementation for calculating sound speed in freshwater 
with dissolved ions based on experimentally retrieved coefficients.

The tool combines:

- Pure water sound speed (Belogolskii et al., 1999)
- Linear temperature dependence of ionic sound speed excess
- Additivity of ionic contributions (validated up to 0.05 mol/L total concentration)

The implementation follows the methodology described in:
[Journal reference here]

# Wbel(T, P) calculates the speed of sound in pure water based on (Belogolskii et al., 1999).

Parameters: T: Temperature in °C, P: Pressure in MPa
Returns: Sound speed in m/s


# List of available ions:
Na+, Cl-, K+, F-, Mg2+, Mn2+, Ca2+, Al3+, (NH4)+, (CO3)-2, (SO4)-2, Fe3+, (HCO3)-, OH-, H+, (NO3)-

Example molarity dictionary:

{
    'Na+': 0.01,
    'Cl-': 0.01,
    'Ca2+': 0.005
}




# Validity Range

- Temperature: 1–30 °C
- Pressure: limnic range (Belogolskii formulation valid 0–40 °C, 0.1–60 MPa)
- Solute concentration: up to 0.05 mol/L (~5 mS/cm)

# 
## Installation

