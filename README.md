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

Parameters: T: Temperature in [°C], P: Pressure in [MPa]
Returns: Sound speed in [m/s]

# calculate_sound_speed(molarities, T, P, save_file=None) calculates total sound speed

Parameters: composition in [mol/L], T: Temperature in [°C], P: Pressure in [MPa], save_file (str, optional): Path to save results as a text file

Example chemical composition dictionary:

{
    'Na+': 0.01,
    'Cl-': 0.01,
    'Ca2+': 0.005
}

Returns: a dictionary with: total_deltaW m/s, Sound speed excess due to dissolved ions
W_belogolskii	m/s	Sound speed of pure water (Belogolskii et al., 1999)
W	m/s	Total sound speed in the sample
sigma0	m/s per (mS/cm)	σ₀ conductivity-normalized sound speed coefficient
sigma1	m/s per (K·mS/cm)	σ₁ temperature derivative of conductivity-normalized coefficient

dict: {
                'temperature': T,
                'pressure': P,
                'total_deltaW': total excess sound speed,
                'W_belogolskii': pure water sound speed,
                'W': total sound speed,
                'sigma0': sigma0,
                'sigma1': sigma1,
                'molarities': molarities
                    }


# List of available ions:
Na+, Cl-, K+, F-, Mg2+, Mn2+, Ca2+, Al3+, (NH4)+, (CO3)-2, (SO4)-2, Fe3+, (HCO3)-, OH-, H+, (NO3)-


# Validity Range

- Temperature: 1–30 °C
- Pressure: limnic range (Belogolskii formulation valid 0–40 °C, 0.1–60 MPa)
- Solute concentration: up to 0.05 mol/L (~5 mS/cm)

# 
## Installation

