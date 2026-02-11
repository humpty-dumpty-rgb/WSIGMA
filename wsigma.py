"""
WSIGMA - Sound Speed Calculator for Freshwater with Solutes
Author: Ivan Selyakov (ivan_sel@yahoo.com), Bertram Boehrer
Helmholtz Centre for Environmental Research, Magdeburg, Germany
Department of Lake Research, Limnophysics and Lake Modelling group
2025

This snippet calculates sound speed in water for a given temperature, pressure,
and ionic composition.
"""

# Constants for pure water sound speed calculation
a0, a1, a2, a3, a4, a5 = 1402.38744, 5.03836171, -5.81172916e-2, 3.34638117e-4, -1.48259672e-6, 3.16585020e-9
a01, a11, a21, a31 = 1.49043589, 1.077850609e-2, -2.232794656e-4, 2.718246452e-6
a02, a12, a22, a32 = 4.31532833e-3, -2.938590293e-4, 6.822485943e-6, -6.674551162e-8
a03, a13, a23, a33 = -1.852993525e-5, 1.481844713e-6, -3.940994021e-8, 3.939902307e-10

# Dictionary of ions with sigma coefficients and electrical conductivity at 25°C
species_properties = {
    'Na+': {'sigma0_bar': 63.8599178743122, 'sigma1_bar': -0.47224026043208, 'C_25': 2.282769525},
    'Cl-': {'sigma0_bar': 0.0, 'sigma1_bar': 0.0, 'C_25': 3.469984222},
    'K+': {'sigma0_bar': 57.947786887924, 'sigma1_bar': -0.494660902938088, 'C_25': 3.350329594},
    'F-': {'sigma0_bar': 38.2247899221143, 'sigma1_bar': 0.489629494611651, 'C_25': 2.494361059},
    'Mg2+': {'sigma0_bar': 125.068115185474, 'sigma1_bar': -0.762273541741449, 'C_25': 3.616122064},
    'Mn2+': {'sigma0_bar': 96.6890194677852, 'sigma1_bar': -0.738806422673466, 'C_25': 3.50528025},
    'Ca2+': {'sigma0_bar': 97.691114421884, 'sigma1_bar': -0.76659002924916, 'C_25': 3.998491692},
    'Al3+': {'sigma0_bar': 167.619452015851, 'sigma1_bar': -1.1234207073516, 'C_25': 4.586898906},
    'NH4+': {'sigma0_bar': 49.906789077711, 'sigma1_bar': -0.206433268618144, 'C_25': 3.325379261},
    'CO3^2-': {'sigma0_bar': 60.2343994848999, 'sigma1_bar': -0.18870088002358, 'C_25': 4.841037404},
    'SO4^2-': {'sigma0_bar': 40.853595284978, 'sigma1_bar': -0.288216480505722, 'C_25': 5.288218303},
    'Fe3+': {'sigma0_bar': 111.3419303, 'sigma1_bar': -1.018031898, 'C_25': 4.980765001},
    'HCO3-': {'sigma0_bar': 30.2002408203893, 'sigma1_bar': -0.34984930394871, 'C_25': 2.088073306},
    'OH-': {'sigma0_bar': 34.5163000359989, 'sigma1_bar': -0.185697034911592, 'C_25': 9.003618345},
    'H+': {'sigma0_bar': 13.7724671112377, 'sigma1_bar': -0.253222454850672, 'C_25': 15.88675774},
    'NO3-': {'sigma0_bar': -19.0535311182391, 'sigma1_bar': -0.162673352869687, 'C_25': 3.239228737},
}

def Wbel(T, P):
    """
    Calculate the speed of sound in pure water.

    Parameters:
        T (float): Temperature in °C
        P (float): Pressure in MPa

    Returns:
        float: Sound speed in m/s
    """
    M1 = a01 + a11 * T + a21 * T**2 + a31 * T**3
    M2 = a02 + a12 * T + a22 * T**2 + a32 * T**3
    M3 = a03 + a13 * T + a23 * T**2 + a33 * T**3
    W0 = a0 + a1*T + a2*T**2 + a3*T**3 + a4*T**4 + a5*T**5
    return W0 + M1*(P - 0.101325) + M2*(P - 0.101325)**2 + M3*(P - 0.101325)**3

def calculate_sound_speed(composition, T, P, save_file=None):
    """
    Calculate total sound speed in a water sample with given ionic composition.

    Parameters:
        composition (dict): ion -> molarity in mol/L (e.g., {'Na+': 0.01, 'Cl-':0.01})
        T (float): Temperature in °C
        P (float): Pressure in MPa
        save_file (str, optional): Path to save results as a text file

    Returns:
        dict: {
            'temperature': T, (°C)
            'pressure': P, (MPa)
            'total_deltaW': total excess sound speed, (m/s)
            'W_belogolskii': pure water sound speed, (m/s)
            'W': total sound speed, (m/s)
            'sigma0': sigma0, (m/s per mS/cm)
            'sigma1': sigma1, (m/s per (K mS/cm))
            'composition': composition dict (mol/L)
        }

    """
    total_deltaW = 0
    total_C25 = 0
    total_deltaW_5 = 0
    total_deltaW_25 = 0

    for ion, properties in species_properties.items():
        m = composition.get(ion, 0)
        C_25 = properties['C_25'] * m / 0.05
        total_C25 += C_25
        if m > 0:
            total_deltaW += m * (properties['sigma0_bar'] + properties['sigma1_bar']*(T-25))
            total_deltaW_5 += m * (properties['sigma0_bar'] + properties['sigma1_bar']*(5-25))
            total_deltaW_25 += m * (properties['sigma0_bar'] + properties['sigma1_bar']*(25-25))

    sigma0 = total_deltaW_25 / total_C25 if total_C25 else 0
    sigma1 = ((total_deltaW_25 - total_deltaW_5)/(25-5))*(1/total_C25) if total_C25 else 0

    W_bel = Wbel(T, P)
    W_total = W_bel + total_deltaW

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

    if save_file:
        with open(save_file, 'w', encoding='utf-8') as f:
            f.write("Sound Speed Calculation Results\n")
            f.write("="*60+"\n")
            for ion, m in composition.items():
                f.write(f"{ion}: {m:.4f} mol/L\n")
            f.write(f"\nTemperature: {T} °C\n")
            f.write(f"Pressure: {P} MPa\n")
            f.write(f"Total sound speed excess: {total_deltaW:.4f} m/s\n")
            f.write(f"Pure water sound speed: {W_bel:.4f} m/s\n")
            f.write(f"Total sound speed: {W_total:.4f} m/s\n")
            f.write(f"σ0: {sigma0:.5f} [m/s]/[mS/cm]\n")
            f.write(f"σ1: {sigma1:.5f} [m/s]/[(K mS/cm)]\n")

    return results

# Example usage
if __name__ == "__main__":
    sample = {'Na+': 0.01, 'Cl-': 0.01, 'Ca2+': 0.005}
    result = calculate_sound_speed(sample, T=20, P=0.101325)
    print(result)
