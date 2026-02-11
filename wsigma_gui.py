"""
ivan selyakov ivan_sel@yahoo.com
bertram boehrer
Helmholtz Centre for Environmental Research, Magdeburg, Germany
Department of Lake Research, Limnophysics and Lake Modelling group
2025
"""

import tkinter as tk

# Constants for pure water sound speed calculation
a0, a1, a2, a3, a4, a5 = 1402.38744, 5.03836171, -5.81172916e-2, 3.34638117e-4, -1.48259672e-6, 3.16585020e-9
a01, a11, a21, a31 = 1.49043589, 1.077850609e-2, -2.232794656e-4, 2.718246452e-6
a02, a12, a22, a32 = 4.31532833e-3, -2.938590293e-4, 6.822485943e-6, -6.674551162e-8
a03, a13, a23, a33 = -1.852993525e-5, 1.481844713e-6, -3.940994021e-8, 3.939902307e-10

# Dictionary of ions with sigma coefficients and el. conductivities at 25 degrees
species = {
    'Na⁺': {'sigma0_bar': 63.8599178743122, 'sigma1_bar': -0.47224026043208, 'C_25': 2.282769525},
    'Cl⁻': {'sigma0_bar': 0.0, 'sigma1_bar': 0.0, 'C_25': 3.469984222},
    'K⁺': {'sigma0_bar': 57.947786887924, 'sigma1_bar': -0.494660902938088, 'C_25': 3.350329594},
    'F⁻': {'sigma0_bar': 38.2247899221143, 'sigma1_bar': 0.489629494611651, 'C_25': 2.494361059},
    'Mg²⁺': {'sigma0_bar': 125.068115185474, 'sigma1_bar': -0.762273541741449, 'C_25': 3.616122064},
    'Mn²⁺': {'sigma0_bar': 96.6890194677852, 'sigma1_bar': -0.738806422673466, 'C_25': 3.50528025},
    'Ca²⁺': {'sigma0_bar': 97.691114421884, 'sigma1_bar': -0.76659002924916, 'C_25': 3.998491692},
    'Al³⁺': {'sigma0_bar': 167.619452015851, 'sigma1_bar': -1.1234207073516, 'C_25': 4.586898906},
    'NH₄⁺': {'sigma0_bar': 49.906789077711, 'sigma1_bar': -0.206433268618144, 'C_25': 3.325379261},
    'CO₃²⁻': {'sigma0_bar': 60.2343994848999, 'sigma1_bar': -0.18870088002358, 'C_25': 4.841037404},
    'SO₄²⁻': {'sigma0_bar': 40.853595284978, 'sigma1_bar': -0.288216480505722, 'C_25': 5.288218303},
    #'Fe³⁺': {'sigma0_bar': 195.90285209595, 'sigma1_bar': -1.0180318979529, 'C_25': 4.980765001},
    'Fe³⁺': {'sigma0_bar': 111.3419303, 'sigma1_bar': -1.018031898, 'C_25': 4.980765001},
    'HCO₃⁻': {'sigma0_bar': 30.2002408203893, 'sigma1_bar': -0.34984930394871, 'C_25': 2.088073306},
    'OH⁻': {'sigma0_bar': 34.5163000359989, 'sigma1_bar': -0.185697034911592, 'C_25': 9.003618345},
    'H⁺': {'sigma0_bar': 13.7724671112377, 'sigma1_bar': -0.253222454850672, 'C_25': 15.88675774},
    'NO₃⁻': {'sigma0_bar': -19.0535311182391, 'sigma1_bar': -0.162673352869687, 'C_25': 3.239228737},
}


def Wbel(t, P):
    
    """
    Calculates the sound speed in pure water based on temperature and pressure.
    
    Reference: Belogolskii, Sekoyan et al 1999
    Range of validity: 0-40 °C, 0.1 - 60 MPa

    Parameters:
    t (float): Temperature in degrees Celsius.
    P (float): Pressure in Pascals.

    Returns:
    float: Sound speed in m/s.
    """
    
    M1 = a01 + a11 * t + a21 * t ** 2 + a31 * t ** 3
    M2 = a02 + a12 * t + a22 * t ** 2 + a32 * t ** 3
    M3 = a03 + a13 * t + a23 * t ** 2 + a33 * t ** 3
    W0 = a0 + a1 * t + a2 * t ** 2 + a3 * t ** 3 + a4 * t ** 4 + a5 * t ** 5
    return W0 + M1 * (P - 0.101325) + M2 * (P - 0.101325) ** 2 + M3 * (P - 0.101325) ** 3

def calculate_sound_speed():
    
    """
    Calculates the sound speed excess based on entered molarities 
    and temperature. Added to Belogol'skii sound speed yeilds sound 
    speed in the sample.
    
    Calculates sound sigma coefficients specific for the sample.
    """
    
    T = float(entry_T.get())
    P = float(entry_P.get())

    molarity_values = {ion: float(entry.get()) for ion, entry in molarity_entries.items()}
    total_deltaW = 0
    total_C_25 = 0
    total_deltaW_5 = 0
    total_deltaW_25 = 0
    
    
    for ion, properties in species.items():
        molarity = molarity_values.get(ion, 0)
        sigma0_bar = properties['sigma0_bar']
        sigma1_bar = properties['sigma1_bar']
        C_25 = properties['C_25'] * molarity / 0.05
        total_C_25 += C_25

        if molarity > 0:
            total_deltaW += molarity * (sigma0_bar + sigma1_bar * (T - 25))
            total_deltaW_5 += molarity * (sigma0_bar + sigma1_bar * (5 - 25))
            total_deltaW_25 += molarity * (sigma0_bar + sigma1_bar * (25 - 25))
            
    sigma0 = total_deltaW_25 / total_C_25 if total_C_25 != 0 else 0
    sigma1 = ((total_deltaW_25 - total_deltaW_5) / (25 - 5)) * (1 / total_C_25) if total_C_25 != 0 else 0

    W_belogolskii = Wbel(T, P)
    W = total_deltaW + W_belogolskii

    # Display results in UI
    result_deltaW.config(text=f"Total sound speed excess: {total_deltaW:.2f} m/s")
    result_Wbel.config(text=f"Pure water sound speed (Belogol'skii 1999): {W_belogolskii:.2f} m/s")
    result_W.config(text=f"Calculated sound speed in sample: {W:.2f} m/s")
    result_sigma0.config(text=f"σ_0: {sigma0:.5f} [m/s] / [mS/cm]")
    result_sigma1.config(text=f"σ_1: {sigma1:.5f} [m/s] / [(K mS/cm)]")

    return {
        "temperature": T,
        "pressure": P,
        "total_deltaW": total_deltaW,
        "W_belogolskii": W_belogolskii,
        "W": W,
        "sigma0": sigma0,
        "sigma1": sigma1,
        "molarities": molarity_values
            }

# Function to set atmospheric pressure
def set_atmospheric_pressure():
    entry_P.delete(0, tk.END)  
    entry_P.insert(0, "0.101325")

# Function to save results to a text file
def save_to_file():
    filename = "results.txt"
    results = calculate_sound_speed()  # Call function to get updated values\
        
    try:
        with open(filename, 'w', encoding='utf-8') as file:
            file.write("Calculation Results\n")
            file.write("=" * 60 + "\n")

            # Write user-input molarities
            file.write("Chemical Composition of the Sample:\n")
            for ion, molarity in results["molarities"].items():
                file.write(f"{ion}: {molarity:.3f} mol/l\n")

            file.write("\n")
            file.write(f"Temperature: {results['temperature']} °C\n")
            file.write(f"Pressure: {results['pressure']} MPa\n")
            file.write(f"Total sound speed excess: {results['total_deltaW']:.2f} m/s\n")
            file.write(f"Pure water sound speed (Belogol'skii 1999): {results['W_belogolskii']:.2f} m/s\n")
            file.write(f"Calculated sound speed in sample: {results['W']:.2f} m/s\n")
            file.write(f"σ_0: {results['sigma0']:.5f} [m/s] / [mS/cm]\n")
            file.write(f"σ_1: {results['sigma1']:.5f} [m/s] / [K mS/cm]\n")

            # Add timestamp
            from datetime import datetime
            file.write("\nCalculation performed on: " + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

        print("Data successfully saved to file!")

    except Exception as e:
        print(f"Error saving file: {e}")

# Function to reset molarity fields
def clear_molarities():
    for entry in molarity_entries.values():
        entry.delete(0, tk.END)
        entry.insert(0, "0.")

# Create GUI window
root = tk.Tk()
root.title("Sound Speed Calculator")
root.geometry("800x750")

# Title label
title_label = tk.Label(root, text="", font=("Arial", 16))
title_label.pack(pady=10)

# Temperature and Pressure input fields
frame_TP = tk.Frame(root)
frame_TP.pack(pady=5)

tk.Label(frame_TP, text="Temperature (°C):").grid(row=0, column=0, padx=5, pady=5)
entry_T = tk.Entry(frame_TP)
entry_T.grid(row=0, column=1, padx=5, pady=5)

tk.Label(frame_TP, text="Pressure (MPa):").grid(row=1, column=0, padx=5, pady=5)
entry_P = tk.Entry(frame_TP)
entry_P.grid(row=1, column=1, padx=5, pady=5)

# Atmospheric pressure button
pressure_button = tk.Button(frame_TP, text="Set Atmospheric Pressure", command=set_atmospheric_pressure)
pressure_button.grid(row=1, column=2, padx=5, pady=5)

# Molarity input fields arranged in a grid layout
frame_molarity = tk.Frame(root)
frame_molarity.pack(pady=10, fill="both", expand=True)

molarity_label = tk.Label(frame_molarity, text="Molarities (mol/l):", font=("Arial", 12))
molarity_label.grid(row=0, column=0, padx=10, pady=5, sticky="w")

molarity_entries = {}  # Store Entry fields separately

row, col = 1, 0
for ion in species.keys():  # Only use ion names
    tk.Label(frame_molarity, text=ion, anchor="w").grid(row=row, column=col, padx=5, pady=5)

    molarity_entry = tk.Entry(frame_molarity, width=8)
    molarity_entry.insert(0, "0.")  # Default molarity value
    molarity_entry.grid(row=row, column=col+1, padx=5, pady=5)

    molarity_entries[ion] = molarity_entry  # Store input fields separately

    col += 2  # Move to the next column (label + input field)

    # Limit number of columns (e.g., max 5 ions per row)
    if col >= 10:
        col = 0
        row += 1

# Button to calculate sound speed
calculate_button = tk.Button(root, text="Calculate Sound Speed", command=calculate_sound_speed)
calculate_button.pack(pady=15)

# Result labels
result_deltaW = tk.Label(root, text="", font=("Arial", 12))
result_deltaW.pack(pady=5)

result_Wbel = tk.Label(root, text="", font=("Arial", 12))
result_Wbel.pack(pady=5)

result_W = tk.Label(root, text="", font=("Arial", 12))
result_W.pack(pady=5)

result_sigma0 = tk.Label(root, text="", font=("Arial", 12))
result_sigma0.pack(pady=5)

result_sigma1 = tk.Label(root, text="", font=("Arial", 12))
result_sigma1.pack(pady=5)

# Save to file button
save_button = tk.Button(root, text="Save Results to File", command=save_to_file)
save_button.pack(pady=15)

# Clear molarity fields button
clear_button = tk.Button(frame_molarity, text="Clear Fields", command=clear_molarities)
clear_button.grid(row=0, column=9, padx=5, pady=15, sticky="w")  # Position it right next to the label

# Start the GUI event loop
root.mainloop()
