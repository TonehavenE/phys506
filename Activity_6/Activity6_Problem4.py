"""
PHYS 506 Activity 6 - Problem 4

For Professor Butbaia
Authors: Luke Abanilla, Eben Quenneville, Augustus Vigorito
Date: 2024-10-20
"""

import numpy as np
import matplotlib.pyplot as plt
import DataAnalysisTools as dat
import math

def read_data(file_name: str) -> tuple[list, list]:
    values = []
    uncertainties = []
    with open(file_name, "r") as f:
        for line in f.readlines():
            value, uncertainty = line.split(",")
            values.append(float(value))
            uncertainties.append(float(uncertainty))
    return (values, uncertainties)

def main():
    resistors, resistor_uncertainties = read_data("resistor_data.txt")
    current, current_uncertainties = read_data("current_data.txt")

    # we measured the internal resistance of the ammmeter to be 11 Ohms when on the 200 mA range
    internal_resistance = 11 # Ohms

    # we offset our measured resistances by this to account for the internal resistance
    resistors = [r + internal_resistance for r in resistors]

    def v_uncertainty(I, R, uncertainty_I, uncertainty_R):
        """Helper function to find the uncertainty in voltage given current and resistance measurements."""
        return math.sqrt(R**2 * uncertainty_I**2 + I**2 * uncertainty_R**2)

    # Calculate the voltage and corresponding uncertainties
    voltages = []
    voltage_uncertainties = []
    for j in range(len(resistors)):
        r, del_r, i, del_i = resistors[j], resistor_uncertainties[j], current[j], current_uncertainties[j] # temp variables
        voltages.append(i * r)
        voltage_uncertainties.append(v_uncertainty(i, r, del_i, del_r))
    
    # Plot Data
    plt.errorbar(resistors, voltages, yerr=voltage_uncertainties)
    plt.xscale("log")
    plt.xlabel("resistance (ohms)")
    plt.ylabel("voltage (V)")

    m1, b1, m1Unc, b1Unc = dat.lineFitWt(resistors, voltages, voltage_uncertainties)
    line1 = np.vectorize(lambda x: m1 * x + b1)
    line_x_coordinates = np.linspace(min(resistors), max(resistors), 100)
    plt.plot(line_x_coordinates, line1(line_x_coordinates), label="Chi Squared Fit")
    print(f"Fit Quality: Q = {dat.fitQuality(resistors, voltages, voltage_uncertainties, m1, b1)}")
    print(f"lineFitWt results:\n* m = ({m1} ± {m1Unc}) V/Ohm\n* b = ({b1} ± {b1Unc}) V")

    plt.legend()
    plt.show()
    plt.savefig("Activity6_Problem4.png")
    
if __name__ == "__main__":
    main()
