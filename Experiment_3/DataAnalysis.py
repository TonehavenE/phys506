from math import log, sqrt
import DataAnalysisTools as dat
import matplotlib.pyplot as plt 
import numpy as np 
import re

def gus_algorithm(voltages, voltage_uncs, currents, current_uncs, dark_v, dark_v_uncs, dark_c, dark_c_uncs):
    # The idea of the following algorithm is to calculate a curve for current vs voltage that is
    # shifted by the piecewise dark characteristic curve, and then find 0s of that corrected curve
    # these should be the intersection values between the two data sets

    def linear(x, x1, y1, x2, y2):
        # linear equation f(x) = m(x - x1) + y1
        m = (y2 - y1)/(x2 - x1)
        return m * (x - x1) + y1
    
    def linear_for_x(y, x1, y1, x2, y2):
        # inverse of the above, x = (y-y1)/m + x1
        m = (y2 - y1)/(x2 - x1)
        return (y - y1)/m + x1

    corrected = []
    for i in range(len(voltages)): # for each voltage
        for beta in range(len(dark_v)):
            if dark_v[beta] > voltages[i]: # iterate till we find a dark characteristic voltage data point bigger than the voltage
                # we interpolate between the first point we find that is bigger than our voltage and the previous point (assumed smaller)
                # corrected is the current at this voltage minus the interpolated value of the characteristic curve at this voltage
                corrected.append(currents[i] - linear(voltages[i], dark_v[beta], dark_c[beta], dark_v[beta-1], dark_c[beta-1]))
                break

    corrected.sort()

    # now we need to find where the shifted current values intersect y = 0
    vl1, vl1_unc, vl2, vl2_unc, cl1, cl1_unc, cl2, cl2_unc, vd1, vd1_unc, vd2, vd2_unc, cd1, cd1_unc, cd2, cd2_unc = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 # get rid of error messages
    for i in range(len(corrected)):
        if corrected[i] >= 0:
            # select data points
            vl1 = voltages[i - 1]
            vl1_unc = voltage_uncs[i - 1]
            vl2 = voltages[i]
            vl2_unc = voltage_uncs[i]
            cl1 = currents[i - 1]
            cl1_unc = current_uncs[i - 1]
            cl2 = currents[i]
            cl2 = current_uncs[i]

            vd1 = dark_v[i - 1]
            vd1_unc = dark_v_uncs[i - 1]
            vd2 = dark_v[i]
            vd2_unc = dark_v_uncs[i]
            cd1 = dark_c[i - 1]
            cd1_unc = dark_c_uncs[i - 1]
            cd2 = dark_c[i]
            cd2_unc = dark_c_uncs[i]
            # v_s = linear_for_x(0, voltages[i], corrected[i], voltages[i-1], corrected[i-1])
            break
    
    # calculate uncertainty
    md = (cd1 - cd2)/(vd1 - vd2)
    ml = (cl1 - cl2)/(vl1 - vl2)
    alpha = md * vd1 - ml * vl1 + cl1 - cd1 # just a shorthand
    beta = (md - ml)**2 # just a shorthand
    v_s = alpha/(md - ml)

    # Partial Derivatives
    # for each of the following, x_by_y means the partial derivative of x with respect to y (∂x/∂y)
    ## Dark Slopes
    md_by_vd1 = - (cd1-cd2)/(vd1 - vd2)**2
    md_by_vd2 = (cd1 - cd2)/(vd1 - vd2)**2

    md_by_cd1 = 1/(vd1 - vd2)
    md_by_cd2 = -1/(vd1 - vd2)

    ## Light slopes
    ml_by_vl1 = - (cl1 - cl2)/(vl1 - vl2)**2
    ml_by_vl2 = (cl1 - cl2)/(vl1 - vl2)**2

    ml_by_cl1 = 1/(vl1 - vl2)
    ml_by_cl2 = -1/(vl1 - vl2)
    
    ## Voltage by each variable
    vs_by_vd1 = ( (md_by_vd1 * vd1 + md) * (md - ml) - alpha * md_by_vd1)/beta
    vs_by_vd2 = ( (md_by_vd2 * vd1) * (md - ml) - alpha * md_by_vd2 )/beta
    vs_by_vl1 = ( - (ml_by_vl1 * vl1 + ml ) * (md - ml) + alpha * ml_by_vl1 ) / beta
    vs_by_vl2 = ( - (ml_by_vl2 * vl1) * (md - ml) + alpha * ml_by_vl2 )/beta

    vs_by_cd1 = ( (md_by_cd1 * vd1 - 1) * (md - ml) - alpha * md_by_cd1)/beta
    vs_by_cd2 = ( (md_by_cd2 * vd1) * (md - ml) - alpha * md_by_cd2)/beta
    vs_by_cl1 = ( (- ml_by_cl1 * vl1 + 1) * (md - ml) + alpha * ml_by_cl1)/beta
    vs_by_cl2 = ( (- ml_by_cl2 * vl1) * (md - ml) + alpha * ml_by_cl2 )/beta


    # total uncertainty
    v_s_unc = sqrt(
            (
                (vs_by_vd1 * vd1_unc)**2 +
                (vs_by_vd2 * vd2_unc)**2 +
                (vs_by_cd1 * cd1_unc)**2 +
                (vs_by_cd2 * cd2_unc)**2 +
                (vs_by_vl1 * vl1_unc)**2 +
                (vs_by_vl2 * vl2_unc)**2 +
                (vs_by_cl1 * cl1_unc)**2 +
                (vs_by_cl2 * cl2_unc)**2 
            ))

    return (v_s, v_s_unc)

def digital_data_plot():
    print("\n--- DIGITAL DATA ---\n")
    # files = ["Analog3650", "Analog4047", "Analog4358", "Analog5461", "Analog5770"]
    files = ["RedLED_611_nm", "YellowLED_588_nm", "GreenLED_525_nm", "CyanLED_505_nm", "BlueLED_472_nm"]
    stopping_voltages = []
    # files = ["RedLED_611_nm", "YellowLED_588_nm"]

    fig, axis = plt.subplots(1, len(files), sharey=True)
    fig.suptitle("Digital Experimental Data")
    fig.supxlabel("Accelerating Voltage (V)")
    fig.supylabel("Photocurrent (A)")

    for index, file in enumerate(files):
        wavelength = float(file.split("_")[1])
        wavelength_unc = 0.05 * wavelength
        voltages = np.loadtxt(f"data/{file}.csv", delimiter=',', usecols=0, skiprows=1)
        voltage_unc = np.loadtxt(f"data/{file}.csv", delimiter=',', usecols=1, skiprows=1)
        currents = np.loadtxt(f"data/{file}.csv", delimiter=',', usecols=2, skiprows=1)
        currents *= 10**(-9)
        current_unc = np.loadtxt(f"data/{file}.csv", delimiter=',', usecols=3, skiprows=1)
        current_unc *= 10**(-9)

        # manufacture dark characteristic data to match as 0 - this works for the digital measurement
        # we found experimentally that the dark current is 0 everywhere for this instrument
        dark_v = voltages
        dark_v_unc = [0 for _ in voltage_unc]
        dark_c = [0 for _ in currents]
        dark_c_unc = [0 for _ in current_unc]

        axis[index].errorbar(voltages, currents, xerr=voltage_unc, yerr=current_unc)
        axis[index].errorbar(dark_v, dark_c, xerr=dark_v_unc, yerr=dark_c_unc)
        axis[index].set_title(f"{file}")

        stopping_voltage, stopping_voltage_unc = gus_algorithm(voltages, voltage_unc, currents, current_unc, dark_v, dark_v_unc, dark_c, dark_c_unc)
        # print(f"Gus's Algorithm: {stopping_voltage} ± {stopping_voltage_unc}")
        axis[index].errorbar([stopping_voltage, stopping_voltage], [min(currents), max(currents)], xerr=[stopping_voltage_unc, stopping_voltage_unc]) 
        stopping_voltages.append((wavelength, wavelength_unc, stopping_voltage, stopping_voltage_unc))

        print(f"({wavelength:.4f} ± {wavelength_unc:.4f}) nm: ({stopping_voltage:.4f} ± {stopping_voltage_unc:.4f}) v")
        # print("")
    
    return stopping_voltages

def analog_data_plot():
    print("\n--- ANALOG DATA ---\n")

    files = ["Analog3650A", "Analog4047A", "Analog4358A", "Analog5461A", "Analog5770A"][::-1]

    data = []
    

    dark_v = np.loadtxt(f"data/AnalogDark.csv", delimiter=',', usecols=0, skiprows=1)
    dark_v_unc = np.abs(np.loadtxt(f"data/AnalogDark.csv", delimiter=',', usecols=0, skiprows=1))
    dark_c = np.loadtxt(f"data/AnalogDark.csv", delimiter=',', usecols=0, skiprows=1)
    dark_c *= 10**(-14)
    dark_c_unc = np.abs(np.loadtxt(f"data/AnalogDark.csv", delimiter=',', usecols=0, skiprows=1))
    dark_c_unc *= 10**(-14)

    fig, axis = plt.subplots(1, len(files), sharey=True)
    fig.suptitle("Analog Experimental Data")
    fig.supxlabel("Accelerating Voltage (V)")
    fig.supylabel("Photocurrent (A)")

    for index, file in enumerate(files):
        wavelength = float(re.sub('[A-Za-z]', '', file)) / 10
        wavelength_unc = 0.05 * wavelength
        voltages = np.loadtxt(f"data/{file}.csv", delimiter=',', usecols=0, skiprows=1)
        voltage_unc = np.loadtxt(f"data/{file}.csv", delimiter=',', usecols=1, skiprows=1)
        currents = np.loadtxt(f"data/{file}.csv", delimiter=',', usecols=2, skiprows=1)
        currents *= 10**(-14)
        current_unc = np.loadtxt(f"data/{file}.csv", delimiter=',', usecols=3, skiprows=1)
        current_unc *= 10**(-14)


        axis[index].errorbar(voltages, currents, xerr=voltage_unc, yerr=current_unc)
        axis[index].errorbar(dark_v, dark_c, xerr=dark_v_unc, yerr=dark_c_unc)
        axis[index].set_title(f"{file}")

        stopping_voltage, stopping_voltage_unc = gus_algorithm(voltages, voltage_unc, currents, current_unc, dark_v, dark_v_unc, dark_c, dark_c_unc)
        axis[index].errorbar([stopping_voltage, stopping_voltage], [min(currents), max(currents)], xerr=[stopping_voltage_unc, stopping_voltage_unc]) 
        data.append((wavelength, wavelength_unc, stopping_voltage, stopping_voltage_unc))
        print(f"({wavelength:.4f} ± {wavelength_unc:.4f}) nm: ({stopping_voltage:.4f} ± {stopping_voltage_unc:.4f}) v")
        # print("")
    

    return data

def stopping_potential_vs_frequency(data: list[tuple[float, float, float, float]]):
    """
    Takes in a tuple of data points, (λ, δλ, v, δv)
    Converts λ to ν using c = λν and then plots voltage vs frequency
    Also plots a line of best fit
    """
    c = 2.998 * 10**8
    voltages = []
    wavelengths = []
    wavelength_uncs = []
    voltage_uncs = []
    for wavelength, wavelength_unc, voltage, voltage_unc in data:
        wavelengths.append(wavelength)
        wavelength_uncs.append(wavelength_unc)
        voltages.append(voltage)
        voltage_uncs.append(voltage_unc)

    def freq_from_wavelength(wavelength: float) -> float:
        """
        Takes a wavelength in nanometers and converts it to frequency.
        """
        return c / (wavelength * 10**(-9))
    
    fig, axis = plt.subplots()
    fig.suptitle("Stopping Potential vs. Frequency")
    fig.supylabel("Voltage (V)")
    fig.supxlabel("Frequency (Hz)")

    frequencies = [freq_from_wavelength(i) for i in wavelengths]
    # for calculating frequency uncertainties, use δν = √[(∂ν/∂λ * δλ)²] = |c/(λ²) δ λ|
    # we assume that there is no uncertainty in the measurement of `c`, the speed of light
    frequency_uncs = [sqrt(((-c / (wavelengths[j]**2)) * wavelength_uncs[j])**2) for j in range(len(wavelengths))]

    # plot raw data
    axis.errorbar(frequencies, voltages, xerr=frequency_uncs, yerr=voltage_uncs)

    # calculate line of best fit
    m, b, m_unc, b_unc = dat.lineFitWt(frequencies, voltages, voltage_uncs)
    line = np.vectorize(lambda x: m * x + b)
    # plot line of best fit
    line_x_coordinates = np.linspace(min(frequencies), max(frequencies), 100) # generate 100 evenly spaced points
    axis.plot(line_x_coordinates, line(line_x_coordinates), label="Chi Squared Fit")
    fig.legend()

    print(f"Fit Quality: Q = {dat.fitQuality(frequencies, voltages, voltage_uncs, m, b)}")
    print(f"lineFitWt results:\n* m = ({m} ± {m_unc}) V/Hz\n* b = ({b} ± {b_unc}) V")

    elementary_charge = 1.602176634 * 10**(-19) # NIST Value: https://physics.nist.gov/cgi-bin/cuu/Value?e
    elementary_charge_unc = 0
    h = m * elementary_charge
    h_unc = sqrt( (elementary_charge * m_unc)**2 + (m * elementary_charge_unc)**2 )
    print(f"h = {h} ± {h_unc}")
    

# for data in (digital_data_plot() + analog_data_plot()):
    # print(data)
data = digital_data_plot()
data = analog_data_plot()
stopping_potential_vs_frequency(data)
# plt.yscale('log')
plt.show()
