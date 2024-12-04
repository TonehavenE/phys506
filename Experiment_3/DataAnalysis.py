from math import log
import DataAnalysisTools as dat
import matplotlib.pyplot as plt 
import numpy as np 

def gus_algorithm(voltages, currents, dark_v, dark_c):
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
        for d in range(len(dark_v)):
            if dark_v[d] > voltages[i]: # iterate till we find a dark characteristic voltage data point bigger than the voltage
                # we interpolate between the first point we find that is bigger than our voltage and the previous point (assumed smaller)
                # corrected is the current at this voltage minus the interpolated value of the characteristic curve at this voltage
                corrected.append(currents[i] - linear(voltages[i], dark_v[d], dark_c[d], dark_v[d-1], dark_c[d-1]))
                break

    corrected.sort()

    # now we need to find where the shifted current values intersect y = 0
    for i in range(len(corrected)):
        if corrected[i] == 0:
            x_value = voltages[i]
            break
        elif corrected[i] > 0:
            x_value = linear_for_x(0, voltages[i], corrected[i], voltages[i-1], corrected[i-1])
            break

    return x_value

def digital_data_plot():
    # files = ["Analog3650", "Analog4047", "Analog4358", "Analog5461", "Analog5770"]
    files = ["RedLED_611_nm", "YellowLED_588_nm", "GreenLED_525_nm", "CyanLED_505_nm", "BlueLED_472_nm"]
    # files = ["RedLED_611_nm", "YellowLED_588_nm"]

    fig, axis = plt.subplots(1, len(files), sharey=True)
    fig.suptitle("Digital Experimental Data")
    fig.supxlabel("Accelerating Voltage (V)")
    fig.supylabel("Photocurrent (nA)")

    for index, file in enumerate(files):
        print(f"{file}")

        voltages = np.loadtxt(f"data/{file}.csv", delimiter=',', usecols=0, skiprows=1)
        voltage_unc = np.loadtxt(f"data/{file}.csv", delimiter=',', usecols=1, skiprows=1)
        currents = np.loadtxt(f"data/{file}.csv", delimiter=',', usecols=2, skiprows=1)
        current_unc = np.loadtxt(f"data/{file}.csv", delimiter=',', usecols=3, skiprows=1)

        # manufacture dark characteristic data to match as 0 - this works for the digital measurement
        # we found experimentally that the dark current is 0 everywhere for this instrument
        dark_v = voltages
        dark_v_unc = [0 for _ in voltage_unc]
        dark_c = [0 for _ in currents]
        dark_c_unc = [0 for _ in current_unc]

        axis[index].plot(voltages, currents)
        axis[index].plot(dark_v, dark_c)
        axis[index].set_title(f"{file}")

        x_value = gus_algorithm(voltages, currents, dark_v, dark_c)
        print("Gus's Algorithm:", x_value)
        axis[index].plot([x_value, x_value], [-2, 0.5 ]) 
        print("")

def analog_data_plot():

    files = ["Analog3650A", "Analog4047A", "Analog4358A", "Analog5461A", "Analog5770A"]
    

    dark_v = np.loadtxt(f"data/AnalogDark.csv", delimiter=',', usecols=0, skiprows=1)
    dark_v_unc = np.loadtxt(f"data/AnalogDark.csv", delimiter=',', usecols=0, skiprows=1)
    dark_c = np.loadtxt(f"data/AnalogDark.csv", delimiter=',', usecols=0, skiprows=1)
    dark_c_unc = np.loadtxt(f"data/AnalogDark.csv", delimiter=',', usecols=0, skiprows=1)

    fig, axis = plt.subplots(1, len(files), sharey=True)
    fig.suptitle("Analog Experimental Data")
    fig.supxlabel("Accelerating Voltage (V)")
    fig.supylabel("Photocurrent (10^-14)")

    for index, file in enumerate(files):
        voltages = np.loadtxt(f"data/{file}.csv", delimiter=',', usecols=0, skiprows=1)
        voltage_unc = np.loadtxt(f"data/{file}.csv", delimiter=',', usecols=1, skiprows=1)
        currents = np.loadtxt(f"data/{file}.csv", delimiter=',', usecols=2, skiprows=1)
        current_unc = np.loadtxt(f"data/{file}.csv", delimiter=',', usecols=3, skiprows=1)


        axis[index].plot(voltages, currents)
        axis[index].plot(dark_v, dark_c)
        axis[index].set_title(f"{file}")

        x_value = gus_algorithm(voltages, currents, dark_v, dark_c)
        print("Gus's Algorithm:", x_value)
        axis[index].plot([x_value, x_value], [min(currents), max(currents)]) 
        print("")

# digital_data_plot()
analog_data_plot()
plt.yscale('log')
plt.show()
