import numpy as np
import matplotlib.pyplot as plt
import DataAnalysisTools as dat

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
    resistor, resistor_uncertainties = read_data("resistor_data.txt")
    current, current_uncertainties = read_data("current_data.txt")

    plt.errorbar(time, velocity, yerr=uncertainty)
    plt.xlabel("time (s)")
    plt.ylabel("velocity (m/s)")

    m1, b1 = dat.lineFit(time, velocity)
    line1 = np.vectorize(lambda x: m1 * x + b1)
    print(f"Least Squares: v = {m1}t + {b1}")
    print(f"Least Squares Fit Quality: Q = {dat.fitQuality(time, velocity, uncertainty, m1, b1)}")
    print()

    m2, b2, m2Unc, b2Unc = dat.lineFitWt(time, velocity, uncertainty)
    line2 = np.vectorize(lambda x: m2 * x + b2)
    print(f"Chi Squared: v = {m2}t + {b2}")
    print(f"Chi Squared Fit Quality: Q = {dat.fitQuality(time, velocity, uncertainty, m2, b2)}")
    print()
    
    print(f"lineFitWt results:\n* m = ({m2} ± {m2Unc}) m/s^2\n* b = ({b2} ± {b2Unc}) m/s")

    line_x = np.linspace(min(time), max(time), 100)

    plt.plot(line_x, line1(line_x), label="Least Squares")
    plt.plot(line_x, line2(line_x), label="Chi Squared")

    plt.legend()
    # plt.show()
    plt.savefig("Activity6_Problem3.png")
    

if __name__ == "__main__":
    main()
