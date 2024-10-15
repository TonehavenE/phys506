import numpy as np
import matplotlib.pyplot as plt
import DataAnalysisTools as dat

def main():
    time, velocity, uncertainty = read_velocity_data()

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


def read_velocity_data() -> tuple[list, list, list]:
    time, velocity, uncertainty = [], [], []
    with open("velocity_data.txt", "r") as f:
        lines = f.readlines()
        for line in lines:
            t, v, u = line.split(" ")
            time.append(float(t))
            velocity.append(float(v))
            uncertainty.append(float(u))

    return (time, velocity, uncertainty)

if __name__ == "__main__":
    main()
