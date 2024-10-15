import numpy as np
import matplotlib.pyplot as plt
import DataAnalysisTools as dat

def main():

    # Get data
    x, y = read_millikan()
    m, b = dat.lineFit(x, y)
    line = np.vectorize(lambda x: m*x + b)

    line_x = np.linspace(min(x), max(x), 100)
    line_y = line(line_x) 

    plt.plot(x, y)
    plt.xlabel("$\nu$")
    plt.ylabel("$V$")

    plt.plot(line_x, line_y)
    plt.show()

    # our slope, m, should be h/e. So m = h/e -> h = me
    e = 1.602176634 * 10**(-19)
    h_measured = m * e
    print(f"Our experimental value for h is: {h_measured = }")

    relative_h_error = dat.relativeError(h_measured, 6.62607015*10**(-34))
    print(f"The relative error is: {relative_h_error * 100}%")

def read_millikan():
    x, y = [], []

    with open("millikan.txt", "r") as f: 
        lines = f.readlines()
        for line in lines: 
            i, j = line.split(" ")
            x.append(float(i))
            y.append(float(j))

    return x, y

if __name__ == "__main__":
    main()
