import math
from numpy import mean
import DataAnalysisTools as dat

magnetic_constant = 1.25663706 * 10**(-6) # m kg / (s^2 A^2)
number_turns = 130

def calculate_coil_radius(internal_ds: list, external_ds: list, type_b: float) -> tuple[float, float]:
    """
    Returns the best coil radius and total uncertainty
    """
    assert len(internal_ds) == len(external_ds)
    average_diameters = [(internal_ds[i] + external_ds[i])/2 for i in range(len(external_ds))]
    best_diameter = float(mean(average_diameters))
    variance = mean([(d - best_diameter)**2 for d in average_diameters])
    type_a = math.sqrt(variance)
    print(type_a * 100)

    uncertainty = 2 * math.sqrt(type_a**2 + type_b**2)

    return (best_diameter/2, uncertainty/2)

def read_diameter_data(file_name: str):
    internal_ds = []
    external_ds = []
    with open(file_name) as f:
        for line in f:
            internal, external = line.split(",")
            internal_ds.append(float(internal)/100) # divide by 100 to convert cm to m
            external_ds.append(float(external)/100)

    return internal_ds, external_ds



def read_voltage_data(file_name: str):
    """
    returns tuple containing:
    1. voltage: int
    2. beam_loop_diameters: list[int]
    3. currents: list[int]
    """
    voltage = int(file_name.split("v.")[0])
    currents = []
    beam_loop_diameters = []
    with open(file_name) as f:
        for line in f:
            current, beam_loop_diameter = line.split(",")
            currents.append(float(current))
            beam_loop_diameters.append(float(beam_loop_diameter) / 100) # divide by 100 to get meters

    return (voltage, beam_loop_diameters, currents)

def magnetic_uncertainty(current: float, current_unc: float, radius: float, radius_unc: float) -> float:
    """
    Returns the uncertainty in the magnetic field given the measurements of the coil radius and the current.
    """
    B_unc = math.sqrt(
            ((4/5)**3 * (magnetic_constant*number_turns/radius)**2 * (current_unc)**2) + ( (4/5)**3 * (-magnetic_constant * number_turns * current / (radius**2) )**2 * (radius_unc)**2 )
            )
    return B_unc

def magnetic_field(current: float, coil_radius: float) -> float:
    """
    Returns the magnetic field given the measurements of the current and coil radius.
    Assumes all measurements given in SI units, and returns Teslas.
    """
    B = (4/5)**(3/2) * magnetic_constant * number_turns * current / coil_radius
    return B

def charge_to_mass(voltage: float, radius: float, magnetic: float) -> float:
    """
    Returns the ratio of charge to mass.
    """
    return 2*voltage/(radius**2 * magnetic**2)

def charge_to_mass_unc(voltage: float, voltage_unc: float, radius: float, radius_unc: float, magnetic: float, magnetic_unc: float) -> float:
    """
    Returns the uncertainty in the charge to mass ratio.
    """
    return (2/(radius ** 2 * magnetic ** 2)) * math.sqrt( voltage_unc**2 + (4/(radius**2)) * radius_unc**2 + (4/(magnetic**2) * magnetic_unc **2))

def main():
    internal_ds, external_ds = read_diameter_data("diameter.txt")
    type_b_diameter = 0.05 / 100 # 0.05 cm to meters, constant based on ruler
    coil_radius, coil_radius_unc = calculate_coil_radius(internal_ds, external_ds, type_b_diameter)
    # print(f"{(coil_radius * 100):.3f}, unc = {(coil_radius_unc * 100):.3f}")

    total_charge_to_mass = []
    total_charge_to_mass_unc = []
    for x in ["200v.txt", "300v.txt", "400v.txt", "500v.txt"]:
        voltage, beam_loop_diameters, currents = read_voltage_data(x)

        voltage_unc = 0.05 * voltage # 5% of voltage?

        # current uncertainty is 5% of current
        B_data = [
                (magnetic_field(current, coil_radius), magnetic_uncertainty(current, 0.05 * current, coil_radius, coil_radius_unc)) for current in currents
        ]

        charge_to_mass_data = []
        for i in range(len(beam_loop_diameters)):
            r = beam_loop_diameters[i] / 2
            r_unc = 0.3 / 100 # 0.3 cm -> meters
            B, B_unc = B_data[i]
            c_to_m = charge_to_mass(voltage, r, B)
            c_to_m_unc = charge_to_mass_unc(voltage, voltage_unc, r, r_unc, B, B_unc)
            total_charge_to_mass.append(c_to_m)
            total_charge_to_mass_unc.append(c_to_m_unc)
            charge_to_mass_data.append((c_to_m, c_to_m_unc))

        print(f"MAGNETIC FIELD DATA")
        for b, b_unc in B_data:
            print(f"{b:.2e} ± {b_unc:.0e}")

        print("\n\n")

        print("CHARGE TO MASS DATA")
        for c_to_m, c_to_m_unc in charge_to_mass_data:
            print(f"{c_to_m:.2e} ± {c_to_m_unc:.0e}")

    # now calculate stuff
    weights = [1/unc**2 for unc in total_charge_to_mass_unc]
    charge_to_mass_total = sum([weights[i] * total_charge_to_mass[i] for i in range(len(weights))]) / sum(weights)
    charge_to_mass_unc_total = 1 / math.sqrt(sum(weights))
    
    charge_to_mass_total *= -1 # make negative
    nist_value = -1.75882001076 * 10**11
    nist_unc = 0.00000000053 * 10**11

    difference = abs(nist_value - charge_to_mass_total)
    total_unc = abs(charge_to_mass_unc_total + nist_unc)
    print(f"{difference} ?< {total_unc}")
    if difference < total_unc:
        print(f"Agreement")

    else:
        print(f"Disagreement")

    print(f"{charge_to_mass_total:.2e}")
    print(f"{charge_to_mass_unc_total:.0e}")

    # print(f"{magnetic_field(2.3, 14.8/100)}")
    # print(f"{charge_to_mass(200, 5.0/200, magnetic_field(2.3, 14.8/100))}")

if __name__ == "__main__":
    main()
