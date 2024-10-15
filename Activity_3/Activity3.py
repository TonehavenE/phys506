import DataAnalysisTools as dat
import numpy as np

# r_330k = [340 * 1000, 340 * 1000, 340 * 1000, 340 * 1000, 340 * 1000, 340 * 1000]
# r_56 = [59, 59, 59, 59, 59]
# r_1 = [4, 4, 4, 4]
# r_390 = [400, 400, 400, 400]
#
# def print_result(data: list, percent: float) -> None:
#     print(f"{dat.totalUncertainty(data, dat.typeBFromPercent(percent, np.mean(data)))}")
#
# print_result(r_330k, 1.5/100)
# print_result(r_56, 1.5/100)
# print_result(r_1, 1.5/100)
# print_result(r_390, 1.5/100)

r_330k_predicted = 330 * 1000 
r_330k_predicted_uncertainty = 2 * 3300

r_330k_actual = 340 * 1000
r_330k_actual_uncertainty = 10200

r_56_predicted = 56
r_56_predicted_uncertainty = 2 * 0.56

r_56_actual = 56
r_56_actual_uncertainty = 2

r_1_predicted = 1
r_1_predicted_uncertainty = 2 * 0.01
r_1_actual = 1
r_1_actual_uncertainty = 0.03

r_390_predicted = 390
r_390_predicted_uncertainty = 2 * 3.9
r_390_actual = 400
r_390_actual_uncertainty = 12


print(f"330k Resistor: {dat.discrepancy(r_330k_predicted, r_330k_predicted_uncertainty, r_330k_actual, r_330k_actual_uncertainty)}")
print(f"56 Resistor: {dat.discrepancy(r_56_predicted, r_56_predicted_uncertainty, r_56_actual, r_56_actual_uncertainty)}")
print(f"1 Resistor: {dat.discrepancy(r_1_predicted, r_1_predicted_uncertainty, r_1_actual, r_1_actual_uncertainty)}")
print(f"390 Resistor: {dat.discrepancy(r_390_predicted, r_390_predicted_uncertainty, r_390_actual, r_390_actual_uncertainty)}")

r_theory = 57
r_theory_uncertainty = 0.56

r_actual = 57
r_actual_uncertainty = 2

print(f"Series: {dat.discrepancy(r_theory, r_theory_uncertainty, r_actual, r_actual_uncertainty)}")

parallel_theory = 49
parallel_theory_uncertainty = 1

parallel_actual = 49
parallel_actual_uncertainty = 2

print(f"Parallel: {dat.discrepancy(parallel_theory, parallel_theory_uncertainty, parallel_actual, parallel_actual_uncertainty)}")

voltage_theory = 5.0
voltage_theory_uncertainty = 0.05

voltage_actual = 4.8
voltage_actual_uncertainty = 0.1

print(f"Voltage: {dat.discrepancy(voltage_theory, voltage_theory_uncertainty, voltage_actual, voltage_actual_uncertainty)}")

