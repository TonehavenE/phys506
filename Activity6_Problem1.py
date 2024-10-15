import DataAnalysisTools as dat

def main():
    # ds = dataset
    ds_1 = [72.2, 77.6, 82.4, 86.3, 88.9] # cm
    ds_1_type_b = 0.3 # cm

    ds_2 = [80.10, 81.45, 81.50, 81.34, 82.01] # cm
    ds_2_type_b = 0.05 # cm

    # Print results
    # i. 
    print("(i.)")
    xbest1, dx1 = dat.reportMeasurement(ds_1, ds_1_type_b)
    xbest2, dx2 = dat.reportMeasurement(ds_2, ds_2_type_b)

    print(f"The first measurement is: {round(xbest1, 1)} \xB1 {round(dx1, 1)}")
    print(f"The second measurement is: {round(xbest2, 1)} \xB1 {round(dx2, 1)}")

    # print(f"The first measurement is: {xbest1} +- {dx1}")
    # print(f"The second measurement is: {xbest2} +- {dx2}")

    # ii.
    print("\n(ii.)")
    relativeUncertainty1 = dat.relativeUncertainty(ds_1, ds_1_type_b)
    relativeUncertainty2 = dat.relativeUncertainty(ds_2, ds_2_type_b)

    print(f"The relative uncertainty of the first measurement is: {relativeUncertainty1}")
    print(f"The relative uncertainty of the second measurement is: {relativeUncertainty2}")
    print(f"The percent uncertainty of the first measurement is: {relativeUncertainty1 * 100}")
    print(f"The percent uncertainty of the second measurement is: {relativeUncertainty2 * 100}")

    # iii.
    print("\n(iii.)")
    discrepancy, significance, agreement = dat.discrepancy(xbest1, dx1, xbest2, dx2)
    if agreement:
        print("The two measurements are in agreement.")
    else:
        print("The two measurements are not in agreement.")

    # iv.
    print("\n(iv.)")
    combined_best, combined_uncertainty = dat.combineMeasurements([xbest1, xbest2], [dx1, dx2])
    print(f"The combined measurement is: {combined_best} \xB1 {combined_uncertainty}")

if __name__ == "__main__":
    main() 


# int(np.ceil((np.log(measurement) / np.log(10)))) # number of decimal places
