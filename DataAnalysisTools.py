import math
import numpy as np

def statisticalUncertainty(xdata: list[float]) -> float:
    """
    Returns the standard deviation of the mean of the data set `xdata`.

    Inputs:
    - xdata: list, repeated measurements

    Returns:
    - float: the standard deviation of `xdata`
    """
    N = len(xdata)
    mean = np.mean(xdata)
    series_sum = sum([(datapoint - mean)**2 for datapoint in xdata])

    return math.sqrt((1/(N * (N-1))) * series_sum)
    

def totalUncertainty(xdata: list[float], typeBUnc: float = 0.0) -> float:
    """
    Returns the total uncertainty associated with a measurement.
    
    Inputs:
    - xdata: list, repeated measurements
    - typeBUnc: float, the Type B uncertainty associated with the measurements

    Returns:
    - float: the total uncertainty, combining Type A and Type B
    """
    typeAUnc = statisticalUncertainty(xdata)
    rawUnc = 2 * math.sqrt(typeAUnc**2 + typeBUnc**2) # 2 is the coverage factor
    return rawUnc # should I filter to 1 sig fig?

def typeBFromPercent(measurement: float, percent: float) -> float:
    return percent * measurement


def reportMeasurement(xdata: list[float], typeBUnc: float = 0.0) -> tuple[float, float]:
    """
    Calculates and returns the best estimate of the measurement from the `xdata` set and `typeBUnc` uncertainty.
    
    Inputs:
    - xdata: list, repeated measurements
    - typeBUnc: float, the Type B uncertainty associated with the measurements

    Returns:
    - tuple[float, float]:
        - the best guess for the value
        - the total uncertainty
    """
    xBest = np.mean(xdata)
    delX = totalUncertainty(xdata, typeBUnc)
    return (xBest, delX)

def relativeUncertainty(xdata: list[float], typeBUnc: float = 0.0) -> float:
    """
    Calculates the relative uncertainty of a set of measurements `xdata` taken with Type B uncertainty given by `typeBUnc`

    Inputs:
    - xdata: list, repeated measurements
    - typeBUnc: float, the Type B uncertainty associated with the measurements

    Returns:
    - float: the relative uncertainty of the measurement
    """
    xBest, delX = reportMeasurement(xdata, typeBUnc)
    return delX/abs(xBest)

def discrepancy(xbest1: float, dx1: float, xbest2: float, dx2: float) -> tuple[float, float, bool]:
    """
    Takes two measurements and computes the discrepancy.
    - xbest1: float, the best estimate for the first measurement
    - dx1: float, the uncertainty of the first measurement
    - xbest2: float, the best estimate for the second measurement
    - dx2: float, the uncertainty of the second measurement
    Returns:
    - tuple[float, float, bool]:
        - |xbest1 - xbest2|
        - |dx1 - dx2|
        - determination of agreement
    """
    discrepancy = abs(xbest1 - xbest2)
    significanceCriterion = abs(dx1 - dx2)

    # determination of agreement?
    return (discrepancy, significanceCriterion, discrepancy < significanceCriterion)

def combineMeasurements(xresults: list[float], dxresults: list[float]) -> tuple[float, float]:
    """
    Takes a set of best estimates and uncertainties from a series of independent trials and returns the combined best measurement result.

    Inputs:
    - xresults: list[float], the best estimates of each trial
    - dxresults: list[float], the uncertainties of each trial

    Returns:
    - tuple[float, float]:
        - the combined best guess for x
        - the combined uncertainty
    """
    assert len(xresults) == len(dxresults)

    def w(i):
        # returns the weight
        return 1/(dxresults[i]**2)
    
    # calculate x_best
    numerator = 0
    denominator = 0
    for i in range(len(xresults)):
        numerator += w(i)* xresults[i]
        denominator += w(i)

    xbest = numerator/denominator

    dxbest = 1/math.sqrt(denominator)

    return (xbest, dxbest)

def lineFit(x: list[float], y: list[float]) -> tuple[float, float]:
    """
    Calculates the line of best fit for a set of data.

    Inputs:
    - x: list, the input x-coordinates
    - y: list, the input y-coordinates

    Returns:
    - tuple[float, float]:
        - float, m: the slope of the line of best fit
        - float, b: the y-intercept of the line of best fit
    """
    assert len(x) == len(y) # the data set should be symmetric

    N = len(x)

    E_x = 1/N * sum(x)
    E_y = 1/N * sum(y)
    E_xx = 1/N * sum([i**2 for i in x])
    E_xy = 1/N * sum([x[j] * y[j] for j in range(N)])

    m = (E_xy - E_x*E_y)/(E_xx - E_x**2)
    b = (E_xx * E_y - E_x * E_xy)/(E_xx - E_x**2)

    return (m, b)

def relativeError(xmeasured: float, xaccepted: float) -> float:
    """
    Returns the realtive error for a measurement.

    Inputs:
    - xmeasured: float, the measured value for x
    - xaccepted: float, the accepted value for x

    Returns:
    - float, the relative error as |x_measured - x_accepted|/x_accepted
    """
    return abs((xmeasured - xaccepted)/xaccepted)


def lineFitWt(x: list[float], y: list[float], dy: list[float]) -> tuple[float, float, float, float]:
    """
    Returns the line of best fit parameters using weighted uncertainty.

    Inputs:
    - x: list, the input x-coordinates
    - y: list, the input y-coordinates
    - dy: list, the uncertainty of y measurements

    Returns:
    - tuple[float, float]:
        - float, m: the slope of the line of best fit
        - float, b: the y-intercept of the line of best fit
        - float: mUnc: the uncertainty of the slope
        - float: bUnc: the uncertainty of the y-intercept
    """
    assert len(x) == len(y)
    assert len(y) == len(dy)

    w = [1/(i**2) for i in dy]
    N = len(x)

    S_wxx = sum([w[i] * x[i]**2 for i in range(N)])
    S_wy = sum([w[i] * y[i] for i in range(N)])
    S_wx = sum([w[i] * x[i] for i in range(N)])
    S_wxy = sum([w[i] * x[i] * y[i] for i in range(N)])
    S_w = sum(w)

    m_w = (S_w * S_wxy - S_wx * S_wy)/(S_w * S_wxx - (S_wx)**2)
    b_w = (S_wxx * S_wy - S_wx * S_wxy)/(S_w * S_wxx - (S_wx)**2)
    
    mUnc = math.sqrt( S_w / (S_w * S_wxx - S_wx**2) )
    bUnc = math.sqrt( S_wxx / (S_w * S_wxx - S_wx**2) )

    return (m_w, b_w, mUnc, bUnc)

def fitQuality(x: list[float], y: list[float], dy: list[float], m: float, b: float) -> float:
    """
    Determines how well a line of best fit fits to a data set by using the chi-squared parameter Q.

    Inputs:
    - x: list, the input x-coordinates
    - y: list, the input y-coordinates
    - dy: list, the uncertainty of y measurements
    - m: float, the slope of the line of best fit
    - b: float, the y-intercept of the line of best fit

    Returns:
    - float: Q, the chi-squared value for the fit
    """
    assert len(x) == len(y)
    assert len(x) == len(dy)

    N = len(x)
    
    total = 0
    for i in range(N):
        total += ((y[i] - m*x[i] - b)/dy[i])**2

    return total/(N-2)




 



    

    
