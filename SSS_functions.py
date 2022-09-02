import numpy as np

'''
SSS functions
'''


def bottom_track(swath, mean, N):
    first_value, new_port_chan = np.zeros(len(swath)), np.zeros(len(swath))
    start = len(swath) // 2
    finish = len(swath) - 1
    n = N // 2
    average = mean / 2
    for i in range(start, len(swath)):
        try:
            if swath[i] >= average and swath[i]-swath[i-1]<=0.2 and abs(swath[i]-swath[i+1])>=0.2:
                first_value[i] = 1 # Set row value equal to 1 to identify ocean floor
        except IndexError:
            pass
    while finish >= start:
        if first_value[finish] == 1:
            value = swath[finish] # First Backscatter value of ocean floor
            break
        else:
            finish -= 1
    return first_value, value, finish

def TVG_factors(ts):
    factors = np.ones_like(ts)
    factor = [[1,0]]
    i=0
    mean = np.mean(ts)
    while i < len(ts):
        try:
            if ts[i+2] == 0:
                pass
            else:
                factor.append([mean / (ts[i+2]), i+2])
                factors[i+2] = mean / (ts[i+2])
            i += 4
        except IndexError:
            factor.append([1,len(ts)-1])
            break
    
    for i, (value, index) in enumerate(factor):
        try:
            a1 = value
            a2 = factor[i+1][0]
            b1 = index
            b2 = factor[i+1][1]
            values = [a1, a2]
            times = [b1, b2]
            interpolator = CubicSpline(times, values)
            for a in range(b1+1, b2):
                factors[a] = interpolator(a)
        except IndexError:
            pass
    return factors
