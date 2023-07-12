import numpy as np

def lonlat2xy(lon, lat, lon0, lat0):
    lon = np.deg2rad(lon)
    lat = np.deg2rad(lat)
    lon0 = np.deg2rad(lon0)
    lat0 = np.deg2rad(lat0)
    a = 6378137
    F = 298.257222101
    m0 = 0.9999
    #rhopp = 180/np.pi
    n = 1/(2*F-1)
    alpha = np.zeros(5)
    alpha[0] = n/2 - 2*(n**2)/3 + 5*(n**3)/16 + 41*(n**4)/180 - 127*(n**5)/288
    alpha[1] = 13*(n**2)/48 - 3*(n**3)/5 + 557*(n**4)/1440 + 281*(n**5)/630
    alpha[2] = 61*(n**3)/240 - 103*(n**4)/140 + 15061*(n**5)/26880
    alpha[3] = 49561*(n**4)/161280 - 179*(n**5)/168
    alpha[4] = 34729*(n**5)/80640
    t = np.sinh(np.arctanh(np.sin(lat))-2*np.sqrt(n)*np.arctanh(2*np.sqrt(n)*np.sin(lat)/(1+n))/(1+n))
    tb = np.sqrt(1+t*t)
    lambda_c = np.cos(lon-lon0)
    lambda_s = np.sin(lon-lon0)
    xi = np.arctan(t/lambda_c)
    eta = np.arctanh(lambda_s/tb)
    A = np.zeros(6)
    A[0] = 1 + (n**2)/4 + (n**4)/64
    A[1] = -3*(n-(n**3)/8-(n**5)/64)/2
    A[2] = 15*(n**2-(n**4)/4)/16
    A[3] = -35*(n**3-5*(n**5)/16)/48
    A[4] = 315*(n**4)/512
    A[5] = -693*(n**5)/1280
    Ab = m0*a*A[0]/(1+n)
    sum = 0
    for i in range(5):
        sum += A[i+1]*np.sin(2*(i+1)*lat0)
    Sphi0 = m0*a*(A[0]*lat0 + sum)/(1+n)
    sumx = 0
    sumy = 0
    for i in range(5):
        sumx += alpha[i]*np.sin(2*(i+1)*xi)*np.cosh(2*(i+1)*eta)
        sumy += alpha[i]*np.cos(2*(i+1)*xi)*np.sinh(2*(i+1)*eta)
    x = Ab*(xi+sumx) - Sphi0
    y = Ab*(eta+sumy)

    xy = np.array([y, x]) # x: North, y: East
    
    return xy
