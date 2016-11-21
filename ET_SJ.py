import numpy as np

def ET_SJ(Tmax,Tmin,RHmax,RHmin,Rs,u2,J,lat,Elev, alpha, z0):
    #Taken and simplified from package Evapotranspiration >> Danlu Guo <danlu.guo@adelaide.edu.au>
    #Daily Actual Evapotranspiration after Szilagyi, J. 2007, doi:10.1029/2006GL028708
    #Time Series as daily values zoo objects: T[degC], RH[%], Rs[MJ/m2], u2[km/h], J[julian day], 
    #lat[deg], Elev[m]
    #alpha->penman landuse, z0->reference evap surface level[m]
    #Example: ET_SJ(Tmax,Tmin,RHmax,RHmin,Rs,u2,J,lat=49.8,Elev=440.0,alpha = 0.23, z0 = 0.02)
    #Result in [mm/day]
    lat_rad=lat*np.pi/180.
    alphaPT = 1.31
    sigma=4.903e-09
    Gsc=0.082
    lambdax=2.45
    Ta = (Tmax + Tmin)/2
    vs_Tmax = 0.6108 * np.exp(17.27 * Tmax/(Tmax + 237.3))
    vs_Tmin = 0.6108 * np.exp(17.27 * Tmin/(Tmin + 237.3))
    vas = (vs_Tmax + vs_Tmin)/2.
    vabar = (vs_Tmin * RHmax/100. + vs_Tmax * RHmin/100.)/2.
    P = 101.3 * ((293. - 0.0065 * Elev)/293)**5.26
    delta = 4098. * (0.6108 * np.exp((17.27 * Ta)/(Ta + 237.3)))/((Ta + 237.3)**2.)
    gamma = 0.00163 * P/lambdax
    d_r2 = 1. + 0.033 * np.cos(2. * np.pi/365. * J)
    delta2 = 0.409 * np.sin(2. * np.pi/365. * J - 1.39)
    w_s = np.arccos(-np.tan(lat_rad) * np.tan(delta2))
    N = 24./np.pi * w_s
    R_a = (1440./np.pi) * d_r2 * Gsc * (w_s * np.sin(lat_rad) * np.sin(delta2) + np.cos(lat_rad) * np.cos(delta2) * np.sin(w_s))
    R_so = (0.75 + (2. * 10.**-5.) * Elev) * R_a
    
    R_nl = sigma * (0.34 - 0.14 * np.sqrt(vabar)) * ((Tmax + 273.2)**4. + (Tmin + 273.2)**4.)/2. * (1.35 * Rs/R_so - 0.35)
    R_nsg = (1. - alpha) * Rs
    R_ng = R_nsg - R_nl
    
    f_u = 1.313 + 1.381 * u2 #version 1956
    #f_u = 2.626 + 1.381 * u2 #version 1948
    
    Ea = f_u * (vas - vabar)
    Epenman_Daily = delta/(delta + gamma) * (R_ng/lambdax) + gamma/(delta + gamma) * Ea
    
    T_e = Ta
    for i in np.arange(99999):
        v_e = 0.6108 * np.exp(17.27 * T_e/(T_e + 237.3))
        T_enew = Ta - 1./gamma * (1. - R_ng/(lambdax * Epenman_Daily)) * (v_e - vabar)
        deltaT_e = T_enew - T_e
        maxdeltaT_e = np.abs(np.max(deltaT_e[-np.isnan(deltaT_e)]))
        T_e = T_enew
        if (maxdeltaT_e < 0.01):
            break
    
    deltaT_e = 4098. * (0.6108 * np.exp((17.27 * T_e)/(T_e + 237.3)))/((T_e + 237.3)**2.)
    E_PT_T_e = alphaPT * (deltaT_e/(deltaT_e + gamma) * R_ng/lambdax)
    E_SJ_Act_Daily = 2. * E_PT_T_e - Epenman_Daily
    
    return(E_SJ_Act_Daily)