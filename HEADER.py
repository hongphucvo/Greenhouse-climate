/*IMPORT */
import pandas as pd
import numpy as np
import math

/*CLEAN DATA*/
data = pd.read_csv("Greenhouse_climate.csv")
data1 = pd.read_csv("meteo.csv")
data = data.dropna()
data1 = data1.dropna()

/*FORMULAR*/
Iglob = data1['Iglob'].values
vwind = data1['Windsp'].values
Tout = data1['Tout'].values
CO2_Air = data['CO2air'].values
CO2_Air = 0.0409 * CO2_Air * 44
TAir = data['Tair'].values
RHair = data['RHair'].values
VP_Air = RHair * 25 
VP_Out = data1['Rhout'].values * 21 

data1.head()
data.head()

/*CONSTANT*/
c_evap1 = 4.3
c_evap2 = 0.54
c_evap3day = 6.1 * 10 ** (-7)
c_evap3night = 1.1 * 10 ** (-11)
c_evap4day = 4.3 * 10 ** (-6)
c_evap4night = 5.2 * 10 ** (-6)
cPAir = 10 ** 3 # J K-1 kg-1 Table 8.1
CBufMax = 20 * 10 ** 3
deltaH = 2.45 * 10 ** 6 # J kg-1 Table 8.1
gamma = 65.8 #Pa K-1
g = 9.81 
M_CH2O = 30 * 10 ** (-3) # mg {CH2O}μmol-1 {CH2O}
M_water = 18  #kg/kmol
M_air = 28.96 #kg/kmol
R = 8.314 * 10 ** 3 # Hằng số lý tưởng J kmol-1 K-1
rb = 275 #s m-1
rs_min = 82 #s m-1
H = 22 * 10 ** 4 #deactivation energy
S = 710

/*VARIABLE*/
# Chọn vùng Texas
AFlr = 78000 #m2 Table 8.2
ACov = 9 * 10 ** 4 #m2 Table 8.2
alpha = 0.385 #μmol {e-} μmol-1 {photons} Table 9.1
Cgh_d = 0.65 # Table 8.2
Cgh_w = 0.09 # Table 8.2
#cLeaf = 0.28 
cHECin = 1.86 #W/m2K Table 8.2
cLeakage = 10 ** (-4) # Table 8.2
c_Gamma = 1.7 #μmol {CO2} mol-1 {air} K-1 Table 9.1
CO2Out = 668 #mg/m^3 2.4.2
COPMechCool = 0 #Hệ thống không có MechCool
C_Buf = 15 * 10 ** 3 #C_Buf < C_MaxBuf
C_MaxBuf = 20 * 10 ** 3 #mg {CH2O} m-2 Table 9.1
E_j = 37 * 10 ** 3 #Jmol-1 Table 9.1
gamma_Gh = 0.78 # Table 9.1
hVent = 0.97 #m 4.6.3
hElevation = 1470 #m 4.6.3
hSideRoof = 3 #m 4.6.3 0?
hRoof = 4 * 10 ** (-3) #m Table 8.2
hAir = 4.7 #m Table 8.2
hMean = 5.1 #m Table 8.2
J_Max25Leaf = 210 #μmol {e-} m-2 {leaf} s-1 Table 9.1
KThScr = 0.25 * 10 ** (-3) #m3 m-2 K-0.66 s-1 Table 8.2
K1 = 0.7 # table 8.1
K2 = 0.7 # table 8.1
LAI = 3.5 # = LAI_max
etaCO2Air_Stom = 0.67 #μmol {CO2} mol-1 {air} Table 9.1
etaGlobAir = 0.1
etaGlobPar = 0.5 #Table 8.1
etaGlobNir = 0.5 #Table 8.1
etamg_ppm = 0.554
etaHeatVap = 4.43 * 10 ** (-8) #kg vapour J-1 Table 8.1
etaHeatCO2 = 0.057 #mg CO2 J-1 Table 8.1
etaSide = 0.7
etaRoof = 0.3
etaSide_Thr = 0.9 #Table 8.1
etaRoof_Thr = 0.9 #Table 8.1
etaPad = 0 # Hệ thống không có Pad
phiFog = 0 # Hệ thống không có Fog System
phiExtCO2 = 4.3 * 10 ** 5 #mg s-1 Table 8.2
phiVentForced = 0 # Hệ thống không có VentForced
phiPad = 0 # Hệ thống không có Pad
PBlow = 0.50 * 10 ** 6 # W 5.3.3.1
PMechCool = 2000
rhoAir0 = 1.20 # kg/m3 Table 8.1
rhoCan = 0.07
rhoFlr = 0.5
RCan_SP = 5
stigmaInsScr = 1
UBlow = 0 #0.5 #[0:1]
UPad = 0 # Hệ thống không có Pad
URoof = 0.279 # (VentLee + VentWind) / 2
UFog = 0 # Hệ thống không có Fog System
USide = 0.4
UExtCO2 = 0.75
UThScr = 0 # EnScr = 0
UVentForced = 0 # Hệ thống không có VentForced
UMechCool = 0.5 # Hệ thống không có MechCool
tauCovPar = 0.177  # tính bằng công thức tau
tauCovNir = 0.177 # tính bằng công thức tau
T_CanK = 292.35 
T_25K = 298.15 #25oC
theta = 0.7 # Table 9.1
xPad = 0 #no data
xOut = 0 #no data

/*RELATED VARIABLE*/
ARoof = AFlr * 0.18
ASide = AFlr * 0.09
capCO2Air = hAir
capCO2Top = hMean - hAir
Cd = Cgh_d
Cw = Cgh_w
rhoAir = rhoAir0 * pow(np.exp(1), (g * M_air * hElevation / (293.15 * R)))
rhoTop = rhoAir
rhoAirMean = (rhoAir + rhoTop) / 2

/*GRAPH*/  
import plotly.express as px
import plotly.graph_objects as go
fig = go.Figure()

from sklearn.metrics import mean_squared_error 
