
/*FORMULAR*/

# (3)
def MCBlowAir(UBlow, PBlow, AFlr):
    return (etaHeatCO2 * UBlow * PBlow) / AFlr
# (4)
def MCExtAir(UExtCO2, phiExtCO2, AFlr):
  return (UExtCO2 * phiExtCO2) / AFlr
# (5)
def MCPadAir(UPad, phiPad, AFlr, CO2Out, CO2Air):
  return (UPad * phiPad * (CO2Out - CO2Air)) / AFlr
# (7) 
# rhoAir và rhoTop có hàm để tính, rhoAir là 8.24
def fThScr(UThScr, KThScr, TAir, TTop, rhoAirMean, rhoAir, rhoTop):
  return UThScr * KThScr * pow(abs(TAir - TTop), 2/3) + (1 - UThScr) * pow((g * (1 - UThScr) * abs(rhoAir - rhoTop)) / (2 * rhoAirMean), 1/2)
# (6)
def MCAirTop(UThScr, KThScr, TAir, TTop, rhoAirMean, rhoAir, rhoTop, CO2Top, CO2Air):
  return fThScr(UThScr, KThScr, TAir, TTop, rhoAirMean, rhoAir, rhoTop) * (CO2Air - CO2Top)# (10)
def fVentRoofSide(Cd, Cw, vWind, URoof, USide, ARoof, ASide, hSideRoof, TAir, TOut, TAirMean, AFlr):
  roof = URoof * ARoof
  side = USide * ASide
  c = (roof * roof * side * side) * 2 * g * hSideRoof * (TAir - TOut)
  d = (roof * roof + side * side) * TAirMean
  a = c / d
  b = pow((roof + side) / 2, 2) * Cw * vWind * vWind
  return (Cd * (a + b) ** (1/2)) / AFlr
# (11)
def etaIS(stigmaInsScr):
  return stigmaInsScr * (2 - stigmaInsScr)
# (12)
def fLeakage(vWind, cLeakage):
  if vWind < 0.25:
    return cLeakage * 0.25
  else:
    return vWind * cLeakage
# (13) dung cho airout; fVentSide_0 la f''ventSide tai Aroof = 0
def fVentSideAO(etaSide, etaInsScr, UThScr, fVentRoofSide, fVentSide_0, vWind, cLeakage):
  if etaSide >= etaSide_Thr:
    return etaInsScr * fVentSide_0 + 0.5 * fLeakage(vWind, cLeakage)
  else:
     return etaInsScr * (UThScr * fVentSide_0 + (1 - UThScr) * fVentRoofSide * etaSide) + 0.5 * fLeakage(vWind, cLeakage)
# (14)
def fVentForcedAO(etaInsScr, UVentForced, phiVentForced, AFlr):
  return (etaInsScr * UVentForced * phiVentForced) / AFlr
# (9)
def MCAirOut(CO2Out, CO2Air, stigmaInsScr, Cd, Cw, vWind, URoof, USide, ARoof, ASide, hSideRoof, TAir, TOut, TAirMean, AFlr, cLeakage, etaSide, UThScr, UVentForced, phiVentForced):
  etaInsScr = etaIS(stigmaInsScr)
  fVentSide_0 = (Cd * USide * ASide * vWind * (Cw ** (1/2))) / (2*AFlr)
  fVentSide_ARoof = fVentRoofSide(Cd, Cw, vWind, URoof, USide, ARoof, ASide, hSideRoof, TAir, TOut, TAirMean, AFlr)
  fVentSide = fVentSideAO(etaSide, etaInsScr, UThScr, fVentSide_ARoof, fVentSide_0, vWind, cLeakage)
  fVentForced = fVentForcedAO(etaInsScr, UVentForced, phiVentForced, AFlr)
  return (fVentSide + fVentForced) * (CO2Air - CO2Out)
  # (16)
def fVentRoofTO(etaRoof, etaInsScr, UThScr, fVentRoofSide, fVentRoof_1, fVentRoof_2, vWind, cLeakage):
  if etaRoof >= etaRoof_Thr:
    return etaInsScr * fVentRoof_1 + 0.5 * fLeakage(vWind, cLeakage)
  else: 
    return etaInsScr * (UThScr * fVentRoof_2 + (1 - UThScr) * fVentRoofSide * etaSide) + 0.5 * fLeakage(vWind, cLeakage)

# Hàm tính fVentRoof'' (17)
def fVentRoof_TO(Cd, URoof, ARoof, hVent, TAir, TOut, TAirMean, Cw, vWind, AFlr):
  return Cd * URoof * ARoof * pow((g * hVent * (TAir - TOut)) / (2 * TAirMean) + Cw * vWind * vWind, 1/2) / (2 * AFlr)
# (15)
def MCTopOut(CO2Top, CO2Out, stigmaInsScr, Cd, Cw, vWind, URoof, USide, ARoof, ASide, hSideRoof, hVent, TAir, TOut, TAirMean, AFlr, cLeakage, etaRoof, UThScr):
  etaInsScr = stigmaInsScr * (2 - stigmaInsScr)
  ffVentRoofSide  = fVentRoofSide(Cd, Cw, vWind, URoof, USide, ARoof, ASide, hSideRoof, TAir, TOut, TAirMean, AFlr)
  ffVentRoof_1 = fVentRoof_TO(Cd, URoof, ARoof, hVent, TAir, TOut, TAirMean, Cw, vWind, AFlr)
  ffVentRoof_2 = (Cd * USide * ASide * vWind * (Cw ** (1/2))) / (2*AFlr)
  ffVentRoof = fVentRoofTO(etaRoof, etaInsScr, UThScr, ffVentRoofSide, ffVentRoof_1, ffVentRoof_2, vWind, cLeakage)
  return ffVentRoof * (CO2Top - CO2Out)def P_(J, CO2_Stom, Gamma): 
  return J * (CO2_Stom - Gamma)/(4 * (CO2_Stom + 2 * Gamma))

def Jpot_(J_Max25Can, E_j, T_CanK, T_25K, S, H): 
  return J_Max25Can * math.exp(E_j*((T_CanK - T_25K)/(R*10**(-3)* T_CanK* T_25K)))*((1+math.exp((S*T_25K - H)/(R*10**(-3)*T_25K)))/(1+math.exp((S*T_CanK - H)/(R*T_CanK))))

def PARCan_(gamma_Gh, etaGlobPar, IGlob, rhoCan, K1, LAI, rhoFlr, K2): 
  PAR_Gh = gamma_Gh * etaGlobPar * IGlob
  PAR_GhCan = PAR_Gh * (1 - rhoCan)*(1 - math.exp(-1 * K1 * LAI))
  PAR_FlrCan = rhoFlr * PAR_Gh * (1 - rhoCan) * math.exp(-1 * K1 * LAI)* (1- math.exp(-1 * K2 * LAI))
  return PAR_GhCan + PAR_FlrCan

def J_(Jpot, alpha, PARCan, theta): 
  return (Jpot + alpha * PARCan - math.sqrt(pow(Jpot + alpha * PARCan, 2) - 4 * theta * Jpot * alpha * PARCan))/(2*theta)

def J_Max25Can_(LAI, J_Max25Leaf):
  return LAI * J_Max25Leaf

def Gamma_(J_Max25Leaf, J_Max25Can, c_Gamma, TCan):
  return (J_Max25Leaf/J_Max25Can) * c_Gamma * (TCan + 273.15) + 20 * c_Gamma * (1 - J_Max25Leaf/J_Max25Can)

def h_cBuf_(C_Buf, C_MaxBuf):
  if C_Buf > C_MaxBuf :
    return 0
  else:
    return 1

def R__(P, Gamma, CO2_Stom):
  return P * (Gamma / CO2_Stom)

def MCAirCan(CO2Air,M_CH2O, C_Buf, C_MaxBuf, LAI, etaCO2Air_Stom, alpha, theta, gamma_Gh , etaGlobPar, IGlob, rhoCan, rhoFlr, K1, K2, J_Max25Leaf, E_j, H, S, T_CanK, T_25K, c_Gamma, TCan):
  J_Max25Can = J_Max25Can_(LAI, J_Max25Leaf)
  Jpot = Jpot_(J_Max25Can, E_j, T_CanK, T_25K, S, H)
  PARCan = PARCan_(gamma_Gh, etaGlobPar, IGlob, rhoCan, K1, LAI, rhoFlr, K2)
  J = J_(Jpot, alpha, PARCan, theta)
  Gamma = Gamma_(J_Max25Leaf, J_Max25Can, c_Gamma, TCan)
  CO2_Stom = etaCO2Air_Stom * CO2Air
  P = P_(J, CO2_Stom, Gamma)
  h_cBuf = h_cBuf_(C_Buf, C_MaxBuf)
  R_ = R__(P, Gamma, CO2_Stom)
  return M_CH2O * h_cBuf* (P - R_)
def dx(CO2Air, CO2Top, TAir, TOut, IGlob, vWind):
  TTop = TAir - 1
  TOut = TAir - 1
  TCan = TAir - 1
  TMechCool = TAir - 1
  TThScr = TAir + 1
  TCov_in = TAir + 1
  TAirMean = (TAir + TOut) / 2
  mcBlowAir = MCBlowAir(UBlow, PBlow, AFlr)
  mcExtAir = MCExtAir(UExtCO2, phiExtCO2, AFlr)
  mcPadAir = MCPadAir(UPad, phiPad, AFlr, CO2Out, CO2Air)
  mcAirTop = MCAirTop(UThScr, KThScr, TAir, TTop, rhoAirMean, rhoAir, rhoTop, CO2Out, CO2Air)
  mcAirOut = MCAirOut(CO2Out, CO2Air, stigmaInsScr, Cd, Cw, vWind, URoof, USide, ARoof, ASide, hSideRoof, TAir, TOut, TAirMean, AFlr, cLeakage, etaSide, UThScr, UVentForced, phiVentForced)
  mcTopOut = MCTopOut(CO2Top, CO2Out, stigmaInsScr, Cd, Cw, vWind, URoof, USide, ARoof, ASide, hSideRoof, hRoof, TAir, TOut, TAirMean, AFlr, cLeakage, etaRoof, UThScr)
  mcAirCan = MCAirCan(CO2Air,M_CH2O, C_Buf, C_MaxBuf, LAI, etaCO2Air_Stom, alpha, theta, gamma_Gh , etaGlobPar, IGlob, rhoCan, rhoFlr, K1, K2, J_Max25Leaf, E_j, H, S, T_CanK, T_25K, c_Gamma, TCan)
  V_CO2Air = (mcBlowAir + mcExtAir + mcPadAir - mcAirCan - mcAirTop - mcAirOut) / capCO2Air
  V_CO2Top = (mcAirTop - mcTopOut) / capCO2Top
  return V_CO2Air, V_CO2Top
  
def dxAir(CO2Air, CO2Top, TAir, TOut, IGlob, vWind):
  TTop = TAir - 1
  TOut = TAir - 1
  TCan = TAir - 1
  TMechCool = TAir - 1
  TThScr = TAir + 1
  TCov_in = TAir + 1
  TAirMean = (TAir + TOut) / 2
  mcBlowAir = MCBlowAir(UBlow, PBlow, AFlr)
  mcExtAir = MCExtAir(UExtCO2, phiExtCO2, AFlr)
  mcPadAir = MCPadAir(UPad, phiPad, AFlr, CO2Out, CO2Air)
  mcAirTop = MCAirTop(UThScr, KThScr, TAir, TTop, rhoAirMean, rhoAir, rhoTop, CO2Out, CO2Air)
  mcAirOut = MCAirOut(CO2Out, CO2Air, stigmaInsScr, Cd, Cw, vWind, URoof, USide, ARoof, ASide, hSideRoof, TAir, TOut, TAirMean, AFlr, cLeakage, etaSide, UThScr, UVentForced, phiVentForced)
  mcAirCan = MCAirCan(CO2Air,M_CH2O, C_Buf, C_MaxBuf, LAI, etaCO2Air_Stom, alpha, theta, gamma_Gh , etaGlobPar, IGlob, rhoCan, rhoFlr, K1, K2, J_Max25Leaf, E_j, H, S, T_CanK, T_25K, c_Gamma, TCan)
  V_CO2Air = (mcBlowAir + mcExtAir + mcPadAir - mcAirCan - mcAirTop - mcAirOut) / capCO2Air
  return V_CO2Air
  
def dxTop(CO2Air, CO2Top, TAir, TOut, IGlob, vWind):
  TTop = TAir - 1
  TOut = TAir - 1
  TCan = TAir - 1
  TMechCool = TAir - 1
  TThScr = TAir + 1
  TCov_in = TAir + 1
  TAirMean = (TAir + TOut) / 2
  mcAirTop = MCAirTop(UThScr, KThScr, TAir, TTop, rhoAirMean, rhoAir, rhoTop, CO2Out, CO2Air)
  mcTopOut = MCTopOut(CO2Top, CO2Out, stigmaInsScr, Cd, Cw, vWind, URoof, USide, ARoof, ASide, hSideRoof, hRoof, TAir, TOut, TAirMean, AFlr, cLeakage, etaRoof, UThScr)
  V_CO2Top = (mcAirTop - mcTopOut) / capCO2Top
  return V_CO2Top


/*METHOD*/
def euler(x_euler, x0, y_euler, y0, h, n, x_e5, y_e5):
x_euler.append(x0)
y_euler.append(y0)
x_e5.append(x0)
y_e5.append(y0)
for i in range(1, n):
x_euler.append(x_euler[i-1] + h * dxAir(x_euler[i-1],y_euler[i-1],TAir[int((i-1) * h/300)], Tout[int((i-1) * h/300)], Iglob[int((i-1) * h/300)], vwind[int((i-1) * h/300)]))
y_euler.append(y_euler[i-1] + h * dxTop(x_euler[i-1],y_euler[i-1],TAir[int((i-1) * h/300)], Tout[int((i-1) * h/300)], Iglob[int((i-1) * h/300)], vwind[int((i-1) * h/300)]))
if (i * h) % 300 == 0:
    x_e5.append(x_euler[i-1] + h * dxAir(x_euler[i-1],y_euler[i-1],TAir[int((i-1) * h/300)], Tout[int((i-1) * h/300)], Iglob[int((i-1) * h/300)], vwind[int((i-1) * h/300)]))
    y_e5.append(y_euler[i-1] + h * dxTop(x_euler[i-1],y_euler[i-1],TAir[int((i-1) * h/300)], Tout[int((i-1) * h/300)], Iglob[int((i-1) * h/300)], vwind[int((i-1) * h/300)]))
return 
  
def rk4(x_rk4, x0, y_rk4, y0, h, n, x_r5, y_r5):
  x_rk4.append(x0)
  y_rk4.append(y0)
  x_r5.append(x0)
  y_r5.append(y0)
  for i in range(1,n - 1):
    k1 = h*dxAir(x_rk4[i-1], y_rk4[i-1], TAir[int((i-1) * h/300)], Tout[int((i-1) * h/300)], Iglob[int((i-1) * h/300)], vwind[int((i-1) * h/300)])
    d1 = h*dxTop(x_rk4[i-1], y_rk4[i-1], TAir[int((i-1) * h/300)], Tout[int((i-1) * h/300)], Iglob[int((i-1) * h/300)], vwind[int((i-1) * h/300)])
    k2 = h*dxAir(x_rk4[i-1] + k1 / 2, y_rk4[i-1] + d1 / 2, TAir[int((i-1) * h/300)], Tout[int((i-1) * h/300)], Iglob[int((i-1) * h/300)], vwind[int((i-1) * h/300)])
    d2 = h*dxTop(x_rk4[i-1] + k1 / 2, y_rk4[i-1] + d1 / 2, TAir[int((i-1) * h/300)], Tout[int((i-1) * h/300)], Iglob[int((i-1) * h/300)], vwind[int((i-1) * h/300)])
    k3 = h*dxAir(x_rk4[i-1] + k2 / 2, y_rk4[i-1] + d2 / 2, TAir[int((i-1) * h/300)], Tout[int((i-1) * h/300)], Iglob[int((i-1) * h/300)], vwind[int((i-1) * h/300)])
    d3 = h*dxTop(x_rk4[i-1] + k2 / 2, y_rk4[i-1] + d2 / 2, TAir[int((i-1) * h/300)], Tout[int((i-1) * h/300)], Iglob[int((i-1) * h/300)], vwind[int((i-1) * h/300)])
    k4 = h*dxAir(x_rk4[i-1] + k3, y_rk4[i-1] + d3, TAir[int((i-1) * h/300)], Tout[int((i-1) * h/300)], Iglob[int((i-1) * h/300)], vwind[int((i-1) * h/300)])
    d4 = h*dxTop(x_rk4[i-1] + k3, y_rk4[i-1] + d3, TAir[int((i-1) * h/300)], Tout[int((i-1) * h/300)], Iglob[int((i-1) * h/300)], vwind[int((i-1) * h/300)])
    x_rk4.append(x_rk4[i-1] + k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6)
    y_rk4.append(y_rk4[i-1] + d1 / 6 + d2 / 3 + d3 / 3 + d4 / 6)
    if (i * h) % 300 == 0:
      x_r5.append(x_rk4[i-1] + k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6)
      y_r5.append(y_rk4[i-1] + d1 / 6 + d2 / 3 + d3 / 3 + d4 / 6)
  return
  
/*RUN METHOD*/
x_euler = []
y_euler = []
x_e5 = []
y_e5 = []
euler(x_euler, CO2_Air[0], y_euler, CO2_Air[0], 10, 34560, x_e5, y_e5)

x_rk4 = []
y_rk4 = []
x_r5 = []
y_r5 = []
rk4(x_rk4, CO2_Air[0], y_rk4, CO2_Air[0], 10, 34560, x_r5, y_r5)
  
fig.add_trace(go.Scatter(y=CO2_Air[0:1150],mode='lines',name='CO2Air'))
fig.add_trace(go.Scatter(y=x_e5,mode='lines',name='Euler'))
fig.add_trace(go.Scatter(y=x_r5,mode='lines',name='RK4'))
fig.update_layout(title='Biểu đồ trực quan giữa kết quả CO2Air của hai giải thuật Euler và RK4 so với thực tế',
                   xaxis_title='Thời gian (5 mins)',
                   yaxis_title='CO2Air (mgm^(-3))')
fig.show()
fig.add_trace(go.Scatter(y=x_e5,mode='lines',name='Euler'))
fig.add_trace(go.Scatter(y=x_r5,mode='lines',name='RK4'))
fig.update_layout(title='Biểu đồ trực quan giữa kết quả CO2Air của hai giải thuật Euler và RK4',
                   xaxis_title='Thời gian (5 mins)',
                   yaxis_title='CO2Air (mgm^(-3))')
fig.show()

rmse = mean_squared_error(CO2_Air[0:1152], x_r5, squared = False)
rmse_ = mean_squared_error(CO2_Air[0:1152], x_e5, squared = False)
print("Root mean squared error của Euler: ", rmse_)
print("Root mean squared error của RK4: ", rmse)