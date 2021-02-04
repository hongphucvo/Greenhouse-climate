/*FORMULAR*/
def VECCanAir_(rhoAir, cPAir, LAI, deltaH, gamma, rb, rs):
  return (2 * rhoAir * cPAir * LAI) / (deltaH * gamma * (rb + rs))
  
def MVCanAir(rhoAir, cPAir, LAI, deltaH, gamma, VPCan, VPAir):
  rs = rs_min
  VECCanAir = VECCanAir_(rhoAir, cPAir, LAI, deltaH, gamma, rb, rs)
  return VECCanAir * (VPCan - VPAir)

def MVBlowAir(etaHeatVap, UBlow, PBlow):
  # HBlowAir = UBlow * PBlow / AFlr
  return etaHeatVap * UBlow * PBlow / AFlr 

def MVPadAir(rhoAir, UPad, phiPad, AFlr, etaPad, xPad, xOut):
  # fPad = phiPad * UPad / AFlr
  return rhoAir * phiPad * UPad * (etaPad * (xPad - xOut) + xOut) / AFlr

def MVAirOut_Pad(UPad, phiPad, AFlr, VPAir, TAir):
  return M_water * VPAir * phiPad * UPad / (R * (TAir + 273.15) * AFlr)

def MVFogAir(UFog, phiFog, AFlr):
  return UFog * phiFog / AFlr

def HECAirThScr_(UThScr, TAir, TThScr):
  return 1.7 * UThScr * pow(abs(TAir - TThScr), 0.33)

def MVAirThScr(VPAir, VPThScr, UThScr, TAir, TThScr):
  HECAirThScr = HECAirThScr_(UThScr, TAir, TThScr)
  if VPAir <= VPThScr:
    return 0
  else:
    return 6.4 * 10 ** (-9) * HECAirThScr * (VPAir - VPThScr)

def HECTopCov_in_(cHECin, ACov, AFlr, TTop, TCov_in): 
  return cHECin * abs(TTop - TCov_in) ** 0.33 * ACov / AFlr

def MVTopCov_in(VPTop, VPCov_in, cHECin, ACov, AFlr, TTop, TCov_in):
  HECTopCov_in = HECTopCov_in_(cHECin, ACov, AFlr, TTop, TCov_in)
  if VPTop <= VPCov_in:
    return 0
  else:
    return 6.4 * 10 ** (-9) * HECTopCov_in * (VPTop - VPCov_in)

def HECMechAir_(UMechCool, COPMechCool, PMechCool, TAir, TMechCool, deltaH, VPAir, VPMechCool, AFlr):
  a = UMechCool * COPMechCool * PMechCool / AFlr
  b = TAir - TMechCool + 6.4 * 10 ** (-9) * deltaH * (VPAir - VPMechCool)
  return a / b

def MVAirMech(UMechCool, COPMechCool, PMechCool, TAir, TMechCool, deltaH, VPAir, VPMechCool, AFlr):
  HECMechAir = HECMechAir_(UMechCool, COPMechCool, PMechCool, TAir, TMechCool, deltaH, VPAir, VPMechCool, AFlr)
  if VPAir <= VPMechCool:
    return 0
  else:
    return 6.4 * 10 ** (-9) * HECMechAir * (VPAir - VPMechCool)

def MVAirTop(UThScr, KThScr, TAir, TTop, rhoAirMean, rhoAir, rhoTop, VPAir, VPTop):
  return M_water * fThScr(UThScr, KThScr, TAir, TTop, rhoAirMean, rhoAir, rhoTop) * (VPAir / (TAir + 273.15) - VPTop / (TTop + 273.15)) / R

def MVAirOut(Cd, Cw, vWind, URoof, USide, ARoof, ASide, hSideRoof, TAir, TOut, TAirMean, AFlr, etaSide, stigmaInsScr, UThScr, cLeakage, UVentForced, phiVentForced, VPAir, VPOut):
  etaInsScr = stigmaInsScr * (2 - stigmaInsScr)
  fVentSide_0 = (Cd * USide * ASide * vWind * (Cw ** (1/2))) / (2*AFlr)
  fVentSide_ARoof = fVentRoofSide(Cd, Cw, vWind, URoof, USide, ARoof, ASide, hSideRoof, TAir, TOut, TAirMean, AFlr)
  ffVentSide = fVentSideAO(etaSide, etaInsScr, UThScr, fVentSide_ARoof, fVentSide_0, vWind, cLeakage)
  ffVentForced = fVentForcedAO(etaInsScr, UVentForced, phiVentForced, AFlr)
  return M_water * (ffVentForced + ffVentSide) * (VPAir / (TAir + 273.15) - VPOut / (TOut + 273.15)) / R

def MVTopOut(Cd, Cw, vWind, URoof, USide, ARoof, ASide, hSideRoof, TAir, TOut, TAirMean, AFlr, hRoof, cLeakage, UThScr, etaRoof, VPTop, VPOut, TTop):
  etaInsScr = stigmaInsScr * (2 - stigmaInsScr)
  ffVentRoofSide  = fVentRoofSide(Cd, Cw, vWind, URoof, USide, ARoof, ASide, hSideRoof, TAir, TOut, TAirMean, AFlr)
  ffVentRoof_1 = fVentRoof_TO(Cd, URoof, ARoof, hRoof, TAir, TOut, TAirMean, Cw, vWind, AFlr) 
  ffVentRoof_2 = (Cd * USide * ASide * vWind * (Cw ** (1/2))) / (2*AFlr)
  ffVentRoof = fVentRoofTO(etaRoof, etaInsScr, UThScr, ffVentRoofSide, ffVentRoof_1, ffVentRoof_2, vWind, cLeakage)
  return M_water * ffVentRoof * (VPTop / (TTop + 273.15) - VPOut / (TOut + 273.15)) / R

def dxVP(VPAir, VPTop, TAir, CO2Air, TOut, IGlob, vWind, VPOut):
  TTop = TAir - 1
  TCan = TAir - 1
  TMechCool = TAir - 1
  TThScr = TAir + 1
  TCov_in = TAir - 1
  TAirMean = (TAir + TOut) / 2
  capVPAir = M_water * hAir / (R * (TAir + 273.15))
  capVPTop = M_water * hAir / (R * (TTop + 273.15))
  RCan = (1 - etaGlobAir) * IGlob * (etaGlobPar * tauCovPar + etaGlobNir * tauCovNir)
  VPCan = 610.78 * pow(np.exp(1), (TCan * 17.2694) / (TCan + 238.3))
  VPThScr = 610.78 * pow(np.exp(1), (TThScr * 17.2694) / (TThScr + 238.3))
  VPMechCool = 610.78 * pow(np.exp(1), (TMechCool * 17.2694) / (TMechCool + 238.3))
  VPCov_in = 610.78 * pow(np.exp(1), (TCov_in * 17.2694) / (TCov_in + 238.3))

  mvCanAir = MVCanAir(rhoAir, cPAir, LAI, deltaH, gamma, VPCan, VPAir)
  mvPadAir = MVPadAir(rhoAir, UPad, phiPad, AFlr, etaPad, xPad, xOut)
  mvFogAir = MVFogAir(UFog, phiFog, AFlr)
  mvBlowAir = MVBlowAir(etaHeatVap, UBlow, PBlow)
  mvAirThScr = MVAirThScr(VPAir, VPThScr, UThScr, TAir, TThScr)
  mvAirTop = MVAirTop(UThScr, KThScr, TAir, TTop, rhoAirMean, rhoAir, rhoTop, VPAir, VPTop)
  mvAirOut = MVAirOut(Cd, Cw, vWind, URoof, USide, ARoof, ASide, hSideRoof, TAir, TOut, TAirMean, AFlr, etaSide, stigmaInsScr, UThScr, cLeakage, UVentForced, phiVentForced, VPAir, VPOut)
  mvAirOut_Pad = MVAirOut_Pad(UPad, phiPad, AFlr, VPAir, TAir)
  mvAirMech = MVAirMech(UMechCool, COPMechCool, PMechCool, TAir, TMechCool, deltaH, VPAir, VPMechCool, AFlr)
  mvTopCov_in = MVTopCov_in(VPTop, VPCov_in, cHECin, ACov, AFlr, TTop, TCov_in)
  mvTopOut = MVTopOut(Cd, Cw, vWind, URoof, USide, ARoof, ASide, hSideRoof, TAir, TOut, TAirMean, AFlr, hRoof, cLeakage, UThScr, etaRoof, VPTop, VPOut, TTop)
  V_VPAir = (mvCanAir + mvPadAir + mvFogAir + mvBlowAir - mvAirThScr - mvAirTop - mvAirOut - mvAirOut_Pad + mvAirMech) / capVPAir
  V_VPTop = (mvAirTop - mvTopCov_in - mvTopOut) / capVPTop
  return V_VPAir, V_VPTop

def dxVPAir(VPAir, VPTop, TAir, CO2Air, TOut, IGlob, vWind, VPOut):
  TTop = TAir - 1
  TCan = TAir - 1
  TMechCool = TAir - 1
  TThScr = TAir + 1
  TCov_in = TAir - 1
  TAirMean = (TAir + TOut) / 2 
  capVPAir = M_water * hAir / (R * (TAir + 273.15))
  capVPTop = M_water * hAir / (R * (TTop + 273.15))
  RCan = (1 - etaGlobAir) * IGlob * (etaGlobPar * tauCovPar + etaGlobNir * tauCovNir)
  VPCan = 610.78 * pow(np.exp(1), (TCan * 17.2694) / (TCan + 238.3))
  VPThScr = 610.78 * pow(np.exp(1), (TThScr * 17.2694) / (TThScr + 238.3))
  VPMechCool = 610.78 * pow(np.exp(1), (TMechCool * 17.2694) / (TMechCool + 238.3))
  VPCov_in = 610.78 * pow(np.exp(1), (TCov_in * 17.2694) / (TCov_in + 238.3))

  mvCanAir = MVCanAir(rhoAir, cPAir, LAI, deltaH, gamma, VPCan, VPAir)#
  mvPadAir = MVPadAir(rhoAir, UPad, phiPad, AFlr, etaPad, xPad, xOut)
  mvFogAir = MVFogAir(UFog, phiFog, AFlr)
  mvBlowAir = MVBlowAir(etaHeatVap, UBlow, PBlow)#
  mvAirThScr = MVAirThScr(VPAir, VPThScr, UThScr, TAir, TThScr)
  mvAirTop = MVAirTop(UThScr, KThScr, TAir, TTop, rhoAirMean, rhoAir, rhoTop, VPAir, VPTop)
  mvAirOut = MVAirOut(Cd, Cw, vWind, URoof, USide, ARoof, ASide, hSideRoof, TAir, TOut, TAirMean, AFlr, etaSide, stigmaInsScr, UThScr, cLeakage, UVentForced, phiVentForced, VPAir, VPOut)#
  mvAirOut_Pad = MVAirOut_Pad(UPad, phiPad, AFlr, VPAir, TAir)
  mvAirMech = MVAirMech(UMechCool, COPMechCool, PMechCool, TAir, TMechCool, deltaH, VPAir, VPMechCool, AFlr)
  V_VPAir = (mvCanAir + mvPadAir + mvFogAir + mvBlowAir - mvAirThScr - mvAirTop - mvAirOut - mvAirOut_Pad + mvAirMech) / capVPAir
  return V_VPAir

def dxVPTop(VPAir, VPTop, TAir, CO2Air, TOut, IGlob, vWind, VPOut):
  TTop = TAir - 1
  TCan = TAir - 1
  TMechCool = TAir - 1
  TThScr = TAir + 1
  TCov_in = TAir - 1
  TAirMean = (TAir + TOut) / 2 
  capVPAir = M_water * hAir / (R * (TAir + 273.15))
  capVPTop = M_water * hAir / (R * (TTop + 273.15))
  RCan = (1 - etaGlobAir) * IGlob * (etaGlobPar * tauCovPar + etaGlobNir * tauCovNir)
  VPCan = 610.78 * pow(np.exp(1), (TCan * 17.2694) / (TCan + 238.3))
  VPThScr = 610.78 * pow(np.exp(1), (TThScr * 17.2694) / (TThScr + 238.3))
  VPMechCool = 610.78 * pow(np.exp(1), (TMechCool * 17.2694) / (TMechCool + 238.3))
  VPCov_in = 610.78 * pow(np.exp(1), (TCov_in * 17.2694) / (TCov_in + 238.3))

  mvAirTop = MVAirTop(UThScr, KThScr, TAir, TTop, rhoAirMean, rhoAir, rhoTop, VPAir, VPTop)
  mvTopCov_in = MVTopCov_in(VPTop, VPCov_in, cHECin, ACov, AFlr, TTop, TCov_in)
  mvTopOut = MVTopOut(Cd, Cw, vWind, URoof, USide, ARoof, ASide, hSideRoof, TAir, TOut, TAirMean, AFlr, hRoof, cLeakage, UThScr, etaRoof, VPTop, VPOut, TTop)
  V_VPTop = (mvAirTop - mvTopCov_in - mvTopOut) / capVPTop
  return V_VPTop






/*METHOD*/
def eulerMV(x_euler, x0, y_euler, y0, h, n, x_e5, y_e5):
  x_euler.append(x0)
  y_euler.append(y0)
  x_e5.append(x0)
  y_e5.append(y0)
  for i in range(1, n):
    x_euler.append(x_euler[i-1] + h * dxVPAir(x_euler[i-1],y_euler[i-1],TAir[int((i-1) * h/300)], CO2_Air[int((i-1) * h/300)], Tout[int((i-1) * h/300)], Iglob[int((i-1) * h/300)], vwind[int((i-1) * h/300)], VP_Out[int((i-1) * h/300)]))
    y_euler.append(y_euler[i-1] + h * dxVPTop(x_euler[i-1],y_euler[i-1],TAir[int((i-1) * h/300)], CO2_Air[int((i-1) * h/300)], Tout[int((i-1) * h/300)], Iglob[int((i-1) * h/300)], vwind[int((i-1) * h/300)], VP_Out[int((i-1) * h/300)]))
    if (i * h) % 300 == 0:
      x_e5.append(x_euler[i-1] + h * dxVPAir(x_euler[i-1],y_euler[i-1],TAir[int((i-1) * h/300)], CO2_Air[int((i-1) * h/300)], Tout[int((i-1) * h/300)], Iglob[int((i-1) * h/300)], vwind[int((i-1) * h/300)], VP_Out[int((i-1) * h/300)]))
      y_e5.append(y_euler[i-1] + h * dxVPTop(x_euler[i-1],y_euler[i-1],TAir[int((i-1) * h/300)], CO2_Air[int((i-1) * h/300)], Tout[int((i-1) * h/300)], Iglob[int((i-1) * h/300)], vwind[int((i-1) * h/300)], VP_Out[int((i-1) * h/300)]))
  return

def rk4MV(x_rk4, x0, y_rk4, y0, h, n, x_r5, y_r5):
  x_rk4.append(x0)
  y_rk4.append(y0)
  x_r5.append(x0)
  y_r5.append(y0)
  for i in range(1,n - 1):
    k1 = h*dxVPAir(x_rk4[i-1], y_rk4[i-1], TAir[int((i-1) * h/300)], CO2_Air[int((i-1) * h/300)], Tout[int((i-1) * h/300)], Iglob[int((i-1) * h/300)], vwind[int((i-1) * h/300)], VP_Out[int((i-1) * h/300)])
    d1 = h*dxVPTop(x_rk4[i-1], y_rk4[i-1], TAir[int((i-1) * h/300)], CO2_Air[int((i-1)/300)], Tout[int((i-1) * h/300)], Iglob[int((i-1) * h/300)], vwind[int((i-1) * h/300)], VP_Out[int((i-1) * h/300)])
    k2 = h*dxVPAir(x_rk4[i-1] + k1 / 2, y_rk4[i-1] + d1 / 2, TAir[int((i-1) * h/300)], CO2_Air[int((i-1) * h/300)], Tout[int((i-1) * h/300)], Iglob[int((i-1) * h/300)], vwind[int((i-1) * h/300)], VP_Out[int((i-1) * h/300)])
    d2 = h*dxVPTop(x_rk4[i-1] + k1 / 2, y_rk4[i-1] + d1 / 2, TAir[int((i-1) * h/300)], CO2_Air[int((i-1) * h/300)], Tout[int((i-1) * h/300)], Iglob[int((i-1) * h/300)], vwind[int((i-1) * h/300)], VP_Out[int((i-1) * h/300)])
    k3 = h*dxVPAir(x_rk4[i-1] + k2 / 2, y_rk4[i-1] + d2 / 2, TAir[int((i-1) * h/300)], CO2_Air[int((i-1) * h/300)], Tout[int((i-1) * h/300)], Iglob[int((i-1) * h/300)], vwind[int((i-1) * h/300)], VP_Out[int((i-1) * h/300)])
    d3 = h*dxVPTop(x_rk4[i-1] + k2 / 2, y_rk4[i-1] + d2 / 2, TAir[int((i-1) * h/300)], CO2_Air[int((i-1)/300)], Tout[int((i-1) * h/300)], Iglob[int((i-1) * h/300)], vwind[int((i-1) * h/300)], VP_Out[int((i-1) * h/300)])
    k4 = h*dxVPAir(x_rk4[i-1] + k3, y_rk4[i-1] + d3, TAir[int((i-1) * h/300)], CO2_Air[int((i-1) * h/300)], Tout[int((i-1) * h/300)], Iglob[int((i-1) * h/300)], vwind[int((i-1) * h/300)], VP_Out[int((i-1) * h/300)])
    d4 = h*dxVPTop(x_rk4[i-1] + k3, y_rk4[i-1] + d3, TAir[int((i-1) * h/300)], CO2_Air[int((i-1) * h/300)], Tout[int((i-1) * h/300)], Iglob[int((i-1) * h/300)], vwind[int((i-1) * h/300)], VP_Out[int((i-1) * h/300)])
    x_rk4.append(x_rk4[i-1] + k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6)
    y_rk4.append(y_rk4[i-1] + d1 / 6 + d2 / 3 + d3 / 3 + d4 / 6)
    if (i * h) % 300 == 0:
      x_r5.append(x_rk4[i-1] + k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6)
      y_r5.append(y_rk4[i-1] + d1 / 6 + d2 / 3 + d3 / 3 + d4 / 6)
  return
  

/*RUN METHOD*/
dxVP(VP_Air[0], VP_Air[0], TAir[0], CO2_Air[0], Tout[0], Iglob[0], vwind[0], VP_Out[0])
x_euler = []
y_euler = []
x_e5 = []
y_e5 = []
eulerMV(x_euler, VP_Air[0], y_euler, VP_Air[0], 50, 6912, x_e5, y_e5)
x_rk4 = []
y_rk4 = []
x_r5 = []
y_r5 = []
rk4MV(x_rk4, VP_Air[0], y_rk4, VP_Air[0], 50, 6912, x_r5, y_r5)

/*COMPARE*/
fig.add_trace(go.Scatter(y=VP_Air[0:1152],
                    mode='lines',
                    name='VPAir'))
fig.add_trace(go.Scatter(y=x_e5,
                    mode='lines',
                    name='Euler'))
fig.add_trace(go.Scatter(y=x_r5,
                    mode='lines',
                    name='RK4'))
fig.update_layout(title='Biểu đồ trực quan giữa kết quả VPAir của hai giải thuật Euler và RK4 so với thực tế',
                   xaxis_title='Thời gian (5 mins)',
                   yaxis_title='VPAir (Pa)')
fig.show()

fig.add_trace(go.Scatter(y=x_e5,
                    mode='lines',
                    name='Euler'))
fig.add_trace(go.Scatter(y=x_r5,
                    mode='lines',
                    name='RK4'))
fig.update_layout(title='Biểu đồ trực quan giữa kết quả VPAir của hai giải thuật Euler và RK4',
                   xaxis_title='Thời gian (5 mins)',
                   yaxis_title='VPAir (Pa)')
fig.show()

rmse = mean_squared_error(VP_Air[0:1152], x_r5, squared = False)
rmse_ = mean_squared_error(VP_Air[0:1152], x_e5, squared = False)
print("Root mean squared error của Euler: ", rmse_)
print("Root mean squared error của RK4: ", rmse)