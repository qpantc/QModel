## Core equations

# Ac = Vc * (ci - cp)/(ci+ kco)
# Aj = E * alpha_l * I * (ci - cp)/(ci+2*cp)

# cs = ca - An * 1.4 * 1000/gb
# gs = g1*An*RH*1000/cs
# ci = cs - An*1.6*1000/gs

# An = (Ca - ci*g1*RH*1000/(ci*g1*RH*1000 -1600))*gb/1400



import numpy as np
import math



# physical constant
class Constants:
    def __init__(self):
        self.R = 8.31*1e-3 # Gas constant, kJ/mol/K
        self.NA = 6.02e23 # Avogadro's constant, /mol
        self.hc = 2e-25 # Planck constant times light speed, J*s times m/s
        self.wavelen = 500e-9 # wavelength of light, m
        self.Ephoton = self.hc/self.wavelen
        self.ca = 400 # Atmospheric CO2 concentration, ppm
        self.Cpmol = 1005*28.97*1e-3 # J/kg/K*kg/mol -> J/mol/K
        self.lambda0 = 2.26*10**6 # heat of vaporization of water, J/kg
        self.gammaV = 100*1005/(self.lambda0*0.622) #in kPa/K, constant in PM equation
        self.a0 = 1.6 # relative diffusivity of h2o to co2 through stomata
        self.U3 = 273.15
        
        self.rgas = 8.3144598        # Perfect gaz constant R in J mol-1 K-1
        self.tfrz = 273.15           # Converstion temperature from K to C, in K      

# leaf parameters
class Leaf:
    def __init__(self):
        
        self.gb = 300             # Leaf boundary layer conductance to water vapour in mmol (H2O) m-2 s-1  
        self.g1 = 13.70     # Stomatal sensitivity coefficient (dimmensionless)
        
        # leaf parameters for carbon assimilation

        self.kc25 = 404.9                              # Michaelis-Menten constant  for CO2 at 25°C in µmol m-2 s-1
        self.ko25 = 278.4                              # Michaelis-Menten constant  for O2 at 25°C in mmol m-2 s-1
        self.cp25 = 42.75                              # Gamma* constant  at 25°C in µmol mol-1 µmol m-2 s-1
        
        # temperature response parameters
        self.kcha = 79430                              # Activation energy dHa for Kc in J mol-1
        self.koha = 36380                              # Activation energy dHa for Ko in J mol-1
        self.cpha = 37830                              # Activation energy dHa for Gamma* in J mol-1
        self.rdha = 46390                              # Activation energy dHa for Rd in J mol-1
        self.vcmaxha = 65330                           # Activation energy dHa for Vcmax in J mol-1
        self.jmaxha  = 43540                           # Activation energy dHa for Jmax in J mol-1
        
        self.hd = 150000                               # Desactivation energy dHd in J mol-1
  
        self.se = 490                                  # Entropy term dS in J mol-1 K-1
        
        # Light response parameters
        self.phi_psii = 0.85                           # Quantuum yield of PSII in µmol m-2 s-1
        self.theta_j = 0.90                            # Curvature of the light response curve (0-1)
        self.alpha_l = 0.88                            # Leaf absorptance (0-1)
        
        # stomatal conductance
        self.gmin=76.20                                # Minimum stomatal conductance to water vapour in umol m-2 s-1 

        self.Ifac=0.5                                  # Coefficient representing the extent to which Rdark is inhibited in the light (dimmensionless)

def co2_to_ca(co2, patm):
    ca = (1.0e-3) * co2 * patm * 10
    return ca

def ft(tl, ha, physcon = Constants()):
    return math.exp(ha / (physcon.rgas * (physcon.tfrz + 25)) * (1 - (physcon.tfrz + 25) / tl))

def fth(tl, hd, se, physcon = Constants()):
    fc25 = 1 + math.exp((-hd + se * (physcon.tfrz + 25)) / (physcon.rgas * (physcon.tfrz + 25)))
    fc = 1 + math.exp((-hd + se * tl) / (physcon.rgas * tl))
    fh = fc25 / fc
    return fh


def saturation_vapor_pressure(T):
    """计算温度下的饱和水蒸气压(kPa)"""
    A = 6.112  # hPa
    B = 17.67
    C = 243.5  # °C
    return A * np.exp(B * T / (C + T)) / 10  # 将结果转换为 kPa

def actual_water_vapor_pressure(H2O_concentration, P, constant):
    """计算实际水蒸气压(kPa)"""
    return H2O_concentration * constant * P / 1000  # 将水汽浓度转换为 kPa

def relative_humidity(T, H2O_concentration, P):
    """计算相对湿度(%)"""
    P_actual = actual_water_vapor_pressure(H2O_concentration, P, 0.018015)  # 水的摩尔质量为 18.015 g/mol,将其转换为 kPa
    P_sat = saturation_vapor_pressure(T)
    RH = (P_actual / P_sat) * 100
    return RH

def vapor_pressure_deficit(T, H2O_concentration, P):
    """计算蒸汽压差(VPD,kPa)"""
    P_actual = actual_water_vapor_pressure(H2O_concentration, P, 0.0987)
    P_sat = saturation_vapor_pressure(T)
    VPD = P_sat - P_actual
    return VPD


def Q_model(leaf = Leaf(),tleaf= 27,co2 = 400,RH = 0.71,PAR = 1000,patm = 101.82,cost = [0.11,0.8,0.00001]):
    
    # Vc_curve = fth (tl = tleaf + 273.15 , hd = leaf.vcmaxha, se = leaf.se) * ft(tl = tleaf + 273.15, ha = leaf.vcmaxha)
    a = cost[0]*(1+(abs(tleaf-30)/100)**0.5) # cost for unit vcmas
    b = cost[1]*1000/PAR# cost for unit E
    c = cost[2]/RH  # cost for unit gs
    
    
    
    # ENV
    Iabs= PAR*0.5
    Oi = 209 # we fic Oi to atmospheric concentrations

    constant = Constants()
    

    Ca=co2_to_ca(co2 = co2, patm = patm)
    
    # We compute the adjusted parameter values for leaf temperature
    # for Kc
    kc = leaf.kc25 #* ft(tleaf + 273.15 , leaf.kcha)
    # for Ko
    ko = leaf.ko25 #* ft(tleaf+ 273.15 , leaf.koha)
    # for Gamma*
    cp = leaf.cp25 #* ft(tleaf+ 273.15 , leaf.cpha)
    # cp = min(60,cp)

    kco = kc*(1+Oi/ko)


    # An = (Ca - ci*g1*RH*1000/(g1*RH*1000 -1600))*gb/1400
    p = Ca*leaf.gb/1400
    q = leaf.g1*leaf.gb*RH*1000/(leaf.g1*RH*1000 -1600)/1400
    
    M = c*leaf.g1*Ca*(leaf.gb**2) * RH * 1000 * q


    A2 = (a*q + b*q/(leaf.alpha_l*Iabs) - q)
    A1 = 2*q*cp + a*q*kco - a*q/(leaf.alpha_l*Iabs)*cp + b*q*cp - a*q*cp - a*q*kco - 3*b*q/(leaf.alpha_l*Iabs)*cp
    A0 = q*cp**2 -a*q*kco*cp - 2*b*q/(leaf.alpha_l*Iabs)*cp**2 + p*(a*(cp+kco)+3*b/(leaf.alpha_l*Iabs)*cp)

    N1 = 1400*q
    N0 = (Ca*leaf.gb - 1400*p)


    p_poly = (A1*N1 + A2*N0)/(A2*N1)
    q_poly = (A0*N1 + A1*N0)/(A2*N1)
    r_poly = (A0*N0 + M)/(A2*N1)

    Q_poly = (p_poly**2 - 3*q_poly)/9
    U_poly = (2*p_poly**3-9*p_poly*q_poly+27*r_poly)/54

    if Q_poly>0 and (U_poly/Q_poly**(3/2) > -1 and U_poly/Q_poly**(3/2) < 1):
        Psii_poly = math.acos(U_poly/Q_poly**(3/2))


        ci1 = 0-2*Q_poly**0.5*np.cos(Psii_poly/3) - p_poly/3
        ci2 = 0-2*Q_poly**0.5*np.cos((Psii_poly + 2*3.141592653)/3) - p_poly/3
        ci3 = 0-2*Q_poly**0.5*np.cos((Psii_poly + 4*3.141592653)/3) - p_poly/3
        
        
        An = p - q * ci2
        Vcmax = An/(ci2-cp)*(ci2+kco)
        E = An/(ci2-cp)*(ci2+2*cp)/(leaf.alpha_l*Iabs)
        cs = Ca - An*1400/leaf.gb
        gs = An*1600/(cs-ci2)
        
        An = An - Vcmax*a - E * b - gs * c

        if(An > 0 and Vcmax > 0):

            return({
                'An': An,
                'cp':cp,
                'Ci': ci2,
                'Vc':Vcmax,
                'E':E,
                'gs': gs
            })
        else:
            return({
                'An': 0,
                'cp':cp,
                'Ci': Ca,
                'Vc':0,
                'E':0,
                'gs': 0
            })
    else:
        # print('parameters error')
        return({
                'An': 0,
                'cp':cp,
                'Ci': Ca,
                'Vc':0,
                'E':0,
                'gs': 0
            })
