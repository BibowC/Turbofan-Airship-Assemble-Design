# -*- coding: utf-8 -*-
"""
Created on Sun Dec 18 20:55:15 2022

@author: Carlos Eduardo Bibow Corrêa
Não autorizo o uso desse código para provas e trabalhos. 
"""
import numpy as np
import sympy as smp
import matplotlib.pyplot as plt
import random 
#calculo de turbofan

K = 273.15 #temp kelvin
Mo = np.linspace(0.1, 1, 10) #mach de voo

pif = 1.1 #taxa de compressão LPC
pic = 32  #taxa de compressão HPC
T04 = 1800 #Temp maxima camara combustão
PCI = 42E6 #poder calorifico
cp = 1004 
gamma = 1.4
ma = 115 # massa de ar
Po = 40.E3 #pressão na entrada
To = -30 + K #Temperatura em K na entrada 
B = 4 #Bypass
R = 287 #constante R
PSL= 101325 #Pressão em nível do mar

#-----------------------
#Intake

P01 = Po * (1+((gamma-1)/2)*Mo**2)**(gamma/(gamma-1))
T01 = To * (1+((gamma-1)/2)*Mo**2)

#-----------------------
#fan LPC

P02f = P01
T02f = T01

P03f = pif * P02f
tauf = pif**((gamma-1)/gamma)
T03f = tauf * T02f

#------------------------
#compressor HPC

P03 = P03f * pic
tauc = pic**((gamma-1)/gamma)
T03 = T03f * tauc

#-----------------------
#câmara de combustão
P04 = P03
f = ((T04/T03)-1)/((PCI/(cp*T03))-(T04/T03))

#-----------------------
#Conectar HPT-HPC
T05 = T04 - ((T03 - T03f)/(1+f))
P05 = P04 * ((T05/T04)**((gamma)/(gamma-1)))
#Conectar LPT-LPC
T06 = T05 - ((1+B)*(T03f-T02f))/(1+f)
P06 = P05 * (T06/T05)**((gamma)/(gamma-1))

#-------------------------
#Bocal 3' - 7'
P07f = P03f

NPRf = P07f/Po
Pef = P07f/1.892
# Bocal 6-7
P07 = P06

NPR = P07/Po
Pe = P07/1.892
#------------------------
Tef = T03f/(1+((gamma-1)/2))
Te = T06/(1+((gamma-1)/2))

Uef = (gamma*R*(Tef))**(1/2)
Ue = (gamma*R*(Te))**(1/2)
Uo = Mo*((gamma*R*To)**(1/2))

rhof = Pef / (R*T03f)
rho = Pe / (R*T06)
#--------------------------------------------
#Cálculo das vazões mássicas e área de saída dos bocais

a = np.array([[5, -1], [1, 1]])
b = np.array([0, 115])
x = np.linalg.solve(a, b)

mac = x[1]
mah = x[0]

Aef = mac/(rhof*Uef)
Ae = mah/(rho*Ue)

#----------------------------------
#Cálculo do impulso

fnh = (mah*((1+f)*(Ue-Uo))) + ((Pe - Po)*Ae)
fnc = (B * mac * (Uef-Uo)) + ((Pef - Po)*Aef )
FN = 2*(fnc+fnh)

plt.figure(1)
plt.legend('Impulso x Mach')
plt.xlabel('Mo')
plt.ylabel('FN')
plt.plot(Mo, FN, 'r')

#--------#--------#-------#----------------------------
#ACOPLAMENTO MOTOR-AERONAVE


Wcf = 2000 #Peso de combustível final desejado
delta1 = Po/101325 #Pressão adimensional
Wempty = 27000 # Peso vazio sem motor
Clmax = 2.1 #CLMAX Admitido
Aw = 116.3 #Área das asas
rhocomb = 0.81 #Densidade combustível (Querosene)
qtdcomb = 20000 #Litros
Wci = rhocomb * qtdcomb #Peso inicial combustível
Wmotor = 2366 #Peso motor
Wpas = 146*70 #peso total de passageiros com peso médio de 70kg
W = (27000 + 2*(Wmotor) + Wci + Wpas)*9.81 #Peso total inicial
Cl = (2*W) / (delta1*101325*(Mo**2)*Aw) #coeficiente de sustentação
Cd =  0.0181 + 0.0362*((Cl)**2) #Polar de arrasto


FNdelta1 = FN/delta1

WW = np.linspace(Wci, Wcf, 10) #Variação do peso inicial até o peso final
newcl = (2*WW) / (delta1*gamma*101325*(Mo**2)*Aw)
cd2 = 0.0181 + 0.0362*((newcl)**2) #Polar de arrasto
FNdelta2 = (gamma/2) * 101325 * cd2 * (Mo**2) * Aw
plt.figure(2)
plt.plot(Mo, FNdelta1, 'r')
plt.plot(Mo, FNdelta2, 'b')



rhoarSL = (101325/(R*298)) #Densidade do ar em SL
Vstall = (W/(Clmax*(1/2)*rhoarSL*Aw))**(1/2) #Velocidade de stall em SL
Vlo = (1.2)*Vstall #Velocidade de liftoff simplificada
FNDEC = 0.707*Vlo #FN Decolagem simplificada
amed = ((FNDEC - (0.02*W))/(W)) * 9.81 #Aceleração média
machfnedc = FNDEC/347
SG = (Vlo**2)/(2*amed)

#Consumo e tempo para decolagem
mf = f[2] * ma
TSFC = mf/FN[2]
tdec = Vlo / amed
consumodec = TSFC*FN[2]*tdec

#-----------------------------
#Atualizando dados de peso/ Subida + Descida

W2 = W-consumodec
D = Cd * 0.5*rhoarSL*Vlo*Aw
dhdt = ((FN[2]-D[2])*Vlo)/W
tsubida = (11000/dhdt) # tempo de subida ()
consumosubida = TSFC*FN[2]*tsubida

tdescida = tsubida
consumodescida = consumosubida

#-------------------------------
#Voo em cruzeiro
W3 = Wci - (consumodec+consumosubida)
rhocruz = Po / (R*To)
Vcruz = Mo * (gamma*R*To)**(1/2)
S = ((2/(9.81*TSFC)) * ((2/(Aw*rhocruz))**(1/2)) * (Cl[1] / Cd[1]**2) * ((W2**(1/2)) - (W3**(1/2))))
tcruzeiro = S / Vcruz
consumocruz = TSFC*FN[9] * tcruzeiro




