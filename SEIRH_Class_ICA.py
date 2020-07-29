#!/usr/bin/python
# -*- coding:utf8 -*-

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import odeint
BBDD_Dir='BBDD/'
"""
En este script se generaliza la implementación del modelo desarrollado
por la Secretaría de Salud al uso de clases.

update:
	26/07/2020: se agrega compartimento adicional con nivel de aislamiento sobre
				IC antes de ir a hospitalización (ICA), se define sigma_CA
				tal que 1/sigma_C + 1/sigma_CA = 7.1
				
	23/07/2020: sigma_C=7.1, escenario R0=3 A=0.6
				Índice de Movilidad hasta 19 de mayo (día 80)
				Guardar vectores con Info Usando Pandas
"""

########################################################################
#################            SEIRHT CLASS              #################
########################################################################
class SEIRH():
	def __init__(self):
		self.dummy=True
	
	def var_trans(self,beta_0,beta_1,beta_H,r,A,Movilidad):
		self.beta_0=beta_0
		self.beta_1=beta_1
		self.beta_H=beta_H
		self.r=r
		self.A=A
		self.Movilidad=Movilidad
	
	def beta(self,t):
		if t<92:
			if self.A==0:
				At=self.A
				return (1-At)*self.beta_0 + At*self.beta_1		
			else:
				At=self.Movilidad[int(t)]
				return (1-At)*self.beta_0 + At*self.beta_1	
		if t>=92:
			At=self.A
			return (1-At)*self.beta_0 + At*self.beta_1		
	
	def var_t_estadia(self, omega, gamma_M, sigma_C, sigma_CA, gamma_HR, 
		nu, gamma_R, sigma_HD, sigma_UD):
		self.omega=omega		#T. prom. latencia
		self.gamma_M=gamma_M	#T. prom de recuperación para IM
		self.sigma_C=sigma_C	#T. inf. hosp
		self.sigma_CA=sigma_CA	#T. inf. hosp
		self.gamma_HR=gamma_HR	#T. prom HR->R
		self.nu=nu				#T. prom HU->IR
		self.gamma_R=gamma_R	#T. prom IR->R
		self.sigma_HD=sigma_HD	#T. prom HD->D
		self.sigma_UD=sigma_UD	#T. prom UD->D
	
	def var_H_UCI(self, delta_M, delta_HR, delta_HD, delta_UR):
		self.delta_M=delta_M
		self.delta_HR=delta_HR
		self.delta_HD=delta_HD
		self.delta_UR=delta_UR
		self.delta_UD=1-delta_HR-delta_HD-delta_UR
	
	def var_ini(self, N0, E0=0, IM0=0, IC0=0, ICA0=0,
		IHR0=0, IUR0=0, IHD0=0, IUD0=0, IR0=0, R0=0, D0=0 ):
		self.N0=N0
		self.E0=E0
		self.IM0=IM0
		self.IC0=IC0
		self.ICA0=ICA0
		self.IHR0=IHR0
		self.IUR0=IUR0
		self.IHD0=IHD0
		self.IUD0=IUD0
		self.IR0=IR0
		self.R0=R0
		self.D0=D0
		self.S0= (self.N0 - self.E0 - self.IM0 - self.IC0 - self.ICA0 
				  - self.IHR0 - self.IUR0 - self.IHD0 - self.IUD0 
				  - self.IR0 - self.R0 - self.D0)
	
	def ODES(self,y,t):
		S, E, IM, IC, ICA, IHR, IUR, IHD, IUD, IR, R, D, N = y
		
		dSdt = -S/float(N)*(self.beta(t)*(IM+IC+self.r*ICA)+self.beta_H*(IHR+IHD+IUD+IUR+IR))
		dEdt =  S/float(N)*(self.beta(t)*(IM+IC+self.r*ICA)+self.beta_H*(IHR+IHD+IUD+IUR+IR)) - self.omega*E
		dIMdt = self.delta_M*self.omega*E - self.gamma_M*IM
		dICdt = (1-self.delta_M)*self.omega*E - self.sigma_C*IC
		dICAdt = self.sigma_C*IC - self.sigma_CA*ICA
		dIHRdt = self.delta_HR*self.sigma_CA*ICA - self.gamma_HR*IHR
		dIURdt = self.delta_UR*self.sigma_CA*ICA - self.nu*IUR
		dIHDdt = self.delta_HD*self.sigma_CA*ICA - self.sigma_HD*IHD
		dIUDdt = self.delta_UD*self.sigma_CA*ICA - self.sigma_UD*IUD
		dIRdt = self.nu*IUR - self.gamma_R*IR
		dRdt = self.gamma_HR*IHR + self.gamma_R*IR + self.gamma_M*IM
		dDdt = self.sigma_HD*IHD + self.sigma_UD*IUD
		dNdt = -self.sigma_HD*IHD - self.sigma_UD*IUD
		
		return [dSdt, dEdt, dIMdt, dICdt, dICAdt, dIHRdt,
				dIURdt, dIHDdt, dIUDdt, dIRdt, dRdt, dDdt, dNdt]
		
	def solve(self,t0,tf):
		self.t0=t0
		self.tf=tf
		y0= [self.S0, self.E0, self.IM0, self.IC0, self.ICA0, 
			 self.IHR0, self.IUR0, self.IHD0, self.IUD0, self.IR0, self.R0, 
			 self.D0, self.N0]
		t_vect= np.linspace(self.t0, self.tf, self.tf*100)
		solution= odeint(self.ODES,y0,t_vect)
		
		self.S_vect=solution.T[0]
		self.E_vect=solution.T[1]
		self.IM_vect=solution.T[2]
		self.IC_vect=solution.T[3]
		self.ICA_vect=solution.T[4]
		self.IHR_vect=solution.T[5]
		self.IUR_vect=solution.T[6]
		self.IHD_vect=solution.T[7]
		self.IUD_vect=solution.T[8]
		self.IR_vect=solution.T[9]
		self.R_vect=solution.T[10]
		self.D_vect=solution.T[11]
		self.N_vect=solution.T[12]
		
########################################################################
#################               PARAMETERS             #################
########################################################################
Indice_Movilidad= pd.read_excel(str(BBDD_Dir)+'A.xlsx')
A_vect=Indice_Movilidad['Movilidad'][:]

model=SEIRH()

#1. Transmisión
model.var_trans(beta_0=1.3750,		#R0=3.0
				#beta_0=1.1439,		#R0=2.5
				beta_1=0.019899665, 
				beta_H=0.01,		#Contacto  HR, UR, HD, UD, R
				r=0.35,
				A=0.6,
				Movilidad= A_vect
				)
				
#2. Tiempos de estadía
model.var_t_estadia(omega=1/4.6,		#T. prom. latencia
					gamma_M=1/2.1,		#T. prom IM
					sigma_C=1/3.,		#T. inf. hosp
					sigma_CA=1/4.1,
					gamma_HR=1/9.5,		#T. prom HR->R
					nu=1/11.3,			#T. prom HU->IR
					gamma_R=1/3.4,		#T. prom IR->R
					sigma_HD=1/7.6,		#T. prom HD->D
					sigma_UD=1/10.1,	#T. prom UD->D
					)

#3. Prob Transición Hospital y UCI
model.var_H_UCI(delta_M = 0.965578477,
				delta_HR = 0.696594546,
				delta_UR = 0.122566499,
				delta_HD = 0.058272457
				)

#4. Condiciones iniciales

#Dist 40% y 60%

model.var_ini(N0= 7592871.0,
			E0= 1018.685,
			IM0= 983.62, 
			IC0= 36.91, 
			IHR0= 107.91, 
			IUR0= 23.5, 
			IHD0= 9.03, 
			IUD0= 23.5,	
			D0  = 49
			)
"""	
#Sin acciones de mitigación
model.var_ini(N0= 7592871.0,
			E0= 243.0,
			IM0= 78.21,
			IC0= 2.79
			)
"""
########################################################################
#################                 PLOT                 #################
########################################################################
T_total=600
T_ini=42
model.solve(t0=T_ini,tf=T_total)

fig = plt.figure()
ax = plt.subplot(111)	
counter=0;
SS_data= pd.read_excel(str(BBDD_Dir)+'RES_SecSalud_Rev.xlsx', sheet_name='A06R3')
SS_time=SS_data['t']

results={'t': np.linspace(T_ini,T_total,len(model.S_vect)),
    'Susceptibles': np.array(model.S_vect),
	'Requieren_HG': np.array(model.IHD_vect)+np.array(model.IR_vect)+np.array(model.IHR_vect),
	'Requieren_UCI': np.array(model.IUR_vect)+np.array(model.IUD_vect),
	'Sintomas_Leves': np.array(model.IM_vect),
	'Fallecidos': np.array(model.D_vect),
	'Recuperados': np.array(model.R_vect)
	}

df=pd.DataFrame(results,columns=['t','Susceptibles','Requieren_HG','Requieren_UCI','Sintomas_Leves','Fallecidos','Recuperados'])
df.to_excel('Data.xlsx')

'''PLOT'''
nombre='Susceptibles'; Data=model.S_vect; SS=SS_data['Susceptibles']
nombre='Requieren_HG'; Data=np.array(model.IHD_vect)+np.array(model.IR_vect)+np.array(model.IHR_vect); SS=SS_data['Requieren_HG']
nombre='Requieren_UCI'; Data=np.array(model.IUR_vect)+np.array(model.IUD_vect); SS=SS_data['Requieren_UCI']
#nombre='Sintomas_Leves'; Data=np.array(model.IM_vect); SS=SS_data['Sintomas_Leves']
#nombre='Fallecidos'; Data=np.array(model.D_vect); SS=SS_data['Fallecidos']
#nombre='Recuperados'; Data=np.array(model.R_vect); SS=SS_data['Recuperados']

t_span=np.linspace(T_ini,T_total,len(Data))



plt.title(nombre)
ax.plot(SS_time,SS, linestyle=' ', marker='*', color='black', label='Sec. Salud')
plt.ylabel(u"Población", weight='bold', fontsize=10)
ax.plot(t_span,Data, linestyle='-', marker=' ', label=r'Aislamiento $I_{CA}$')

ax.ticklabel_format(axis='y',style='sci',scilimits=(0,3))
plt.yticks(size=8)
plt.xticks(size=8)
ax.legend(loc='best')

plt.savefig(str(nombre)+'.png',dpi = 300)
plt.show()


