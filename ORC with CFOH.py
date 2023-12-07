import CoolProp
import CoolProp.CoolProp as CP
print('CoolProp version:', CoolProp.__version__)
from CoolProp.CoolProp import PropsSI
from CoolProp.Plots import PropertyPlot
from CoolProp.Plots import SimpleCompressionCycle
import matplotlib.pyplot as plt
import numpy as np
from sympy import Eq, solve
from CoolProp.CoolProp import PropsSI
from CoolProp.HumidAirProp import HAPropsSI

#ORC with CFOH

#Input Paramters

eta_pump=.7
eta_turbine1=.75
eta_turbine2=.75
T_i_T=429.42
B_P=3
C_o_T=300
Intermediate_P=0.9367

T_H=600
T_L=293
T_O=273

#States

#State 1
x1=0
t1=C_o_T
p1=PropsSI("P","T",t1,"Q",x1,'R245fa')/1000000
h1=PropsSI("H","T",t1,"Q",x1,'R245fa')/1000
s1=PropsSI("S","T",t1,"Q",x1,'R245fa')/1000

#State 2s
p2s=B_P
s2s=s1
t2s=PropsSI("T","P",p2s*1000000,"S",s2s*1000,'R245fa')
h2s=PropsSI("H","P",p2s*1000000,"S",s2s*1000,'R245fa')/1000

#State 2
p2=B_P
h2=h1+(h2s-h1)/eta_pump
t2=PropsSI("T","P",p2*1000000,"H",h2*1000,'R245fa')
s2=PropsSI("S","P",p2*1000000,"H",h2*1000,'R245fa')/1000

#State 8
x8=0
p8=Intermediate_P
t8=PropsSI("T","P",p8*1000000,"Q",x8,'R245fa')
h8=PropsSI("H","P",p8*1000000,"Q",x8,'R245fa')/1000
s8=PropsSI("S","P",p8*1000000,"Q",x8,'R245fa')/1000

#State 9s
p9s=B_P
s9s=s8
t9s=PropsSI("T","P",p9s*1000000,"S",s9s*1000,'R245fa')
h9s=PropsSI("H","P",p9s*1000000,"S",s9s*1000,'R245fa')/1000

#State 9
p9=B_P
h9=h8+(h9s-h8)/eta_pump
t9=PropsSI("T","P",p9*1000000,"H",h9*1000,'R245fa')
s9=PropsSI("S","P",p9*1000000,"H",h9*1000,'R245fa')/1000

#State 5
p5=B_P
t5=T_i_T
h5=PropsSI("H","P",p5*1000000,"T",t5,'R245fa')/1000
s5=PropsSI("S","P",p5*1000000,"T",t5,'R245fa')/1000

#State 6s
p6s=Intermediate_P
s6s=s5
t6s=PropsSI("T","P",p6s*1000000,"S",s6s*1000,'R245fa')
h6s=PropsSI("H","P",p6s*1000000,"S",s6s*1000,'R245fa')/1000

#State 6
p6=Intermediate_P
h6=h5-(h5-h6s)*eta_turbine1
t6=PropsSI("T","P",p6*1000000,"H",h6*1000,'R245fa')
s6=PropsSI("S","P",p6*1000000,"H",h6*1000,'R245fa')/1000
y1=(h3-h2)/((h6-h8)+(h3-h2))

#State 7s
p7s=p1
s7s=s6
t7s=PropsSI("T","P",p7s*1000000,"S",s7s*1000,'R245fa')
h7s=PropsSI("H","P",p7s*1000000,"S",s7s*1000,'R245fa')/1000

#State 7
p7=p1
h7=h6-(h6-h7s)*eta_turbine2
t7=PropsSI("T","P",p7*1000000,"H",h7*1000,'R245fa')
s7=PropsSI("S","P",p7*1000000,"H",h7*1000,'R245fa')/1000
y2=(1-y1)

#State 3
p3=B_P
t3=t8
h3=PropsSI("H","P",p3*1000000,"T",t3,'R245fa')/1000
s3=PropsSI("S","P",p3*1000000,"T",t3,'R245fa')/1000

#State 4
p4=B_P
h4=y1*h9+y2*h3
t4=PropsSI("T","P",p4*1000000,"H",h4*1000,'R245fa')
s4=PropsSI("S","P",p4*1000000,"H",h4*1000,'R245fa')/1000

#Performance Paramters

pump1_e=(1-y1)*(h2-h1)
cfoh_e=(1-y1)*(h3-h2)
pump2_e=y1*(h9-h8)
mixer_e=y1*h9+(1-y1)*h3
evaporator_e=h5-h4
expander_e=(h5-h6s)*eta_turbine1+(1-y1)*(h6-h7s)*eta_turbine2
condenser_e=(1-y1)*(h7-h1)

pump1_i=T_O*(1-y1)*(s2-s1)
cfoh_i=T_O*(y1*(s8-s6)+(1-y1)*(s3-s2))
pump2_i=T_O*y1*(s9-s8)
mixer_i=T_O*(s4-y1*s9-(1-y1)*s3)
evaporator_i=T_O*((s5-s4)-(h5-h4)/T_H)
expander_i=T_O*((s6-s5)+(1-y1)*(s7-s6))
condenser_i=T_O*(1-y1)*((s1-s7)-(h1-h7)/T_L)

eta_thermal=1-(condenser_e/evaporator_e)

#Mass Flow Rate

q_in=50                                   #MW
m=q_in/evaporator_e
w_net=eta_thermal*q_in

#Result

print('E=Work & I=Exergy')
print('Pump1_E =',pump1_e,'kW')
print('CFOH_E =',cfoh_e,'kW')
print('Pump2_E =',pump2_e,'kW')
print('Mixer_E =',mixer_e,'kW')
print('Evaporator_E =',evaporator_e,'kW')
print('Expander_E =',expander_e,'kW')
print('Condenser_E =',condenser_e,'kW')
print('Pump1_I =',pump1_i,'kW')
print('CFOH_I =',cfoh_i,'kW')
print('Pump2_I =',pump2_i,'kW')
print('Mixer_I =',mixer_i,'kW')
print('Evaporator_I =',evaporator_i,'kW')
print('Expander_I =',expander_i,'kW')
print('Condenser_I =',condenser_i,'kW')
print('Thermal_efficiency =',eta_thermal)

print('Properties of state 1: p1 =',p1, 'MPa , t1 = ',t1, ' K, h1=',h1,'kJ/kh2s , s1=',s1,'kJ/kh2sK')
print('Properties of state 2: p2 =',p2, 'MPa , t2 = ',t2, ' K, h2=',h2,'kJ/kh2s , s2=',s2,'kJ/kh2sK')
print('Properties of state 3: p3 =',p3, 'MPa , t3 = ',t3, ' K, h3=',h3,'kJ/kh2s , s3=',s3,'kJ/kh2sK')
print('Properties of state 4: p4 =',p4, 'MPa , t4 = ',t4, ' K, h4=',h4,'kJ/kh2s , s4=',s4,'kJ/kh2sK')
print('Properties of state 5: p5 =',p5, 'MPa , t5 = ',t5, ' K, h5=',h5,'kJ/kh2s , s5=',s5,'kJ/kh2sK')
print('Properties of state 6: p6 =',p6, 'MPa , t6 = ',t6, ' K, h6=',h6,'kJ/kh2s , s6=',s6,'kJ/kh2sK')
print('Properties of state 7: p7 =',p7, 'MPa , t7 = ',t7, ' K, h7=',h7,'kJ/kh2s , s7=',s7,'kJ/kh2sK')
print('Properties of state 8: p8 =',p8, 'MPa , t8 = ',t8, ' K, h8=',h8,'kJ/kh2s , s8=',s8,'kJ/kh2sK')
print('Properties of state 9: p9 =',p9, 'MPa , t9 = ',t9, ' K, h9=',h9,'kJ/kh2s , s9=',s9,'kJ/kh2sK')

print('W_net =',w_net,'MW')