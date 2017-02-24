# Optimize a gasoline-powered general-aviation aircraft, representative
# of the Cessna 402

import math
import numpy as np
from matplotlib import pyplot as plt
from gpkit import Variable, VectorVariable, Model
import pint

# Set up units
ureg = pint.UnitRegistry()

print "This is a test of the Git branching feature (branch = object_oriented)."

# Design variables

W_TO = Variable("W_TO","N","Takeoff weight")
W_fuel = Variable("W_fuel","N","Fuel weight")
W_ZF = Variable("W_ZF","N","Zero-fuel weight")
W_wing = Variable("W_wing","N","Wing weight")
S = Variable("S","m^2","Wing reference area")
b = Variable("b","m","Wingspan")
AR = Variable("AR","-","Aspect ratio")
S_HT = Variable("S_HT","m^2","Tailplane reference area")
S_VT = Variable("S_VT","m^2","Fin reference area")
CD_cruise = Variable("CD_cruise","-","Cruise drag coefficient")
CL_cruise = Variable("CL_cruise","-","Cruise lift coefficient")
CD_TO = Variable("CD_TO","-","Drag coefficient at takeoff")

S_wet = Variable("S_wet","m^2","Aircraft wetted area")
Swet_wing = Variable("Swet_wing","m^2","Wing wetted area")
Swet_HT = Variable("Swet_HT","m^2","Tailplane wetted area")
Swet_VT = Variable("Swet_VT","m^2","Fin wetted area")

y = Variable("y","-","Takeoff-distance parameter 1")
zeta = Variable("zeta","-","Takeoff-distance parameter 2")

N_ult = Variable("N_ult","-","Ultimate load factor")


# Constants (both models)

mission_range = Variable("Range",1.852*np.linspace(1000,2000,10),"km","Mission range (km)")#sweep
s_TO = Variable("s_TO",1100,"ft","Takeoff distance")

V_cruise = Variable("V_cruise",180,"knots","Cruise speed")
h_cruise = Variable("h_cruise",8000,"ft","Cruise altitude")
rho_cruise = Variable("rho_cruise",0.96287,"kg/m^3","Air density at altitude")
rho_SL = Variable("rho_SL",1.225,"kg/m^3","Air density at sea level")

num_crew = Variable("num_crew",2,"-","Number of crew members")
crew_mass = Variable("crew_mass",200,"lbs","Mass per crew member")
num_passengers = Variable("num_passengers",8,"-","Number of passengers")
passenger_mass = Variable("passenger_mass",180+40,"lbs","Mass per passenger (including baggage)")

c_HT = Variable("c_HT",0.8,"-","Tailplane volume coefficient")
c_VT = Variable("c_VT",0.07,"-","Fin volume coefficient")
l_HT = Variable("l_HT",5.00,"m","Tailplane moment arm")
l_VT = Variable("l_VT",6.79,"m","Fin moment arm")
AR_HT = Variable("AR_HT",4.20,"-","Tailplane aspect ratio")
AR_VT = Variable("AR_VT",1.16,"-","Fin aspect ratio")
t_c = Variable("t_c",0.15,"-","Thickness-to-chord ratio") # used for wing, HT, and VT

Swet_fuse = Variable("Swet_fuse",43.6,"m^2","Fuselage wetted area")
Swet_nacelle = Variable("Swet_nacelle",8.172,"m^2","Nacelle wetted area")
wf_wing = Variable("wf_wing",1.6+2,"m","Fuselage width at wing") #Nacelle effect accounted for
wf_tail = Variable("wf_tail",0.7,"m","Fuselage width at tail")
CDA_gear = Variable("CDA_gear",1.149,"m^2","Landing-gear drag area")

c_fe = Variable("c_fe",0.0045,"-","Equivalent skin-friction coefficient")
e = Variable("e",0.8,"-","Oswald efficiency factor")
CLmax_cruise = Variable("CLmax_cruise",1.6,"-","Maximum lift coefficient (cruise)")
CLmax_TO = Variable("CLmax_TO",2.0,"-","Maximum lift coefficient (takeoff)")
CLmax_LDG = Variable("CLmax_LDG",2.6,"-","Maximum lift coefficient (landing)")

N_lim = Variable("N_lim",3.5,"-","Limit load factor")
rho_AL = Variable("rho_AL",2810,"kg/m^3","Density of 7075-T6 aluminum")
sigma_AL = Variable("sigma_AL",503e6,"N/m^2","Allowable stress of 7075-T6 aluminum")

wing_SA_multiplier = Variable("wing_SA_multiplier",3.60,"kg/m^2","Wing surface-area weight multiplier")
weight_multiplier_tail = Variable("weight_multiplier_tail",10,"kg/m^2","Tail weight multiplier")
weight_multiplier_fuse = Variable("weight_multiplier_fuse",7,"kg/m^2","Fuselage weight multiplier")

g = Variable("g",9.807,"m/s**2","Gravitational constant")

# Engine-related variables & constants
P_fixedEngine = Variable("P_fixedEngine",242,"kW","Engine power (fixed engine)")
m_engine = Variable("m_engine",1.4*227,"kg","Installed engine mass (1 engine)")
num_engines = Variable("num_engines",2,"-","Number of engines")
eta_prop = Variable("eta_prop",0.8,'-',"Propeller efficiency")
SFC = Variable("SFC",6.94e-8,"s**2/m**2",\
        "Engine power-specific fuel consumption") # careful with unit conversion

#Rubber-engine parameters
P_rubberEngine = Variable("P_rubberEngine","kW","Engine power (rubber engine)")
K_p = Variable("K_p",1.4*9.2e-3,"N/W","Engine-weight scale factor (Newtons per watt)")


# Other constants
pi = math.pi
lambda_wing = 0.6 #Included here to avoid posynomial-compatibility issues
fuel_allowance = 0.03+0.015+0.005+0.06 #allowance for takeoff, climb, landing, reserves, and trapped fuel

# Initialize constraints
constraints = [] #constraints that apply to both aircraft
constraints_fixedEngine = [] 
constraints_rubberEngine = []

# Geometric model and constraints (including wetted-area constraints)

constraints += [AR == b**2/S,
                Swet_wing >= 1.6*S*(1+0.25*t_c),
                Swet_HT >= 1.6*S_HT*(1+0.25*t_c),
                Swet_VT >= 1.6*S_VT*(1+0.25*t_c),
                S_wet >= Swet_fuse + num_engines*Swet_nacelle + Swet_wing + Swet_HT + Swet_VT]

# Aerodynamic model and constraints
CD0_cruise = c_fe * S_wet / S
CDi_cruise = (CL_cruise**2) / (pi*e*AR)
constraints += [CD_cruise >= CD0_cruise + CDi_cruise]

# Stability & control model and constraints
MAC = (S/AR)**(1./2) #mean aerodynamic chord
b = (AR * S)**(1./2) #wingspan
constraints += [S_VT == c_VT*b*S / l_VT,
                S_HT == c_HT*MAC*S / l_HT]

# Load-factor model and constraints (FAR 23.337)

constraints += [ N_ult == N_lim*1.5]

# Structure/weight model and constraints
lambda_term = (1+2*lambda_wing) / (1+lambda_wing)
WTO_term = (b**3 * N_ult * W_TO) / (8 * S * t_c)

constraints += [W_wing >= g*((4./5)*(0.56*rho_AL/sigma_AL + rho_AL/sigma_AL) \
                *WTO_term*lambda_term + wing_SA_multiplier*Swet_wing)]

W_HT = g * weight_multiplier_tail * S_HT
W_VT = g * weight_multiplier_tail * S_VT
W_fuse = g * weight_multiplier_fuse * Swet_fuse
W_gear = 0.057 * W_TO #landing-gear weight
W_all_else = 0.1 * W_TO
W_crew = g * num_crew * crew_mass
W_payload = g * num_passengers * passenger_mass

W_fixedEngines = g * num_engines * m_engine #installed engine weight (fixed)
W_rubberEngines = num_engines * K_p * P_rubberEngine #installed engine weight (rubber)


#note the 5% margin on empty weight
constraints_fixedEngine += [W_ZF >= 1.05 * (W_wing + W_HT + W_VT + W_fuse + W_gear
                + W_fixedEngines + W_all_else) + W_crew + W_payload,
                W_TO >= W_ZF + W_fuel]
constraints_rubberEngine += [W_ZF >= 1.05 * (W_wing + W_HT + W_VT + W_fuse + W_gear
                + W_rubberEngines + W_all_else) + W_crew + W_payload,
                W_TO >= W_ZF + W_fuel]                

#Performance model (cruise) and constraints
z = (g * mission_range * SFC * CD_cruise) / (eta_prop * CL_cruise)

constraints += [W_TO == 0.5 * rho_cruise * V_cruise**2 * S * CL_cruise, #steady level flight
                W_fuel >= (1+fuel_allowance)*(z + z**2/np.math.factorial(2)
			+ z**3/np.math.factorial(3)
                        + z**4/np.math.factorial(4.)) * W_ZF] #Breguet range (Taylor)


#Performance model (takeoff) and constraints
V_TO = 1.1*((2*W_TO) / (rho_SL*S*CLmax_TO))**(1./2)#1.1*V_stall
D_TO = 0.5 * rho_SL * V_TO**2 * S * CD_TO
T_TO_fixedEngine = eta_prop * num_engines * P_fixedEngine / V_TO
T_TO_rubberEngine = eta_prop * num_engines * P_rubberEngine / V_TO

constraints_fixedEngine += [CD_TO >= CD0_cruise + CDA_gear/S,
		zeta == D_TO / T_TO_fixedEngine,
                2*g*s_TO*T_TO_fixedEngine / (W_TO*V_TO**2) >= 1 + y,
                1 >= 0.0464*zeta**2.73/y**2.88 + 1.044*zeta**0.296/y**0.049]
constraints_rubberEngine += [CD_TO >= CD0_cruise + CDA_gear/S,
                zeta == D_TO / T_TO_rubberEngine,
                2*g*s_TO*T_TO_rubberEngine / (W_TO*V_TO**2) >= 1 + y,
                1 >= 0.0464*zeta**2.73/y**2.88 + 1.044*zeta**0.296/y**0.049]


# Assemble and run the models
model_fixedEngine = Model(W_TO,constraints + constraints_fixedEngine)
sol_fixedEngine = model_fixedEngine.solve(verbosity=0)

model_rubberEngine = Model(W_TO,constraints + constraints_rubberEngine)
sol_rubberEngine = model_rubberEngine.solve(verbosity=0)


# Plotting commands
plt.ion()

# Cessna 402 data
cessna_402_data = {'W_TO':7210*ureg.lbf,
	"W_fuel":1278*ureg.lbf,
        "W_wing":860*ureg.lbf,#Cessna 404 (not 402)
	"b":13.44*ureg.m,
	"AR":8.6,
	"range":1234*ureg.nautical_mile,
        "P":P_fixedEngine.value}

fig1 = plt.figure(figsize=(16, 12), dpi=80)
plt.show()


#MTOW vs. Mission Range
plt.subplot(2,3,1)
plt.plot(sol_fixedEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude,
	sol_fixedEngine["freevariables"][W_TO].to(ureg.lbf).magnitude,
        color="black", linewidth=1.5, linestyle="-", marker='s', markersize=8,
        label='GPKit model (fixed engine)')
plt.plot(sol_rubberEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude,
        sol_rubberEngine["freevariables"][W_TO].to(ureg.lbf).magnitude,
        color="black", linewidth=1.5, linestyle="-", marker='^', markersize=8,
        label='GPKit model (rubber engine)')
plt.plot(cessna_402_data["range"].to(ureg.nautical_mile).magnitude,
	cessna_402_data["W_TO"].to(ureg.lbf).magnitude, 
	marker='o',color="black", markersize=12, label="Cessna 402")		 
plt.grid()
plt.xlabel('Mission range (nm)', fontsize = 16)
plt.ylabel('MTOW (lbs)', fontsize = 16)
plt.title("Maximum Takeoff Weight",fontsize = 20)
plt.legend(numpoints = 1,loc='lower left', fontsize = 12)
plt.xlim(xmin = np.min(sol_fixedEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude),
         xmax = np.max(sol_fixedEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude))
plt.ylim(ymin = 0,ymax = 1.1*np.max(sol_fixedEngine["freevariables"][W_TO].to(ureg.lbf).magnitude))


#Fuel Weight vs. Mission Range
plt.subplot(2,3,2)
plt.plot(sol_fixedEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude,
	sol_fixedEngine["freevariables"][W_fuel].to(ureg.lbf).magnitude,
         color="black", linewidth=1.5, linestyle="-", marker='s', markersize=8,
         label='GPKit model (fixed engine)')
plt.plot(sol_rubberEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude,
        sol_rubberEngine["freevariables"][W_fuel].to(ureg.lbf).magnitude,
         color="black", linewidth=1.5, linestyle="-", marker='^', markersize=8,
         label='GPKit model (rubber engine)')
plt.plot(cessna_402_data["range"].to(ureg.nautical_mile).magnitude,
	cessna_402_data["W_fuel"].to(ureg.lbf).magnitude, 
	marker='o',color="black", markersize=12, label="Cessna 402")
plt.grid()
plt.xlabel('Mission range (nm)', fontsize = 16)
plt.ylabel('Fuel weight (lbs)', fontsize = 16)
plt.title("Fuel Weight",fontsize = 20)
plt.legend(numpoints = 1, loc='lower right', fontsize = 12)
plt.xlim(xmin = np.min(sol_fixedEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude),
         xmax = np.max(sol_fixedEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude))

#Wing Weight vs. Mission Range
plt.subplot(2,3,3)
plt.plot(sol_fixedEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude,
        sol_fixedEngine["freevariables"][W_wing].to(ureg.lbf).magnitude,
         color="black", linewidth=1.5, linestyle="-", marker='s', markersize=8,
         label='GPKit model (fixed engine)')
plt.plot(sol_rubberEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude,
        sol_rubberEngine["freevariables"][W_wing].to(ureg.lbf).magnitude,
         color="black", linewidth=1.5, linestyle="-", marker='^', markersize=8,
         label='GPKit model (rubber engine)')
plt.plot(cessna_402_data["range"].to(ureg.nautical_mile).magnitude,
        cessna_402_data["W_wing"].to(ureg.lbf).magnitude, 
        marker='o',color="black", markersize=12, label="Cessna 404 (not 402)")
plt.grid()
plt.xlabel('Mission range (nm)', fontsize = 16)
plt.ylabel('Wing weight (lbs)', fontsize = 16)
plt.title("Wing Weight",fontsize = 20)
plt.legend(numpoints = 1, loc='lower right', fontsize = 12)
plt.xlim(xmin = np.min(sol_fixedEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude),
         xmax = np.max(sol_fixedEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude))

#Wingspan vs. Mission Range
plt.subplot(2,3,4)
plt.plot(sol_fixedEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude,
	sol_fixedEngine["freevariables"]["b"].to(ureg.ft).magnitude,
         color="black", linewidth=1.5, linestyle="-", marker='s', markersize=8,
         label='GPKit model (fixed engine)')
plt.plot(sol_rubberEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude,
        sol_rubberEngine["freevariables"]["b"].to(ureg.ft).magnitude,
         color="black", linewidth=1.5, linestyle="-", marker='^', markersize=8,
         label='GPKit model (rubber engine)')
plt.plot(cessna_402_data["range"].to(ureg.nautical_mile).magnitude,
	cessna_402_data["b"].to(ureg.ft).magnitude, 
	marker='o',color="black", markersize=12, label="Cessna 402")
plt.grid()
plt.xlabel('Mission range (nm)', fontsize = 16)
plt.ylabel('Wingspan (ft)', fontsize = 16)
plt.title("Wingspan",fontsize = 20)
plt.legend(numpoints = 1, loc='lower right', fontsize = 12)
plt.xlim(xmin = np.min(sol_fixedEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude),
         xmax = np.max(sol_fixedEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude))

#Aspect Ratio vs. Mission Range
plt.subplot(2,3,5)
plt.plot(sol_fixedEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude,
	sol_fixedEngine["freevariables"][AR].magnitude,
         color="black", linewidth=1.5, linestyle="-", marker='s', markersize=8,
         label='GPKit model (fixed engine)') 
plt.plot(sol_rubberEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude,
        sol_rubberEngine["freevariables"][AR].magnitude,
         color="black", linewidth=1.5, linestyle="-", marker='^', markersize=8,
         label='GPKit model (rubber engine)')
plt.plot(cessna_402_data["range"].to(ureg.nautical_mile).magnitude,
	cessna_402_data["AR"], 
	marker='o',color="black", markersize=12, label="Cessna 402")		 	 
plt.grid()
plt.xlabel('Mission range (nm)', fontsize = 16)
plt.ylabel('Aspect ratio (dimensionless)', fontsize = 16)
plt.title("Aspect Ratio",fontsize = 20)
plt.legend(numpoints=1, loc='upper right', fontsize = 12)
plt.xlim(xmin = np.min(sol_fixedEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude),
         xmax = np.max(sol_fixedEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude))

#Engine Power vs. Mission Range
plt.subplot(2,3,6)
plt.plot(sol_rubberEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude,
        sol_rubberEngine["freevariables"]["P_rubberEngine"].to(ureg.horsepower).magnitude,
         color="black", linewidth=1.5, linestyle="-", marker='^', markersize=8,
         label='GPKit model (rubber engine)')
plt.plot(cessna_402_data["range"].to(ureg.nautical_mile).magnitude,
        cessna_402_data["P"].to(ureg.horsepower).magnitude, 
        marker='o',color="black", markersize=12, label="Cessna 402")                     
plt.grid()
plt.xlabel('Mission range (nm)', fontsize = 16)
plt.ylabel('Engine power (HP)', fontsize = 16)
plt.title("Engine Power",fontsize = 20)
plt.legend(numpoints=1, loc='lower right', fontsize = 12)
plt.xlim(xmin = np.min(sol_fixedEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude),
         xmax = np.max(sol_fixedEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude))

title_str = "Key Design Variables vs. Mission Range for a Twin-Engined GA Aircraft"
plt.suptitle(title_str,fontsize = 24)

plt.tight_layout()#makes sure subplots are spaced neatly
plt.subplots_adjust(left=0.125,right=0.9,bottom=0.1,top = 0.9)#adds space at the top for the title
