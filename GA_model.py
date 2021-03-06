# Optimize a gasoline-powered general-aviation aircraft, representative
# of the Cessna 402

import math
import numpy as np
from matplotlib import pyplot as plt
from gpkit import Variable, VectorVariable, Model
import pint

# Set up units
ureg = pint.UnitRegistry()

# Design variables

class Aircraft(Model): #Overall vehicle model
      
    def setup(self, engineType="fixed"):
        W_ZF = Variable("W_ZF","N","Zero-fuel weight")
        W_fuel = Variable("W_fuel","N","Fuel weight")
        W_TO = Variable("W_TO","N","Takeoff weight")
        self.Swet = Variable("Swet","m^2","Aircraft wetted area")
        
        self.W_TO = W_TO

        self.fuse = Fuselage()
        self.wing = Wing(self)
        self.tailplane = Tailplane()
        self.fin = VerticalFin()
        self.engines = Engines(engineType)
        self.nacelles = Nacelles()
        self.landing_gear = LandingGear(self)
        self.crew = Crew()
        self.passengers = Passengers()
        self.allOtherWeights = AllOtherWeights(self)

        self.components = [self.fuse, self.wing, self.tailplane,self.fin,
            self.engines, self.nacelles,self.landing_gear,self.crew,
            self.passengers,self.allOtherWeights]

        constraints = []
        constraints += [self.components,#all constraints implemented at component level
            W_TO >= W_ZF + W_fuel, #top-level weight constraint
            W_ZF >= 1.05*sum(c["W"] for c in self.components), #zero-fuel weight constraint (with 5% margin)
            self.Swet >= sum(c["Swet"] for c in self.components)] #wetted-area constraint

        # Stability & control constraints (volume coefficients)
        constraints += [self.fin["S_VT"] == self.fin["c_VT"]*self.wing["b"]*self.wing["S"]/self.fin["l_VT"],
            self.tailplane["S_HT"] == self.tailplane["c_HT"]*self.wing["MAC"]*self.wing["S"]/self.tailplane["l_HT"]]


        return constraints

class Fuselage(Model):
    def setup(self):
        g = Variable("g",9.807,"m/s**2","Gravitational constant")
        weight_multiplier_fuse = Variable("weight_multiplier_fuse",7,
            "kg/m^2","Fuselage weight multiplier")
        Swet = Variable("Swet",43.6,"m^2","Fuselage wetted area")
        W = Variable("W","N","Fuselage weight")

        return [W == g * weight_multiplier_fuse * Swet]

class Wing(Model):
    def setup(self,aircraft):
        W = Variable("W","N","Wing weight")
        S = Variable("S","m^2","Wing reference area")
        b = Variable("b","m","Wingspan")
        AR = Variable("AR","-","Aspect ratio")
        MAC = Variable("MAC","m","Mean aerodynamic chord")
        Swet = Variable("Swet","m^2","Wing wetted area")
        t_c = Variable("t_c",0.15,"-","Thickness-to-chord ratio")
        N_lim = Variable("N_lim",3.5,"-","Limit load factor")
        N_ult = Variable("N_ult","-","Ultimate load factor")
        rho_AL = Variable("rho_AL",2810,"kg/m^3","Density of 7075-T6 aluminum")
        sigma_AL = Variable("sigma_AL",503e6,"N/m^2",
            "Allowable stress of 7075-T6 aluminum")
        wing_SA_multiplier = Variable("wing_SA_multiplier",
            3.60,"kg/m^2","Wing surface-area weight multiplier")
        g = Variable("g",9.807,"m/s**2","Gravitational constant")
        WTO_term = Variable("WTO_term","N*m","Takeoff-weight term")

        lambda_wing = 0.6 #Included here to avoid posynomial-compatibility issues
        lambda_term = (1+2*lambda_wing) / (1+lambda_wing)


        return [AR == b**2/S, 
            Swet >= 1.6*S*(1+0.25*t_c),
            MAC == (S/AR)**(1./2), #mean aerodynamic chord
            b == (AR*S)**(1./2), #wingspan
            N_ult == N_lim*1.5, 
            WTO_term == (b**3 * N_ult * aircraft.W_TO) / (8 * S * t_c),
            W >= g*((4./5)*(0.56*rho_AL/sigma_AL + rho_AL/sigma_AL) \
                *WTO_term*lambda_term + wing_SA_multiplier*Swet)]


class Tailplane(Model):
    def setup(self):
        S_HT = Variable("S_HT","m^2","Tailplane reference area")
        Swet = Variable("Swet","m^2","Tailplane wetted area")
        c_HT = Variable("c_HT",0.8,"-","Tailplane volume coefficient")
        l_HT = Variable("l_HT",5.00,"m","Tailplane moment arm")
        AR_HT = Variable("AR_HT",4.20,"-","Tailplane aspect ratio")
        t_c = Variable("t_c",0.15,"-","Thickness-to-chord ratio")
        weight_multiplier_tail = Variable("weight_multiplier_tail",10,
            "kg/m^2","Tail weight multiplier")
        g = Variable("g",9.807,"m/s**2","Gravitational constant")
        W = Variable("W","N","Tailplane weight")

        return [Swet >= 1.6*S_HT*(1+0.25*t_c),
            W == g * weight_multiplier_tail * S_HT]
        

class VerticalFin(Model):
    def setup(self):
        S_VT = Variable("S_VT","m^2","Fin reference area")
        Swet = Variable("Swet","m^2","Fin wetted area")
        c_VT = Variable("c_VT",0.07,"-","Fin volume coefficient")
        l_VT = Variable("l_VT",6.79,"m","Fin moment arm")
        AR_VT = Variable("AR_VT",1.16,"-","Fin aspect ratio")
        t_c = Variable("t_c",0.15,"-","Thickness-to-chord ratio")
        weight_multiplier_tail = Variable("weight_multiplier_tail",10,
            "kg/m^2","Tail weight multiplier")
        g = Variable("g",9.807,"m/s**2","Gravitational constant")
        W = Variable("W","N","Vertical-fin weight")

        return [Swet >= 1.6*S_VT*(1+0.25*t_c),
            W == g * weight_multiplier_tail * S_VT]

class Engines(Model): #Fixed engine only for now
    def setup(self,engineType):
        # Engine-related variables & constants
        num_engines = Variable("num_engines",2,"-","Number of engines")
        eta_prop = Variable("eta_prop",0.8,'-',"Propeller efficiency")
        SFC = Variable("SFC",6.94e-8,"s**2/m**2",\
            "Engine power-specific fuel consumption") # careful with unit conversion
        g = Variable("g",9.807,"m/s**2","Gravitational constant")
        W = Variable("W","N","Engine(s) weight")
        Swet = Variable("Swet",0,"m^2","Wetted area of engines")
        P = Variable("P","kW","Engine power")
        
        constraints = []

        #Fixed-engine parameters
        if engineType == "fixed":
            P_fixedEngine = Variable("P_fixedEngine",242,"kW","Engine power (fixed engine)")
            m_engine = Variable("m_engine",1.4*227,"kg","Installed engine mass (1 engine)")
            constraints += [W == g * num_engines * m_engine,
                P == P_fixedEngine]
        #Rubber-engine parameters
        elif engineType == "rubber":
            K_p = Variable("K_p",1.4*9.2e-3,"N/W","Engine-weight scale factor (Newtons per watt)")
            constraints += [W == num_engines * K_p * P]

        return constraints


class Nacelles (Model):
    def setup(self):
        Swet_oneNacelle = Variable("Swet_oneNacelle",8.172,"m^2","Wetted area of 1 nacelle")
        num_nacelles = Variable("num_nacelles",2,"-","Number of nacelles")
        Swet = Variable("Swet","m^2","Wetted area of nacelles")
        W = Variable("W",0,"N","Nacelle weight")
        return [Swet == num_nacelles * Swet_oneNacelle]


class LandingGear(Model):
    def setup(self,aircraft):
        CDA_gear = Variable("CDA_gear",1.149,"m^2","Landing-gear drag area")
        W = Variable("W","N","Landing-gear weight")
        Swet = Variable("Swet",0,"m^2","Wetted area of landing gear")#not included in cruise drag estimate

        return [W == 0.057*aircraft.W_TO]           


class Crew(Model):
    def setup(self):
        num_crew = Variable("num_crew",2,"-","Number of crew members")
        crew_mass = Variable("crew_mass",200,"lbs","Mass per crew member")
        g = Variable("g",9.807,"m/s**2","Gravitational constant")
        W = Variable("W","N","Crew weight")
        Swet = Variable("Swet",0,"m^2","Wetted area of crew")
        
        return [W == g * num_crew * crew_mass]


class Passengers(Model):
    def setup(self):
        num_passengers = Variable("num_passengers",8,"-","Number of passengers")
        passenger_mass = Variable("passenger_mass",180+40,"lbs",
            "Mass per passenger (including baggage)")
        g = Variable("g",9.807,"m/s**2","Gravitational constant")
        W = Variable("W","N","Passenger weight")
        Swet = Variable("Swet",0,"m^2","Wetted area of passengers")
        
        return [W == g * num_passengers * passenger_mass]


class AllOtherWeights(Model):
    def setup(self,aircraft):
        W = Variable("W","N","Landing-gear weight")
        Swet = Variable("Swet",0,"m^2","Wetted area of everything else")
        return [W == 0.1*aircraft.W_TO]


class Mission(Model):
    def setup(self,aircraft):
        
        self.aircraft = aircraft
        W_TO = aircraft["W_TO"]
        g = Variable("g",9.807,"m/s**2","Gravitational constant")

        constraints = []

        #Aerodynamics model
        c_fe = Variable("c_fe",0.0045,"-","Equivalent skin-friction coefficient")
        e = Variable("e",0.8,"-","Oswald efficiency factor")
        CLmax_cruise = Variable("CLmax_cruise",1.6,"-","Maximum lift coefficient (cruise)")
        CLmax_TO = Variable("CLmax_TO",2.0,"-","Maximum lift coefficient (takeoff)")
        CLmax_LDG = Variable("CLmax_LDG",2.6,"-","Maximum lift coefficient (landing)")
        CD_cruise = Variable("CD_cruise","-","Cruise drag coefficient")
        CL_cruise = Variable("CL_cruise","-","Cruise lift coefficient")
        CD_TO = Variable("CD_TO","-","Drag coefficient at takeoff")

        # Aerodynamics constraints (cruise)
        CD0_cruise = c_fe * aircraft.Swet / aircraft.wing["S"]
        CDi_cruise = (CL_cruise**2) / (np.pi*e*aircraft.wing["AR"])
        constraints += [CD_cruise >= CD0_cruise + CDi_cruise]

        # Aerodynamics constraints (takeoff)
        constraints += [CD_TO >= CD0_cruise + aircraft.landing_gear["CDA_gear"]/aircraft.wing["S"]]

        #Performance model (cruise)
        V_cruise = Variable("V_cruise",180,"knots","Cruise speed")
        h_cruise = Variable("h_cruise",8000,"ft","Cruise altitude")
        rho_cruise = Variable("rho_cruise",0.96287,"kg/m^3","Air density at altitude")
        mission_range = Variable("Range",np.linspace(1000,2000,10),"nautical_mile","Mission range")#sweep
        fuel_allowance = 0.03+0.015+0.005+0.06 #allowance for takeoff, climb, landing, reserves, and trapped fuel

        #Performance constraints (cruise)
        z = (g * mission_range * aircraft.engines["SFC"] * CD_cruise) / (aircraft.engines["eta_prop"] * CL_cruise)
        constraints += [aircraft.W_TO == 0.5*rho_cruise*V_cruise**2 * aircraft.wing["S"] * CL_cruise] #steady level flight
        constraints += [aircraft["W_fuel"] >= (1+fuel_allowance)*(z + z**2/np.math.factorial(2)
            + z**3/np.math.factorial(3)
            + z**4/np.math.factorial(4.)) * aircraft["W_ZF"]] #Breguet range (Taylor)

        #Performance model (takeoff)
        y = Variable("y","-","Takeoff-distance parameter 1")
        zeta = Variable("zeta","-","Takeoff-distance parameter 2")
        s_TO = Variable("s_TO",1100,"ft","Takeoff distance")
        rho_SL = Variable("rho_SL",1.225,"kg/m^3","Air density at sea level")

        #Performance constraints (takeoff)
        V_TO = 1.1*((2*aircraft.W_TO) / (rho_SL*aircraft.wing["S"]*CLmax_TO))**(1./2)#1.1*V_stall
        D_TO = 0.5 * rho_SL * V_TO**2 * aircraft.wing["S"] * CD_TO
        T_TO = aircraft.engines["eta_prop"]*aircraft.engines["num_engines"]*aircraft.engines["P"]/V_TO

        constraints += [zeta == D_TO / T_TO,
            2*g*s_TO*T_TO / (aircraft.W_TO*V_TO**2) >= 1 + y,
            1 >= 0.0464*zeta**2.73/y**2.88 + 1.044*zeta**0.296/y**0.049]
        
        return constraints

GA_aircraft_fixedEngine = Aircraft(engineType="fixed")
GA_aircraft_rubberEngine = Aircraft(engineType="rubber")

GA_mission_fixedEngine = Mission(GA_aircraft_fixedEngine)
GA_mission_rubberEngine = Mission(GA_aircraft_rubberEngine)

GA_model_fixedEngine = Model(GA_aircraft_fixedEngine.W_TO,[GA_aircraft_fixedEngine, GA_mission_fixedEngine])
GA_model_rubberEngine = Model(GA_aircraft_rubberEngine.W_TO,[GA_aircraft_rubberEngine, GA_mission_rubberEngine])

GA_solution_fixedEngine = GA_model_fixedEngine.solve(verbosity=0)
GA_solution_rubberEngine = GA_model_rubberEngine.solve(verbosity=0)

#print GA_solution_fixedEngine.summary()
#GA_mission.debug()

# Plotting commands
plt.ion()

# Cessna 402 data
cessna_402_data = {'W_TO':7210*ureg.lbf,
	"W_fuel":1278*ureg.lbf,
    "W_wing":860*ureg.lbf,#Cessna 404 (not 402)
	"b":13.44*ureg.m,
	"AR":8.6,
	"range":1234*ureg.nautical_mile,
    "P":GA_aircraft_fixedEngine.engines["P_fixedEngine"].value}

fig1 = plt.figure(figsize=(16, 12), dpi=80)
plt.show()


#MTOW vs. Mission Range
plt.subplot(2,3,1)
plt.plot(GA_solution_fixedEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude,
	GA_solution_fixedEngine["freevariables"]["W_TO"].to(ureg.lbf).magnitude,
    color="black", linewidth=1.5, linestyle="-", marker='s', markersize=8,
    label='GPKit model (fixed engine)')
plt.plot(GA_solution_rubberEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude,
    GA_solution_rubberEngine["freevariables"]["W_TO"].to(ureg.lbf).magnitude,
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
plt.xlim(xmin = np.min(GA_solution_fixedEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude),
     xmax = np.max(GA_solution_fixedEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude))
plt.ylim(ymin = 0,ymax = 1.1*np.max(GA_solution_fixedEngine["freevariables"]["W_TO"].to(ureg.lbf).magnitude))


#Fuel Weight vs. Mission Range
plt.subplot(2,3,2)
plt.plot(GA_solution_fixedEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude,
	GA_solution_fixedEngine["freevariables"]["W_fuel"].to(ureg.lbf).magnitude,
     color="black", linewidth=1.5, linestyle="-", marker='s', markersize=8,
     label='GPKit model (fixed engine)')
plt.plot(GA_solution_rubberEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude,
    GA_solution_rubberEngine["freevariables"]["W_fuel"].to(ureg.lbf).magnitude,
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
plt.xlim(xmin = np.min(GA_solution_fixedEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude),
     xmax = np.max(GA_solution_fixedEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude))
plt.ylim(ymin = 0)


#Wing Weight vs. Mission Range
plt.subplot(2,3,3)
plt.plot(GA_solution_fixedEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude,
    GA_solution_fixedEngine["freevariables"]["W_Aircraft, Wing"].to(ureg.lbf).magnitude,
     color="black", linewidth=1.5, linestyle="-", marker='s', markersize=8,
     label='GPKit model (fixed engine)')
plt.plot(GA_solution_rubberEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude,
    GA_solution_rubberEngine["freevariables"]["W_Aircraft, Wing"].to(ureg.lbf).magnitude,
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
plt.xlim(xmin = np.min(GA_solution_fixedEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude),
     xmax = np.max(GA_solution_fixedEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude))
plt.ylim(ymin = 0)


#Wingspan vs. Mission Range
plt.subplot(2,3,4)
plt.plot(GA_solution_fixedEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude,
	GA_solution_fixedEngine["freevariables"]["b_Aircraft, Wing"].to(ureg.ft).magnitude,
     color="black", linewidth=1.5, linestyle="-", marker='s', markersize=8,
     label='GPKit model (fixed engine)')
plt.plot(GA_solution_rubberEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude,
    GA_solution_rubberEngine["freevariables"]["b_Aircraft, Wing"].to(ureg.ft).magnitude,
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
plt.xlim(xmin = np.min(GA_solution_fixedEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude),
     xmax = np.max(GA_solution_fixedEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude))
plt.ylim(ymin = 0)


#Aspect Ratio vs. Mission Range
plt.subplot(2,3,5)
plt.plot(GA_solution_fixedEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude,
	GA_solution_fixedEngine["freevariables"]["AR_Aircraft, Wing"],
     color="black", linewidth=1.5, linestyle="-", marker='s', markersize=8,
     label='GPKit model (fixed engine)') 
plt.plot(GA_solution_rubberEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude,
    GA_solution_rubberEngine["freevariables"]["AR_Aircraft, Wing"],
     color="black", linewidth=1.5, linestyle="-", marker='^', markersize=8,
     label='GPKit model (rubber engine)')
plt.plot(cessna_402_data["range"].to(ureg.nautical_mile).magnitude,
	cessna_402_data["AR"], 
	marker='o',color="black", markersize=12, label="Cessna 402")		 	 
plt.grid()
plt.xlabel('Mission range (nm)', fontsize = 16)
plt.ylabel('Aspect ratio (dimensionless)', fontsize = 16)
plt.title("Aspect Ratio",fontsize = 20)
plt.legend(numpoints=1, loc='lower right', fontsize = 12)
plt.xlim(xmin = np.min(GA_solution_fixedEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude),
     xmax = np.max(GA_solution_fixedEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude))
plt.ylim(ymin = 0)


#Engine Power vs. Mission Range
plt.subplot(2,3,6)
plt.plot(GA_solution_rubberEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude,
    GA_solution_rubberEngine["freevariables"]["P_Aircraft, Engines"].to(ureg.horsepower).magnitude,
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
plt.xlim(xmin = np.min(GA_solution_rubberEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude),
     xmax = np.max(GA_solution_rubberEngine["sweepvariables"]["Range"].to(ureg.nautical_mile).magnitude))
plt.ylim(ymin = 0)


title_str = "Key Design Variables vs. Mission Range for a Twin-Engined GA Aircraft"
plt.suptitle(title_str,fontsize = 24)

plt.tight_layout()#makes sure subplots are spaced neatly
plt.subplots_adjust(left=0.125,right=0.9,bottom=0.1,top = 0.9)#adds space at the top for the title