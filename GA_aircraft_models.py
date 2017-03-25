#Component models for a twin-engined general-aviation aircraft
#Representative of the Cessna 402


import math
import numpy as np
from gpkit import Variable, VectorVariable, Model

class Aircraft(Model): #Overall vehicle model
          
    def dynamic(self,state,segmentRange):
        return AircraftPerformance(self,state,segmentRange)

    def setup(self, engineType="fixed"):
        W_ZF = Variable("W_ZF","lbf","Zero-fuel weight")
        W_fuel = Variable("W_fuel","lbf","Fuel weight")
        W_TO = Variable("W_TO","lbf","Takeoff weight")
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
        
        mission_range = Variable("Range",np.linspace(1000,2000,10),"nautical_mile","Mission range")#sweep
        fuel_allowance = 0.03+0.015+0.005+0.06 #allowance for takeoff, climb, landing, reserves, and trapped fuel

        fs = FlightSegment(aircraft,segmentRange = mission_range)
        Wburn = fs.aircraftp["W_{burn}"]
        Wfuel = fs.aircraftp["W_{fuel}"]

        constraints += [Wfuel >= (1+fuel_allowance)*Wburn,
                        self.aircraft["W_fuel"] >= Wfuel,
                        fs]

        return constraints

class FlightSegment(Model):
    "Combines a flight state and an aircraft"

    def setup(self,aircraft,segmentRange):

        constraints = []

        self.flightstate = FlightState()
        self.aircraftp = aircraft.dynamic(self.flightstate,segmentRange)

        constraints += [self.flightstate, self.aircraftp]
        return constraints


class AircraftPerformance(Model):
    "Aircraft performance model. Currently cruise only."

    def setup(self,aircraft,state,segmentRange):
        self.aircraft = aircraft
        self.aerodynamics = Aerodynamics(self.aircraft) 

        Wfuel = Variable("W_{fuel}","lbf","Fuel weight (at beginning of segment)")
        Wburn = Variable("W_{burn}","lbf","Fuel burned during segment")
        CL = self.aerodynamics["C_L"]
        CD = self.aerodynamics["C_D"]
        
        g = Variable("g",9.807,"m/s**2","Gravitational constant")

        constraints = []
        constraints += [self.aerodynamics]

        #lift >= weight
        constraints += [0.5*state["rho"]*state["V"]**2 * aircraft.wing["S"]*CL >= aircraft["W_ZF"] + Wfuel]

        #Breguet range (Taylor)
        z = (g * segmentRange * aircraft.engines["SFC"] * CD) / (aircraft.engines["eta_prop"] * CL)
        constraints += [Wburn >= (z + z**2/np.math.factorial(2)
            + z**3/np.math.factorial(3)
            + z**4/np.math.factorial(4)) * (aircraft["W_ZF"]+Wfuel)] 
        
        constraints += [self.aerodynamics]
        return constraints

class Aerodynamics(Model):
    "Aerodynamic variables and constraints"

    def setup(self,aircraft,stateType="cruise"):

        self.aircraft = aircraft

        c_fe = Variable("c_fe",0.0045,"-","Equivalent skin-friction coefficient (twin-engined GA)")
        e = Variable("e",0.8,"-","Oswald efficiency factor")
        CLmax_cruise = Variable("CLmax_cruise",1.6,"-","Maximum lift coefficient (cruise)")
        CLmax_TO = Variable("CLmax_TO",2.0,"-","Maximum lift coefficient (takeoff)")
        CLmax_LDG = Variable("CLmax_LDG",2.6,"-","Maximum lift coefficient (landing)")

        CL = Variable("C_L","-","Lift coefficient")
        CD = Variable("C_D","-","Drag coefficient")
        CD0 = Variable ("C_{D_0}","-","Zero-lift drag coefficient")
        CDi = Variable ("C_{D_i}","-","Induced drag coefficient")
        
        constraints = []
        constraints += [CDi == (CL**2) / (np.pi*e*aircraft.wing["AR"])]

        if stateType == "cruise":
            #lift
            constraints += [CL <= CLmax_cruise]

            #drag
            constraints += [CD0 == c_fe * aircraft.Swet / aircraft.wing["S"],
                            CD >= CD0 + CDi]
        
        if stateType == "takeoff":
            #lift
            constraints += [CL == CLmax_TO] #not actually used

            #drag
            constraints += [CD0 >= c_fe * aircraft.Swet / aircraft.wing["S"] + aircraft.landing_gear["CDA_gear"]/aircraft.wing["S"],
                            CD >= CD0]#induced drag not counted; lift assumed 0 prior to rotation
        
        return constraints


class FlightState(Model):
    
    "Context for flight physics"
    def setup(self,stateType="cruise"):
        
        constraints = []

        g = Variable("g",9.807,"m/s**2","Gravitational constant")
        constraints += [g == g]

        if stateType == "cruise":
            V = Variable("V",180,"knots","Airspeed")
            h = Variable("h",8000,"ft","Altitude")
            rho = Variable("rho",0.96287,"kg/m^3","Air density")
            constraints += [V == V, h == h, rho == rho]

        if stateType == "takeoff":
            h = Variable("h",0,"ft","Altitude")
            rho = Variable("rho",1.225,"kg/m^3","Air density")
            constraints += [h == h, rho == rho]

        return constraints 

class TakeoffConstraint(Model):

    def setup(self,aircraft,state=FlightState(stateType="takeoff")):
        self.aircraft = aircraft
        self.state = state
        self.aerodynamics = Aerodynamics(aircraft,stateType="takeoff")
        
        CLmax_TO = self.aerodynamics["CLmax_TO"]
        CD = self.aerodynamics["C_D"]
        g = self.state["g"]

        y = Variable("y","-","Takeoff-distance parameter 1")
        zeta = Variable("zeta","-","Takeoff-distance parameter 2")
        s_TO = Variable("s_TO",1763,"ft","Takeoff distance")

        V_TO = 1.1*((2*aircraft.W_TO) / (state["rho"]*aircraft.wing["S"]*CLmax_TO))**(1./2)#1.1*V_stall
        D_TO = 0.5 * state["rho"] * V_TO**2 * aircraft.wing["S"] * CD
        T_TO = aircraft.engines["eta_prop"]*aircraft.engines["num_engines"]*aircraft.engines["P"]/V_TO

        constraints = []

        constraints += [zeta == D_TO / T_TO,
            2*g*s_TO*T_TO / (aircraft.W_TO*V_TO**2) >= 1 + y,
            1 >= 0.0464*zeta**2.73/y**2.88 + 1.044*zeta**0.296/y**0.049]
        
        constraints += [self.state, self.aerodynamics]

        return constraints