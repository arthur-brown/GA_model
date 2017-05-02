# Optimize a gasoline-powered general-aviation aircraft, representative
# of the Cessna 402. Produces a plot of 

import math
import numpy as np
from gpkit import Variable, VectorVariable, Model, Vectorize
from matplotlib import pyplot as plt
import pint
from GA_aircraft_models import *

# Set up units
ureg = pint.UnitRegistry()

takeoff_distance_ft = 1763
range_nm = 957
N = 5 #number of cruise segments
stall_speed_kts = 78
cruise_speed_kts = np.linspace(130,180,10)

GA_aircraft_fixedEngine = Aircraft(engineType="fixed")
GA_aircraft_rubberEngine = Aircraft(engineType="rubber")

W_TO_fixedEngine = np.zeros(np.size(cruise_speed_kts))
W_TO_rubberEngine = np.zeros(np.size(cruise_speed_kts))

W_fuel_fixedEngine = np.zeros(np.size(cruise_speed_kts))
W_fuel_rubberEngine = np.zeros(np.size(cruise_speed_kts))

for i,V_cruise in enumerate(cruise_speed_kts):

    GA_mission_fixedEngine = Mission(GA_aircraft_fixedEngine,missionLength_nm=range_nm,
        numCruiseSegments=N,cruiseSpeed_kts=V_cruise)
    GA_mission_rubberEngine = Mission(GA_aircraft_rubberEngine,missionLength_nm=range_nm,
        numCruiseSegments=N,cruiseSpeed_kts=V_cruise)

    GA_takeoffConstraint_fixedEngine = TakeoffDistance(GA_aircraft_fixedEngine,
        s_TO_ft=takeoff_distance_ft)
    GA_takeoffConstraint_rubberEngine = TakeoffDistance(GA_aircraft_rubberEngine,
        s_TO_ft=takeoff_distance_ft)

    stallSpeed_fuelFraction = 0.5 #50% cruise
    GA_stallSpeedConstraint_fixedEngine = StallSpeed(GA_aircraft_fixedEngine,
        Vstall_kts=stall_speed_kts,fuel_fraction=stallSpeed_fuelFraction)
    GA_stallSpeedConstraint_rubberEngine = StallSpeed(GA_aircraft_rubberEngine,
        Vstall_kts=stall_speed_kts,fuel_fraction=stallSpeed_fuelFraction)

    gamma_req = 0.01 #regulations say "demonstratably positive"
    GA_OEI_ClimbConstraint_fixedEngine = ClimbConstraint(GA_aircraft_fixedEngine,
        constraintType="ClimbGradient",velocityType="1.3Vstall",flightStateType="sea_level",
        aeroStateType="OEI_climb",gamma_req=gamma_req,fuel_fraction=1)
    GA_OEI_ClimbConstraint_rubberEngine = ClimbConstraint(GA_aircraft_rubberEngine,
        constraintType="ClimbGradient",velocityType="1.3Vstall",flightStateType="sea_level",
        aeroStateType="OEI_climb",gamma_req=gamma_req,fuel_fraction=1)

    V_v_req_serviceCeiling_fpm = 100 #definition of service ceiling
    GA_serviceCeilingConstraint_fixedEngine = ClimbConstraint(GA_aircraft_fixedEngine,
        constraintType="ClimbRate",velocityType="unconstrained",flightStateType="service_ceiling",
        aeroStateType="cruise",V_v_req_fpm=V_v_req_serviceCeiling_fpm,fuel_fraction=1)
    GA_serviceCeilingConstraint_rubberEngine = ClimbConstraint(GA_aircraft_rubberEngine,
        constraintType="ClimbRate",velocityType="unconstrained",flightStateType="service_ceiling",
        aeroStateType="cruise",V_v_req_fpm=V_v_req_serviceCeiling_fpm,fuel_fraction=1)

    V_v_req_seaLevel_fpm = 1450 
    GA_seaLevelClimbRateConstraint_fixedEngine = ClimbConstraint(GA_aircraft_fixedEngine,
        constraintType="ClimbRate",velocityType="unconstrained",flightStateType="sea_level",
        aeroStateType="cruise",V_v_req_fpm=V_v_req_seaLevel_fpm,fuel_fraction=1)
    GA_seaLevelClimbRateConstraint_rubberEngine = ClimbConstraint(GA_aircraft_rubberEngine,
        constraintType="ClimbRate",velocityType="unconstrained",flightStateType="sea_level",
        aeroStateType="cruise",V_v_req_fpm=V_v_req_seaLevel_fpm,fuel_fraction=1)

    constraints_fixedEngine = [GA_aircraft_fixedEngine, GA_mission_fixedEngine,
                            GA_takeoffConstraint_fixedEngine, 
                            GA_stallSpeedConstraint_fixedEngine,
                            GA_OEI_ClimbConstraint_fixedEngine,
                            GA_serviceCeilingConstraint_fixedEngine]

    constraints_rubberEngine = [GA_aircraft_rubberEngine, GA_mission_rubberEngine,
                            GA_takeoffConstraint_rubberEngine, 
                            GA_stallSpeedConstraint_rubberEngine,
                            GA_OEI_ClimbConstraint_rubberEngine,
                            GA_serviceCeilingConstraint_rubberEngine]

    GA_model_fixedEngine = Model(GA_aircraft_fixedEngine.W_TO,constraints_fixedEngine)
    GA_model_rubberEngine = Model(GA_aircraft_rubberEngine.W_TO,constraints_rubberEngine)

    GA_solution_fixedEngine = GA_model_fixedEngine.solve(verbosity=0)
    GA_solution_rubberEngine = GA_model_rubberEngine.solve(verbosity=0)

    W_TO_fixedEngine[i] = GA_solution_fixedEngine["variables"]["W_TO"].magnitude
    W_TO_rubberEngine[i] = GA_solution_rubberEngine["variables"]["W_TO"].magnitude
    
    W_fuel_fixedEngine[i] = GA_solution_fixedEngine["variables"]["W_fuel"].magnitude
    W_fuel_rubberEngine[i] = GA_solution_rubberEngine["variables"]["W_fuel"].magnitude


#print GA_solution_fixedEngine.summary()
#print GA_solution_rubberEngine.summary()


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

fig1 = plt.figure(figsize=(14, 8), dpi=80)
plt.show()


#MTOW vs. Cruising Speed
plt.subplot(1,2,1)
plt.plot(cruise_speed_kts,W_TO_fixedEngine,
    color="black", linewidth=1.5, linestyle="-", marker='s', markersize=8,
    label='GPKit model (fixed engine)')
plt.plot(cruise_speed_kts,W_TO_rubberEngine,
     color="black", linewidth=1.5, linestyle="-", marker='^', markersize=8,
     label='GPKit model (rubber engine)')
#plt.plot(cessna_402_data["range"].to(ureg.nautical_mile).magnitude,
#    cessna_402_data["W_TO"].to(ureg.lbf).magnitude, 
#    marker='o',color="black", markersize=12, label="Cessna 402")         
plt.grid()
plt.xlabel('Cruising speed (knots)', fontsize = 16)
plt.ylabel('MTOW (lbs)', fontsize = 16)
plt.title("Maximum Takeoff Weight vs. Cruising Speed",fontsize = 20)
plt.legend(numpoints = 1,loc='lower left', fontsize = 12)
#plt.xlim()
plt.ylim(ymin = 0, ymax = 9000)

#Fuel Weight vs. Cruising Speed
plt.subplot(1,2,2)
plt.plot(cruise_speed_kts,W_fuel_fixedEngine,
    color="black", linewidth=1.5, linestyle="-", marker='s', markersize=8,
    label='GPKit model (fixed engine)')
plt.plot(cruise_speed_kts,W_fuel_rubberEngine,
     color="black", linewidth=1.5, linestyle="-", marker='^', markersize=8,
     label='GPKit model (rubber engine)')
#plt.plot(cessna_402_data["range"].to(ureg.nautical_mile).magnitude,
#    cessna_402_data["W_TO"].to(ureg.lbf).magnitude, 
#    marker='o',color="black", markersize=12, label="Cessna 402")         
plt.grid()
plt.xlabel('Cruising speed (knots)', fontsize = 16)
plt.ylabel('Fuel weight (lbs)', fontsize = 16)
plt.title("Fuel Weight vs. Cruising Speed",fontsize = 20)
plt.legend(numpoints = 1,loc='lower left', fontsize = 12)
#plt.xlim(xmin = np.min())
plt.ylim(ymin = 0)


'''
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
'''