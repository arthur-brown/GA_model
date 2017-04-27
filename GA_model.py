# Optimize a gasoline-powered general-aviation aircraft, representative
# of the Cessna 402

import math
import numpy as np
from gpkit import Variable, VectorVariable, Model, Vectorize
from matplotlib import pyplot as plt
import pint
from GA_aircraft_models import *

# Set up units
ureg = pint.UnitRegistry()

takeoff_distance_ft = 1763
range_nm = np.linspace(1000,2000,10)
N = 5 #number of cruise segments
stall_speed_kts = 78
cruise_speed_kts = 190

GA_aircraft_fixedEngine = Aircraft(engineType="fixed")
GA_aircraft_rubberEngine = Aircraft(engineType="rubber")

GA_mission_fixedEngine = Mission(GA_aircraft_fixedEngine,missionLength_nm=range_nm,numCruiseSegments=N,cruiseSpeed_kts=cruise_speed_kts)
GA_mission_rubberEngine = Mission(GA_aircraft_rubberEngine,missionLength_nm=range_nm,numCruiseSegments=N,cruiseSpeed_kts=cruise_speed_kts)

GA_takeoffConstraint_fixedEngine = TakeoffDistance(GA_aircraft_fixedEngine,
    s_TO_ft=takeoff_distance_ft)
GA_takeoffConstraint_rubberEngine = TakeoffDistance(GA_aircraft_rubberEngine,
    s_TO_ft=takeoff_distance_ft)

GA_stallSpeedConstraint_fixedEngine = StallSpeed(GA_aircraft_fixedEngine,Vstall_kts=stall_speed_kts)
GA_stallSpeedConstraint_rubberEngine = StallSpeed(GA_aircraft_rubberEngine,Vstall_kts=stall_speed_kts)

GA_OEI_ClimbConstraint_fixedEngine = OEI_ClimbConstraint(GA_aircraft_fixedEngine)
GA_OEI_ClimbConstraint_rubberEngine = OEI_ClimbConstraint(GA_aircraft_rubberEngine)

constraints_fixedEngine = [GA_aircraft_fixedEngine, GA_mission_fixedEngine,
                            GA_takeoffConstraint_fixedEngine, 
                            GA_stallSpeedConstraint_fixedEngine,
                            GA_OEI_ClimbConstraint_fixedEngine]

constraints_rubberEngine = [GA_aircraft_rubberEngine, GA_mission_rubberEngine,
                            GA_takeoffConstraint_rubberEngine, 
                            GA_stallSpeedConstraint_rubberEngine,
                            GA_OEI_ClimbConstraint_rubberEngine]

GA_model_fixedEngine = Model(GA_aircraft_fixedEngine.W_TO,constraints_fixedEngine)
GA_model_rubberEngine = Model(GA_aircraft_rubberEngine.W_TO,constraints_rubberEngine)

GA_solution_fixedEngine = GA_model_fixedEngine.solve(verbosity=0)
GA_solution_rubberEngine = GA_model_rubberEngine.solve(verbosity=0)

#print GA_solution_fixedEngine.summary()
#print GA_solution_rubberEngine.summary()


#GA_mission.debug()

# Plotting commands
plt.ion()

# Cessna 402 data

P_available = GA_aircraft_fixedEngine.engines["P_fixedEngine"].value \
     * GA_aircraft_fixedEngine.engines["num_engines"].value

cessna_402_data = {'W_TO':7210*ureg.lbf,
	"W_fuel":1278*ureg.lbf,
    "W_wing":860*ureg.lbf,#Cessna 404 (not 402)
	"b":13.44*ureg.m,
	"AR":8.6,
	"range":1234*ureg.nautical_mile,
    "P":P_available}

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
plt.ylim(ymin = 0, ymax = 800)


title_str = "Key Design Variables vs. Mission Range for a Twin-Engined GA Aircraft"
plt.suptitle(title_str,fontsize = 24)

plt.tight_layout()#makes sure subplots are spaced neatly
plt.subplots_adjust(left=0.125,right=0.9,bottom=0.1,top = 0.9)#adds space at the top for the title
