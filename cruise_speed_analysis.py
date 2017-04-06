# Cruise speed analysis for the Cessna 310 and Embraer 190

import pint
import math
import numpy as np
from matplotlib import pyplot as plt
import pint

# Set up units
ureg = pint.UnitRegistry()

def max_range(W,rho,S,CD0,K,aircraftType="jet"):
	#Computes the cruising speed and lift coefficient of an aircraft for maximum range
	
	#If aircraft is a jet, assume TSFC is proportional to sqrt(M)
		#Means that Range is proportional to sqrt(V)*L/D
	
	#If aircraft is prop driven, assume that PSFC is constant
		#Means that range is proportional to L/D

	if aircraftType == "jet":
		CL = np.sqrt((3*CD0)/(5*K))

	if aircraftType == "prop":
		CL = np.sqrt(CD0/K)

	V = np.sqrt((2*W)/(rho*S*CL))

	return V, CL

if __name__ == "__main__":
	
	#Atmospheric data
	ureg.define('slug = 14.5939 * kg')
	slug_per_ft3 = ureg.slug / (ureg.ft)**3
	ft_per_s = ureg.ft/ureg.s
	
	rho = {"h = 0":0.00237*slug_per_ft3,
		"h = 7500 ft":0.00189*slug_per_ft3,
		"h = 10000 ft":0.00175*slug_per_ft3,
		"h = 36000 ft":0.00070*slug_per_ft3}

	a = {"h = 0":1116.449*ft_per_s,
		"h = 7500 ft":1087.282*ft_per_s,
		"h = 10000 ft":1077.384*ft_per_s,
		"h = 36000 ft":968.4699*ft_per_s}

	#Aircraft data
	square_foot = (ureg.ft)**2
	square_meter = (ureg.m)**2

	cessna_310 = {"W_DG":4830*ureg.lbf,
		"S":175*square_foot,
		"e":0.817*ureg.dimensionless,
		"AR":7.0*ureg.dimensionless,
		"CD0":0.03*ureg.dimensionless}
	cessna_310["K"] = 1/(math.pi * cessna_310["e"] * cessna_310["AR"])

	E190 = {"W_DG":105822*ureg.lbf,
		"S":102.25*square_meter,
		"e":0.8547*ureg.dimensionless,
		"AR":7.6839*ureg.dimensionless,
		"CL":0.675*ureg.dimensionless,#lift coefficient at which CD_wing and CD_nonwing were computed
		"CD_wing":0.0295*ureg.dimensionless,#includes induced drag
		"CD_nonwing":0.0125*ureg.dimensionless}
	E190["S"] = E190["S"].to(square_foot)

	#separate the E190 wing drag into zero-lift and induced components
	E190["K"] = 1/(math.pi * E190["e"] * E190["AR"])
	E190["CDi"] = E190["K"] * E190["CL"]**2
	E190["CD0_wing"] = E190["CD_wing"] - E190["CDi"]
	E190["CD0"] = E190["CD0_wing"] + E190["CD_nonwing"]

	#Computations for the Embraer 190 (h = 36,000 ft)
	E190["V_maxRange"], E190["CL_maxRange"] = max_range(W=E190["W_DG"],\
		rho=rho["h = 36000 ft"],S=E190["S"],CD0=E190["CD0"],K=E190["K"],aircraftType="jet")
	E190["V_maxRange"] = E190["V_maxRange"].to(ft_per_s)
	E190["M_maxRange"] = E190["V_maxRange"] / a["h = 36000 ft"]

	#Computations for the Cessna 310 (h = 7,500 ft)
	cessna_310["V_maxRange"], cessna_310["CL_maxRange"] = max_range(W=cessna_310["W_DG"],\
		rho=rho["h = 7500 ft"],S=cessna_310["S"],CD0=cessna_310["CD0"],K=cessna_310["K"],aircraftType="prop")
	cessna_310["V_maxRange"] = cessna_310["V_maxRange"].to(ft_per_s)
	cessna_310["M_maxRange"] = cessna_310["V_maxRange"] / a["h = 7500 ft"]

	
	#Print results to command line
	print ""
	print "Analysis of the Embraer 190 (h = 36,000 ft)"
	print ""
	print "W_DG: %0.0f lbs" % E190["W_DG"].to(ureg.lbf).magnitude
	print "S: %0.0f ft^2" % E190["S"].to(square_foot).magnitude
	print "CD0: %0.4f" % E190["CD0"]
	print "K: %0.4f" % E190["K"]
	print "CL for maximum range: %0.4f" % E190["CL_maxRange"]
	print "V for maximum range: %0.0f ft/s (%0.0f knots)" % \
		(E190["V_maxRange"].to(ft_per_s).magnitude, E190["V_maxRange"].to(ureg.knot).magnitude)
	print "M for maximum range: %0.4f" % E190["M_maxRange"].magnitude

	print ""
	print "Analysis of the Cessna 310 (h = 7,500 ft)"
	print ""
	print "W_DG: %0.0f lbs" % cessna_310["W_DG"].to(ureg.lbf).magnitude
	print "S: %0.0f ft^2" % cessna_310["S"].to(square_foot).magnitude
	print "CD0: %0.4f" % cessna_310["CD0"]
	print "K: %0.4f" % cessna_310["K"]
	print "CL for maximum range: %0.4f" % cessna_310["CL_maxRange"]
	print "V for maximum range: %0.1f ft/s (%0.1f knots)" % \
		(cessna_310["V_maxRange"].to(ft_per_s).magnitude, cessna_310["V_maxRange"].to(ureg.knot).magnitude)
	print "M for maximum range: %0.4f" % cessna_310["M_maxRange"].magnitude
	

	

	#Save results to file
	output_data = open("cruise_speed_analysis_data.txt","w")

	output_data.write("Analysis of the Embraer 190 (h = 36,000 ft)\n\n")
	for key in E190.keys():
		output_data.write("\t%s: %0.4f %s\n" % (key, E190[key].magnitude, E190[key].units))

	output_data.write("\nAnalysis of the Cessna 310 (h = 7,500 ft)\n\n")
	for key in cessna_310.keys():
		output_data.write("\t%s: %0.4f %s\n" % (key, cessna_310[key].magnitude, cessna_310[key].units))

	output_data.close()
	
