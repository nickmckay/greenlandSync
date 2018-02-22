#Load relevant libraries
library(lipdR)
library(geoChronR)
library(here)

#and functions
source(here("syncFunctions.R"))


#Load in data
D = readLipd(here("data"))
TS = extractTs(D)


# Sea ice records:
#   CC04.Gibb.2015: sea ice in the labrador sea. 14C age model. Standard marine reservoir correction (Delta R is zero)
# CC70.Gibb.2015: sea ice in southern Baffin Bay. 14C age model. Standard marine reservoir correction (Delta R is zero)
# LS009.Ledu.2010: sea ice in Lancaster Sound.  14C age model. Delta R is 400. 14C only covers the early Holocene (so we’re fine for this study). I couldn’t figure out how to add extra “calibrated age” columns, which will eventually include the paleomagnetic tie points. I didn’t bother since these don’t affect the Early Holocene. Let me know if you’d like that info.
# It’d be great to compare each of the sea ice records with our aquatic leaf wax d2H, epsilon, and temperature records.
# 
# 




# Ice core records:
#   GISP2.Kobashi.2017: temp, age. Maximum counting error for the Early Holocene is 2% (Vintner et al. 2006, JGR)
# GISP2.Stuiver.1995: d18O, age, depth. Would be nice to figure out how to put this on the same age model as Kobashi, although Kobashi doesn’t give us depth, so I’m not totally sure how to approach that…maybe there’s a version of this record in Iso2k on the GICC05 age model, which is what Kobashi used?
# NEEM.Rasmussen.2013: accumulation rate, age & depth (plus accumulation uncertainty)
# It’d be nice to compare the GISP2 temperature data to our temperature record.
# Also would be nice to compare NEEM accumulation to our epsilon record
# All three of these ice core records are on (or in the case of Stuiver 1995 should be on) the GICC05 scale. For this age model, maximum counting error for the Early Holocene is 2% (Vinther et al. 2006, JGR)