# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 13:24:13 2021

@author: JDiaz
"""

######################################################################################################################
######################################################################################################################
######################################################################################################################
########## This code is a working draft, not for publication!!########################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################


#Let's delete variables so we can free up some RAM
#for element in dir():
#    if element[0:2] != "__":
#        del globals()[element]

import os
import numpy as np
import pandas as pd
import datetime
import glob
#import arcpy
#arcpy.CheckOutExtension("Spatial")
#from arcpy.sa import *

#arcpy.env.addOutputsToMap = 0

#arcpy.env.overwriteOutput = True



###Let's set the working directory
os.chdir(r"\\hydro-nas\Team\Projects\5635_SFEI Carbon and GHG\Accretion\20220408")

#if not arcpy.Exists(os.path.join(os.getcwd(),"SEDCALC.gdb")):
    #Let's create geodatabase
    #arcpy.CreateFileGDB_management(os.getcwd(), "SEDCALC.gdb")

#arcpy.env.workspace = os.path.join(os.getcwd(),"SEDCALC.gdb")

###################################################SEDCALC###################
#############################################################################
#############################################################################



def rtemin(y,
           tdrnge=300,
           mw=0):
    
    # *************************************************************
    # *                       RTEMIN                              *
    # * Subroutine for determining the rate of mineral sed input  *
    # *************************************************************
    #
    #
    # this section calculates the amount of mineral sediment 
    # input each year - based on the relative elevation of the 
    # marsh surface.  Y is the relative elevation of 
    # the core at at given time.
    # tdrnge is the tidal range in meters
    # mhw is the relative elevation of mhw in meters (0 in this case)
    #
   tdhght = (y-mw)/(tdrnge/2.0) 
   if tdhght<=0:
       rtemin = 1.0
   else:
      rtemin = 1 - (min(tdhght, 1.0)) 
   return rtemin


def poresp(c,k1):
    # *****************************************************
    # *            PORESP                                 *
    # * Subroutine for determining changes in pore space  *
    # *****************************************************
    #
    # the next section calculates the pore space for each section.
    # This is where changes due to compaction occur. I am assuming that
    # all of the pore spaces are filled with water, and that any compaction 
    # is due to the loss of water and decrease in pore space volume.
    # Pore space is assumed to be a function of the amount of material 
    # (both organic and mineral) that is in a given section, as well
    # as the mass above that particular section.
    #
    # a = totorg(t2)
    # b = min(t2)
    # c = densabv(t2)
    # d = oldorg(t2)
    # e = min(t2)
    #
    # k1 is a constant that affects the curve for compaction
    # k2 affects the relative importance of organic versus mineral
    # matter in determining pore space.
    # K2 > 1 - organic matter more important
    # k2 < 1 - organic matter less important
    # k2 = 1 - organic and mineral matter the same
    #
    #
    # p1 & p2 are just temporary variable to make calculations easier.
    #
    poresp=1-(c/(k1+c))
    #	k2 = 0.1
    #	poresp = (1/(1+(k1*c)))
    #
    # everything below here has been commented out in order to  
    # "simplify" the calculation of pore space.
    #
    #  	write(6,2005) a,b,c,d,e
    #  	if (((k2*d)+e).le.0) then
    #		p2 = 1
    #	else
    #		p2 = sqrt(((k2*a)+b)/((k2*d)+e))
    #	end if
    #	poresp = p1*p2	     

    #     Last change:  SJD  23 Aug 2007    4:35 pm
    return poresp

def rtprod(d):
    #c **********************************************
    #c *  RTPROD                                    *
    #c *  Root production subroutine                *
    #c **********************************************
    #c
    #c
    #c  What follows is a subroutine for determining organic
    #c  production at various time / depths
    #c  1/11/94 - I am changing this so root production decreases \
    #c  exponentially with depth.
    #c  sed2.for has the old version of root production.
    #c
    #c
    #c parameters:
        #c depth(t2) is the only parameter - it is passed 
        #c as the single variable "d"
        #c 
        #c variables:
            #c	undpro - total underground production (g/cm^3)
            #c	kdist - controls the decay of the root production function
            #c    
    undpro = 0.06
    kdist = 0.40
    rtprod =np.exp(-kdist*d)*kdist*undpro
    
    #c
    #c
    #c this section calculates root production at a particular depth
    #c based on the parameters/curve that are designated above.
    #c
    #c
    

    # this is sed5.for  -  sediment accretion program for fortran
    # this version uses an exponentially decreasing underground production function
    # this change is in the subroutine - rtprod
    #
    # additionally - this version uses a new decomposition model - still 3
    # different rates - but the rates for each year class are from an
    # exponential decay curve.
    #
    # this version (4/20/2007) reads organic and mineral inputs from a file
    # declaring variables
    #
    #      real*8 org, min, orgden, minden, h2oden, pore, minin, orgin
    #      real*8 orgbd, minbd, bulkd, porg, depth, massabv, intelv, relelv
    #      real*8 slr, subsid, totorg, totvol, orgvol, minvol, densabv
    #      real*8 dt, mindev, porelim, refrac,h2oin
    #      real*8 pctpore, tpore
    #      real*8 acc1000, org1000, min1000, acc4000, org4000, min4000
    #      real*8 acc4900, org4900, min4900, finelv
    #      integer time, t2, endtim
    #      dimension org(0:7000,4), min(0:7000), pore(0:7000), totvol(7000)
    #      dimension orgbd(7000), minbd(7000), porg(7000), minvol(7000)
    #      dimension depth(0:7000),  massabv(0:7000), densabv(0:7000)
    #      dimension bulkd(7000),relelv(0:7000),totorg(0:7000),orgvol(7000)
    #      dimension pctpore(0:7000)
    #      dimension tpore(0:7000)
    #      dimension minin(7000)
    #      dimension orgin(7000)
    #      dimension h2oin(7000)
    #      dimension porelim(7000)
    #c
    return rtprod


def decmp1(d,kdec1):
    #c
    #c *************************************************
    #c *           DECMP1                              *
    #c * Decomposition of "youngest" organic material  *
    #c ************************************************* 
    #c
    #c
    #c the next section is the decomposition subroutine FOR 1st year org matter
    #c it is a function that gives a decomposition rate (from 0 to 1)
    #c based on the depth of each section. 
    #c
    #c decomp is a RATE (units g lost/g present) so it has 
    #c to be multiplied by the organic mass (org) of each section 
    #c
    #c as with the production function the only thing that determines 
    #c this is the depth. 
    #c
    #c Variables:
        #c 	mx1 - maximum rate of decay for this age class
        #c	kdec1 - k for exponential decay curve for this 
        #c 		age class decompostion curve
    mx1 = 0.92

    #	decmp1 = (exp(-kdec1*d))*mx1
    decmp1=np.exp(-kdec1*d)*mx1
    return decmp1
    

def decmp2(d):
    #c
    #c
    #c
    #c *************************************************
    #c *           DECMP2                              *
    #c * Decomposition of "medium" organic material    *
    #c ************************************************* 
    #c
    #c
    #c
    #c the next section is the decomposition subroutine FOR 2nd year org matter
    #c it is a function that gives a decomposition rate (from 0 to 1)
    #c based on the depth of each section. 
    #c
    #c decomp is a RATE (units g lost/g present) so it has 
    #c to be multiplied by the organic mass (org) of each section 
    #c
    #	real*8 function decmp2(d)
    #	real*8 mx2, kdec2, d
    #c
    #c as with the production function the only thing that determines 
    #c this is the depth. 
    #c
    #c Variables:
        #c 	mx2 - maximum rate of decay for this age class
        #c	kdec2 - k for exponential decay curve for this 
        #c 		age class decompostion curve
    mx2 = 0.37
    kdec2 = 0.57
    decmp2 = (np.exp(-kdec2*d))*mx2
    return decmp2
 
def	decmp3(d):
    #c
    #c
    #c
    #c *************************************************
    #c *           DECMP3                              *
    #c * Decomposition of "oldest" organic material    *
    #c ************************************************* 
    #c
    #c
    #c
    #c the next section is the decomposition subroutine FOR old org matter
    #c it is a function that gives a decomposition rate (from 0 to 1)
    #c based on the depth of each section. 
    #c
    #c The rates are LOWEST for this group of organic material
    #c
    #c decomp is a RATE (units g lost/g present) so it has 
    #c to be multiplied by the organic mass (org) of each section 
    #c
    #	real*8 mx3, kdec3, d
    #c
    #c as with the production function the only thing that determines 
    #c this is the depth. 
    #c
    #c Variables:
        #c 	mx3 - maximum rate of decay for this age class
        #c	kdec3 - k for exponential decay curve for this 
        #c 		age class decompostion curve
    mx3 = 0.16
    kdec3 = 0.1
    #c
    decmp3 = (np.exp(-kdec3*d))*mx3
    #	end
    return decmp3
    

def sedcalc(endtim = 84,
            h2oden = 1.00,
            h2oin=pd.read_csv("h2oin.csv"), #Initial pore space (fraction)
            #################################This file needs to be in the working
            #################################directory
            intelv = -26,
            k1=2.5,        #Consolidation constant
            mindev = 0.0,
            minden = 2.61, #Mineral particle density (g cm-2)
            minin=pd.read_csv("minin.csv"), ##Surface mineral matter deposition (g/cm2)
            #################################This file needs to be in the working
            #################################directory
            orgden = 1.14, #Organic particle density (g cm-2)
            orgin=pd.read_csv("orgin.csv"),#Surface organic matter deposition (g/cm2)
            #################################This file needs to be in the working
            #################################directory
            porelim=pd.read_csv("porelim.csv"), #Final pore space (fraction)
            #################################This file needs to be in the working
            #################################directory            
            refrac = 0.4,
            slr = 0.0, #Sea level rise (cm yr -1)
            strint = 5.0,
            strdev = 0.3,
            subsid = 0.0, #Subsidence (cm yr -1) from gas withdrawal
            kdec1=0.41
            ):
    
    # h20in is a percent.  It has to be converted to a volume
    # to be useful for calculations.  The conversion from % to volume is:
    # porespace volume = ((%)/(1-%))*(minvol + orgvol)
    #Let's initialize bulkd
    bulkd=np.zeros(endtim+1)
    #Let's initialize densabv
    densabv=np.zeros(endtim+1)

    #Let's initialize depth
    depth=np.zeros(endtim+1)

    #Let's initialize massabv
    massabv=np.zeros(endtim+1)

    #Let's initialize minbd
    minbd=np.zeros(endtim+1)

    #Let's initialize mini (min is a reserved name in Python)
    mini=np.zeros(endtim+1)

    #Let's initialize minvol
    minvol=np.zeros(endtim+1)

    #Let's initialize org
    org=np.zeros((endtim+1,5))

    #Let's initialize orgbd
    orgbd=np.zeros(endtim+1)

    #Let's initialize orgvol
    orgvol=np.zeros(endtim+1)

    #Let's initialize pctpore
    pctpore=np.zeros(endtim+1)

    #Let's initialize pore
    pore=np.zeros(endtim+1)

    #Let's initialize porg
    porg=np.zeros(endtim+1)

    #Let's initialize relelv
    relelv=np.zeros(endtim+1)

    #Let's initialize totorg
    totorg=np.zeros(endtim+1)

    #Let's initialize totvol
    totvol=np.zeros(endtim+1)

    #Let's initialize tpore
    tpore=np.zeros(endtim+1)




    #************************************************
    # this is the beginning of the main control loop *
    # ************************************************

    for time in range(1,endtim+1):

        #c this section moves all values down one section
        #c before the next round of growth, new input and decomposition
        #c

        for t2 in range(time-1,-1,-1):
            org[t2+1,4]=org[t2,4]
            org[t2+1,3]=org[t2,3]+org[t2,2]
            org[t2+1,2]=org[t2,1]
            org[t2+1,1]=0
            mini[t2+1]=mini[t2]
            pctpore[t2+1]=pctpore[t2]
            tpore[t2+1]=tpore[t2]

    #c ************************************************************
    #c * these are the new inputs of material onto the surface of *
    #c * the marsh (into the first position in the array).        *
    #c ************************************************************

        org[1,1]=float(orgin.loc[time-1])*(1-refrac)
        org[1,4]=float(orgin.loc[time-1])*(refrac)
        mini[1]=float(minin.loc[time-1])*rtemin(relelv[time-1])
        tpore[1]=1
        pctpore[1]=float(h2oin.loc[0])
        pore[1]=((float(h2oin.loc[time-1])/(1-float(h2oin.loc[time-1])))*((org[1,1]/orgden)+(mini[1]/minden)))
    
    # ******************************************************************
    # * the following section is where the "yearly" calculations       *
    # * take place.  It combines all of the other calculation sections *
    # * from earlier versions of the model (7/10/93).                  *
    # ******************************************************************
    #
    # this section calculates the volume of each section
    # based on the mass of organic matter, mineral matter, and
    # water.  It will also use a compaction subfunction in
    # the future. Compaction will be a function of the
    # mass that is on top of the current section.
    #
    #
    # this is where the new roots and rhizomes are put into
    # the sediment.  rtprod is a subroutine/function
    # that will determine root production based on depth/time
    #
    # this is also the decomposition section.  Again decomp is
    # a subroutine based on depth/time.
    #

        for t2 in range(1,time+1):

            istep = 10
            dt=1/istep
            for ie in range(1,istep+1):           
                totorg[t2]=org[t2,1]+org[t2,2]+org[t2,3]+org[t2,4]
                massabv[t2]=massabv[t2-1]+totorg[t2-1]+mini[t2-1]+pore[t2-1]   
                if depth[t2-1]==0:
                    densabv[t2]=0
                else:
                    densabv[t2]=massabv[t2]/depth[t2-1]
                orgvol[t2]=totorg[t2]/orgden
                minvol[t2]=mini[t2]/minden
                if t2<=1:
                    pctpore[t2]=pctpore[t2]
                else:
                    dum1=float(h2oin.loc[time-1])
                    dum2=float(porelim.loc[time-1])
                    tpore[t2]=tpore[t2]-(tpore[t2]-tpore[t2]*poresp(densabv[t2],k1))*dt
                    pctpore[t2]=dum2+(dum1-dum2)*tpore[t2]

                pore[t2]=((pctpore[t2]/(1-pctpore[t2]))*(orgvol[t2] + minvol[t2]))

                #c the line above is for running the model without compaction
                #c pore space is constant for all sections.
                #c
                totvol[t2]=orgvol[t2]+minvol[t2]+ pore[t2]
                depth[t2]=depth[t2-1]+totvol[t2]
                porg[t2]=totorg[t2]/(totorg[t2] + mini[t2])
                bulkd[t2]= (totorg[t2]+mini[t2])/totvol[t2]           
                orgbd[t2]=totorg[t2]/totvol[t2]
                minbd[t2] = mini[t2]/totvol[t2]
                org[t2,1]=org[t2,1]+((dt*(rtprod(depth[t2])*totvol[t2]))*(1-refrac)) \
                    -((dt*(((decmp1(depth[t2],kdec1))*org[t2,1]))))
                org[t2,2]=org[t2,2]-\
                    (dt*(((decmp2(depth[t2]))*org[t2,2])))
                org[t2,3]=org[t2,3]\
                    -(dt*((decmp3(depth[t2]))*org[t2,3]))
                org[t2,4]=org[t2,4]+((dt*(rtprod(depth[t2])*totvol[t2]))*refrac)

#c
#c commenting out the 3 "&" lines above, cuts out decompostion
#c
#c           write(29,*) decmp3(depth(t2)), depth(t2)
#c

#c
#c **********************************************************************
#c * this section calculates the relative elevation of the marsh at the *
#c * end of the year                                                    *
#c **********************************************************************
#c

        relelv[time]=intelv+depth[time]-(slr*time)-(subsid*time)
    


    sksxx=pd.DataFrame({"totorg":totorg[1:],
                        "min":mini[1:],
                        "pore": pore[1:],
                        "totvol":totvol[1:],
                        "orgvol":orgvol[1:],
                        "minvol": minvol[1:],
                        "porg": porg[1:],
                        "bulkd": bulkd[1:],
                        "depth":depth[1:],
                        "massabv": massabv[1:],
                        "densabv":densabv[1:],
                        "relelv": np.flip(relelv,0)[1:],
                        "time":np.array(range(1,endtim+1))})
    sksxx = sksxx.reindex(
        columns=['totorg', "min", "pore", "totvol", "orgvol", "minvol", "porg", "bulkd", "depth", "massabv", "densabv",
                 "relelv", "time"])


    sksxx.to_csv("sksxx.csv",index=False)
    
    return sksxx


#############################################################################
##########################CODE###############################################
#############################################################################

#Last version
###Raster template. All the numpy matrices must be created using rasters with the same resolution, domain and reference
## as the template
template_name="Template.tif"
template_path=os.path.join(os.getcwd(),template_name)

###We use this part of the code to call SEDCALC

t0=datetime.datetime.now()
#####Parameters###############

###Base tidal datum numpy array name. We have to move this array manually to the 
###np_path folder we'll set up later in the code
BTD_nam="MTL.npy"

#MTL_cutoff is the elevation in cm above the tidal datum (MTL) where accretion
#will stop.

MTL_cutoff=0

#First year of simulation
year0=2017

yearf=2101
#Number of years of simulation
sim_years=yearf-year0

#SEDCALC scenarios: We use this list to set the names of the different SEDCALC scenarios we want to run.
#These are different than the SLR scenarios

SEDCALC_scenarios=["SEDCALC"]

#SLR scenarios
SLR_scenarios=["1_1_ft","2_6_ft"]

#SLR masks

masks={"1_1_ft":np.load('SLR_1_1.npy'),
       "2_6_ft":np.load('SLR_2_6.npy')}

Grid_Codes=[1,9]

#CSV file with SLR time series. This file must be in the working directory
SLR_nam="SLR_ts.csv"

#Dictionary of k1 parameters for the SEDCALC scenarios defined in SEDCALC_scenarios
k1={"SEDCALC":2.5}

#Name of the numpy array where we have our DEM
DEM_name="TopoBathy_10m.npy"

#Directory where we'll store numpy arrays

np_path=os.path.join(os.getcwd(),"np_arrays")


if not os.path.exists(np_path):
    os.makedirs(np_path)

#Directory where we'll store rasters
rasters_path=os.path.join(os.getcwd(),"rasters")

if not os.path.exists(rasters_path):
    os.makedirs(rasters_path)

#Base tidal datum
MTL_0=np.load(os.path.join(os.getcwd(),BTD_nam))

#Let's convert to cm
#MTL_0=MTL_0*30.48


#Let's create sea level rise time series






#Let's import SLR time series
SLR_ts=pd.read_csv(SLR_nam,usecols=["Year",
                                    "1_1_ft_Gridcode_9",
                                    "1_1_ft_Gridcode_1",
                                    "2_6_ft_Gridcode_9",
                                    "2_6_ft_Gridcode_1"],
                   index_col=["Year"])

##Let's use the SLR time series to create tidal datums for every simulation year

for slr in SLR_scenarios:
    #Let's get numpy indices of each gridcode
    ind_dum= {}
    for grid_code in Grid_Codes:
        ind_dum[grid_code] = masks[slr] == grid_code

    #Let's loop through simulation years
    for t in range(sim_years):
        year=year0+t
        MTL_dum = np.empty(MTL_0.shape)
        MTL_dum[:]= np.NaN
        #Let's check if tidal datum array already exists for that year and scenario
        if not os.path.exists(os.path.join(np_path,"MTL_"+str(slr)+"_"+str(year)+".npy")):
            #Let's loop through grid codes
            for grid_code in Grid_Codes:
                print(SLR_ts.loc[year,slr+'_Gridcode_'+str(grid_code)])
                MTL_dum[ind_dum[grid_code]]=MTL_0[ind_dum[grid_code]]+SLR_ts.loc[year,slr+'_Gridcode_'+str(grid_code)]
            #We create a dummy array for the year
  #          MTL_dum=MTL_0+SLR_ts.loc[t,slr]
            #Let's set save the array to a numpy file
            np.save(os.path.join(np_path,"MTL_"+str(slr)+"_"+str(year)+".npy"),MTL_dum)

#Let's delete the dummy array to release some RAM
if 'MTL_dum' in locals():
    del MTL_dum




#Let's import LiDAR 2017 array

LiDAR_2017=np.load(os.path.join(os.getcwd(),DEM_name))

LiDAR_2017[np.isnan(MTL_0)]=np.nan
LiDAR_2017[LiDAR_2017>MTL_0+MTL_cutoff]=np.nan

#Number of rows
nrows=LiDAR_2017.shape[0]

#Number of columns
ncols=LiDAR_2017.shape[1]

#Dictionary where we will store elevations for the different SLR scenarios
Elevations_all={}


#Empty data frame where we will save total carbon
#totC=pd.DataFrame({"Year":range(year0+1,year0+sim_years)})

#Let's loop through sensitivity scenarios
for sensitivity in SEDCALC_scenarios:
    
    #Let's run SEDCALC
    sedcalc_dum=sedcalc(endtim=sim_years,k1=k1[sensitivity])
    #We export SEDCALC output table to a csv
    sedcalc_dum.to_csv("SEDCALC_out_"+sensitivity+".csv",index=False)
    relelv_dum=sedcalc_dum["relelv"].values
    depth_dum=sedcalc_dum["depth"].values
    porg_dum=sedcalc_dum["porg"].values
    bulkd_dum=sedcalc_dum["bulkd"].values
    #Let's create data frame with accretion rates
    acc_rate_dum=pd.DataFrame(columns=["Year","Yearly Accretion (cm)","Yearly Accretion (ft)"])
#    totorg_dum=pd.DataFrame(columns=["Year","Accreted Organic matter (g/cm3)"])
    totC_dum=pd.DataFrame(columns=["Year","Yearly Accretion (cm)","gC/cm3","gC/cm2"])

    #Let's get accretion rates and organic matter
    year=year0+1
    for i in range(sim_years-1):
        #Index for the backwards sum
        i_bw=sim_years-i-1
                                           
#        else:
        acc_rate_dum_dum=pd.DataFrame({"Year":[year],
                                        "Yearly Accretion (cm)":[(depth_dum[i_bw]-depth_dum[i_bw-1])],
                                        "Yearly Accretion (ft)":[(depth_dum[i_bw]-depth_dum[i_bw-1])*0.0328084]})
        totC_dum_dum=pd.DataFrame({"Year":[year],
                                   "Yearly Accretion (cm)":[(depth_dum[i_bw]-depth_dum[i_bw-1])],
                                   "gC/cm3":[porg_dum[i_bw]*bulkd_dum[i_bw]/2],
                                    "gC/cm2":[(depth_dum[i_bw]-depth_dum[i_bw-1])*porg_dum[i_bw]*bulkd_dum[i_bw]/2]})



        acc_rate_dum=acc_rate_dum.append(acc_rate_dum_dum,ignore_index=True)
        totC_dum=totC_dum.append(totC_dum_dum,ignore_index=True)
        year=year+1
    #totC[sensitivity]=totC_dum["gC/cm2"]
    #Let's export the total carbon table as a csv file
    totC_dum.to_csv(sensitivity+"_totC_output.csv",index=False)

    #Let's export the SEDCALC accretion rates
    acc_rate_dum.to_csv(sensitivity+"_SEDCALC_Accretion.csv",index=False)

        
        #Let's loop through scenarios
    for scenario in SLR_scenarios:
        #Years to MTL +10cm (0.328084 ft)
        years_0= np.zeros((nrows,ncols))
        years_0[np.isnan(MTL_0)]=np.nan
        #Let's loop through years now
        #First year is 2017
        elevations_0=LiDAR_2017


        carbon_cum = np.zeros((nrows, ncols))
        carbon_cum[np.isnan(MTL_0)]=np.nan

        #The following matrix keeps track of whether the manage wetland surface
        #elevation has caught up with the MTL+MTL_cutoff, therefore converting the wetland
        #into tidal. 1 is for managed wetland and 0 is for tidal wetland. It starts as
        #a matrix of ones. If accretion caughts up with the tidal datum, it will turn into a cero for the
        #rest of the simulation
        #MWL=np.ones((nrows, ncols))
        MWL = np.zeros((nrows, ncols))
        MWL[:]=np.nan
        MWL[~np.isnan(MTL_0)]=1

        #let's create directory for numpy arrays
        dir_dum=os.path.join(np_path,sensitivity+"_"+scenario)

        if not os.path.exists(dir_dum):
               os.makedirs(dir_dum)

        for year in range(year0+1, year0+sim_years):
            #MTL at the beginning of the year
            MTL_dum=np.load(os.path.join(np_path,"MTL_"+str(scenario)+"_"+str(year)+".npy"))
            #Elevation accreted during the year
            yearly_acc_rate=float(acc_rate_dum.loc[acc_rate_dum["Year"]==year,"Yearly Accretion (cm)"])
            print(sensitivity+" SLR: "+scenario+ " year: "+str(year)+" accretion: "+str(round(yearly_acc_rate,1))+" cm")
            #Let's calculate elevation at the beginning of year
            elevations_1=np.minimum(elevations_0+yearly_acc_rate,MTL_dum+MTL_cutoff)

            #Let's calculate cummulative accretion
            cum_ac=np.maximum(elevations_1-LiDAR_2017,0)

            MWL=MWL*(elevations_1<MTL_dum+MTL_cutoff).astype(int)
            #Let's add 1 to the cells which haven't reached MTL + MTL_cutoff elevation yet
            years_0=years_0+MWL*(elevations_1<MTL_dum+MTL_cutoff).astype(int)
            np.save(os.path.join(dir_dum, sensitivity+"_"+str(scenario)+"_"+str(year)+
                                 "_years_to_MTL_"+str(int(MTL_cutoff))+".npy"),years_0)
            #Let's add updated array to dictionary
#            Elevations_dum[year]=elevation_new_dum
            #Let's export matrix of elevations
            np.save(os.path.join(dir_dum, sensitivity+"_"+str(scenario)+"_"+str(year)+".npy"),elevations_1)

            # Let's export cumulative accretion
            np.save(os.path.join(dir_dum, sensitivity + "_" + str(scenario) + "_" + str(year) + "cum_ac_cm.npy"), cum_ac)

            #Let's do carbon now
            carbon_dum=np.zeros((nrows, ncols))
            carbon_dum[:]=np.nan

            #Cells where full accretion occured
            ind_dum=np.where(elevations_1<MTL_dum+MTL_cutoff)
            carbon_dum[ind_dum]=float(totC_dum.loc[totC_dum.Year==year,"gC/cm2"])

            #Cells with partial accretion
            ind_dum=np.where(elevations_0+yearly_acc_rate>MTL_dum+MTL_cutoff)
            #We add carbon proportionally to the accretion rate
            carbon_dum[ind_dum]=np.maximum(0,MTL_dum[ind_dum]+MTL_cutoff-elevations_0[ind_dum])/yearly_acc_rate*float(totC_dum.loc[totC_dum.Year==year,"gC/cm2"])


            #Let's export yearly rates
            np.save(os.path.join(dir_dum, sensitivity + "_" + str(scenario) + "_" + str(year) + "carbon_rate_gC_cm2_yr.npy"), carbon_dum)

            #Let's update cumulative carbon
            carbon_cum=carbon_cum+carbon_dum

            # Let's export cumulative carbon
            np.save(os.path.join(dir_dum, sensitivity + "_" + str(scenario) + "_" + str(year) + "carbon_cum_gC_cm2.npy"),
                carbon_cum)

            #Update elevations t-1 at the end of year t
            elevations_0=elevations_1


t1=datetime.datetime.now()
print("Elapsed time:",t1-t0)

########################################################################################################################
#######POST PROCESSING###################################################################################################
#######################################################################################################################

#In this part we convert the numpy matrices we generated into rasters

###We extract the geospatial properties of the raster template

template=arcpy.Raster(template_path)

arcpy.env.outputCoordinateSystem=template_path

inRas = arcpy.Raster(template_path)

lowerLeft = arcpy.Point(inRas.extent.XMin,inRas.extent.YMin)

cellSize = inRas.meanCellWidth

#Let's get names of subdirectories within numpy arrays folder
Scenarios=[x[0][x[0].rfind("\\")+1:] for x in os.walk(np_path)][1:]

for scenario in Scenarios:
   # let's create directory for rasters
    dir_dum = os.path.join(rasters_path, scenario)

    if not os.path.exists(dir_dum):
        os.makedirs(dir_dum)
    for path in glob.glob(os.path.join(np_path,scenario)+"/*"):
        name=path[path.rfind("\\")+1:-4]
        arr = np.load(path)
        arr=arr.astype('float32')
        NewRaster = arcpy.NumPyArrayToRaster(arr, lowerLeft, cellSize, value_to_nodata=np.NAN)
        outpath = os.path.join(dir_dum,name+".tif")
        #This line creates a tif file
        arcpy.CopyRaster_management(NewRaster, outpath)
        #This line saves the raster in the geodatabase
        arcpy.CopyRaster_management(NewRaster, name)

    ## clip by extent
    # runspath = r"P:\Projects\5635_SFEI Carbon and GHG\GIS\delta_wide\SUBCALC\runs"
    # extent = r"P:\Projects\5635_SFEI Carbon and GHG\GIS\vector\peat_extent_1.shp"
    # clippath = os.path.join(runspath,os.path.splitext(os.path.basename(outpath))[0] + "_extent.tif")
    # arcpy.Clip_management(outpath,"#",clippath,extent,"0","ClippingGeometry")
