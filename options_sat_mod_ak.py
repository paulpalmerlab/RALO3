#!/usr/bin/env python

"""options_sat_mod_ak.py : User options for sat_mod_ak"""

from datetime import datetime as dt
from datetime import timedelta as td
from datetime import date as date
import os

def replaceYYYYMMDD(string,d):
    YYYYMMDD = d.strftime("%Y%m%d")
    YYYYMM    = d.strftime("%Y%m")
    YYYY      = d.strftime("%Y")
    
    string = string.replace("YYYYMMDD",YYYYMMDD)
    string = string.replace("YYYYMM"  ,YYYYMM  )
    string = string.replace("YYYY"    ,YYYY    )
    
    return string
    
class opts:
    """A class to hold all of the options in one object"""
    def __init__(self,
           start_date,end_date,
           cycle_type,
           domain,
           sat_instrument,
           ret_version,
           sat_time_res,
           do_gsc,do_MACC,do_MACC_wAK,do_prior,
           geos_species,
           mask_as_all,
           do_prodloss,prodloss_path,
           do_geos_nox,do_geos_o3_wAK,
           do_3D,
           mod_data_path,sat_data_path,
           sav_pre
           ):
        self.start_date      = start_date
        self.end_date        = end_date
        self.cycle_type      = cycle_type
        self.domain          = domain
        self.sat_instrument  = sat_instrument
        self.ret_version     = ret_version
        self.sat_time_res    = sat_time_res
        self.mask_as_all     = mask_as_all
        self.do_gsc          = do_gsc
        self.do_MACC         = do_MACC
        self.do_MACC_wAK     = do_MACC_wAK
        self.do_prior        = do_prior
        self.geos_species    = geos_species
        self.do_prodloss     = do_prodloss
        self.prodloss_path   = prodloss_path
        self.do_geos_nox     = do_geos_nox
        self.do_geos_o3_wAK  = do_geos_o3_wAK
        self.do_3D           = do_3D
        self.mod_data_path   = mod_data_path
        self.sat_data_path   = sat_data_path
        self.sav_pre         = sav_pre

def get_options():

    print("Reading options from options_sat_mod_ak.py")
    #==Date and time==

    start_date = dt(2016, 1,   1) #Start date          :: (YYYY,MM,DD)
    end_date   = dt(2016, 1,  31) #End date (inclusive):: (YYYY,MM,DD)

    #==Output options==
    cycle_type = "month"       #Time resolution of output:
                               #"month" -> take averages for individual months
                               #"day"   -> produce output on daily resolution

    #domain
    domain = [-11.,55.,60.,150.]  #Spatial span of output in format [south, north, west, east] degrees 

    #satellite options
    sat_instrument = "omi"     #Choice of satellite instrument: options are "omi" and "g2a" 
    ret_version = "fv0300_xc"  #Retrieval versions:
                               #OMI -> "fv0214", "fv0214_xc" or "fv0300_xc"
                               #G2A -> "fv0300" only 

    sat_time_res = "m"         #Temporal resolution of satellite file input
                               #Options are "m"=monthly, and "d"=daily
                               #Monthly available for all instruments and retrieval versions
                               #Daily only available for OMI fv0214
                               
    mask_as_all = True         #Mask all data as the satellite is masked
                               #For computing averages
                               
    do_gsc = True              #Process the retrieved o3 column
    do_MACC= True              #Process output from the MACC model
    do_MACC_wAK = True         #Process output from the MACC model after AK application
    do_prior = True            #Apply the prior

    #Model options
    
    #List species to read in, using their GEOS-Chem names.
    #Make list empty to do no model 
    geos_species = [
                    "O3"
                   # "NO",
                   # "NO2",
                   # "CO",
                   # "CH2O"
                   ]                 
    
    #Production-loss
    do_prodloss = False
    prodloss_path = '/path_to_gc_prodloss_dir/rate.YYYYMMDD'
    
    #Special species
    do_geos_nox = False           #Output "pure" GEOS-Chem tropospheric NOx (sum of NO and NO2)
    do_geos_o3_wAK = True         #Apply satellite AKs to GEOS-Chem O3                
    
    #Choose whether to return 3D fields of model output and satellite retrieval
    do_3D = False
                            
    #==Data locations==
    #YYYYMMMDD, YYYYMM and YYYYY will be replaced in model processing 
    #Path of model data
    if sat_instrument == "omi":
        mod_data_path = "/home/users/mmarvin/geoschem/nested_as_2016/ND51/ts_satellite_pm.YYYYMMDD.bpch"
    else:
        mod_data_path = "/home/users/mmarvin/geoschem/nested_as_2016/ND51/ts_satellite_am.YYYYMMDD.bpch"

    #Path of satellite data
    if sat_instrument == "omi":
        if ret_version == "fv0214":
            if sat_time_res == "m":
                sat_data_path = '/gws/nopw/j04/nceo_generic/nceo_ral/ozone_omi_fv0214/monthly_lag-90_90_2.5_log-180_180_2.5/YYYY/o3p_bin_vs_model_YYYYMMXX_MACC_mcef0.2_sza80_mcost120_mcb1_lzr_omfra_ak.str'
            elif sat_time_res == "d":
                sat_data_path = '/gws/nopw/j04/nceo_generic/nceo_ral/ozone_omi_fv0214/daily_lag-90_90_2.5_log-180_180_2.5/YYYYMM/o3p_bin_vs_model_YYYYMMDD_MACC_mcef0.2_sza80_mcost120_mcb1_lzr_omfra_ak.str'
        elif ret_version == "fv0214_xc":
            sat_data_path = '/gws/nopw/j04/nceo_generic/nceo_ral/ozone_omi_lv2fv0214_lv3fv0100.xc01/monthly_lag-90_90_1.5_log-180_180_1.5/YYYY/o3p_bvm_YYYYMMXX_omi_MACC_l2fv0214_l3fv0100_mcef0.2_sza80_mcost120_lzr_omfra2_xtxscp_ak_xtcor1.str'
        else:
            sat_data_path = '/gws/nopw/j04/nceo_generic/nceo_ral/ozone_omi_lv2fv0300_lv3fv0100.xc01/monthly_lag-90_90_1.5_log-180_180_1.5/YYYY/o3p_bvm_YYYYMMXX_omi_MACC_l2fv0300_l3fv0100_mcef0.2_sza80_mcost120_lzr_omfra2_xtxscp_ak_xtcor1.str'
    else:
        sat_data_path = '/gws/nopw/j04/nceo_generic/nceo_ral/ozone_gome2a_fv0300/monthly_lag-90_90_2.5_log-180_180_2.5/gome2_metopa/latest_processing/gome2_metopa-rat16-cclev2-saap-lin-lip-cat4-af0-kfk-cmacc-nkd-pcbira-salbv2-salbrat1-cxt2-ecm-eci-sfg2s-rsfw8-ng37-pol7-sc10-noxfr-ch2o15-tglin-i0off-retpoly-s5xs-sckxs-ring4-smcl2-mplap-g2lk-nomr-cci-newfm2b-b2leak-dxs2-0p006-cgo2-superg/YYYY/o3p_bvm_YYYYMMXX_g2a_MACC_mcef0.2_sza80_mcost120_ak.str'

    #==save location, and file prefix==
    #Set save file prefix. The rest of the filename will be the date span (from start_date to end_date).
    sav_pre = "./gc_"+sat_instrument+"_"+ret_version
    
    return(
           start_date,end_date,
           cycle_type,
           domain,
           sat_instrument,
           ret_version,
           sat_time_res,
           do_gsc,do_MACC,do_MACC_wAK,do_prior,
           geos_species,
           mask_as_all,
           do_prodloss,prodloss_path,
           do_geos_nox,do_geos_o3_wAK,
           do_3D,
           mod_data_path,sat_data_path,
           sav_pre
           )

def check_options_valid(options):
   
   #expand out the options 
    (
    start_date,end_date,
    cycle_type,
    domain,
    sat_instrument,
    ret_version,
    sat_time_res,
    do_gsc,do_MACC,do_MACC_wAK,do_prior,
    geos_species,
    mask_as_all,
    do_prodloss,prodloss_path,
    do_geos_nox,do_geos_o3_wAK,
    do_3D,
    mod_data_path,sat_data_path,
    sav_pre
    ) = (
    options.start_date,options.end_date,
    options.cycle_type,
    options.domain,
    options.sat_instrument,
    options.ret_version,
    options.sat_time_res,
    options.do_gsc,options.do_MACC,options.do_MACC_wAK,options.do_prior,
    options.geos_species,
    options.mask_as_all,
    options.do_prodloss,options.prodloss_path,
    options.do_geos_nox,options.do_geos_o3_wAK,
    options.do_3D,
    options.mod_data_path,options.sat_data_path,
    options.sav_pre
    )
    
    major_error = 0 #this flag will increase by 1 for all fundamentally invalid options. The script will not continue if this is not 0 after the end of this module
    minor_error = 0 #this flag will increase by 1 for all options which appear invalid, but an apparent workaround is possible.
    
    #checks that the options make sense. Reports all errors.
    
    #Dates.
    if end_date < start_date: #end_date must be after or same as start_date
        print("OPTIONS ERROR - FATAL: start_date (%s) is after end_date (%s)" %(start_date.strftime("%Y%m%d"),end_date.strftime("%Y%m%d")))
        major_error += 1
    if cycle_type == "season":
        if start_date.strftime("%m%d") not in ["0101","0401","0701","1001"]: #if season, start_date MMDD should be either 0101, 0401, 0701 or 1001
            print("OPTIONS ERROR - MINOR: Seasonal output is requested, but start_date (%s) is not at the start of a season (01-Jan, 01-Apr, 01-Jul, 01-Oct). Output for the first season will be the average of a shorter time frame" %(start_date.strftime("%Y%m%d")))
            minor_error += 1
        if end_date.strftime("%m%d") not in ["0331","0630","0930","1231"]: #if season, end_date MMDD must be either 0331, 0630,0930, or 1231
            print("OPTIONS ERROR - MINOR: Seasonal output is requested, but end_date (%s) is not at the end of a season (31-Mar, 30-Jun, 30-Sep, 31-Dec). Output for the last season will be the average of a shorter time frame" %(end_date.strftime("%Y%m%d")))
            minor_error += 1
    elif cycle_type == "month":
        if start_date.day != 1: #if monthly, we should start on the 1st day of a month.
            print("OPTIONS ERROR - MINOR: Monthly output is requested, but start_date (%s) is not at the start of a month. Output for the first month will be the average of a shorter time frame" %(start_date.strftime("%Y%m%d")))
            minor_error += 1
        if (end_date+td(days=1)).day != 1: #if monthly, we should end on the last day of a month.
            print("OPTIONS ERROR - MINOR: Montly output is requested, but end_date (%s) is not at the end of a month. Output for the last month will be the average of a shorter time frame" %(end_date.strftime("%Y%m%d")))
            minor_error += 1
            
    #domain.
    [south, north, west, east] = domain
    if south > north or west > east:
        print("OPTIONS ERROR - FATAL: domain [ south = %g, north = %g, west = %g, east = %g ] is invalid"%(south,north,west,east))
        major_error += 1
    
    #prodloss needs model
    if do_prodloss and geos_species == []:
        print("OPTIONS ERROR - MINOR: Can't do prodloss without having at least one GEOS-Chem species")
        minor_error += 1
    
    #mask_as_all resquires GSC
    if mask_as_all and not do_gsc:
        print("OPTIONS ERROR - FATAL: mask_as_all=True requires do_gsc=True")
        major_error += 1
     
    #Satellite data
    #Are we using satellite data?
    use_sat_data = do_gsc or do_MACC or do_MACC_wAK or do_prior or do_geos_o3_wAK

    if use_sat_data: #following checks only relevant *if* satellite data is used

        if sat_instrument != "omi" and sat_instrument != "g2a":
            print("OPTIONS ERROR - FATAL: Selected satellite instrument is not currently supported.")
            major_error += 1
        if ret_version != "fv0214" and ret_version != "fv0214_xc" and ret_version != "fv0300" and ret_version != "fv0300_xc":
            print("OPTIONS ERROR - FATAL: Selected retrieval version is not currently supported.")
            major_error += 1
        if (sat_instrument == "omi" and ret_version == "fv0300") or (sat_instrument == "g2a" and ret_version != "fv0300"):
            print("OPTIONS ERROR - FATAL: Selected retrieval version is not compatible with selected satellite instrument.")
            major_error += 1

        #time resolution of satellite input data vs. start date and end date
        if sat_time_res == "d" and (sat_instrument != "omi" or ret_version != "fv0214"):
            print("OPTIONS ERROR - FATAL: Daily input is not available for selected satellite instrument and retrieval version.")
            major_error += 1
        if sat_time_res == "m" and start_date.day != 1:
            print("OPTIONS ERROR - FATAL: Monthly input is used, but start_date (%s) is not at the start of a month." %(start_date.strftime("%Y%m%d")))
            major_error += 1
        if sat_time_res == "m" and (end_date+td(days=1)).day != 1:
            print("OPTIONS ERROR - MINOR: Monthly input is used, but end_date (%s) is not at the end of a month. Be aware that satellite data used for the final month's output will be the average of all of %s" %(end_date.strftime("%Y%m%d"),end_date.strftime("%Y%m")))
            minor_error += 1

        #check data exists
        date_loop = start_date
        while date_loop <= end_date:
            #only check if daily data used or monthly data used AND it's the first day of the month
            if sat_time_res == "d" or (sat_time_res == "m" and date_loop.day == 1):
                file_to_check = replaceYYYYMMDD(sat_data_path,date_loop)
                if not os.path.isfile(file_to_check):
                    print("OPTIONS ERROR - MINOR: File %s is requested but does not exist, time period will be skipped"%file_to_check)
                    minor_error += 1
            date_loop += td(days=1)
    
    #GEOS-Chem data
    #Are we using GEOS-Chem data?
    
    if geos_species != []: #following check only relevant *if* GEOS-Chem data is used
        
        #check data exists
        date_loop = start_date
        while date_loop <= end_date:
            #check each day for a valid GEOS-Chem output file
            file_to_check = replaceYYYYMMDD(mod_data_path,date_loop)
            if not os.path.isfile(file_to_check):
                print("OPTIONS ERROR - FATAL: File %s is requested but does not exist"%file_to_check)
                major_error += 1
            date_loop += td(days=1)
    
    if do_prodloss:        
        #check data exists
        date_loop = start_date
        while date_loop <= end_date:
            #check each day for a valid GEOS-Chem output file
            file_to_check = replaceYYYYMMDD(prodloss_path,date_loop)
            if not os.path.isfile(file_to_check):
                print("OPTIONS ERROR - FATAL: File %s is requested but does not exist"%file_to_check)
                major_error += 1
            date_loop += td(days=1)
        
    #make sure special options are compatible with normal options.
    if do_geos_nox and ("NO" not in geos_species or "NO2" not in geos_species): #NO and NO2 must be extracted to calculate NOx
        print("OPTIONS ERROR - FATAL: do_geos_nox = True is not compatible with geos_species list that does not contain NO and NO2")
        major_error += 1
    if do_geos_nox and ("O3" not in geos_species): #O3 must be extracted to calculate O3 with AK
        print("OPTIONS ERROR - FATAL: do_o3_wAK = True is not compatible with geos_species list that does not contain O3")
        major_error += 1    
    
    #check species are in files?
            
    #Report results
    print("Options have been checked. Found %i fatal errors and %i minor errors"%(major_error,minor_error))
    
    if major_error != 0:
        raise IOError("At least one fatal error in input options")        
             
