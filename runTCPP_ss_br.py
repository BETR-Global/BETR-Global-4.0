## This is a run for BETR-Reserach with annually changing env.-parameters ####
import os
import time
import BETRS
reload(BETRS)
from BETRS import *
import pdb

t_s = time.time() # Start time

"""
SET OPTIONS FOR THE FAST SOLVER AND THE FLUX_INTEGRATION AND
PRIMARY/SECONDARY_EMISSION MODE
"""

use_odespy = False       #SWITCH ON/OFF FAST SOLVER
track_fluxes = False    #SWITCH ON/OFF FLUX INTEGRATION
track_se = False         #SWICH ON/OFF TRACKING OF SECONDARY EMISSIONS
use_correction = False   #SWICH ON/OFF flow with correction

## options
"""
Change to current Directory to ensure that relative paths are set correctly

Tends to cause problem under Windows otherwise

"""

abspath=os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

runID = ['TCPP_ss_br'] # output names 
years = [range(1)]*len(runID)    # range of modeling run (years)

#emisdir = ['Emission_BDE209_10_comparts']  # emission inventory ('Emissions/annual/')
emisfile = ['Emission_TCPP_10_comparts.txt']*len(runID) # emission inventory ('Emissions/')
seasparfile = ['seasonal_parameters_20C3M.a2.comp.ivm.moh_add2mix_esrl_br.txt']*len(runID)#  seasonally varying parameters ('Environment/)
constparfile = ['const_parameters_20C3M.comp.ivm_corrdiffstrato_7scenarios_br.txt']*len(runID)  # seasonally constant parameters ('Environment/')
flowdirectory = ['OcMix.At100_20C3M.a2.comp.ivm_br']*len(runID)  # flows in the atmosphere, ocean and fresh water ('Flows/)

chemdata = ['chemicals.txt']*len(runID)  # chemical properties ('Chemicals/')
chemnr = [3]*len(runID)  # selection of chemical from chemical properties files
compfile = ['compartments_4T.txt']*len(runID)   # compartments used in the model ('Environment/') 

procfile = ['processes_separate_particles_from_gas.txt']*len(runID)   # processes used in the model  ('Processes/')
contfile = ['control_default.txt']*len(runID)      # some options ('Control/')
solvfile = ['solvparams_default.txt']*len(runID)    # options for ODE solver  ('Solver/')
mkendfile = False


for v in [emisfile, seasparfile, chemdata, chemnr, constparfile, compfile, flowdirectory, procfile, contfile, solvfile]:
    if len(v) != len(runID):
        sys.exit('Warning: one of your input lists is not of same length as number of runs specified')        
    

## now run the model

for i in range(0, len(runID)):
    print('\n\nStarting run ' + runID[i])
    ## model first year and write temporary result to text file
    print('\n\nBETR run ' + runID[i] + ' for year ' + str(years[i][0]))
    m=Model(chemical = chemnr[i],
            run = runID[i],
            chemdb = chemdata[i],
            #seasonalparfile = os.path.join('annual', seasdir[i], str(years[i][0]) + '.txt'),
            seasonalparfile = seasparfile[i],    
            constantparfile = constparfile[i], 
            compartmentfile = compfile[i], 
            #flowdir = os.path.join('annual', flowdirectory[i], str(years[i][0])),
            flowdir = os.path.join(flowdirectory[i]),
            processfile = procfile[i], 
            controlfile=contfile[i],
            use_correction = use_correction, 
            track_flows = (track_fluxes or track_se)
            )
    m.update_emissions(emisfile[i])
    m.solve_ss()
    if mkendfile == True:
        m.output_end_txt(filename = 'endstate.txt')
    result = m.ss_res
    m.ss_res = result
    m.output_ss(filename = 'ss', units = ['mol_per_m3', 'Pa','ng_per_m3'], netcdf = True, cpk = True)
    m.output_ss_txt(filename = 'ss.txt')

    if track_fluxes:
        m.flux_res = flux_rkg # experimenting to get the whole fluxes written, RKG, 04.07.2014
        m.output_fluxes(netcdf=True)
    if track_se:
        m.output_se(cpk=False, netcdf=True,scenario="air")

    t_e = time.time() # End time

    print 'Simulation completed!'
    print 'TOTAL SIMULATION TIME: %f minutes' % ((t_e - t_s) / 60.0)

    
