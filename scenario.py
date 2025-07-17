# %% Import necessary modules
import time
import math
import pandas as pd
import numpy as np

import biosteam as bst 
import thermosteam as tmo
from thermo import SRK
#  import compounds and set thermo  
from _compounds import *
from _Hydrocracking_Unit import *
from _Hydrogen_production import *



from _Grinder import *
from _CHScreen import *
from _RYield import *
from _Cyclone import *
from _Sand_Furnace import *
from _UtilityAgents import *    # Created heat utilities that can heat to high temperatures and cool to sub zero temperatures 
from _process_yields import *
from _Compressor import *
from _feed_handling import *
from _teapyrolysis import *
from _tea_wax_mfsp import *
from _pass_unit import *
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 14})

bst.nbtutorial() # Light-mode html diagrams and filter warnings

# Set CEPCI to year of analysis to be 2020
bst.settings.CEPCI = 596.2 # CEPCI 2021 = 708, 2020 = 596.2, 2019	= 607.5

for c in compounds:
    c.default()

# compounds["CO2"].Cn.l.add_method(71, Tmin=-700, Tmax=120)

bst.HeatUtility().cooling_agents.append(NH3_utility)
bst.HeatUtility().heating_agents.append(Liq_utility)
bst.HeatUtility().heating_agents.append(Gas_utility)


# PREFERENCES
bst.preferences.update(flow='kg/hr', T='degK', P='Pa', N=100, composition=True)
 
#  Prices 
actual_prices = {
"HDPE": 25.05/1000, # $/kg from 22 per MT, because feed is defined per MT
"Ethylene": 0.86, # Gracida Al
"Propylene":0.8, # $/kg
"Butene": 1.27,# Butene 1.27 $/kg from Yadav et al. 2023
"Naphtha": 0.65,  # Gracida Alvarez et al
"Diesel": 0.64,  # Gracida Alvarez et al
"Wax": 0.3, #   1989 USD/MT   source: https://www.chemanalyst.com/Pricing-data/paraffin-wax-1205
"NG": 7.40 * 1000 * 1.525/28316.8, 
"Hydrocracking catalyst": 15.5 * 2.20462262,      #15.5 $/lb 2.20462262 is factor to convert to $/kg from Jones et al. 2009 PNNL report SRI international 2007
"Hydrotreating catalyst": 15.5 * 2.20462262,      #15.5 $/lb from Jones et al. 2009 PNNL report SRI international 2007
"Hydrogen plant catalyst": 3.6/885.7,      #3.6 /1000scf/28.317m3/885.71kg 2007 quote from Jones et al. 2009 PNNL report SRI international 2007
"Hydrogen": 2.83      #2.83 USD/kg from Gracida Alvarez et al. 2.1 from Borja Hernandez et al. 2019
}
bst.settings.electricity_price = 0.065 # Gracida Alvarez et al

#  Create dictionaries to store various technologies information
#  NAtural gas required for furnace activity = 2.2 GJ/tonne of HDPE from Gracida Alvarez
#  44.1 MJ/kg of natural gas; 1 GJ =  MJ; 1 tonne = 1000 kg;2.2 GJ/49.887 kg of natural gas
# Jessica report %wt closure 92.7% for CPY, 97% for TOD.
# Jessica assumes equivalent ratio of 7% Oxygen by volume, Means mass closure for TOD in products would be97.231*93/100 assuming oxygen is fully used up and ends up in pyrolysis products 
cpy = {"Technology":"CPY","Hydrocracking": "No", "Yield": cpy_comp_yield,"Reactor size" : 1,"wt_closure": 92.5, "NG_req":2.0,"residence_time" :"low"}
cpy_hc = {"Technology":"CPY","Hydrocracking": "Yes", "Yield": cpy_comp_yield,"Reactor size" : 1,"wt_closure": 92.5, "NG_req":2.0,"residence_time" : "low"}
tod = {"Technology":"TOD","Hydrocracking": "No", "Yield": tod_comp_yield, "Reactor size" : 0.01,"wt_closure":90.4, "NG_req":0.0,"residence_time" :"low"}
tod_hc = {"Technology":"TOD","Hydrocracking": "Yes", "Yield": tod_comp_yield,"Reactor size" : 0.4,"wt_closure":90.4, "NG_req":0.0,"residence_time" :"low"}
hrt = {"Technology":"CPY","Hydrocracking": "No", "Yield": hrt_comp_yield, "Reactor size" : 1, "wt_closure": 92.5, "NG_req":2.0,"residence_time" :"high"}
hrt_hc = {"Technology":"CPY","Hydrocracking": "Yes", "Yield": hrt_comp_yield,"Reactor size" : 1,"wt_closure": 92.5, "NG_req":2.0,"residence_time" :"high"}

scenarios = [cpy,cpy_hc,tod,tod_hc,hrt,hrt_hc]
scenarios_labels = ["CPY","CPY-HC","TOD","TOD-HC","HRT","HRT-HC"]
plant_capacity = 250 # tonnes per day
capacity = plant_capacity
scenario = scenarios[4]
irr = 0.1


def run_scenario(scenario = cpy,capacity=plant_capacity,prices=actual_prices):
    bst.main_flowsheet.set_flowsheet("Plastic Pyrolysis" + time.strftime("%Y%m%d-%H%M%S"))
    # Feed stream, HDPE plastic 
    feed = bst.Stream('HDPE_Feed',
                        HDPE=1,
                        units='kg/hr',
                        T = 298,
                            price = prices['HDPE']/1000  # 22 $/MT; divide by 1000 to get $/kg 
                            )

    feed_mass = capacity             # 250 tonnes per dat    test was 143.435 # kg/hr
    feed.set_total_flow(feed_mass, 'tonnes/day')
    feed.price = prices['HDPE']  # 22 $/MT; divide by 1000 to get $/kg

    # Natural gas and water for hydrogen production and furnace 
    sand_stream = bst.Stream('sand', Sand=100, T=25 + 273.15, P=101325, phase='s')
    natural_gas = bst.Stream('natural_gas', CH4=100, T=25 + 273.15, P=101325, phase='g')
    comb_nat_gas = bst.Stream('comb_nat_gas', CH4=100, T=25 + 273.15, P=101325, phase='g')
    natural_gas.price = prices["NG"]
    comb_nat_gas.price = prices["NG"]
    water = bst.Stream('water',H2O = 100, T=25 + 273.15, P=101325, phase='l')

    # Oxygen for autothermal pyrolysis at 7% equivalence ratio
    pyrolysis_oxygen = bst.Stream('pyrolysis_oxygen',O2=1,units='kg/hr',T=298,price=0.000)
    oxygen_mass = 0.07 * feed_mass * 100/93   # 7% equivalence ratio from Polin et al. 2019
    pyrolysis_oxygen.set_total_flow(oxygen_mass, 'kg/hr')

    # fluidizing gas for the reactor
    fluidizing_gas = bst.Stream('fluidizing_gas',N2=1,units='kg/hr',T=298,price=0.000)
    fluidizing_gas_mass = 15   # fluidizing gas is 20kg/hr for now
    fluidizing_gas.set_total_flow(fluidizing_gas_mass, 'kg/hr')



    recycle = bst.Stream('recycle')
    char_sand = bst.Stream('S108')
    CRHDPE = bst.Stream('S104')
    rec_NCG = bst.Stream('S235')
    hydrocracked = bst.Stream('Cracked_HC')

    ref_Methane = bst.Stream('Methane_ref',Methane=1,units='kg/hr',T= 298,price=0.000)
    ref_Methane.set_total_flow(911*2,units="kg/hr")
    ref_Methane.T = 273.15 - 90

    ref_Methane2 = bst.Stream('Methane_ref2',Methane=1,units='kg/hr',T= 298,price=0.000)
    ref_Methane2.set_total_flow(911*2,units="kg/hr")
    ref_Methane2.T = 273.15 - 90

    ref_ethane = bst.Stream('ethane_ref',Ethane=1,units='kg/hr',T= 298,price=0.000)
    ref_ethane.set_total_flow(911,units="kg/hr")
    ref_ethane.T = 273.15 - 50



    ref_Tetrafluoroethane = bst.Stream('Tetrafluoroethane_ref',Tetrafluoroethane=1,units='kg/hr',T= 298,price=0.000)
    ref_Tetrafluoroethane.set_total_flow(250,units="kg/hr")
    ref_Tetrafluoroethane.T = 273.15 - 50

    ref_Propane = bst.Stream('Propane_ref',Propane=1,units='kg/hr',price=0.000)
    ref_Propane.set_total_flow(250,units="kg/hr")
    ref_Propane.T = 273.15 - 50

    ref_Propene = bst.Stream('Propene_ref',Propene=1,units='kg/hr',T= 298,price=0.000)
    ref_Propene.set_total_flow(1402,units="kg/hr")
    ref_Propene.T = 273.15 - 50
    ref_Propene.P = 1 * 101325
    ref_Propene.phase = 'g'


    HC_hydrogen = bst.Stream('Hydrogen',H2=1,units='kg/hr',T= 298,price=prices["Hydrogen"])


# 2.2 GJ of natural gas per tonne of HDPE needed for combustion
# 49.88 kg of NG need for 1 GJ of energy assuming 44.1 MJ/kg of NG

    #--------------------------------------------------------------------------------------------------------------
    # Pretreatment & Pyrolysis
    #--------------------------------------------------------------------------------------------------------------
    with bst.System('sys_pretreatment') as sys_pretreatment:
        # Pretreatment
        handling = Feed_handling('Handling',ins=feed,outs=("S102"))
        M1 = bst.units.Mixer('Mixer',ins = [handling-0,recycle])     # Mix for feed and recycled NCG stream
        grinder = Grinder('Grinder',ins=[M1-0],outs="S103")       #crush the feed
        CHscreen = Screen("CHScreen",ins=grinder-0,outs=[CRHDPE,recycle]) # screen the feed
        if scenario['Technology'] == "CPY":

            comb_nat_gas.set_total_flow(scenario["NG_req"] * 26.68 * feed.get_total_flow("tonne/hr"),units="m3/hr")
            # comb_nat_gas.set_total_flow(scenario["NG_req"] *49.88* feed.get_total_flow("tonne/hr"),units="kg/hr")
            furnace = Combustor("furnace", ins=[comb_nat_gas,sand_stream,rec_NCG,char_sand],outs=('S105'))
            M2 = bst.units.Mixer('Mixer2',ins = [CRHDPE,pyrolysis_oxygen,fluidizing_gas,furnace-0]) #mix oxygen, fluidizing gas and feed

            reactor = RYield('CFB_Reactor',ins=M2-0,outs=("S106"),yields=scenario["Yield"],factor= scenario["Reactor size"],wt_closure=scenario["wt_closure"])
            # separate the gas and solids in products stream
            Cyclone1 = Cyclone('Cyclone1',ins= reactor-0,outs=['S107',char_sand],efficiency=0.99)
        elif scenario["Technology"] == "TOD":
            M2 = bst.units.Mixer('Mixer2',ins = [CRHDPE,rec_NCG]) #mix oxygen, fluidizing gas and feed

            reactor = RYield('CFB_Reactor',ins=M2-0,outs=("S106"),yields=scenario["Yield"],factor= scenario["Reactor size"],wt_closure=scenario["wt_closure"])
            # separate the gas and solids in products stream
            Cyclone1 = Cyclone('Cyclone1',ins= reactor-0,outs=['S107',"S108"],efficiency=0.99)
        cooler1 = bst.units.HXutility('cooler', ins=Cyclone1-0, outs='S109', T=273.15 +10, rigorous=False) # rigorous = False ignores VLE calculations


    #--------------------------------------------------------------------------------------------------------------
    # product fractionation 
    #--------------------------------------------------------------------------------------------------------------

    with bst.System('sys_Product_Fractionation') as sys_product_fractionation:
        F1 = bst.units.Flash('Condenser', ins=cooler1-0, outs=('S201','S239'), P=101325, T = (cooler1-0).T)     

        H7 = bst.HXutility('Heater7',ins = F1-1, outs=("S232"),T = 273.15 + 150, rigorous=False)
        F3 = bst.units.Flash('FlashSeparator', ins= H7-0, outs=("S233","S234"), P= 1.01*101.325 ,T=273.15) # T = (heater4-0).T)
        K1 = bst.units.IsentropicCompressor('Compressor1',ins=F1-0,outs=("S202"),P = 2 * 101325, eta=0.8)

#         # # Reduce temperature of gas stream to 30C and then use refrigeration cycle to reduce to -40 C
        H2 = bst.units.HXutility('Heater2',ins=K1-0,outs=("S203"),T=273.15 + 30, rigorous=False)
 
        H3 = bst.HXprocess('evaporator_ref', ins = (ref_Propane,H2-0),outs=("","S204"), U=1000, phase0='g')
        H3_K = bst.units.IsentropicCompressor('compressor_ref',ins = H3-0,P=2 * 101325)
        H3_Con = bst.units.HXutility('condenser_ref', H3_K-0,T=273.15 - 50, V=1)
        H3_Exp = bst.units.IsenthalpicValve('expansion_device', H3_Con-0,outs=ref_Propane,P=1 * 101325)

#         # Compress the gaseous stream to 7 bars
        K2 = bst.units.IsentropicCompressor('Compressor2',ins=H3-1,outs=("S205"),P = 7 * 101325, eta=0.8) # originally 7 

        H4 = bst.HXprocess('evaporator_ref2', ins = (ref_ethane,K2-0),outs = ("","S206"), U=1000, phase0='g',T_lim1=273.15-50)
        # H3 = bst.HXprocess('evaporator_ref', ins = (ref_Propene,olu), U=1000, phase0='g',T_lim0=273.15-40)
        H4_K = bst.units.IsentropicCompressor('compressor_ref2',ins = H4-0,P=2 * 101325)
        H4_Con = bst.units.HXutility('condenser_ref2', H4_K-0,T=273.15 - 50, V=1)
        H4_Exp = bst.units.IsenthalpicValve('expansion_device2', H4_Con-0,outs=ref_ethane,P=1 * 101325) 


        H5 = bst.HXprocess('evaporator_ref3', ins = (ref_Methane,H4-1),outs=("","S207"), U=1000, phase0='g',T_lim1=273.15-80)
        # # H3 = bst.HXprocess('evaporator_ref', ins = (ref_Propene,olu), U=1000, phase0='g',T_lim0=273.15-40)
        H5_K = bst.units.IsentropicCompressor('compressor_ref3',ins = H5-0,P=2 * 101325)
        H5_Con = bst.units.HXutility('condenser_ref3', H5_K-0,T=273.15 - 50, V=1)
        H5_Exp = bst.units.IsenthalpicValve('expansion_device3', H5_Con-0,outs=ref_Methane,P=1 * 101325) 

        H5_2 = bst.HXprocess('evaporator_ref4', ins = (ref_Methane2,H5-1),outs=("","S208"), U=1000, phase0='g',T_lim1=273.15-90)
        H5_K2 = bst.units.IsentropicCompressor('compressor_ref3',ins = H5_2-0,P=2 * 101325)
        H5_Con2 = bst.units.HXutility('condenser_ref3', H5_K2-0,T=273.15 - 50, V=1)
        H5_Exp2 = bst.units.IsenthalpicValve('expansion_device3', H5_Con2-0,outs=ref_Methane2,P=1 * 101325) 


        F2 = bst.units.Flash('Condenser2', ins=H5_2-1, outs=("S210","S209"), P= (H5_2-1).P,T=273.15 - 110) # T = (heater4-0).T)

        P1 = bst.units.Pump('Pump',ins=F2-1,outs=("S211"),P = 25 * 101325) # 25 bars
        H6 = bst.HXutility('Heater6',ins = P1-0, outs=("S212"),T = 273.15 +2, rigorous=False)


        D1 = bst.units.BinaryDistillation('De_ethanizer', ins=H6-0,
                                outs=('S213',"S214"),   # ethylene
                                LHK=('C2H4', 'C3H8'),
                                y_top=0.99, x_bot=0.01, k=2,
                                is_divided=True)     #  (97.2% purity)
        D1.check_LHK = False

        D1_spl = bst.units.Splitter("EthyleneFractionator",ins=(D1-0),outs=("S215","S216"), split={'C2H4':0.99,'CO2':0.10,'C3H8':0.05,'O2':1,'CO':1,'H2':1})
        D1_spllMx = bst.Mixer('D1_spplMX', ins = D1_spl-0,outs= ("Ethylene"))
        ethyleneOut = (D1_spl-0)


        ethyleneOut = (D1-0)
        H8 = bst.HXutility('Heater8',ins = D1-1, outs=("S217"),T = 273.15 +100, rigorous=False)
        D2 = bst.units.BinaryDistillation('Depropanizer', ins=H8-0,
                                outs=('S218','S219'),   # propylene
                                LHK=('C3H8', 'C4H8'),
                                y_top=0.99, x_bot=0.01, k=2,
                                is_divided=True)
        H9 = bst.HXutility('Heater9',ins = D2-0, outs=("S220"),T = 273.15 +70, rigorous=False)
        KD2 = bst.units.IsentropicCompressor('CompressorD2',ins=H9-0,outs=("S221"),P = 22 * 101325, eta=0.8) # 25 bars
        D2_spl = bst.units.Splitter("PropyleneFractionator",ins=(KD2-0),outs=("S222","S223"), split={'C3H8':0.99,'C2H4':1,'C3H8':0.05,'O2':1,'CO':1,'H2':1})
        D2_spllMx = bst.Mixer('D2_spplMX', ins = D2_spl-0,outs= ("Propylene"))
        propyleneOut = (D2_spl-0)  

        M3 = bst.Mixer('Mixer3',ins = [F3-0,D2-1],outs=("S224"))    #,D2_spl-1,D1_spl-1])
        Mrec = bst.Mixer('Mixer_rec',ins = [D2_spl-1,F2-0,D1_spl-1],outs=rec_NCG)    #,D2_spl-1,D1_spl-1])
        D3 = bst.units.BinaryDistillation('Debutanizer', ins=M3-0,
                                outs=('S225','S226'),
                                LHK =('C4H8','C10H22'),      #=('C10H22', 'C14H30'),
                                y_top=0.99, x_bot=0.01, k=2,
                                is_divided=True)
        D3_mx = bst.Mixer('D3_MX',ins = D3-0,outs=("Butene"))
        buteneOut = (D3-0)


        if scenario['Hydrocracking'] == "No":          
            if scenario['residence_time'] == "high":
                D4 = bst.units.BinaryDistillation('NaphthaSplitter', ins=D3-1,
                                        outs=('S227','S228'),
                                        LHK =('C10H22','C14H30'),      #=('C10H22', 'C14H30'),
                                        y_top=0.99, x_bot=0.01, k=2,
                                        is_divided=True)
                D4_mx = bst.Mixer('D4_MX',ins = D4-0,outs=("Naphtha"))
                Mdiesel = bst.Mixer('MixerDiesel',ins = D4-1,outs=("S229"))     #,D2_spl-1,D1_spl-1]) 
                Md_mx = bst.Mixer('MixerDiesel_mx',ins = Mdiesel-0,outs=("Diesel"))   
                Mwax = bst.Mixer('Mwax',ins = F3-1, outs = ("S236"))
                Mwax_mx = bst.Mixer('Mwax_mx',ins = Mwax-0, outs = ("Wax"))

                naphthaOut1 = (D4-0)                  
                dieselOut = (Mdiesel-0)
                waxOut = Mwax-0
            else:   
                D4 = bst.units.BinaryDistillation('NaphthaSplitter', ins=D3-1,
                                        outs=('S227',"S228"),
                                        LHK =('C10H22','C14H30'),      #=('C10H22', 'C14H30'),
                                        y_top=0.99, x_bot=0.01, k=2,
                                        is_divided=True)
                D4_mx = bst.Mixer('D4_MX',ins = D4-0,outs=("Naphtha"))
                naphthaOut1 = (D4-0)   

                D5 = bst.units.BinaryDistillation('DieselSplitter', ins=D4-1,
                                    outs=('S229','S230'),
                                    LHK =('C14H30','C24H50'),      #=('C10H22', 'C14H30'),
                                    y_top=0.99, x_bot=0.01, k=2,
                                    is_divided=True)
                D5_mx = bst.Mixer('D5_MX',ins = D5-0,outs=("Diesel"))
                dieselOut1 = (D5-0)
                M4 = bst.Mixer('Mixer4',ins = [F3-1,D5-1], outs = ("S236"))
                M4_mx = bst.Mixer('Mixer4_mx',ins = M4-0, outs = ("Wax"))
                waxOut = (M4-0)

        elif scenario['Hydrocracking'] == "Yes":
            D4 = bst.units.BinaryDistillation('NaphthaSplitter', ins=D3-1,
                        outs=('S227',"S228"),
                        LHK =('C10H22','C14H30'),      #=('C10H22', 'C14H30'),
                        y_top=0.99, x_bot=0.01, k=2,
                        is_divided=True)  
            if scenario['residence_time'] == "low":
                D5 = bst.units.BinaryDistillation('DieselSplitter', ins=D4-1,
                                        outs=('S229','S230'),
                                        LHK =('C14H30','C24H50'),      #=('C10H22', 'C14H30'),
                                        y_top=0.99, x_bot=0.01, k=2,
                                        is_divided=True)
                M4 = bst.Mixer('Mixer4',ins = [F3-1,D5-1],outs=())

    #--------------------------------------------------------------------------------------------------------------
    # Hydrogen production 
    #--------------------------------------------------------------------------------------------------------------
    if scenario['Hydrocracking'] == "Yes":
        # with bst.System('sys_Hydrogen_Production') as sys_hydrogen_production:
        #     smr_catalyst = bst.Stream('H2_production_catalyst',Nickel_catalyst = 1 , units='kg/hr')
        #     smr_catalyst.price = prices["Hydrogen plant catalyst"]

        #     smr_pre = SteamReformer("Steam_Reformer", ins = (natural_gas,water,smr_catalyst),outs=())
        #     M5 = bst.units.Mixer('Mixer5',ins=(smr_pre-0,smr_pre-1),outs=())
        #     smr = bst.units.MixTank('Hydrogen_production', ins=(M5-0),outs=('H2andCO2'))
        #     hydrogenPSA = bst.units.Splitter("PSA",ins=(smr-0),outs=("","Other_Gases"), split={'H2':0.99,})  
    #     ["","Other gases"]
    #--------------------------------------------------------------------------------------------------------------
    # Hydrocracking Unit 
    #--------------------------------------------------------------------------------------------------------------
    # Hydrocracking catalyst calculation (reactor 1  + reactor 2)/(catalyst life year * feed (MT/day) * 365 * stream factor)
    # Hydrocracker details from Dutta et al. 2015. Liquid feed = 16,654 lb/hr;Total H2 feed = 2109 lb/hr; H2 purity = 90%, Make up H2 pure = 647 lb/hr,
    # Power in KW/hr for hydroprocessing 1811, hydrocracking electrucak cpnsumption 369kw, recycle compressor = 31
   
        with bst.System('sys_Hydrocracking') as sys_Hydrocracking:
            K3 = Compressor('Compressor3',ins=HC_hydrogen,outs=("S301"),P = 89.7 * 101325, eta=0.8) # 25 bars            
            hydrocracking_catalyst = bst.Stream('hydrocracking_catalyst',Zeolite = 1 , units='lb/hr')
            hydrocracking_catalyst.set_total_flow(200,"kg/hr")
            hydrocracking_catalyst.price = prices["Hydrocracking catalyst"]            
            
            if scenario['residence_time'] == "low":
                hydro_crack = Hydrocrack("Hydrocracking",ins=(M4-0,K3-0,hydrocracking_catalyst))
            else:
                hydro_crack = Hydrocrack("Hydrocracking",ins=(F3-1,K3-0,hydrocracking_catalyst))

            hydrocrack = bst.units.MixTank('Hydrocracking_Unit', ins=(hydro_crack-0,hydro_crack-1),outs=()) #hydrocracking_catalyst),outs=())
            splitter1 = bst.units.Splitter("H2split",ins=(hydrocrack-0),outs=("ExcessH2"), split={'H2':0.99,})
            H2excess = (splitter1-0)
            # cat_monitor = bst.units.Mixer('cat_monitor',ins=(hydrocracking_catalyst-0),outs=('outlet'))

            D6 = bst.units.BinaryDistillation('NaphthaSplitter2', ins=splitter1-1,
                                    outs=("S302",""),
                                    LHK =('C10H22','C14H30'),      #=('C10H22', 'C14H30'),
                                    y_top=0.99, x_bot=0.01, k=2,
                                    is_divided=True)
        
            D7 = bst.units.BinaryDistillation('DieselSplitter2', ins=D6-1,
                                    outs=("S303","S304"),
                                    LHK =('C14H30','C24H50'),      #=('C10H22', 'C14H30'),
                                    y_top=0.99, x_bot=0.01, k=2,
                                    is_divided=True)
            D7_mx = bst.Mixer('D7_MX',ins = D7-1,outs=("Wax"))

            waxOut = (D7-1)

# Consolidate product streams
        M6 = bst.units.Mixer('mix_Naphtha_out',ins=(D4-0,D6-0),outs=("Naphtha"))
        naphthaOut = (M6-0)

        if scenario['residence_time'] == "low":
            M7 = bst.units.Mixer('mix_Diesel_out',ins=(D5-0,D7-0),outs=("Diesel"))
            dieselOut = (M7-0)
        else:
            M7 = bst.units.Mixer('mix_Diesel_out',ins=(D4-1,D7-0),outs=("Diesel"))
            dieselOut = (M7-0)

    #  Create system
    sys = bst.main_flowsheet.create_system('sys')

    #  Include specifications
    def check_flash1():
        try:
            F1._run()
            # print("flash1 run successful")    
        except:
            # print("flash1 failed")
            pass
        top = F1.outs[0].copy()
        bottom = F1.outs[1].copy()
        flash_in = F1.ins[0].copy()
        for chem in top.available_chemicals:
            if str(chem) in ["O2","CO", "H2","CH4"] and chem.Tb < 273.15 + 10:
                F1.outs[0].imass[chem.ID] = flash_in.imass[chem.ID] * 0.999
                F1.outs[1].imass[chem.ID] = flash_in.imass[chem.ID] * (1-0.999)

            elif str(chem) in ["C2H4"] and chem.Tb < 273.15 + 10:
                F1.outs[0].imass[chem.ID] = flash_in.imass[chem.ID] * 0.989
                F1.outs[1].imass[chem.ID] = flash_in.imass[chem.ID] * (1-0.989)

            elif str(chem) in ["CO2"] and chem.Tb < 273.15 + 10:
                F1.outs[0].imass[chem.ID] = flash_in.imass[chem.ID] * 0.9877
                F1.outs[1].imass[chem.ID] = flash_in.imass[chem.ID] * (1-0.9877)

            elif str(chem) in ["C3H8"] and chem.Tb < 273.15 + 10:
                F1.outs[0].imass[chem.ID] = flash_in.imass[chem.ID] * 0.92
                F1.outs[1].imass[chem.ID] = flash_in.imass[chem.ID] * (1-0.92)
            
            elif str(chem) in ["C4H8"] and chem.Tb < 273.15 + 10:
                F1.outs[0].imass[chem.ID] = flash_in.imass[chem.ID] * 0.69
                F1.outs[1].imass[chem.ID] = flash_in.imass[chem.ID] * (1-0.69)

            elif str(chem) in ["C10H22"] and chem.Tb > 273.15 + 10:
                F1.outs[0].imass[chem.ID] = flash_in.imass[chem.ID] * (1-0.9964)
                F1.outs[1].imass[chem.ID] = flash_in.imass[chem.ID] * 0.9964

            elif str(chem) in ["C14H30","C24H50", "C40H82"] and chem.Tb >= 273.15 + 10:
                F1.outs[0].imass[str(chem)] = flash_in.imass[str(chem)] * 0.00
                F1.outs[1].imass[str(chem)] = flash_in.imass[str(chem)] * 1
            else:
                F1.outs[0].imass[str(chem)] = flash_in.imass[str(chem)] * 0.00
                F1.outs[1].imass[str(chem)] = flash_in.imass[str(chem)] * 1            
        F1.outs[0].P = F1.outs[1].P = flash_in.P
        F1.outs[0].T = F1.outs[1].T = 273.15 + 15     #  flash_in.T
        F1.outs[0].phase = 'g'
        F1.outs[1].phase = 'l'

    F1.add_specification(check_flash1)

    def check_flash2():
        try:
            for x in range(0,5):
                F2._run()    
        except:
            # print("flash2 failed")
            top = F2.outs[0].copy()
            bottom = F2.outs[1].copy()
            flash_in = F2.ins[0].copy()
            for chem in top.chemicals:
                if chem.Tb < (heater4-0).T: # 273.15 - 136:
                    F2.outs[0].imass[str(chem)] = flash_in.imass[str(chem)] * 0.99
                    F2.outs[1].imass[str(chem)] = flash_in.imass[str(chem)] * 0.01
                else:
                    F2.outs[0].imass[str(chem)] = flash_in.imass[str(chem)] * 0.01
                    F2.outs[1].imass[str(chem)] = flash_in.imass[str(chem)] * 0.99
            F2.outs[0].P = F2.outs[1].P = flash_in.P
            F2.outs[0].T = F2.outs[1].T = flash_in.T
            F2.outs[0].phase = 'g'
            F2.outs[1].phase = 'l'
    F2.add_specification(check_flash2)


    def check_flash3():
        try:
            F3._run()
            # print("F3 run successfully")
        except:
            # print ("F3 failed")
            pass
        # print("flash3 fixing outs")
        top = F3.outs[0].copy()
        bottom = F3.outs[1].copy()
        flash_in = F3.ins[0].copy()
        # F3.T = flash_in.T
        for chem in top.available_chemicals:
            if str(chem) in ["O2","CO","CO2","H2","CH4","C2H4", "C3H8","C4H8"] and chem.Tb < 580:
                F3.outs[0].imass[chem.ID] = flash_in.imass[chem.ID] * 1
                F3.outs[1].imass[chem.ID] = flash_in.imass[chem.ID] * 0

            elif str(chem) in ["C10H22","C14H30"] and chem.Tb < 580:
                F3.outs[0].imass[chem.ID] = flash_in.imass[chem.ID] * 0.95
                F3.outs[1].imass[chem.ID] = flash_in.imass[chem.ID] * 0.05

            else:
                F3.outs[0].imass[chem.ID] = flash_in.imass[chem.ID] * 0.05
                F3.outs[1].imass[chem.ID] = flash_in.imass[chem.ID] * 0.95
        F3.outs[0].P = F3.outs[1].P = flash_in.P
        F3.outs[0].T = F3.outs[1].T = flash_in.T
        F3.outs[0].phase = 'g'
        F3.outs[1].phase = 'l'
    F3.add_specification(check_flash3)
    #--------------------------------------------------------------------------------------------------------------
    # Hydrogen production
    #--------------------------------------------------------------------------------------------------------------
    if scenario['Hydrocracking'] == "Yes":

        hydrocracking_reaction = tmo.ParallelReaction([
        tmo.Reaction('C24H50 + 1.4H2 -> 2.4C10H22',  reactant='C24H50',  X=0.49),
        tmo.Reaction('C24H50 + 2.858H2 -> 2.4C14H30',  reactant='C24H50',  X=0.49),
        tmo.Reaction('C40H82 + 3H2 -> 4C10H22',  reactant='C40H82',  X=0.49),
        tmo.Reaction('C40H82 + 1.857H2 -> C14H30',  reactant='C40H82',  X=0.49),  
        ])  
        def hc_reaction():
            rxn=hydrocracking_reaction
            hydrocrack._run()
            rxn(hydrocrack.outs[0])
        hydrocrack.add_specification(hc_reaction)

        # def smr_reaction():
        #     rxn=tmo.Reaction('CH4 + H2O -> CO + 3H2',reactant = 'CH4',X=0.6)
        #     smr.run()
        #     rxn(smr.outs[0])

        # smr.add_specification(smr_reaction)

        def K3_ins():
            H2_req = 2*((0.0252 * hydro_crack.ins[0].get_flow("kg/hr", "C24H50")) +  (0.025 * hydro_crack.ins[0].get_flow("kg/hr", "C40H82")))   # from  Production of gasoline and diesel from biomass via fast pyrolysis, hydrotreating and hydrocracking: a design case... Jones et al 2009
            K3.ins[0].set_flow(H2_req,"kg/hr","H2")
            # water.set_flow(2*3 * H2_req *1/0.6 ,"kg/hr","H2O")    #= water_req multiply by 3 to produce excess H2
            # natural_gas.set_flow(2.67 * H2_req * 1/0.6,"kg/hr","CH4")  # multiply by 2.67 to produce excess H2
            # Ni_catalyst_req = H2_req   # catalyst cost Jones et al. 2009 3.6 $/1000 scf/28.3168466m3/885.71 kg of H2 {1scf = 0.0283168466 m3; 1000scf = 28.3168466 }
            # H2_req = 450 * D5.outs[1].get_total_flow("bbl/hr")* 0.18  #0.18 bbl = 1 scf  this is in barrels/hr
            # water.set_flow(3/3 * H2_req * 3.72965 * 1 ,"kg/hr","H2O")    #= water_req 0.6 representing 60% conversion
            # # natural_gas.set_flow(2.67/3 * H2_req * 3.72965 * 1,"kg/hr","CH4")     
            # Ni_catalyst_req = H2_req   # Confirm catalyst makeup
            # smr_catalyst.set_flow(Ni_catalyst_req,"kg/hr","Nickel_catalyst")
            K3._run()
        K3.add_specification(K3_ins)

    #--------------------------------------------------------------------------------------------------------------
    # Hydrocracking 
    #--------------------------------------------------------------------------------------------------------------
        #calculate the amount of catalyst needed for hydrocracking and hydrogen production

        # (12*Chemical("H2").MW)/(5*Chemical("C24H50").MW) = 0.07100591715976332
        # 1kg of C24H50 requires 0.07100591715976332 of H2
        # (4*Chemical("H2").MW)/(Chemical("C40H82").MW) = 0.014319973858809155
        # 1kg of C40H82 requires 0.014319973858809155 of H2
        #  Excess hydrogen is needed in the hydrocracking unit ina ratio of like 1 mole of HC to about 2-4 mole of H2. I will 
        # assume 3 moles of H2 for this study
        # mass of H2 needed for hydrocracking = 3 * (hydrocrack.ins[0].get_flow("kg/hr","C24H50")*0.07100591715976332 + hydrocrack.ins[0].get_flow("kg/hr","C40H82")*0.014319973858809155) =1025.2906444333576
        # 
        def catalyst_calc():
            catalyst_density = 420 # lb/ft3
            reactor_volume = 0.65 # LHSV
            stream_flow = hydro_crack.ins[0].get_total_flow("gal/hr")
            catalyst_life = 1    # catalyst lifetime years
            stream_factor = 0.9  # stream factor
            reactor_one = catalyst_density/(reactor_volume * 7.4802) * stream_flow       # reactor 1 volume m
            reactor_two = 0
            # print (f"{feed_mass}")
            # hydrocracking_catalyst_req = (reactor_one +reactor_two)/(catalyst_life * feed_mass * 365 * stream_factor) # lb/MT of feed
            # hydrocracking_catalyst.set_flow(hydrocracking_catalyst_req * feed.get_total_flow("tonnes/hr"),"lb/hr","Zeolite")

            hydrocracking_catalyst_req = (reactor_one +reactor_two)/(catalyst_life * feed_mass * 365 * stream_factor) # lb/MT of feed
            # hydrocracking_catalyst_req = 6.7
            hydrocracking_catalyst.set_flow(hydrocracking_catalyst_req * feed.get_total_flow("tonnes/hr"),"lb/hr","Zeolite")
            #             
            # hydrocracking_catalyst.set_flow(25,"lb/hr","Zeolite")
            #    # catalyst requirement in lb/hr
            hydro_crack.run()
        hydro_crack.add_specification(catalyst_calc)            

        # print (f"{hydrocracking_catalyst.get_total_flow('lb/hr')} is the catalyst flow rate")
        #--------------------------------------------------------------------------------------------------------------
 


    for stream in sys.products:
        try:
            stream.price = prices[str(stream)]
        except:
            pass

    sys.simulate()
    #  Set the price of hydrocracking and hydrogen production
    if scenario['Hydrocracking'] == "Yes":

        # D6.purchase_costs = {"NaphthaSplitter2": 0}
        # D6.installed_costs = {"NaphthaSplitter2":0}
        # D7.purchase_costs = {"DieselSplitter2": 0}
        # D7.installed_costs = {"DieselSplitter2":0}

        hydrocrack.purchase_costs = {"Hydrocrackng_Unit_mix":0} 
        hydrocrack.installed_costs = {'Tank':0}

    # Check if the two flash units are unable to set the price of Heat exchanger - Floating head unit in the flash. I am scaling the cost based on the input   
    flash_probs = [F1,F2]
    for unit in flash_probs:
        if math.isnan(unit.purchase_cost):
            unit.purchase_costs['Heat exchanger - Floating head'] = 4698.2 * (unit.ins[0].get_total_flow("kmol/hr")/22)**0.65    
            unit.installed_costs['Heat exchanger - Floating head'] = unit.F_BM['Heat exchanger - Floating head'] * unit.purchase_costs['Heat exchanger - Floating head']
    #   Heat exchanger units also have issues estimating the floating head price so i approximate 
    prob_units = [cooler1, H3, H4, H5]
    for unit in prob_units:
        if math.isnan(unit.purchase_cost):
            unit.purchase_costs['Floating head'] = 4698.2 * (unit.ins[0].get_total_flow("kmol/hr")/22)**0.65    
    if math.isnan(D3.installed_cost):
        D3.installed_costs = {}
        D3.installed_costs = {"Distillation_column":324295 * (D3.ins[0].get_total_flow("kg/hr")/685)**0.65}
        D3.purchase_costs = {}
        D3.purchase_costs = {"Distillation_column":324295 * (D3.ins[0].get_total_flow("kg/hr")/685)**0.65}

    return sys
# %%
