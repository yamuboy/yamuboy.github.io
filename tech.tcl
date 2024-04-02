# A set of layer thicknesses for ADS DemoKit2 that matches the cross
# section used by EM simulation of ADS DemoKit2. 

proc rephbtTech { flow }  { 

    gdSetVar verbosity 1   

    # temperatures are in deg C

    set GaAs_K_T {       0 51.72   10 49.45   20 47.35   30 45.40   40 43.60 
                        50 41.92   60 40.35   70 38.89   80 37.52   90 36.23 
                       100 35.02  110 33.88  120 32.81  130 31.79  140 30.83
                       150 29.93  160 29.06  170 28.25  180 27.47  190 26.73
                       200 26.03  210 25.35  220 24.71  230 24.10  240 23.52 }

    set Si3N4_K_T {      0 17.94   10 18.15   20 18.36   30 18.56   40 18.76
                        50 18.96   60 19.15   70 19.34   80 19.52   90 19.70
                       100 19.88  110 20.06  120 20.23  130 20.40  140 20.56
                       150 20.72  160 20.88  170 21.04  180 21.20  190 21.35
                       200 21.50  210 21.65  220 21.80  230 21.94  240 22.09 }
 
    #---------------------------------
    # TECHNOLOGY: MATERIAL DEFINITIONS 
    #---------------------------------
    # THERMAL CONDUCTIVIY
    # { name    thermal-conductivity in W/(Km) [metal e- mean free path]} 
    set  materialKList "
        { GaAs      {$GaAs_K_T}    }  
        { Si3N4     {$Si3N4_K_T}   }
        { NiCr      11.3 }
        { Au       317.0   4.10e-8 }
		{ Cu       401.0  }
        { Air        0.027         }
		{ PI		0.52		   } 
    "                          ;# using quotes, so $GaAs_K et al are evaluated

    # VOLUMETRIC HEAT CAPACITY (use material names from materialKList)
    # { name   volumetric-heat-capacity in J/(K*m^3) = c*d }
    set materialCvList {   
        { GaAs        1.73e6 }  
        { Si3N4       2.75e6 }  
        { NiCr        3.78e6 }
        { Au          2.49e6 }  
		{ Cu          3.45e6 }  
        { Air         1.18e3 }
		{ PI		  1.64e6 }
    }                          ;# no $variables used, so braces {} are ok

    #---------------------------------------
    # TECHNOLOGY: THERMAL LAYER  DEFINITIONS
    #---------------------------------------
    # {name thickness(m) bkgnd-material {layer1 material1}...{lyrN matN}}
    set  layerList {        
      { substrate   100.000e-6   GaAs                         {BACKVIA Au} }
      { cmesa        0.450e-6   Si3N4 {CMESA GaAs} {CMET Au} }
      { bmesa        1.350e-6   Si3N4 {BMESA GaAs} {CMET Au} {VIA1 Au} }
      { emesa        0.200e-6   Si3N4 {EMESA GaAs} {BMET Au} {VIA1 Au} }
      { emet         0.035e-6   Si3N4              {EMET Au} {VIA1 Au} }
      { via1         0.200e-6   Si3N4                        {VIA1 Au} }
      { tfr          0.030e-6   Si3N4             {TFR NiCr} {VIA1 Au} }
      { met1         1.000e-6   Si3N4                        {VIA2 Au} {MET1 Au} }
      { via2         0.100e-6   Si3N4                        {VIA2 Au}           }
      { cgr          0.300e-6   Si3N4                        {CGR  Au}           }
      { met2         4.000e-6   Si3N4                        {MET2 Au}           }
      { protect      4e-6   	Si3N4                        {PROTECT Air}       }
    }


    # in the technology definition above, each sub-list within layerList
    # defines a horizontal slice called a thermal layer. The syntax of 
    # the sublist is defined above. "name" is the thermal-layer-name.
    # Zero or more mask-layers (= layout-layers) may be associated with
    # the thermal-layer. The "bkgnd-material" will be assumed to exist wherever
    # there are no shapes (geometries) on the mask-layers. The {layer material}
    # pairs define what material each shape on a given layer represents.
    # The first {layer material} pair has the highest precedence. The 
    # bkgnd-material has the lowest precedence.
    # Note that the mask-layer names are those in the stream layer table

    gdCxtInfo "defined K and Cv of materials in stack" ;# context info 
    gdDefineTech flow $materialKList $layerList $materialCvList 
}
