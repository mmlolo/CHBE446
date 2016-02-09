# -*- coding: utf-8 -*-

def powerreq(L,x):
    "Initialize constants"    
    Mw_L = 42.39
    Mw_H = 18.02
    T_F = 313.71
    T_I = 295.37
    
    "Average molecular weight"
    Mw_avg = x*Mw_L + (1-x)*Mw_H
    
    "Converted value of heat capacity to mole basis (KJ/(mol*K)"
    Cp = 3.57/1000 * Mw_avg
    
    "Power in KW"
    
    P = L*Cp*(T_F-T_I)
    return[P]
