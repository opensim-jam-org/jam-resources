# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 13:54:48 2021

@author: csmith
"""
import paraview.simple as pv

import os

from scipy.spatial.transform import Rotation as R

import numpy as np



# Rotations to model reference frame

femur_T=np.array([
    [-0.996637, 0.0819236, 0.0014398, 0.00800132], 
    [0.0819363, 0.996486, 0.0173925, -0.418005],
    [-9.8835e-06, 0.017452, -0.999848, -0.0039984], 
    [0, 0, 0, 1]])

tibia_T=np.array([
    [-0.997565, 0.0696635, 0.00362601, 0.010001 ],
    [0.0697579, 0.996198, 0.0522061, -0.0400074 ],
    [2.46729e-05, 0.0523319, -0.998631, 0.00399948], 
    [0, 0, 0, 1]])

    
femur_r = femur_T[0:3,0:3]
femur_t =femur_T[0:3,3]


    
tibia_r = tibia_T[0:3,0:3]
tibia_t = tibia_T[0:3,3]

femur_R = R.from_matrix(femur_r)
tibia_R = R.from_matrix(tibia_r)
patella_R = R.from_euler('y',180,degrees=True)

femur_angles = femur_R.as_euler('xyz',True)
tibia_angles = tibia_R.as_euler('xyz',True)
patella_angles = patella_R.as_euler('xyz',True)

geometry_folder = "P:\\Projects\\LMB_grand_challenge_database\\read_only\\DataforSixthCompetition\\Geometry Data\\"
out_geometry_folder = "C:\\Users\\csmith\\github\\jam-resources\\models\\knee_tka\\grand_challenge\\DM\\new\\Geometry\\"

femur_files = ["Femur.stl","Femoral Component.stl"]
tibia_files = ["Tibia.stl","Tibial Insert.stl","Tibial Tray.stl","Fibula.stl"]
patella_files = ["Patella.stl","Patellar Button.stl"]

local_files = ["Femoral Component.stl","Tibial Insert.stl","Tibial Tray.stl","Patellar Button.stl"]
num_decimate_iterations = 5
for entry in os.scandir(geometry_folder):
    if entry.name.endswith(".stl"):  
        print(entry.name)
        stl = pv.OpenDataFile(entry.path)
        pv.UpdatePipeline()
        print(stl.GetDataInformation().GetNumberOfCells())
        
        T = pv.Transform(stl)
        T.Transform.Scale = [0.001,0.001,0.001]
        
        if entry.name in femur_files:
            TT = pv.Transform(T)
            TT.Transform.Rotate = femur_angles
            TT.Transform.Translate = femur_t
            
        if entry.name in tibia_files:
            TT = pv.Transform(T)
            TT.Transform.Rotate = tibia_angles
            TT.Transform.Translate = tibia_t
            
        if entry.name in patella_files:
            TT = pv.Transform(T)
            TT.Transform.Rotate = patella_angles
            
        if entry.name in local_files:
            TT = T
                    
        pv.SaveData(out_geometry_folder + entry.name, TT)

