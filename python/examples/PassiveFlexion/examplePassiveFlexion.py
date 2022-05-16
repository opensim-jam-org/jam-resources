# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 09:34:41 2021

@author: csmith
"""
# =============================================================================
# import sys
# osim_path = 'C:\\Users\\csmith\\github\\opensim-core-jam\\source-install\\sdk\\Python'
# sys.path.insert(0,osim_path)
# =============================================================================

import os
#os.add_dll_directory("C:\\github\\opensim-core-jam\\source-install\\bin")

import opensim as osim
import numpy as np

osim.Logger.setLevelString("trace")
useVisualizer=True

model_file = "../../../models/knee_healthy/lenhart2015/lenhart2015.osim";
results_basename = "passive_flexion";
forsim_result_dir = "./results/forsim";
jnt_mech_result_dir = "./results/joint-mechanics";

### Create Input Files

os.makedirs("inputs", exist_ok=True)
os.makedirs("results", exist_ok=True)
os.makedirs("graphics", exist_ok=True)

prescribed_coord_file = "./inputs/prescribed_coordinates.sto"

# Simulation Time
# Simulation consists of two phases:
#   settle : allow unconstraind knee DOFs to settle into equilbrium 
#   flex : prescribe the tibiofemoral flexion
# All time units are in seconds 

time_step = 0.01;

settle_duration = 0.5; 
flex_duration = 2.0; 

settle_time = np.arange(0,settle_duration,time_step)
flex_time = np.arange(settle_duration,flex_duration + settle_duration, time_step)

time = np.concatenate((settle_time, flex_time));

time_points = [0, settle_duration, settle_duration + flex_duration];

nSettleSteps = settle_time.shape[0];
nFlexSteps = flex_time.shape[0];
nTimeSteps = nSettleSteps + nFlexSteps

# Prescribe Knee Flexion
max_knee_flex = 90;

knee_flex = [0,0,max_knee_flex]

#would be better to use smooth spline interpolation here instead of linear
smooth_knee_flex = np.reshape(np.interp(time,time_points,knee_flex),(nTimeSteps,1))


pelvis_tilt = np.ones((nSettleSteps+nFlexSteps,1))*90;

time_std = osim.StdVectorDouble()
for t in np.nditer(time):
    time_std.push_back(float(t))

data_array = np.concatenate((smooth_knee_flex,pelvis_tilt),1)

data_matrix = osim.Matrix.createFromMat(data_array)

labels = osim.StdVectorString()

# labels.append("time")
labels.append("knee_flex_r")
labels.append("pelvis_tilt")

prescribed_coord_table = osim.TimeSeriesTable(time, data_matrix, labels)

sto = osim.STOFileAdapter();
sto.write(prescribed_coord_table,prescribed_coord_file);

# %% Plot Simulation Inputs 
# coord_fig = figure('name','prescribed_coordinates','Position',  [100, 100, 667, 300]);

# subplot(1,2,1);
# plot(time,prescribed_coord_data.pelvis_tilt,'LineWidth',2)
# ylim([0.0 100])
# xlabel('Time [s]')
# ylabel('Angle [^o]')
# title('Pelvis Tilt (pelvis\_tilt\_r)')
# box off

# subplot(1,2,2);
# plot(time,prescribed_coord_data.knee_flex_r,'LineWidth',2)
# ylim([0.0 100])
# xlabel('Time [s]')
# ylabel('Angle [^o]')
# title('Knee Flexion (knee\_flex\_r)')
# box off

# saveas(coord_fig,'./results/graphics/prescribed_coordinates.png')

### Perform Simulation with ForsimTool
forsim = osim.ForsimTool()
forsim.set_model_file(model_file)
forsim.set_results_directory(forsim_result_dir)
forsim.set_results_file_basename(results_basename)
forsim.set_start_time(-1)
forsim.set_stop_time(-1)
forsim.set_integrator_accuracy(1e-2) # Note this should be 1e-6 for research
forsim.set_constant_muscle_control(0.02) # Set all muscles to 2% activation to represent passive state
# forsim.set_ignore_activation_dynamics(True)
# forsim.set_ignore_tendon_compliance(True)
forsim.set_unconstrained_coordinates(0,'/jointset/knee_r/knee_add_r')
forsim.set_unconstrained_coordinates(1,'/jointset/knee_r/knee_rot_r')
forsim.set_unconstrained_coordinates(2,'/jointset/knee_r/knee_tx_r')
forsim.set_unconstrained_coordinates(3,'/jointset/knee_r/knee_ty_r')
forsim.set_unconstrained_coordinates(4,'/jointset/knee_r/knee_tz_r')
forsim.set_unconstrained_coordinates(5,'/jointset/pf_r/pf_flex_r')
forsim.set_unconstrained_coordinates(6,'/jointset/pf_r/pf_rot_r')
forsim.set_unconstrained_coordinates(7,'/jointset/pf_r/pf_tilt_r')
forsim.set_unconstrained_coordinates(8,'/jointset/pf_r/pf_tx_r')
forsim.set_unconstrained_coordinates(9,'/jointset/pf_r/pf_ty_r')
forsim.set_unconstrained_coordinates(10,'/jointset/pf_r/pf_tz_r')
forsim.set_prescribed_coordinates_file(prescribed_coord_file)
forsim.set_use_visualizer(True)
forsim.printToXML('./inputs/forsim_settings.xml')

print('Running Forsim Tool...')
# forsim.run();

### Perform Analysis with JointMechanicsTool

jnt_mech = osim.JointMechanicsTool()
jnt_mech.set_model_file(model_file)
jnt_mech.set_input_states_file(forsim_result_dir + '/' + results_basename + '_states.sto')
jnt_mech.set_results_file_basename(results_basename)
jnt_mech.set_results_directory(jnt_mech_result_dir)
jnt_mech.set_start_time(-1)
jnt_mech.set_stop_time(-1)
jnt_mech.set_normalize_to_cycle(False)
jnt_mech.set_contacts(0,'all')
jnt_mech.set_ligaments(0,'all')
jnt_mech.set_muscles(0,'none')
jnt_mech.set_muscle_outputs(0,'none')
jnt_mech.set_attached_geometry_bodies(0,'/bodyset/femur_distal_r')
jnt_mech.set_attached_geometry_bodies(1,'/bodyset/tibia_proximal_r')
jnt_mech.set_attached_geometry_bodies(2,'/bodyset/patella_r')
jnt_mech.set_output_orientation_frame('ground')
jnt_mech.set_output_position_frame('ground')
jnt_mech.set_write_vtp_files(True)
jnt_mech.set_write_h5_file(True)
jnt_mech.set_h5_kinematics_data(True)
jnt_mech.set_h5_states_data(True)
jnt_mech.set_use_visualizer(useVisualizer)
jnt_mech.printToXML("./inputs/joint_mechanics_settings.xml")

print('Running JointMechanicsTool...')
jnt_mech.run()
