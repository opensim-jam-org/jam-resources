import os
import json

import opensim as osim

matlab_path = os.path.abspath('../../../../matlab/')
models_path = os.path.abspath('../../../../models')


# model_file = "/root/jam-resources/models/knee_healthy/lenhart2015/lenhart2015.osim";
model_file = os.path.join(models_path, 'knee_healthy', 'lenhart2015', 'lenhart2015.osim')
# coordinates_file = "/root/jam-resources/python/examples/Walking/Walking/results/comak-inverse-kinematics/overground_17_ik.mot";
coordinates_file = './results/comak-inverse-kinematics/overground_17_ik.mot'
# external_loads_file = "/root/jam-resources/python/examples/Walking/Walking/overground_17_ext_loads.xml";
external_loads_file = os.path.join(models_path, 'knee_healthy', 'experimental_data', 'motion_analysis', 'overground_17_ext_loads.xml')
# forceset_file =  "/root/jam-resources/models/knee_healthy/lenhart2015/lenhart2015_reserve_actuators.xml";
forceset_file = os.path.join(models_path, 'knee_healthy', 'lenhart2015', 'lenhart2015_reserve_actuators.xml')
# comak_result_dir = "/root/jam-resources/python/examples/Walking/Walking/results/comak";
comak_result_dir = './results/comak'
results_basename = "walking";
# save_xml_path = "/root/jam-resources/python/examples/Walking/Walking/inputs/comak_settings2.xml";
save_xml_path = './inputs/comak_settings.xml'

comak = osim.COMAKTool();
comak.set_model_file(model_file);
comak.set_coordinates_file(coordinates_file);
comak.set_external_loads_file(external_loads_file);
comak.set_results_directory(comak_result_dir);
comak.set_results_prefix(results_basename);
comak.set_replace_force_set(False);
comak.set_force_set_file(forceset_file);
comak.set_start_time(1.16);
comak.set_stop_time(2.32);
comak.set_time_step(0.01);
comak.set_lowpass_filter_frequency(6);
comak.set_print_processed_input_kinematics(False);
comak.set_prescribed_coordinates(0,"/jointset/gnd_pelvis/pelvis_tx");
comak.set_prescribed_coordinates(1,"/jointset/gnd_pelvis/pelvis_ty");
comak.set_prescribed_coordinates(2,"/jointset/gnd_pelvis/pelvis_tz");
comak.set_prescribed_coordinates(3,"/jointset/gnd_pelvis/pelvis_tilt");
comak.set_prescribed_coordinates(4,"/jointset/gnd_pelvis/pelvis_list");
comak.set_prescribed_coordinates(5,"/jointset/gnd_pelvis/pelvis_rot");
comak.set_prescribed_coordinates(6,"/jointset/subtalar_r/subt_angle_r");
comak.set_prescribed_coordinates(7,"/jointset/mtp_r/mtp_angle_r");
comak.set_prescribed_coordinates(8,"/jointset/hip_l/hip_flex_l");
comak.set_prescribed_coordinates(9,"/jointset/hip_l/hip_add_l");
comak.set_prescribed_coordinates(10,"/jointset/hip_l/hip_rot_l");
comak.set_prescribed_coordinates(11,"/jointset/pf_l/pf_l_r3");
comak.set_prescribed_coordinates(12,"/jointset/pf_l/pf_l_tx");
comak.set_prescribed_coordinates(13,"/jointset/pf_l/pf_l_ty");
comak.set_prescribed_coordinates(14,"/jointset/knee_l/knee_flex_l");
comak.set_prescribed_coordinates(15,"/jointset/ankle_l/ankle_flex_l");
comak.set_prescribed_coordinates(16,"/jointset/subtalar_l/subt_angle_l");
comak.set_prescribed_coordinates(17,"/jointset/mtp_l/mtp_angle_l");
comak.set_prescribed_coordinates(18,"/jointset/pelvis_torso/lumbar_ext");
comak.set_prescribed_coordinates(19,"/jointset/pelvis_torso/lumbar_latbend");
comak.set_prescribed_coordinates(20,"/jointset/pelvis_torso/lumbar_rot");
comak.set_prescribed_coordinates(21,"/jointset/torso_neckhead/neck_ext");
comak.set_prescribed_coordinates(22,"/jointset/torso_neckhead/neck_latbend");
comak.set_prescribed_coordinates(23,"/jointset/torso_neckhead/neck_rot");
comak.set_prescribed_coordinates(24,"/jointset/acromial_r/arm_add_r");
comak.set_prescribed_coordinates(25,"/jointset/acromial_r/arm_flex_r");
comak.set_prescribed_coordinates(26,"/jointset/acromial_r/arm_rot_r");
comak.set_prescribed_coordinates(27,"/jointset/elbow_r/elbow_flex_r");
comak.set_prescribed_coordinates(28,"/jointset/radioulnar_r/pro_sup_r");
comak.set_prescribed_coordinates(29,"/jointset/radius_hand_r/wrist_flex_r");
comak.set_prescribed_coordinates(30,"/jointset/acromial_l/arm_add_l");
comak.set_prescribed_coordinates(31,"/jointset/acromial_l/arm_flex_l");
comak.set_prescribed_coordinates(32,"/jointset/acromial_l/arm_rot_l");
comak.set_prescribed_coordinates(33,"/jointset/elbow_l/elbow_flex_l");
comak.set_prescribed_coordinates(34,"/jointset/radioulnar_l/pro_sup_l");
comak.set_prescribed_coordinates(35,"/jointset/radius_hand_l/wrist_flex_l");

comak.set_primary_coordinates(0,"/jointset/hip_r/hip_flex_r");
comak.set_primary_coordinates(1,"/jointset/hip_r/hip_add_r");
comak.set_primary_coordinates(2,"/jointset/hip_r/hip_rot_r");
comak.set_primary_coordinates(3,"/jointset/knee_r/knee_flex_r");
comak.set_primary_coordinates(4,"/jointset/ankle_r/ankle_flex_r");

secondary_coord_set = osim.COMAKSecondaryCoordinateSet();
secondary_coord = osim.COMAKSecondaryCoordinate();

secondary_coord.setName("knee_add_r");
secondary_coord.set_max_change(0.01);
secondary_coord.set_coordinate("/jointset/knee_r/knee_add_r");
secondary_coord_set.cloneAndAppend(secondary_coord);

secondary_coord.setName("knee_rot_r");
secondary_coord.set_max_change(0.01);
secondary_coord.set_coordinate("/jointset/knee_r/knee_rot_r");
secondary_coord_set.cloneAndAppend(secondary_coord);

secondary_coord.setName("knee_tx_r");
secondary_coord.set_max_change(0.05);
secondary_coord.set_coordinate("/jointset/knee_r/knee_tx_r");
secondary_coord_set.cloneAndAppend(secondary_coord);

secondary_coord.setName("knee_ty_r");
secondary_coord.set_max_change(0.05);
secondary_coord.set_coordinate("/jointset/knee_r/knee_ty_r");
secondary_coord_set.cloneAndAppend(secondary_coord);

secondary_coord.setName("knee_tz_r");
secondary_coord.set_max_change(0.05);
secondary_coord.set_coordinate("/jointset/knee_r/knee_tz_r");
secondary_coord_set.cloneAndAppend(secondary_coord);

secondary_coord.setName("pf_flex_r");
secondary_coord.set_max_change(0.01);
secondary_coord.set_coordinate("/jointset/pf_r/pf_flex_r");
secondary_coord_set.cloneAndAppend(secondary_coord);

secondary_coord.setName("pf_rot_r");
secondary_coord.set_max_change(0.01);
secondary_coord.set_coordinate("/jointset/pf_r/pf_rot_r");
secondary_coord_set.cloneAndAppend(secondary_coord);

secondary_coord.setName("pf_tilt_r");
secondary_coord.set_max_change(0.01);
secondary_coord.set_coordinate("/jointset/pf_r/pf_tilt_r");
secondary_coord_set.cloneAndAppend(secondary_coord);

secondary_coord.setName("pf_tx_r");
secondary_coord.set_max_change(0.005);
secondary_coord.set_coordinate("/jointset/pf_r/pf_tx_r");
secondary_coord_set.cloneAndAppend(secondary_coord);

secondary_coord.setName("pf_ty_r");
secondary_coord.set_max_change(0.005);
secondary_coord.set_coordinate("/jointset/pf_r/pf_ty_r");
secondary_coord_set.cloneAndAppend(secondary_coord);

secondary_coord.setName("pf_tz_r");
secondary_coord.set_max_change(0.005);
secondary_coord.set_coordinate("/jointset/pf_r/pf_tz_r");
secondary_coord_set.cloneAndAppend(secondary_coord);

comak.set_COMAKSecondaryCoordinateSet(secondary_coord_set);

comak.set_settle_secondary_coordinates_at_start(True);
comak.set_settle_threshold(1e-3);
comak.set_settle_accuracy(1e-2);
comak.set_settle_internal_step_limit(10000);
comak.set_print_settle_sim_results(True);
comak.set_settle_sim_results_directory(comak_result_dir);
comak.set_settle_sim_results_prefix("walking_settle_sim");
comak.set_max_iterations(25);
comak.set_udot_tolerance(1);
comak.set_udot_worse_case_tolerance(50);
comak.set_unit_udot_epsilon(1e-6);
comak.set_optimization_scale_delta_coord(1);
comak.set_ipopt_diagnostics_level(3);
comak.set_ipopt_max_iterations(500);
comak.set_ipopt_convergence_tolerance(1e-4);
comak.set_ipopt_constraint_tolerance(1e-4);
comak.set_ipopt_limited_memory_history(200);
comak.set_ipopt_nlp_scaling_max_gradient(10000);
comak.set_ipopt_nlp_scaling_min_value(1e-8);
comak.set_ipopt_obj_scaling_factor(1);
comak.set_activation_exponent(2);
comak.set_contact_energy_weight(0);
comak.set_non_muscle_actuator_weight(1000);
comak.set_model_assembly_accuracy(1e-12);
comak.set_use_visualizer(False);
# // comak.set_verbose(2)
comak.setDebugLevel(1);
# comak.print(save_xml_path);
comak.printToXML(save_xml_path);

print("Starting COMAK Tool!")
comak.run()
print('Finished COMAK Tool!')