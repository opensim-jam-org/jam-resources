import os
import json

import opensim as osim

models_path = os.path.abspath('./model')

model_file = os.path.join(models_path, 'scaled_lenhart_model_updated_wrap.osim')
results_basename = "cycling"

external_loads_file = os.path.join(models_path, 'P017', 'external_loads_pedal.xml')

coordinates_file = './results/comak-inverse-kinematics/cycling_ik.mot'

comak_result_dir = './results/comak'
results_basename = "cycling"
save_xml_path = './inputs/comak_settings.xml'

with open('./prescribed_coordinates.json', 'r') as f:
    prescribed_coordinates = json.load(f)

with open('./primary_coordinates.json', 'r') as f:
    primary_coordinates = json.load(f)

with open('./secondary_coordinates.json', 'r') as f:
    secondary_coordinates = json.load(f)


# Settings
start_time = 120.57 #120.60 #120.61
start_pad = 0.
stop_time = 121.26#120.66

comak = osim.COMAKTool();
comak.set_model_file(model_file);
comak.set_coordinates_file(coordinates_file);
comak.set_external_loads_file(external_loads_file);
comak.set_results_directory(comak_result_dir);
comak.set_results_prefix(results_basename);
comak.set_replace_force_set(False);
comak.set_start_time(start_time - start_pad);
comak.set_stop_time(stop_time);
comak.set_time_step(1/112.5);
comak.set_lowpass_filter_frequency(6);
comak.set_print_processed_input_kinematics(False);

for coordinate_number, (key, path) in enumerate(prescribed_coordinates.items()):
    comak.set_prescribed_coordinates(int(coordinate_number), path)

for coord_number, (key, path) in enumerate(primary_coordinates.items()):
    comak.set_primary_coordinates(int(coord_number), path)

secondary_coord_set = osim.COMAKSecondaryCoordinateSet();
secondary_coord = osim.COMAKSecondaryCoordinate();

for coord, dict_ in secondary_coordinates.items():
    secondary_coord.setName(coord)
    secondary_coord.set_max_change(dict_['max_change']);
    secondary_coord.set_coordinate(dict_['coordinate']);
    secondary_coord_set.cloneAndAppend(secondary_coord);

comak.set_COMAKSecondaryCoordinateSet(secondary_coord_set);


comak.set_settle_secondary_coordinates_at_start(True);
comak.set_settle_threshold(1e-3);
comak.set_settle_accuracy(1e-2);
comak.set_settle_internal_step_limit(10000);
comak.set_print_settle_sim_results(True);
comak.set_settle_sim_results_directory(comak_result_dir);
comak.set_settle_sim_results_prefix("cycling_settle_sim");
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
comak.set_non_muscle_actuator_weight(1);  # The higher this is the more costly it is to use actuators
comak.set_model_assembly_accuracy(1e-12);
comak.set_use_visualizer(False);
# // comak.set_verbose(2)
comak.setDebugLevel(1);
# comak.print(save_xml_path);
comak.printToXML(save_xml_path);

print("Starting COMAK Tool!")
comak.run()
print('Finished COMAK Tool!')