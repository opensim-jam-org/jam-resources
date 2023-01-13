import os
import json

import opensim as osim

matlab_path = os.path.abspath('../../../../matlab/')
models_path = os.path.abspath('../../../../models')


model_file = os.path.join(models_path, 'knee_healthy', 'lenhart2015', 'lenhart2015.osim')
coordinates_file = './results/comak-inverse-kinematics/overground_17_ik.mot'
external_loads_file = os.path.join(models_path, 'knee_healthy', 'experimental_data', 'motion_analysis', 'overground_17_ext_loads.xml')
forceset_file = os.path.join(models_path, 'knee_healthy', 'lenhard2015', 'lenhart2015_reserve_actuators.xml')
results_basename = "walking"

markerset_file = os.path.join(models_path, 'knee_healthy/experimental_data/motion_analysis/overground_17.trc')

secondary_constraint_function_file = './results/comak-inverse-kinematics/secondary_coordinate_constraint_functions.xml'
constrained_model_file = './results/comak-inverse-kinematics/ik_constrained_model.osim'

ik_result_dir = './results/comak-inverse-kinematics'

save_xml_path = './inputs/comak_inverse_kinematics_settings.xml'

os.makedirs('./inputs', exist_ok=True)
os.makedirs('./results/graphics', exist_ok=True)

with open('./marker_weights.json', 'r') as f:
    marker_weights = json.load(f)

# with open('./prescribed_coordinates.json', 'r') as f:
#     prescribed_coordinates = json.load(f)

# with open('./primary_coordinates.json', 'r') as f:
#     primary_coordinates = json.load(f)

with open('./secondary_coordinates.json', 'r') as f:
    secondary_coordinates = json.load(f)



# Settings
perform_secondary_constraint_sim = True
secondary_constraint_sim_settle_threshold = 1e-4
secondary_constraint_sim_sweep_time = 3.0
secondary_coupled_coordinate_start_value = 0
secondary_coupled_coordinate_stop_value = 100
secondary_constraint_sim_integrator_accuracy = 1e-2
secondary_constraint_sim_internal_step_limit = 10000
constraint_function_num_interpolation_points = 20
print_secondary_constraint_sim_results = True
perform_inverse_kinematics = True

report_errors = True
report_marker_locations = False
ik_constraint_weight = 100
ik_accuracy = 1e-5
use_visualizer = False

comak_ik = osim.COMAKInverseKinematicsTool();
comak_ik.set_model_file(model_file);
comak_ik.set_results_directory(ik_result_dir);
comak_ik.set_results_prefix(results_basename);
comak_ik.set_perform_secondary_constraint_sim(perform_secondary_constraint_sim);
comak_ik.set_secondary_coordinates(0,'/jointset/knee_r/knee_add_r');
comak_ik.set_secondary_coordinates(1,'/jointset/knee_r/knee_rot_r');
comak_ik.set_secondary_coordinates(2,'/jointset/knee_r/knee_tx_r');
comak_ik.set_secondary_coordinates(3,'/jointset/knee_r/knee_ty_r');
comak_ik.set_secondary_coordinates(4,'/jointset/knee_r/knee_tz_r');
comak_ik.set_secondary_coordinates(5,'/jointset/pf_r/pf_flex_r');
comak_ik.set_secondary_coordinates(6,'/jointset/pf_r/pf_rot_r');
comak_ik.set_secondary_coordinates(7,'/jointset/pf_r/pf_tilt_r');
comak_ik.set_secondary_coordinates(8,'/jointset/pf_r/pf_tx_r');
comak_ik.set_secondary_coordinates(9,'/jointset/pf_r/pf_ty_r');
comak_ik.set_secondary_coordinates(10,'/jointset/pf_r/pf_tz_r');
comak_ik.set_secondary_coupled_coordinate('/jointset/knee_r/knee_flex_r');
comak_ik.set_secondary_constraint_sim_settle_threshold(secondary_constraint_sim_settle_threshold);
comak_ik.set_secondary_constraint_sim_sweep_time(secondary_constraint_sim_sweep_time);
comak_ik.set_secondary_coupled_coordinate_start_value(secondary_coupled_coordinate_start_value);
comak_ik.set_secondary_coupled_coordinate_stop_value(secondary_coupled_coordinate_stop_value);
comak_ik.set_secondary_constraint_sim_integrator_accuracy(secondary_constraint_sim_integrator_accuracy);
comak_ik.set_secondary_constraint_sim_internal_step_limit(secondary_constraint_sim_internal_step_limit);
comak_ik.set_secondary_constraint_function_file(secondary_constraint_function_file);
comak_ik.set_constraint_function_num_interpolation_points(constraint_function_num_interpolation_points);
comak_ik.set_print_secondary_constraint_sim_results(print_secondary_constraint_sim_results);
comak_ik.set_constrained_model_file(constrained_model_file);
comak_ik.set_perform_inverse_kinematics(perform_inverse_kinematics);
comak_ik.set_marker_file(markerset_file);

comak_ik.set_output_motion_file('overground_17_ik.mot');
comak_ik.set_time_range(0, 0);
comak_ik.set_time_range(1, 2.36);
comak_ik.set_report_errors(report_errors);
comak_ik.set_report_marker_locations(report_marker_locations);
comak_ik.set_ik_constraint_weight(ik_constraint_weight);
comak_ik.set_ik_accuracy(ik_accuracy);
comak_ik.set_use_visualizer(use_visualizer);
# comak_ik.set_verbose(10);

ik_task_set = comak_ik.get_IKTaskSet()

ik_task_set = osim.IKTaskSet()
ik_task = osim.IKMarkerTask()

for marker, weight in marker_weights.items():
    ik_task.setName(marker)
    ik_task.setWeight(weight)
    ik_task_set.cloneAndAppend(ik_task)

comak_ik.set_IKTaskSet(ik_task_set)
comak_ik.printToXML(save_xml_path)

print('Running COMAKInverseKinematicsTool...')
comak_ik.run();