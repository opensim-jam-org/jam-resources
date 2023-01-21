import os
import json

import opensim as osim

models_path = os.path.abspath('./model')

model_file = os.path.join(models_path, 'scaled_lenhart_model_updated_wrap.osim')
results_basename = "cycling"

path_ik_settings = os.path.join(models_path, 'P017', 'ik_settings.xml')
markerset_file = os.path.join(models_path, 'P017', 'marker_data.trc')

secondary_constraint_function_file = './results/comak-inverse-kinematics/secondary_coordinate_constraint_functions.xml'
constrained_model_file = './results/comak-inverse-kinematics/ik_constrained_model.osim'

ik_result_dir = './results/comak-inverse-kinematics'

save_xml_path = './inputs/comak_inverse_kinematics_settings.xml'

os.makedirs('./inputs', exist_ok=True)
os.makedirs('./results/graphics', exist_ok=True)


with open('./secondary_coordinates.json', 'r') as f:
    secondary_coordinates = json.load(f)
    


# Settings
start_time = 120.222
start_pad = 0.
stop_time = 122 #122.9

perform_secondary_constraint_sim = True
secondary_constraint_sim_settle_threshold = 1e-4
secondary_constraint_sim_sweep_time = 3.0
secondary_coupled_coordinate_start_value = 0
secondary_coupled_coordinate_stop_value = 110
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

for idx, (coord, dict_) in enumerate(secondary_coordinates.items()):
    comak_ik.set_secondary_coordinates(int(idx), dict_['coordinate'])


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

comak_ik.set_output_motion_file('cycling_ik.mot');
comak_ik.set_time_range(0, start_time-start_pad);
comak_ik.set_time_range(1, stop_time);
comak_ik.set_report_errors(report_errors);
comak_ik.set_report_marker_locations(report_marker_locations);
comak_ik.set_ik_constraint_weight(ik_constraint_weight);
comak_ik.set_ik_accuracy(ik_accuracy);
comak_ik.set_use_visualizer(use_visualizer);


ik_tool = osim.InverseKinematicsTool(path_ik_settings)
ik_task_set = ik_tool.get_IKTaskSet()

idx = 0
while idx < ik_task_set.getSize():
    name  = ik_task_set.get(idx).getName()
    if ('hjc' in name) or ('kjc' in name) or ('ajc' in name) or ('floor' in name) or ('mid' in name):
        ik_task_set.remove(idx)
    else: 
        idx += 1

comak_ik.set_IKTaskSet(ik_task_set)


comak_ik.printToXML(save_xml_path)

print('Running COMAKInverseKinematicsTool...')
comak_ik.run();