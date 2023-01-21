import os
import json

import opensim as osim


models_path = os.path.abspath('./model')

# model_file = os.path.join(models_path, 'scaled_lenhart_model.osim')
model_file = os.path.join(models_path, 'scaled_lenhart_model_complete_wrap_surfaces_new_quads_wrap.osim')

external_loads_file = os.path.join(models_path, 'P017', 'external_loads_pedal.xml')
results_basename = "cycling"

comak_result_dir = './results/comak_new_wrap'
id_result_dir = './results/comak-inverse-dynamics'

save_xml_path = './inputs/comak_inverse_dynamics_settings.xml'

# Settings
start_time = 120.57
start_pad = 0.0
stop_time = 121.26



inverse_dynamics = osim.InverseDynamicsTool()

inverse_dynamics = osim.InverseDynamicsTool()
inverse_dynamics.set_results_directory(id_result_dir)
inverse_dynamics.setModelFileName(model_file)
inverse_dynamics.setStartTime(start_time - start_pad)
inverse_dynamics.setEndTime(stop_time)

exclude_frc = osim.ArrayStr()
exclude_frc.append('ALL')

inverse_dynamics.setExcludedForces(exclude_frc)
inverse_dynamics.setExternalLoadsFileName(external_loads_file)

inverse_dynamics.setCoordinatesFileName(os.path.join(comak_result_dir, results_basename + '_values.sto'))

inverse_dynamics.setLowpassCutoffFrequency(6)
inverse_dynamics.setOutputGenForceFileName(results_basename + '_inverse-dynamics.sto')

inverse_dynamics.printToXML(save_xml_path)

inverse_dynamics.run()