import os
import json

import opensim as osim


matlab_path = os.path.abspath('../../../../matlab/')
models_path = os.path.abspath('../../../../models')


model_file = os.path.join(models_path, 'knee_healthy', 'lenhart2015', 'lenhart2015.osim')
external_loads_file = os.path.join(models_path, 'knee_healthy', 'experimental_data', 'motion_analysis', 'overground_17_ext_loads.xml')
results_basename = "walking"

comak_result_dir = './results/comak'
id_result_dir = './results/comak-inverse-dynamics'

save_xml_path = './inputs/comak_inverse_dynamics_settings.xml'

start_time = 1.16
start_pad = 0.0
stop_time = 2.36



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