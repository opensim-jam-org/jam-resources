
First need to create a model: 
- get_scale_factors_for_each_body.ipynb
- scale_markerset.ipynb
- scale_lenhart_to_subject.ipynb

Upate the model's wrap surfaces: 
- update_the_model_wrapping.ipynb

Run simulation pipelines
- comak_ik_only.py
- comaktool_only.py
- comak_id_only.py
- jntmech_only.py

Then can view all of the results
- plot_comak_simulation_results.ipynb

Finally, can load the meshes from `./results/joint-mechanics` 
into paraview to visualize cartilate pressures. 