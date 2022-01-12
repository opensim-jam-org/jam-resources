import h5py
import os
import numpy as np
import matplotlib.pyplot as plt

def get_h5_output(filepath, outcome):
    with h5py.File(filepath, "r") as f:
        data = list(f[outcome])   
    return data

def get_h5_type(filepath, outcome):
    with h5py.File(filepath, "r") as f:
        if isinstance(f[outcome], h5py.Dataset):
            return "Dataset"
        elif isinstance(f[outcome], h5py.Group):
            return "Group"

def get_h5_groups_datasets(filepath, path, outcomes):
    groups = []
    datasets = []
    for outcome in outcomes:
        if get_h5_type(filepath, path + outcome) == 'Group':
            groups.append(outcome)
        elif get_h5_type(filepath, path + outcome) == 'Dataset':
            datasets.append(outcome)
    return groups, datasets




class JamAnalysis:
    def __init__(
        self,
    ):
        self.h5_file_list = []
        self.num_files = 0
        self.base_name = 'model'
        self.names = []
        self.num_missing_files = 0
        self.missing_files = []

        self.forceset_name = 'forceset'
        self.coordset_name = 'coordinateset'
        self.frametransformsset_name = 'frametransformsset'

        self.SmithArticularContactName = 'Smith2018ArticularContactForce'

        self.forceset = {}
        self.coordinateset = {}
        self.frametransformsset = {}
        self.comak = {}
    
    def process_forceset(self, h5_filepath, h5_file_idx):
        # Get the components & iterate over them 
        components = get_h5_output(h5_filepath, f'/{self.base_name}/{self.forceset_name}')
        for component_idx, component in enumerate(components):
            if component not in self.forceset:
                self.forceset[component] = {}
            # Get the forcesets & iterate over them
            forcesets = get_h5_output(h5_filepath, f'/{self.base_name}/{self.forceset_name}/{component}')
            for forceset_idx, forceset in enumerate(forcesets):
                if forceset not in self.forceset[component]:
                    self.forceset[component][forceset] = {} 
                if component == self.SmithArticularContactName:
                    # Get mesh names & iterate over them
                    mesh_names = get_h5_output(h5_filepath, f'/{self.base_name}/{self.forceset_name}/{component}/{forceset}')
                    for mesh_idx, mesh_name in enumerate(mesh_names):
                        if mesh_name not in self.forceset[component][forceset]:
                            self.forceset[component][forceset][mesh_name] = {x:{} for x in range(6)}
                        # get params from mesh, break into "groups" and "datasets" & iterate over them as appropriate
                        params = get_h5_output(h5_filepath, f'/{self.base_name}/{self.forceset_name}/{component}/{forceset}/{mesh_name}')
                        groups, datasets = get_h5_groups_datasets(
                            h5_filepath, 
                            f'/{self.base_name}/{self.forceset_name}/{component}/{forceset}/{mesh_name}/',
                            params
                        )
                        
                        for group_idx, group in enumerate(groups):
                            # Iterate over "groups" & store results data in dict
                            regions = get_h5_output(h5_filepath, f'/{self.base_name}/{self.forceset_name}/{component}/{forceset}/{mesh_name}/{group}')
                            for region_idx, region in enumerate(regions):
                                data = np.asarray(get_h5_output(h5_filepath, f'/{self.base_name}/{self.forceset_name}/{component}/{forceset}/{mesh_name}/{group}/{region}'))
                                n_rows, n_cols = data.shape
                                if region_idx not in self.forceset[component][forceset][mesh_name]:
                                    self.forceset[component][forceset][mesh_name][region_idx] = {}
                                if group not in self.forceset[component][forceset][mesh_name][region_idx]:
                                    self.forceset[component][forceset][mesh_name][region_idx][group] = np.zeros((n_rows, n_cols, self.num_files))
                                self.forceset[component][forceset][mesh_name][region_idx][group][:, :, h5_file_idx] = data                                                                     
                        
                        for dataset_idx, dataset in enumerate(datasets):
                            data = np.asarray(get_h5_output(h5_filepath, f'/{self.base_name}/{self.forceset_name}/{component}/{forceset}/{mesh_name}/{dataset}'))
                            if 'regional' in dataset:
                                for region_idx in range(6):
                                    if dataset not in self.forceset[component][forceset][mesh_name][region_idx]:
                                        n_rows = data.shape[0]
                                        self.forceset[component][forceset][mesh_name][region_idx][dataset] = np.zeros((n_rows, self.num_files))
                                    self.forceset[component][forceset][mesh_name][region_idx][dataset][:, h5_file_idx] = data[:, region_idx]
                            else:
                                if len(data.shape) == 1:
                                    if dataset not in self.forceset[component][forceset][mesh_name]:
                                        self.forceset[component][forceset][mesh_name][dataset] = np.zeros((data.shape[0], self.num_files))
                                    self.forceset[component][forceset][mesh_name][dataset][:, h5_file_idx] = data
                                else:
                                    if dataset not in self.forceset[component][forceset][mesh_name]:
                                        self.forceset[component][forceset][mesh_name][dataset] = np.zeros((data.shape[0], data.shape[1], self.num_files))
                                    self.forceset[component][forceset][mesh_name][dataset][:, :, h5_file_idx] = data
                
                else:
                    params = get_h5_output(h5_filepath, f'/{self.base_name}/{self.forceset_name}/{component}/{forceset}')
                    groups, datasets = get_h5_groups_datasets(
                        h5_filepath, 
                        f'/{self.base_name}/{self.forceset_name}/{component}/{forceset}/',
                        params
                    )

                    for dataset_idx, dataset in enumerate(datasets):
                        data = np.asarray(get_h5_output(h5_filepath, f'/{self.base_name}/{self.forceset_name}/{component}/{forceset}/{dataset}'))
                        if dataset not in self.forceset[component][forceset]:
                            self.forceset[component][forceset][dataset] = np.zeros((self.num_time_steps, self.num_files))
                        self.forceset[component][forceset][dataset][:, h5_file_idx] = data
    
    def process_coordinateset(self, h5_filepath, h5_file_idx):
        coordinates = get_h5_output(h5_filepath, f'/{self.base_name}/{self.coordset_name}')
        for coord_idx, coord in enumerate(coordinates):
            if coord not in self.coordinateset:
                self.coordinateset[coord] = {}
            params = get_h5_output(h5_filepath, f'/{self.base_name}/{self.coordset_name}/{coord}')
            groups, datasets = get_h5_groups_datasets(
                h5_filepath, 
                f'/{self.base_name}/{self.coordset_name}/{coord}/',
                params
            )
            # Iteratve over the datasets
            for dataset_idx, dataset in enumerate(datasets):
                data = np.asarray(get_h5_output(h5_filepath, f'/{self.base_name}/{self.coordset_name}/{coord}/{dataset}'))
                if h5_file_idx == 0:
                    self.coordinateset[coord][dataset] = np.zeros((self.num_time_steps, self.num_files))
                self.coordinateset[coord][dataset][:, h5_file_idx] = data
    
    def process_frametransformsset(self, h5_filepath, h5_file_idx):
        ##### JAN.11.2022 - AAGATTI
        ##### HAVENT TESTED THE BELOW YET - THE EXAMPLE .h5 I USED DID NOT HAVE THIS DATATYPE
        print('PROCESSING FRAMETRANSFORMSET => THIS HAS NOT BEEN TESTED BEFORE, PLEASE VERIFY IF IT WORKS')
        frames = get_h5_output(h5_filepath, f'/{self.base_name}/{self.frametransformsset_name}')
        for frame_idx, frame in enumerate(frames):
            outcomes = get_h5_output(h5_filepath, f'/{self.base_name}/{self.frametransformsset_name}/{frame}')
            for outcome_idx, outcome in enumerate(outcomes):
                if 'coordinates' in outcome:
                    transform_type = 'coordinates'
                elif 'transformation_matrix' in outcome:
                    transform_type = 'transformation_matrix'

                params = get_h5_output(h5_filepath, f'/{self.base_name}/{self.frametransformsset_name}/{frame}/{outcome}')
                groups, datasets = get_h5_groups_datasets(
                    h5_filepath, 
                    f'/{self.base_name}/{self.frametransformsset_name}/{frame}/{outcome}/',
                    params
                )
                for dataset_idx, dataset in enumerate(datasets):
                    data = np.asarray(get_h5_output(h5_filepath, f'/{self.base_name}/{self.frametransformsset_name}/{frame}/{outcome}/dataset'))
                    if dataset not in self.frametransformsset[transform_type]:
                        self.frametransformsset[transform_type] = {
                            dataset: np.zeros((self.num_time_steps, self.num_files))
                        }
                    self.frametransformsset[transform_type][dataset][:, h5_file_idx] = data



    def jam_analysis(
        self, 
        h5_file_list,
        base_name=None,
        names=None          # REMOVE? => WHAT IS THE "NAMES" variable used for? 
    ):
        if base_name is not None:
            self.base_name = base_name

        if type(h5_file_list) not in (list, tuple):
            raise Exception(f'`h5_file_list` is type: {type(h5_file_list)} and should be type `list` or `tuple`')
        else:
            self.h5_file_list = h5_file_list
        
        self.num_files = len(h5_file_list)

        # I dont know that below is needed - the `names` or `obj.names` variables are not used
        # in the matlab script
        if names is not None:
            self.names = names
        else:
            for idx in range(self.num_files):
                self.names.append(str(idx))
        
        # Read the h5 file: 
        for h5_file_idx, h5_filepath in enumerate(self.h5_file_list):
            # Test to make sure file exists
            if os.path.exists(h5_filepath) is False:
                print(f'File does not exist: {h5_filepath}')
                self.num_missing_files += 1
                self.missing_files.append(
                    {'idx': h5_file_idx,
                     'path': h5_filepath 
                    }
                )
                continue

            # Test to make sure that it is an .h5 file being passed 
            path = os.path.dirname(h5_filepath)
            filename = os.path.basename(h5_filepath)
            name, ext = os.path.splitext(filename)

            if ext != '.h5':
                raise Exception(f'File: {filename} is not `.h5` format!')

            # If its the first file, then get the time information. 
            if h5_file_idx == 0:
                self.time = get_h5_output(h5_filepath, '/time')
                self.num_time_steps = len(self.time)
            
            # Get the data groups in the h5 file
            h5_groups = get_h5_output(h5_filepath, f'/{self.base_name}')

            for group_idx, group in enumerate(h5_groups):
                # Forceset
                if group == self.forceset_name:
                    self.process_forceset(h5_filepath, h5_file_idx)
                elif group == self.coordset_name:
                    self.process_coordinateset(h5_filepath, h5_file_idx)
                elif group == self.frametransformsset_name:
                    self.process_coordinateset(h5_filepath, h5_file_idx)
            
            # See if COMAK data exists. If it does, add its contents
            base_h5_groups = get_h5_output(h5_filepath, '/')
            if 'comak' in base_h5_groups:
                params = get_h5_output(h5_filepath, '/comak')
                groups, datasets = get_h5_groups_datasets(
                    h5_filepath, 
                    f'/comak/',
                    params
                )
                for dataset_idx, dataset in enumerate(datasets):
                    data = np.asarray(get_h5_output(h5_filepath, f'/comak/{dataset}'))
                    if dataset not in self.comak:
                        self.comak[dataset] = np.zeros((data.shape[0], self.num_files))
                    self.comak[dataset][:, h5_file_idx] = data
    
    def plot_muscle_output(self, muscle_name, param_name, fontsize=20, linewidth=2, ax=None, label=None):
        for file_idx in range(self.num_files):
            if ax is None:
                plt.plot(
                    self.forceset['Muscle'][muscle_name][param_name][:,file_idx], 
                    linewidth=linewidth,
                    label=label
                )
                # plt.ylabel(param_name)
            else:
                ax.plot(
                    self.forceset['Muscle'][muscle_name][param_name][:,file_idx], 
                    linewidth=linewidth,
                    label=label
                )
                # ax.set_ylabel(param_name)
                    

