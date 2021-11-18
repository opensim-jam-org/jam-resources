%jam_analysis
%==========================================================================
%Author: Colin Smith
%--------------------------------------------------------------------------
%
%
%
%==========================================================================


classdef jam_analysis   
    properties
        h5_file_list
        num_files
        names = {}
        forceset
        coordinateset
        frametransformsset
%         model_name
        time
        num_time_steps
        comak
        num_missing_files = 0
        missing_files
        %states
    end
    methods
        function obj = jam_analysis(h5_file_list,varargin)
            
            if(~iscell(h5_file_list))
                h5_file_list = {h5_file_list};
            end
                        
            obj.h5_file_list = h5_file_list;
            obj.num_files = length(h5_file_list);
            num_vargin = nargin -1;
            for k = 1:num_vargin
                if(contains(varargin{k},'names'))
                    obj.names = varargin{k+1};
                    
%                 elseif (contains(varargin{k},'model_files'))
%                     if (iscell(model_name))
%                         obj.model_name = model_name;
%                     else
%                         obj.model_name = cell(obj.num_files,1);
%                         obj.model_name(:) = {model_name};
%                     end
                end
            end
            
            if(length(obj.names) == 0)
                names = cell(obj.num_files,1);
                for i = 1:obj.num_files
                    names{i} = int2str(i);
                end
            end                                      
            
            % Read h5 file
            for n = 1:length(h5_file_list)
                h5_file = h5_file_list{n};
                
                forceset_path = ['/model/forceset'];
                coordset_path = ['/model/coordinateset'];
                frametransformsset_path = ['/model/frametransformsset'];
                %Read h5 file
                if (exist(h5_file,'file') ~=2)
                    fprintf('File does not exist:\n')
                    fprintf('%s\n',h5_file);
                    obj.num_missing_files = obj.num_missing_files+1;
                    obj.missing_files(obj.num_missing_files) = n;
                    continue;                
                end

                [pathstr,name,ext] = fileparts(h5_file);
                if (strcmp(ext,'.h5') ~= 1)
                    fprintf('File: %s',h5_file);
                    error('File is not .h5 format!\n')
                end
                
                if n==1
                    obj.time = h5read(h5_file,'/time');
                    obj.num_time_steps = length(obj.time);
                end
            
                h5_info = h5info(h5_file,['/model']);
                h5_grp = {h5_info.Groups.Name};
                n_h5_grp = length(h5_grp);

                for g = 1:n_h5_grp
                    
                    % ForceSet

                    if strcmp(h5_grp{g},forceset_path)
                        info = h5info(h5_file,forceset_path);

                        nComponents = length(info.Groups);

                        for i=1:nComponents
                            label = strsplit(info.Groups(i).Name,'/');
                            comp = label{4};   
                            
                            nForce = length(info.Groups(i).Groups);

                            for j=1:nForce
                                label = strsplit(info.Groups(i).Groups(j).Name,'/');
                                force = label{5}; 
                            
                                if (contains(comp,'Smith2018ArticularContactForce'))
                                    mesh1_label = strsplit(info.Groups(i).Groups(j).Groups(1).Name,'/');
                                    mesh2_label = strsplit(info.Groups(i).Groups(j).Groups(2).Name,'/');
                                    mesh_names{1} = mesh1_label{6};
                                    mesh_names{2} = mesh2_label{6};

                                    for m=1:2
                                        mesh = mesh_names{m}; 
                                        
                                        %Regional Vec3
                                        nRegVec3 = length(info.Groups(i).Groups(j).Groups(m).Groups);
                                        
                                        for v = 1:nRegVec3
                                            param_label = strsplit(info.Groups(i).Groups(j).Groups(m).Groups(v).Name,'/');
                                            param = param_label{7}; 
                                            for r = 1:6
                                                region = info.Groups(i).Groups(j).Groups(m).Groups(v).Datasets(r).Name;
                                                data_set_name = [ info.Groups(i).Groups(j).Groups(m).Groups(v).Name '/' region];
                                                data = h5read(h5_file,data_set_name)';
                                                [nRow, nCol] = size(data);
                                                if(n==1)
                                                    obj.forceset.(comp).(force).(mesh).region(r).(param) = nan(nRow,nCol,obj.num_files);%h5read(h5_file,data_set_name)';
                                                end
                                                obj.forceset.(comp).(force).(mesh).region(r).(param)(:,:,n) = data;
                                                
                                            end
                                        end
                                        

                                        nDataSet = length(info.Groups(i).Groups(j).Groups(m).Datasets);

                                        for k=1:nDataSet
                                            param = info.Groups(i).Groups(j).Groups(m).Datasets(k).Name;
                                            data_set_name = [ info.Groups(i).Groups(j).Groups(m).Name '/' param];                                            
                                            
                                            if(contains(param,'regional'))
                                                for r = 1:6
%                                                     region = info.Groups(i).Groups(j).Groups(m).Groups(v).Datasets(r).Name;
%                                                     data_set_name = [ info.Groups(i).Groups(j).Groups(m).Groups(v).Name '/' region];
                                                    data = h5read(h5_file,data_set_name);
                                                    if(n==1)
                                                        [nCol, nRow] = size(data(r,:));
                                                        obj.forceset.(comp).(force).(mesh).region(r).(param) = nan(nRow,obj.num_files);
                                                    end
                                                    obj.forceset.(comp).(force).(mesh).region(r).(param)(:,n) = data(r,:);
                                                end
                                            else
                                                [nCol, nRow] = size(h5read(h5_file,data_set_name));
                                                if nRow ==1
                                                    if(n==1)                                                    
                                                        obj.forceset.(comp).(force).(mesh).(param) = nan(nCol,obj.num_files);
                                                    end
                                                    obj.forceset.(comp).(force).(mesh).(param)(:,n) = h5read(h5_file,data_set_name)';
                                                else
                                                    if(n==1)                                                    
                                                        obj.forceset.(comp).(force).(mesh).(param) = nan(nRow,nCol,obj.num_files);
                                                    end
                                                    obj.forceset.(comp).(force).(mesh).(param)(:,:,n) = h5read(h5_file,data_set_name)';
                                                end
                                            end
                                        end
                                        
                                    end
                                else % all other forces
                                    nDataSet = length(info.Groups(i).Groups(j).Datasets);

                                    for k=1:nDataSet
                                        param = info.Groups(i).Groups(j).Datasets(k).Name;
                                        data_set_name = [ info.Groups(i).Groups(j).Name '/' param];

                                        if(n==1)
                                            obj.forceset.(comp).(force).(param) = nan(obj.num_time_steps,obj.num_files);
                                        end
                                        obj.forceset.(comp).(force).(param)(:,n) = h5read(h5_file,data_set_name);
                                    end                                    
                                end
                            end
                        end

                    % Coordinate Set                    
                    elseif strcmp(h5_grp{g}, coordset_path)
                        
                        info = h5info(h5_file,coordset_path);

                        nCoords = length(info.Groups);

                        for i=1:nCoords
                            label = strsplit(info.Groups(i).Name,'/');
                            coord = label{4};                


                            nDataSet = length(info.Groups(i).Datasets);

                            for k=1:nDataSet
                                param = info.Groups(i).Datasets(k).Name;
                                data_set_name = [ info.Groups(i).Name '/' param];
                                
                                if(n==1)
                                    obj.coordinateset.(coord).(param) = nan(obj.num_time_steps,obj.num_files);%h5read(h5_file,data_set_name);
                                end
                                
                                obj.coordinateset.(coord).(param)(:,n) = h5read(h5_file,data_set_name);

                            end
                        end

                    % Joint Mechanics Frame Transforms 
                    elseif strcmp(h5_grp{g}, frametransformsset_path)
                        info = h5info(h5_file,frametransformsset_path);
                        
                        nFrameTransforms = length(info.Groups);
                        
                        for i=1:nFrameTransforms
                            %frame_transform_info = h5info(h5_file,[frametransformsset_path '/' info.Groups(i).Name]);
                            frame_transform_info = h5info(h5_file, info.Groups(i).Name); 
                            %frame_transform = length(info.Groups);
                            
                            for j = 1:length(frame_transform_info.Groups)
                                % Coordinates
                                if(contains(frame_transform_info.Groups(i).Name, 'coordinates'))
                                    nDataSet = length(frame_transform_info.Groups(i).Datasets);

                                    for k=1:nDataSet
                                        coord = frame_transform_info.Groups(i).Datasets(k).Name;
                                        data_set_name = [ frame_transform_info.Groups(i).Name '/' coord];

                                        if(n==1)
                                            obj.frametransformsset.coordinates.(coord) = nan(obj.num_time_steps,obj.num_files);%h5read(h5_file,data_set_name);
                                        end

                                        obj.frametransformsset.coordinates.(coord)(:,n) = h5read(h5_file,data_set_name);

                                    end
                                end
                                
                                % Tranformation Matrix 
                                if(contains(frame_transform_info.Groups(i).Name, 'transformation_matrix'))
                                    nDataSet = length(frame_transform_info.Groups(i).Datasets);

                                    for k=1:nDataSet
                                        coord = frame_transform_info.Groups(i).Datasets(k).Name;
                                        data_set_name = [ frame_transform_info.Groups(i).Name '/' coord];

                                        if(n==1)
                                            obj.frametransformsset.transformation_matrix.(coord) = nan(obj.num_time_steps,obj.num_files);%h5read(h5_file,data_set_name);
                                        end

                                        obj.frametransformsset.transformation_matrix.(coord)(:,n) = h5read(h5_file,data_set_name);

                                    end
                                end
                            end
                            
                        end
                    end
                    
                    % COMAK Convergence
                    info = h5info(h5_file);
                    
                    if any(contains({info.Groups.Name}, '/comak'))                        
                        id = find(contains({info.Groups.Name}, '/comak'));
                        nDataSet = length(info.Groups(id).Datasets);
                        for k=1:nDataSet
                            param = info.Groups(id).Datasets(k).Name;
                            data_set_name = [ info.Groups(id).Name '/' param];

                            if(n==1)
                                nSteps = length(h5read(h5_file,data_set_name));
                                obj.comak.(param) = nan(nSteps,obj.num_files);
                            end

                            obj.comak.(param)(:,n) = h5read(h5_file,data_set_name);

                        end
                        
                    end
                end
            end
        end
        
        function plot_muscle_output(obj,h,muscle_name,param_name)
        %==================================================================
        % PLOT_MUSCLE_OUTPUT
        %------------------------------------------------------------------
        %
        %
        %==================================================================      
       

            hold on;
            for i = 1:obj.num_files                   
                plot(h,obj.forceset.Muscle.(muscle_name).(param_name)(:,i),'LineWidth',2);                
            end
            %xlabel('time [s]')
            ylabel(param_name)
            
            %legend(muscle_names)
        end
        

    end   
end