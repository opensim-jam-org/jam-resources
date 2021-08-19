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
        nFiles
        forceset
        coordinateset
        model_name
        time
        comak
        %states
    end
    methods
        function obj = jam_analysis(model_name,h5_file_list)
            obj.model_name = model_name;
            
            obj.h5_file_list = h5_file_list;
            obj.nFiles = length(h5_file_list);
            
            forceset_path = ['/' model_name '/forceset'];
            coordset_path = ['/' model_name '/coordinateset'];
            
            % Read h5 file
            for n = 1:length(h5_file_list)
                h5_file = h5_file_list{n};
                
                %Read h5 file
                if (exist(h5_file,'file') ~=2)
                    fprintf('File: %s',h5_file);
                    error('File does not exist!\n')                
                end

                [pathstr,name,ext] = fileparts(h5_file);
                if (strcmp(ext,'.h5') ~= 1)
                    fprintf('File: %s',h5_file);
                    error('File is not .h5 format!\n')
                end
                
                if(n==1)
                    obj.time = h5read(h5_file,'/time');
                end
            
                h5_info = h5info(h5_file,['/' obj.model_name]);
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


                                nDataSet = length(info.Groups(i).Groups(j).Datasets);

                                for k=1:nDataSet
                                    param = info.Groups(i).Groups(j).Datasets(k).Name;
                                    data_set_name = [ info.Groups(i).Groups(j).Name '/' param];

                                    if(n==1)
                                        obj.forceset.(comp).(force).(param)(:,:,obj.nFiles) = h5read(h5_file,data_set_name);
                                    end
                                    obj.forceset.(comp).(force).(param)(:,:,n) = h5read(h5_file,data_set_name);
%                                     a = h5read(h5_file,data_set_name);
                                
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
                                    obj.coordinateset.(coord).(param)(:,obj.nFiles) = h5read(h5_file,data_set_name);
                                end
                                
                                obj.coordinateset.(coord).(param)(:,n) = h5read(h5_file,data_set_name);

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
                                obj.comak.(param)(:,obj.nFiles) = h5read(h5_file,data_set_name);
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
            for i = 1:obj.nFiles                   
                plot(h,obj.forceset.Muscle.(muscle_name).(param_name)(:,i),'LineWidth',2);                
            end
            %xlabel('time [s]')
            ylabel(param_name)
            
            %legend(muscle_names)
        end
        
        function addExternalLoads(obj,file)
            
        end
    end   
end