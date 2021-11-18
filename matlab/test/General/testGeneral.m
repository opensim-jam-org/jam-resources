%% Test General
%==========================================================================
%
%
%
%
%
%==========================================================================
clear; close all;
import org.opensim.modeling.*
Logger.setLevelString('info');

%% Test OpenSim Spline Fits
%==========================================================================
time.step = 0.01;
time.hold1_duration = 0.1;
time.slide_duration = 0.1;
time.hold2_duration = 0.2;

inputs.start_displacement = 0.1;
inputs.max_displacement = 0.12;

% Simulation Time
time.hold1_time = 0 : time.step : time.hold1_duration;
time.slide_time = time.hold1_duration + time.step : time.step : time.hold1_duration + time.slide_duration;
time.hold2_time = time.hold1_duration + time.slide_duration + time.step : time.step : time.hold1_duration + time.slide_duration + time.hold2_duration;
time.values = [time.hold1_time, time.slide_time, time.hold2_time];

time_points = [0,time.hold1_duration,...
    time.hold1_duration + time.slide_duration,...
    time.hold1_duration + time.slide_duration + time.hold2_duration];

time.num_hold1_steps = length(time.hold1_time);
time.num_slide_steps = length(time.slide_time);
time.num_hold2_steps = length(time.hold2_time);
time.num_steps = length(time.values);


% Prescribed Coordinates

inputs.displacement = [inputs.start_displacement,inputs.start_displacement,inputs.max_displacement ,inputs.max_displacement];

inputs.smooth_displacement = interp1(time_points, inputs.displacement, time.values,'pchip');


coord_data.displacement = inputs.smooth_displacement';
coord_data.time = time.values;

% Fit OpenSim Splines

% Plot Spline Fits
simm = SimmSpline();
simm.setName('simm');

gcv = GCVSpline();
gcv.setName('gcv');
gcv.setDegree(1);

for i = 1:time.num_steps
    simm.addPoint(time.values(i),coord_data.displacement(i));
    gcv.addPoint(time.values(i),coord_data.displacement(i));
end


for i = 1:time.num_steps
    simm_displacement(i) = simm.calcValue(Vector(1,time.values(i)));
    gcv_displacement(i) = gcv.calcValue(Vector(1,time.values(i)));
end

%% Plot Input and Curve Fits
figure('name','TEST: OpenSim Spline Fit')
hold on
plot(coord_data.time,coord_data.displacement);
plot(coord_data.time,simm_displacement,'o');
plot(coord_data.time,gcv_displacement,'s');