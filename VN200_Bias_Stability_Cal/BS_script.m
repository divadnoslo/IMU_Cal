% Bias Stability Script

% Getting Ready
clear;                          % Clear all variables from the workspace
close all;                      % Close all windows
clc;                            % "Clean" the command window
Date_Time = clock;              % Obtain a date/time stamp to name files
addpath(genpath(pwd));          % Add all subfolders the the path (current session)

%% Data Collection Settings

% Number of Runs
num_runs = 20;

% VN200 Settings
Fs = 50;                        % Set Sample frequency <= 100 (Hz)
dT = 1/Fs;                      % Sample interval (sec)
nSec = 10 * 60;                 % Duration of cdata collection (sec)
N = Fs*nSec;                    % Number of samples to collect (dimless)
%fprintf('Collecting data for %2i sec at %d Hz\n\n', nSec, Fs);

% Initialize arrays
for ii = 1 : num_runs
    gyro{ii}    = zeros(3, N);          % gyro data:           3 floats per sample
    accel{ii}   = zeros(3, N);          % accel data:          3 floats per sample
    temp{ii}    = zeros(3,1);           % temperature data:    1 float per sample
end

%% Begin Data Collection

% For Each Run
for ii = 1 : num_runs
    
    [s, SN] = Initialize_VN200_IMU(Fs); % Initialize the VN200
    
    fprintf(['Begin Run #', num2str(ii), '\n'])
    
    % Warm Up the IMU - only record tempurature
    for k = 1:N                        
        
        % Run IMU for warming only, Record Temperature Only
        [~, ~, ~, ~, ~] = Read_VN200_IMU(s); % Get VN200 IMU data
        
        % Print Status
        if (k == 1)
            [~, ~, ~, temp{ii}(1), ~] = Read_VN200_IMU(s); % Get VN200 IMU data
            fprintf(['Begining Warm-Up Cycle: temp = ', num2str(temp{ii}(1)), 'C \n'])
            fprintf(['In American Units: temp = ', num2str((temp{ii}(1)*9/5)+32), 'F \n'])
        end
        
    end
    
    % Record IMU Data at operating tempurature
    for k = 1:N                        
        
        % Print Status
        if (k == 1)
            [~, accel{ii}(:,k), gyro{ii}(:,k), temp{ii}(2), ~] = Read_VN200_IMU(s); % Get VN200 IMU data
            fprintf(['Begining Warm-Up Cycle: temp = ', num2str(temp{ii}(2)), 'C \n'])
            fprintf(['In American Units: temp = ', num2str((temp{ii}(2)*9/5)+32), 'F \n'])
        else
            [~, accel{ii}(:,k), gyro{ii}(:,k), ~, ~] = Read_VN200_IMU(s); % Get VN200 IMU data
        end
        
    end
    
    % Cool down the VN200, no recording (only if not last cycle)
    if (ii ~= num_runs)
        
        % Print Status
        [~, ~, ~, temp{ii}(3), ~] = Read_VN200_IMU(s); % Get VN200 IMU data
        fprintf(['Begining Cool Down Cycle: temp = ', num2str(temp{ii}(3)), 'C \n'])
        fprintf(['In American Units: temp = ', num2str((temp{ii}(3)*9/5)+32), 'F \n\n'])
        
        pause(nSec)
        
    else
        
        % Print Status
        [~, ~, ~, temp{ii}(3), ~] = Read_VN200_IMU(s); % Get VN200 IMU data
        fprintf(['Begining Cool Down Cycle: temp = ', num2str(temp{ii}(3)), 'C \n'])
        fprintf(['In American Units: temp = ', num2str((temp{ii}(3)*9/5)+32), 'F \n\n'])
        
    end
    
    clear s;                        % Clear the serial port object <=> Close serial port
    
end

%% Process Data

% Load in Fixed Biases
load('IMU_Cal_Const_error_Sources.mat')

% Process Data
for ii = 1 : num_runs
    
    b_a_BS_runs(:,ii) = mean(accel{ii}')' - [0; 0; -9.80772] - b_a_FB;
    b_g_BS_runs(:,ii) = mean(gyro{ii}')' - b_g_FB;
    
end

% Average all runs for the average bias stability
b_a_BS = mean(b_a_BS_runs')';
b_g_BS = mean(b_g_BS_runs')';

% Print Results
row_names = {'x-accel (m/s^2)', 'y-accel (m/s^2)', 'z-accel (m/s^2)', ...
             'x-gyro  (deg/s)', 'y-gyro  (deg/s)', 'z-gyro  (deg/s)'};
column_name = {'IMU Bias Stability Values'};
T = array2table([b_a_BS; b_g_BS * 180/pi], 'RowNames', row_names, 'VariableNames', column_name);
disp(T)

%% Plot Data

kk = 0 : 3 : (3*num_runs - 1);

figure
hold on
for ii = 1 : num_runs
    plot(kk(ii)/3, temp{ii}(1), 'r*')
    plot((kk(ii)+1)/3, temp{ii}(2), 'b*')
    plot((kk(ii)+2)/3, temp{ii}(3), 'k*')
end
hold off
title('Warming/Cooling Cycles')
xlabel('Cycle')
ylabel('Temperature (C)')
grid on

figure
subplot(3, 1, 1)
hold on
plot([1:num_runs], b_a_BS_runs(1,:), 'r*')
title('x accel BS values')
grid on
hold off
subplot(3, 1, 2)
hold on
plot([1:num_runs], b_a_BS_runs(2,:), 'g*')
title('y accel BS values')
grid on
hold off
subplot(3, 1, 3)
hold on
plot([1:num_runs], b_a_BS_runs(3,:), 'b*')
title('z accel BS values')
grid on
hold off

figure
subplot(3, 1, 1)
hold on
plot([1:num_runs], b_g_BS_runs(1,:), 'r*')
title('x gyro BS values')
grid on
hold off
subplot(3, 1, 2)
hold on
plot([1:num_runs], b_g_BS_runs(2,:), 'g*')
title('y gyro BS values')
grid on
hold off
subplot(3, 1, 3)
hold on
plot([1:num_runs], b_g_BS_runs(3,:), 'b*')
title('z gyro BS values')
grid on
hold off

%% Save Data

save('IMU_BS.mat', 'b_a_BS', 'b_g_BS', 'SN', 'Date_Time')
fprintf('Data has been saved. Ensure to move data file before running next experiment. \n\n')





