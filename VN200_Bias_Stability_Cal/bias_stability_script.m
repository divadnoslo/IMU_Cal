% Bias Stability Script

% Getting Ready
clear;                          % Clear all variables from the workspace
close all;                      % Close all windows
clc;                            % "Clean" the command window
Date_Time = clock;              % Obtain a date/time stamp to name files
addpath(genpath(pwd));          % Add all subfolders the the path (current session)

%% Data Collection Settings

% Number of Runs
num_runs = 10;

% VN200 Settings
Fs = 50;                        % Set Sample frequency <= 100 (Hz)
dT = 1/Fs;                      % Sample interval (sec)
nSec = 1;                 % Duration of cdata collection (sec)
N = Fs*nSec;                   % Number of samples to collect (dimless)
%fprintf('Collecting data for %2i sec at %d Hz\n\n', nSec, Fs);

% Initialize arrays
for ii = 1 : num_runs
    gyro{ii}    = zeros(3, N);          % gyro data:           3 floats per sample
    accel{ii}   = zeros(3, N);          % accel data:          3 floats per sample
    temp{ii}    = zeros(3, N);          % temperature data:    1 float per sample
    t{ii}       = zeros(1, 3*N);
end

%% Begin Data Collection



% For Each Run
for ii = 1 : num_runs
    
    [s, SN] = Initialize_VN200_IMU(Fs); % Initialize the VN200
    
    fprintf(['Begin Run #', num2str(ii), '\n'])
    
    % Warm Up the IMU - only record tempurature
    for k = 1:N                        
        
        % Run IMU for warming only, Record Temperature Only
        [~, ~, ~, temp{ii}(k), ~] = Read_VN200_IMU(s); % Get VN200 IMU data
        
        % Print Status
        if (k == 1)
            fprintf(['Begining Warm-Up Cycle: temp = ', num2str(temp{ii}(k)), 'C \n'])
            fprintf(['In American Units: temp = ', num2str((temp{ii}(k)*9/5)+32), 'F \n'])
        end
        
        % Create Time Vector for later plotting functions
        if (ii == 1)
            t{ii}(1:N) = 0 : dT : (nSec - dT);
        else
            t{ii}(1:N) = (t{ii-1}(end)+dT) : dT : (t{ii-1}(end)+nSec);
        end
        
    end
    
    % Record IMU Data at operating tempurature
    for k = 1:N                        
        
        % Record IMU Data
        [~, accel{ii}(:,k), gyro{ii}(:,k), temp{ii}(N+k), ~] = Read_VN200_IMU(s); % Get VN200 IMU data
        
        % Print Status
        if (k == 1)
            fprintf(['Begining Data Collection Cycle: temp = ', num2str(temp{ii}(N+k)), 'C \n'])
            fprintf(['In American Units: temp = ', num2str((temp{ii}(N+k)*9/5)+32), 'F \n'])
        end
        
        % Create Time Vector for later plotting functions
        t{ii}(N+k) = t{ii}(N-1+k) + dT;
        
    end
    
    % Cool down the VN200, no recording (only if not last cycle)
    if (ii ~= num_runs)
        
        % Print Status
        fprintf(['Begining Cool Down Cycle: temp = ', num2str(temp{ii}(N+N)), 'C \n'])
        fprintf(['In American Units: temp = ', num2str((temp{ii}(N+N)*9/5)+32), 'F \n\n'])
        
        for k = 1 : N
            pause(dT) % One Time Step
            t{ii}((2*N)+k) = t{ii}((2*N)-1+k) + dT;
        end
        
    else
        
        % Print Status
        fprintf(['Begining Cool Down Cycle: temp = ', num2str(temp{ii}(N+N)), 'C \n'])
        fprintf(['In American Units: temp = ', num2str((temp{ii}(N+N)*9/5)+32), 'F \n\n'])
        fprintf('Warming/Cooling Cycles are complete.  Data has been saved and plotted. \n\n\n')
        
        for k = 1 : N
            t{ii}((2*N)+k) = t{ii}((2*N)-1+k) + dT;
        end
        
        
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
             'x-gyro (rad/s)', 'y-gyro (rad/s)', 'z-gyro (rad/s)'};
column_name = {'IMU Bias Stability Values'};
T = array2table([b_a_BS; b_g_BS], 'RowNames', row_names, 'VariableNames', column_name);
disp(T)

%% Plot Data

clr_str = ['r', 'b', 'k'];
kk = 1 : N;

figure
hold on
for ii = 1 : num_runs
    plot(t{ii}(kk), temp{ii}(1,:), 'r')
    plot(t{ii}(N+kk), temp{ii}(2,:), 'b')
    plot(t{ii}(N+N+kk), temp{ii}(3,:), 'k')
end
hold off
title('Warming/Cooling Cycles')
xlabel('Time(s)')
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

save_data = input('Do you want to save this data? (y/n): ', 's');

if (save_data(1) == 'y' || save_data(1) == 'Y')
    save('IMU_BS.mat', 'b_a_BS', 'b_g_BS', 'SN', 'Date_Time')
end




