% EE 440 Modern Nav
% Test code for the VN200 IMU
%  - Plots Gyro, accel, Barometer, compass data, and IMU temperature
%
% Author: S. Bruder

clear;                          % Clear all variables from the workspace
close all;                      % Close all windows
clc;                            % "Clean" the command window
Date_Time = clock;              % Obtain a date/time stamp to name files
addpath(genpath(pwd));          % Add all subfolders the the path (current session)

Fs = 50;                        % Set Sample frequency <= 100 (Hz)
dT = 1/Fs;                      % Sample interval (sec)
nSec = 30;                      % Duration of cdata collection (sec)
N = Fs*nSec+1;                  % Number of samples to collect (dimless)
fprintf('Collecting data for %2i sec at %d Hz\n\n', nSec, Fs);

% Initialize arrays
compass = zeros(3, N);          % compass data:        3 floats per sample
gyro    = zeros(3, N);          % gyro data:           3 floats per sample
accel   = zeros(3, N);          % accel data:          3 floats per sample
temp    = zeros(1, N);          % temperature data:    1 float per sample
baro    = zeros(1, N);          % Baro data:           1 float per sample

%% Collect VN200 measurements
[s, SN] = Initialize_VN200_IMU(Fs); % Initialize the VN200
 for k = 1:N                        % Retrieve IMU data from VN200
    [compass(:,k), accel(:,k), gyro(:,k), temp(k), baro(k)] = Read_VN200_IMU(s); % Get VN200 IMU data
    if ~mod(N-k+1,Fs)
        fprintf('Please wait %i more secomds!!!\n', round((N-k+1) / Fs));
    end
 end
clear s;                        % Clear the serial port object <=> Close serial port

t = 0:dT:(N-1)*dT;              % Time vector (sec)

% Plot the compass/magnetometer data (in Gauss)
figure,
plot(t, compass(1,:), 'r', t, compass(2,:), 'g', t, compass(3,:), 'b');
title('Plot of VN200 Compass Data', 'FontSize', 12);
xlabel('Time (sec)');
ylabel('Magnetic Field (Gauss)')
legend('m_x', 'm_y', 'm_z')
grid

% Plot the accel data (in m/s^2)
figure,
plot(t, accel(1,:), 'r', t, accel(2,:), 'g', t, accel(3,:), 'b');
title('Plot of VN200 Accel Data', 'FontSize', 12);
xlabel('Time (sec)');
ylabel('Accel (m/s^2)')
legend('a_x', 'a_y', 'a_z')
grid

% Plot the Gyro data (in deg/s)
figure,
plot(t, gyro(1,:)*180/pi, 'r', t, gyro(2,:)*180/pi, 'g', t, gyro(3,:)*180/pi, 'b');
title('Plot of VN200 Gyro Data', 'FontSize', 12);
xlabel('Time (sec)');
ylabel('Angular Rate (°/s)')
legend('\omega_x', '\omega_y', '\omega_z')
grid

% Plot the Temperature data (in deg C)
figure,
plot(t, temp, 'k');
title('Plot of VN200 Temperature Data', 'FontSize', 12);
xlabel('Time (sec)');
ylabel('Temperature (°C)')
grid

% Plot the Brometric pressure data (in kPa)
figure,
plot(t, baro, 'k');
title('Plot of VN200 Barometric Pressure Data', 'FontSize', 12);
xlabel('Time (sec)');
ylabel('Pressure (kPa)')
grid

%% Save the data to a time-stamped *.mat file
% Save_file_name = sprintf('VN200_data_%d_%02.0f_%02.0f_%02.0f_%02.0f_%02.0f',Date_Time(1), Date_Time(2),Date_Time(3),...
%     Date_Time(4),Date_Time(5),round(Date_Time(6)));

in1 = input("Do you want to save this data? (y/n):  ", 's');

if (in1(1) == 'y' || in1 == 'Y')
    
    Save_file_name = input("File Name:  ", 's');
    
    UNITS.compass = 'Gauss';
    UNITS.accel = 'm/s^2';
    UNITS.gyro = 'rad/s';
    UNITS.temp = '°C';
    UNITS.baro = 'kPa';
    UNITS.Fs = 'Hz';
    
    save([Save_file_name, '.mat'], 'compass', 'accel', 'gyro', 'baro', 'temp', 'UNITS', 'Fs', 'SN');
    
end

