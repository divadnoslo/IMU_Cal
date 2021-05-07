clear;
clc
close all;

load('VN200_accel_gyro_8_hr.mat');

[avar_x,tau] = allanvar(accel(:,1),'octave',Fs);
[avar_y,tau] = allanvar(accel(:,2),'octave',Fs);
[avar_z,tau] = allanvar(accel(:,3),'octave',Fs);

figure
loglog(tau/3600,sqrt(avar_x), 'r', tau/3600,sqrt(avar_y), 'g', tau/3600,sqrt(avar_z), 'b')
title('Accel Allan Deviation')
xlabel('\tau (hr)')
ylabel('\sigma m/s^2')
grid on


[avar_x,tau] = allanvar(gyro(:,1),'octave',Fs);
[avar_y,tau] = allanvar(gyro(:,2),'octave',Fs);
[avar_z,tau] = allanvar(gyro(:,3),'octave',Fs);

figure
loglog(tau/3600,sqrt(avar_x), 'r', tau/3600,sqrt(avar_y), 'g', tau/3600,sqrt(avar_z), 'b')
title('Gyro Allan Deviation')
xlabel('\tau (hr)')
ylabel('\sigma m/s^2')
grid on