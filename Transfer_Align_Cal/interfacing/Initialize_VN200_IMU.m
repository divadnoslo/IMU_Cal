function [s, SN] = Initialize_VN200_IMU( Fs )
% Function Description:
%   Initializes the VN200 IMU and configures the serial port to run from
%   the VN200's internal clock (async mode) at a freq of Fs Hz.
%
% INPUTS:
%   Fs = VN200 sample rate (Hz)  Must be <= 100
%   
% OUTPUTS:
%   s = Serial port object
%   SN = Device serial number
%
% NOTES:
%   The default BaudRate of the VN200 is 115200.
%
% Reference: VN-200 User Manual
%   https://www.vectornav.com/support/documentation
%
% Author: S. Bruder

if Fs > 100     % Limit the sample rate to <= 100 Hz
   error('The sample freq must be <= 100 Hz!!'); 
end

% Define a few VN200 commands
VN200_cmds.Diasble_Async = 'VNWRG,06,0';        % Disable asynchronous data output
VN200_cmds.Get_SN        = 'VNRRG,03';          % Request the IMU's serial number
VN200_cmds.Set_Fs        = ['VNWRG,07,',num2str(Fs)]; % Set the asynchronous data output freq to Fs Hz
VN200_cmds.Enable_IMU    = 'VNWRG,06,19';       % Enable async IMU Measurements on VN200
 
% Configure serial port
BaudRate = 115200;                              % BaudRate: 115200, 128000, 230400, 460800, or 921600
ComPorts = serialportlist("available");         % Determine available serial ports

s = serialport(ComPorts(end), BaudRate);        % Create a serial port object
configureTerminator(s,"CR/LF");                 % Set both the read and write terminators to "CR/LF"

%% Configure the VN200 to output IMU data at Fs Hz
Send_VN200_cmd(s, VN200_cmds.Diasble_Async);    % Disable asynchronous data output
pause(0.1);                                     % A "small" delay
[ret, N] = Send_VN200_cmd(s, VN200_cmds.Diasble_Async);    % Disable asynchronous data output
pause(0.1);                                     % A "small" delay
if ~strcmp(ret(2:N+1), VN200_cmds.Diasble_Async)
   error('Failed to obtain Disable Async data output. Please try again.'); 
end

[ret, N] = Send_VN200_cmd(s, VN200_cmds.Get_SN);    % Request the IMU's serial number 
if ~strcmp(ret(2:N+1), VN200_cmds.Get_SN)
   error('Failed to obtain valid SN. Please try again.'); 
end
SN = ret(11:end-3);                                 % Extract the Serial Number

[ret, N] = Send_VN200_cmd(s, VN200_cmds.Set_Fs);    % Set the asynchronous data output freq to Fs Hz
if ~strcmp(ret(2:N+1), VN200_cmds.Set_Fs)
   error('Failed to set the sample Freq to Fs. Please try again.'); 
end

[ret, N] = Send_VN200_cmd(s, VN200_cmds.Enable_IMU);% Enable async IMU Measurements on VN200
if ~strcmp(ret(2:N+1), VN200_cmds.Enable_IMU)
   error('Failed to enable async IMU Measurements. Please try again.'); 
end
 
end     % End of function "Initialize_VN200_IMU"

