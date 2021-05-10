function [ret, N] = Send_VN200_cmd(s, cmd)
% Description:
%   A function to issue commands to the VN200
%
% Inputs:
%   s       The serial port object
%   cmd     The command as a char array
%
%  Output:
%   ret     The confirmed response returned from the VN200
%
% Reference: VN-200 User Manual
%   https://www.vectornav.com/support/documentation
%
% Author: S. Bruder

N = length(cmd);                        % Length of the cmd string

%   Compute the 8-bit checksum (i.e., XOR) of the command bytes
checksum = uint8(cmd(1));               % Convert to type unsigned 8-bit integer
for i = 2:length(cmd)
    checksum = bitxor(checksum, uint8(cmd(i)), 'uint8');
end
checksum = dec2hex(checksum, 2);        % Convert to type ASCII - Must have 2-bytes

flush(s);                               % Flush beffers in the serial port
writeline(s, ['$',cmd,'*',checksum]);   % Issue the cmd with checksum
ret = char(readline(s));                % Read "echo" of command
    
end     % End function "send_VN200_cmd"
