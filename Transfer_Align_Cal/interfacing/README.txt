This version of the code uses the new MATLAB serial interface functions:

s = serialport(ComPorts(end), BaudRate);        % Create a serial port object

vs older MATLAB serial port access functions:

s = serial(comPort(end), 'BaudRate', BaudRate); % Create a serial port object