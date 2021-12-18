clc, clear all, close all;
%% Load the log file and extract vectors which match certain strings
Log = fileread('log');
time = str2double(regexp(Log, '(?<=Time =[^0-9])\-?[0-9]\.?[0-9]+\e?\-?\+?[0-9]+','match'));
conc = str2double(regexp(Log, '(?<=Flux out =[^0-9])\-?[0-9]\.?[0-9]+\e?\-?\+?[0-9]+','match'));
%% Keep the significant data and rescale
c = conc(conc>1e-12)/max(conc); % Keep significant concentrations and rescale
Tadv = 2/4e-6; % Average travel time: domain length / advective velocity 
t = time(c>1e-12)/Tadv; % Re-scaling time axis on the average travel time for simple advection
%% Plot the curve with the name of its simulation
path = pwd; % Path of the folder that FINISHES WITH THE SIMULATION NUMBER
figure(1)
plot(t, c, 'LineWidth', 3, 'DisplayName', sprintf('TO%d',str2double(path(end))))
xlabel('t/t* [-]'); ylabel('c/c_{max} [-]');
legend
figure(2)
plot(t, 1-c, 'LineWidth', 3, 'DisplayName', sprintf('TO%d',str2double(path(end))))
xlabel('t/t* [-]'); ylabel('1-c/c_{max} [-]');
legend
figure(3)
semilogy(t, 1-c, 'LineWidth', 3, 'DisplayName', sprintf('TO%d',str2double(path(end))))
xlabel('t/t* [-]'); ylabel('c/c_{max} [-]');
legend
figure(4)
loglog(t, 1-c, 'LineWidth', 3, 'DisplayName', sprintf('TO%d',str2double(path(end))))
xlabel('t/t* [-]'); ylabel('1-c/c_{max} [-]');
legend