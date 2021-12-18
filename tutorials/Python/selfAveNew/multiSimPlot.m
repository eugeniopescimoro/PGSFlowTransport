clc, clear all, close all;

for i = 1:3
    cd(sprintf('TK%d',i));
    %% Load the log file and extract vectors which match certain strings (UNSTABLE)
    % Log = fileread('log');
    % time = str2double(regexp(Log, '(?<=Time =[^0-9])\-?[0-9]\.?[0-9]+\e?\-?\+?[0-9]+','match'));
    % conc = str2double(regexp(Log, '(?<=Flux out =[^0-9])\-?[0-9]\.?[0-9]+\e?\-?\+?[0-9]+','match'));
    %% Using the "grep" bash command and redirecting output to matlab variable
    system('cat log | grep "Total mass =" > mass');
    massID = fopen('mass');
    massFormat = 'Total mass = %f';
    mass = cell2mat(textscan(massID, massFormat));
    fclose(massID);
    system('cat log | grep "Mean vel =" > mvel');
    mvelID = fopen('mvel');
    mvelFormat = 'Mean vel = (%f %f %f)';
    mvel = cell2mat(textscan(mvelID, mvelFormat));
    fclose(mvelID);
    system('cat log | grep "Flux out =" > conc');
    concID = fopen('conc');
    concFormat = 'Flux out = %f';
    conc = cell2mat(textscan(concID, concFormat));
    fclose(concID);
    system('cat log | grep "Time =" > time');
    timeID = fopen('time');    
    timeFormat = 'Time = %f';    
    time = cell2mat(textscan(timeID, timeFormat));
    fclose(timeID);
    %% Keep the significant data and rescale
    c = conc(conc>1e-12); % Keep significant concentrations
    Tadv = 2/mvel(1, 1); % Average travel time: domain length / mean velocity
    t = time(conc>1e-12)/Tadv; % Re-scaling time axis on the average travel time for simple advection
    %% Plot the curve with the name of its simulation
    path = pwd; % Path of the folder that FINISHES WITH THE SIMULATION NUMBER
    CL = ["Kv/Kh=1", "Kv/Kh=0.1", "Kv/Kh=0.01"];
    subplot(2, 2, 1)
    plot(t, c, 'LineWidth', 3)
    xlabel('t/t* [-]'); ylabel('c/c_{max} [-]');
    hold on
    subplot(2, 2, 2)
    plot(t, 1-c, 'LineWidth', 3)
    xlabel('t/t* [-]'); ylabel('1-c/c_{max} [-]');
    hold on
    subplot(2, 2, 3)
    semilogy(t, 1-c, 'LineWidth', 3)
    xlabel('t/t* [-]'); ylabel('log(1-c/c_{max}) [-]');
    hold on
    subplot(2, 2, 4)
    loglog(t, 1-c, 'LineWidth', 3, 'DisplayName', sprintf('TK%d: %s Pé=%.1f',str2double(path(end)), CL{i}, 2*mvel(1,1)/1e-9))
    xlabel('log(t/t*) [-]'); ylabel('log(1-c/c_{max}) [-]');
    legend
    hold on
    cd(sprintf('..'));
end
