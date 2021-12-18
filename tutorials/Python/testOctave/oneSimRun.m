clear all;
close all;

load('./../corrLengths.mat'); % Load the correlation lengths table
load('./../verticalKratios.mat'); % Load the vertical permeabilities table
% Lcorr is a matrix that needs to be parsed by 'Allproperties' one row at time
% Kz is a vector which provides the vertical/horizontal permeability ratio
path = pwd; % Path of the folder that FINISHES WITH THE SIMULATION NUMBER
key = '/TS';
index = strfind(path, key);
simNum = str2num(path(index+length(key):end));
system(sprintf('./Allproperties %f %f %f %f', Lcorr(simNum, :), Kz(simNum)));
% system('./Allrun');
