% Internship Eugenio Pescimoro
%
% Compute transiograms and transition matrices from the most probable
% facies map of the Herten data case.
%
% Dimitri D'Or - Ephesia Consult - 08/09/2020
% Eugenio Pescimoro - University of Nottingham - 13/01/2022

clear;clc
close all

% appinfo = matlab.apputil.install('./TransitionMatrix.mlappinstall');
% cd (appinfo.location);
% matlab.apputil.run('TransitionMatrixAPP')

%% Transform OpenFOAM permeability matrix K into categories
x = 66;
y = 50;
z = 50;
skip = 22;
fileID = fopen('K_scat6_HC_proc0', 'r');
formatSpec = '%s %f %f %f %f %s';
K = textscan(fileID, formatSpec, x*y*z, 'HeaderLines', skip);
fclose(fileID);

fileID = fopen('Kcategory','w');
Kcateg = zeros(x*y*z, 1);
for i=1:1:x*y*z
    if K{1, 4}(i) == 1e-09
        Kcateg(i) = 1.000000;
    end
    if K{1, 4}(i) == 1e-11
        Kcateg(i) = 2.000000;
    end
    if K{1, 4}(i) == 1e-13
        Kcateg(i) = 3.000000;
    end
    if K{1, 4}(i) == 1e-15
        Kcateg(i) = 4.000000;
    end
    fprintf(fileID, '%f\n', Kcateg(i));
end
fclose(fileID);

%% Parameters

datafilename='Kcategory';  % Name of the file containing the map
% fileID = fopen('best4PropMapKhomo','r');
% formatSpec = '%d';
% datafilename = cell2mat(textscan(fileID, formatSpec));

%% Computing the transition matrix

[Z,transitionmatrix]=TransitionMatrix(datafilename, 7);

%% Viewing the data set

view3d