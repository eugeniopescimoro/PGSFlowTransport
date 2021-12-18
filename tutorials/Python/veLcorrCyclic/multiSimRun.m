clear all;
close all;

n = 30; % Number of simulations to be run

system('rm -r TS*'); % Deletes all the simulations called TO+Number present in the folder

for i=1:n
  system(sprintf('cp -r testOctave TS%d',i)); % It copies the baseCase folder "n" times 
endfor

Lcorr = [0.1 0.05 0.05; % Correlation lengths table
         0.2 0.05 0.05;
         0.3 0.05 0.05;
         0.4 0.05 0.05;
         0.5 0.05 0.05;
         0.6 0.05 0.05;
         0.7 0.05 0.05;
         0.8 0.05 0.05;
         0.9 0.05 0.05;
         1.0 0.05 0.05;
         1.0 0.05 0.05;
         1.0 0.10 0.05;
         1.0 0.15 0.05;
         1.0 0.20 0.05;
         1.0 0.25 0.05;
         1.0 0.30 0.05;
         1.0 0.35 0.05;
         1.0 0.40 0.05;
         1.0 0.45 0.05;
         1.0 0.50 0.05;
         1.0 0.50 0.05;
         1.0 0.50 0.10;
         1.0 0.50 0.15;
         1.0 0.50 0.20;
         1.0 0.50 0.25;
         1.0 0.50 0.30;
         1.0 0.50 0.35;
         1.0 0.50 0.40;
         1.0 0.50 0.45;
         1.0 0.50 0.50];
save corrLengths.mat Lcorr;

Kz = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; % Vertical permeability ratios
save verticalKratios.mat Kz;
     
for i=1:n
  cd(sprintf('TS%d',i));
  system(sprintf('octave oneSimRun.m')); % It runs the Octave single simulation scripts inside every folder
  cd(sprintf('..'));
endfor
