clear all;
close all;

n = 10; % Number of simulations to be run

system('rm -r TS*'); % Deletes all the simulations called TO+Number present in the folder

for i=1:n
  system(sprintf('cp -r testOctave TS%d',i)); % It copies the baseCase folder "n" times 
endfor

Lcorr = [5e-1 5e-1 25e-2; % Correlation lengths table
         5e-1 5e-1 25e-2;
         5e-1 5e-1 25e-2;
         5e-1 5e-1 25e-2;
         5e-1 5e-1 25e-2;
         5e-1 5e-1 25e-2;
         5e-1 5e-1 25e-2;
         5e-1 5e-1 25e-2;
         5e-1 5e-1 25e-2;
         5e-1 5e-1 25e-2]; 
save corrLengths.mat Lcorr;

Kz = [1 1 1 1 1 1 1 1 1 1]; % Vertical permeability ratios
save verticalKratios.mat Kz;
     
for i=1:n
  cd(sprintf('TS%d',i));
  system(sprintf('octave oneSimRun.m')); % It runs the Octave single simulation scripts inside every folder
  cd(sprintf('..'));
endfor
