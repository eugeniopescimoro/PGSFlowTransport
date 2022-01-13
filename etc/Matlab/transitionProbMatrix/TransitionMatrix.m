% Internship Eugenio Pescimoro
%
% Compute transiograms and transition matrices from the most probable
% facies map of the Herten data case.
%
% INPUTS:
%   datafilename         name of the data file
%   icode                set 7 for conditional transition probabilities or
%                        5 for bivariate probabilities
%
% Dimitri D'Or - Ephesia Consult - 08/09/2020

function [Z,transitionmatrix]=TransitionMatrix(datafilename, icode)

%% Loading the data set

mat = load(datafilename);
%Z=mat;
Z=reshape(mat,66,50,50);

origin=[0 0 0];
step=[1 1 1];
nb_steps=size(Z);
c=creategrid(origin,step,nb_steps);

%% Putting the facies in the right order (as in the Herten description document Bayer2015)

ncat=length(unique(Z(:))); %Number of facies

Z=Z+1;
Z(Z==(ncat+1))=1;

%% Computing the proportions
   
zc=Z(:);
prop_temp=tabulateFSS(zc);
prop=prop_temp(:,3)/100


%% Computing the omnidirectional transiograms

[Xmesh,Ymesh,Zmesh,gh,nh]=varioFFT3D(c,zc,icode,1,1);

%% Extracting omnidirectional transition matrix (as the mean of the transitions in all directions)

s=1;
display=0;
[Xsub,Ysub,Zsub,ghe]=modelextracting(Xmesh,Ymesh,Zmesh,gh,s,display);

transitionmatrix=NaN*ones(ncat,ncat);
for i=1:ncat
    for j=1:ncat
        d=0;
        if(i==j)
            if(icode==7)
                d=1;
            else
                d=prop(i);
            end
        end
        transitionmatrix(i,j)=(sum(ghe{i,j}(:))-d)/(length(ghe{i,j}(:))-1);
    end
end

transitionmatrix

