function [Xsub,Ysub,Zsub,ghe]=modelextracting(Xmesh,Ymesh,Zmesh,gh,s,display)

% function [Xsub,Ysub,Zsub,ghe]=modelextracting(Xmesh,Ymesh,Zmesh,gh,s,display); 
%
% In case a transiogram or bivariate probability model is computed by FFT 
% on a grid containing NaN values, it is necessary to select only the central
% portion of the maps. Too far from the center, the computed probabilities
% are unaccurate and thus useless. The purpose of this function is to
% extract the central portion of the maps.
%
% INPUT : 
%
% Xmesh    nx by ny matrix of X coordinates
% Ymesh    nx by ny matrix of Y coordinates
% Zmesh    nx by ny matrix of Z coordinates
% gh       nvar by nvar cell array of nx by ny (direct- and cross-) maps for
%                     variables i and j.
% s        scalar     size of the half of the isotropic window. Choose it not
%                     too large in order to avoid considering inaccurate 
%                     values located far from the center.
% display  boolean    tells if a plot must be displayed (set to 1) or not
%                     (set to 0, default)
%                     
% OUTPUT :
%
% Xsub      nx by ny matrix of X coordinates for the extracted portion
% Ysub      nx by ny matrix of Y coordinates for the extracted portion
% Zsub      nx by ny matrix of Z coordinates for the extracted portion
% ghe       nvar by nvar cell array of nx by ny by nz extracted (direct- and cross-) 
%                     maps for variables i and j depending on icode.
%
% Author: Dimitri D'Or (Ephesia Consult)
%
% 2015/01/23

testiscell=1;
if ~iscell(gh),
    testiscell=0;
    ghtemp=gh;
    clear gh
    gh{1}=ghtemp;
    clear ghtemp;
end;

ncat=size(gh,1);

xrange=[-s s];
yrange=[-s s];
zrange=[-s s];

for i=1:length(xrange),
    ix(i)=find(Xmesh==xrange(i));
    iy(i)=find(Ymesh==yrange(i));
    iz(i)=find(Zmesh==zrange(i));
end;

Xsub=Xmesh(ix(1):ix(2));
Ysub=Ymesh(iy(1):iy(2));
Zsub=Zmesh(iz(1):iz(2));

for i=1:ncat,
    for j=1:ncat,
        if ~isempty(gh{i,j}),
            ghe{i,j}=gh{i,j}(ix(1):ix(2),iy(1):iy(2),iz(1):iz(2));
        end;
    end;
end;

if testiscell==0,
    ghetemp=ghe;
    clear ghe;
    ghe=ghetemp{1};
    clear ghetemp;
end;

% Display graph

if display && (length(size(gh{1,1}))==2) % if display and we are in 2D
    figure
    for i=1:ncat,
        for j=1:ncat,
            if isempty(gh{i,j}),
                continue;
            end;
            subplot(ncat,ncat,j+ncat*(i-1));
            pcolor(Xsub,Ysub,ghe{i,j});
            title(['Var. ',num2str(i),' vs. Var. ',num2str(j)]);
            axis equal
            axis([min(Xsub) max(Xsub) min(Ysub) max(Ysub)]);
            shading flat;
            caxis([0 1]);
            colorbar
        end;
    end;
end;

if display && (length(size(gh{1,1}))~=2)
    warning('Plots cannot be drawn if dimension is not 2D');
end;