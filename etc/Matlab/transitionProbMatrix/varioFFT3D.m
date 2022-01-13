function [Xmesh,Ymesh,Zmesh,gh,nht,mem]=varioFFT3D(c,z,icode,categ,display,hmax)
%
% function [Xmesh,Ymesh,Zmesh,gh,nht,mem]=varioFFT3D(c,z,icode,categ,display,hmax);
%
% function to compute variograms, cross-variograms, covariograms,
% cross-covariograms and pseudo-cross-variograms in 1D, 2D or 3D for up to 3 variables.
% the data are on a (possibly incomplete) regular grid
% the program computes variograms in the frequency domain using 2D-FFT.
%
% INPUT :
%
% c        n by d     matrix of coordinates (regular grid)
% z        n by nvar  matrix of values for the variables. Each line is
%                     associated with the corresponding line vector of
%                     coordinates in the c matrix, and each column corresponds
%                     to a different variable.
%                     Missing values are indicated by NaN
%                     IMPORTANT: At a point, either all variables are observed or all are missing
% icode               a code to indicate which function to compute
%                      =1 : variograms and cross-variograms;
%                      =2 : covariograms and cross-covariograms
%                      =3 : variograms and pseudo-cross-variograms
%                      =4 : a mean is computed for the whole field instead of according to the lags
%                      =5 : covariance non-centr�e (bivariate probabilities, for categorical data)
%                      =6 : transiograms (for categorical data)
%                      =7 : non-ergodic transiograms
%                      =8 : local marginal proportions
% categ    boolean    tells if the variables are categorical (set to 1) or not
%                     (set to 0, default)
% display  boolean    tells if a plot must be displayed (set to 1) or not
%                     (set to 0, default)%
% hmax     1 x d      optional maximum distances desired along the dimensions
%
% OUTPUT :
%
% Xmesh    nx by ny matrix of X coordinates
% Ymesh    nx by ny matrix of Y coordinates
% Zmesh    nx by ny matrix of Z coordinates
% gh       nvar by nvar cell array of nx by ny (direct- and cross-) maps for
%                     variables i and j depending on icode.
% nht      same structure as gh for continuous data. For categorical data, a single cell
% mem      scalar   used memory volume (kb)
%
% This program uses the functions FFT2, IFFT2, FFT2SHIFT and CONJ which are
% standard MATLAB functions.
%
% Author: Dimitri D'Or - Ephesia Consult - 2014/11/17
% Modified from variof2.m written by D. Marcotte, dmarcotte@mail.polymtl.ca
%
% Reference :
%
% Marcotte D., 1996. Fast Variogram computation with FFT. Computers & Geoscience, 22, 10, 1175-1186.

%space1=memory;

if nargin<5,
    display=0;
end;
if nargin<4,
    categ=0;
end;

[n_nodes,nvar]=size(z);
if categ,
    nvar=length(unique(z(~isnan(z))));
end;

%%% Finding the parameters of the grid

minc=min(c);
maxc=max(c);

for i=1:size(c,2),
    dum1=diff(unique(c(:,i)));
    if isempty(dum1)
        dc(i)=1;
    else
        dc(i)=dum1(1);
    end
end;
nc=((maxc-minc)./dc)+1;

%%% Reformatting the data

if categ,
    idnan=isnan(z);
    zc=zeros(n_nodes,nvar);
    for i=1:nvar,
        zc(z==i,i)=1;
        zc(idnan,i)=nan;       % to keep the nan in the computation
    end;
    z=zc;
    clear zc;
end;

Z=cell(nvar,1);
Zid=cell(nvar,1);

if length(nc)==1
    nc=[nc 1];
end
for i=1:nvar,
    Z{i}=reshape(z(:,i),nc);
end;

[n,p,q]=size(Z{1});          % dimensions of data matrix

if nargin<6
    hmax=[n-1,p-1,q-1];
else
    hmax=ceil(hmax./dc);
    hmax=min([hmax;nc-1]);
    if length(hmax)<3;
        hmax(3)=0;
    end
end

% Initialization of the graphics elements

Xmesh=NaN;
Ymesh=NaN;
Zmesh=NaN;

if size(c,2)>=1,
    Xmesh=unique(c(:,1));
    Xmesh=[-flipud(Xmesh(2:end)) ;Xmesh(1:end)];     % symmetrize the vector around zero
    r=[-hmax(1) hmax(1)]*dc(1);
    for i=1:length(r),
        ix(i)=find(Xmesh==r(i));
    end;
    Xmesh=Xmesh(ix(1):ix(2));
end;
if size(c,2)>=2,
    Ymesh=unique(c(:,2));
    Ymesh=[-flipud(Ymesh(2:end)) ;Ymesh(1:end)];     % symmetrize the vector around zero
    r=[-hmax(2) hmax(2)]*dc(2);
    for i=1:length(r),
        iy(i)=find(Ymesh==r(i));
    end;
    Ymesh=Ymesh(iy(1):iy(2));
end;
if size(c,2)==3,
    Zmesh=unique(c(:,3));
    Zmesh=[-flipud(Zmesh(2:end)) ;Zmesh(1:end)];     % symmetrize the vector around zero
    r=[-hmax(3) hmax(3)]*dc(3);
    for i=1:length(r),
        iz(i)=find(Zmesh==r(i));
    end;
    Zmesh=Zmesh(iz(1):iz(2));
end;
%%% Initialization

if icode==6,
    for i=1:nvar,
        idxnotnan=~isnan(z(:,i));
        prop(i)=sum(z(idxnotnan,i))/sum(idxnotnan);
    end;
end;

clear z

gh=cell(nvar,nvar);

% find the closest multiple of 8 to obtain a good compromise between
% speed (a power of 2) and memory required

nrows=2*n-1;
ncols=2*p-1;
nz=2*q-1;
nr2=ceil(nrows/8)*8;
nc2=ceil(ncols/8)*8;
if nz>8
    nz2=ceil(nz/8)*8;
else
    nz2=nz;
end
nv=[nr2,nc2,nz2];

% form an indicator  matrix:  1's for all data values
%                             0's for missing values
% in data matrix, replace missing values by 0;

for i=1:nvar,
    Zid{i}=~isnan(Z{i});                        % 1 for a data value; 0 for missing
    Z{i}(~Zid{i})=0;                            % missing replaced by 0
end;

% Preparation

if icode==4;
    for i=1:nvar,
        m(i)=sum(sum(Z{i}(Zid{i})))/sum(sum(Zid{i}));
        Z{i}(Zid{i})=Z{i}(Zid{i})-m(i);
    end
end;

% Compute number of pairs

if categ,
    fxid=fftn(Zid{1},nv);
    nh=round(real(ifftn(conj(fxid).*fxid)));
end

if icode==7 || icode==8,
    Zall=Z{1};
    for i=2:nvar
        Zall=Zall+Z{i}(:,:,:);
    end
    fx_all=fftn(Zall,nv);                      % fourier transform of Z{i}
end

for i=1:nvar,
    fxi=fftn(Z{i},nv);                         % fourier transform of Z{i}
    fxid=fftn(Zid{i},nv);
    fxid_ii=fftn(Zid{i}.*Zid{i},nv);
    
    if icode==7 || icode==8,
        propi=round(real(ifftn(conj(fxi).*fx_all)));
    end
    
    if icode==1 || icode==3, % These functions are symmetric
        jdeb=i;
    else
        jdeb=1;
    end
    
    for j=jdeb:nvar,
        
        if icode==8 && i~=j,
            continue;
        end;
        
        fxj=fftn(Z{j},nv);                         % fourier transform of Z{i}
        fxid_ij=fftn(Zid{i}.*Zid{j},nv);
        fxid_jj=fftn(Zid{j}.*Zid{j},nv);
        
        if ~categ
            nh=round(real(ifftn(conj(fxid).*fxid_ij)));            % number of pairs
        end
        
        switch icode
            case 1                    % variograms and cross-variograms are computed
                t1=fftn(Z{i}.*Zid{j},nv);
                t2=fftn(Z{j}.*Zid{i},nv);
                fx2=fftn(Z{i}.*Z{j},nv);
                gh{i,j}=real(ifftn(conj(fxid_ij).*fx2+conj(fx2).*fxid_ij-conj(t1).*t2-t1.*conj(t2)))./max(nh,1)/2;
                
            case 2                    % covariograms and cross-covariograms are computed
                m_tail=real(ifftn(conj(fxi).*fxid_jj))./max(nh,1); % computes the tail means
                m_head=real(ifftn(conj(fxid_ii).*fxj))./max(nh,1); % computes the head means
                gh{i,j}=real((ifftn(conj(fxi).*fxj))./max(nh,1)-m_tail.*m_head);
                
            case 3                         % variograms and pseudo-cross-variograms are computed
                fx2ii=fftn(Z{i}.*Z{i},nv);
                fx2jj=fftn(Z{j}.*Z{j},nv);
                gh{i,j}=real(ifftn(fxid_jj.*conj(fx2ii)+conj(fxid_ii).*fx2jj-2*conj(fxi).*fxj))./max(nh,1)/2;
                
            case {4,5}
                gh{i,j}=real((ifftn(conj(fxi).*fxj))./max(nh,1));
                
            case 6                        % Transiograms are computed
                gh{i,j}=(real((ifftn(conj(fxi).*fxj))./max(nh,1)))/prop(i);
                
            case 7                        % non ergodic transiograms
                gh{i,j}=real((ifftn(conj(fxi).*fxj)))./max(propi,1);
            case 8                        % local marginal proportions
                gh{i,i}=max(propi,1)./nh;
        end
        
        % reduce matrices to required size
        
        t=floor(nv/2)+1;
        if ~categ
            nh=fftshift(nh);
            if nv(3)==1,
                nht{i,j}=nh(t(1)-hmax(1):t(1)+hmax(1),t(2)-hmax(2):t(2)+hmax(2));
            else
                nht{i,j}=nh(t(1)-hmax(1):t(1)+hmax(1),t(2)-hmax(2):t(2)+hmax(2),t(3)-hmax(3):t(3)+hmax(3));
            end;
        end
        gh{i,j}=fftshift(gh{i,j});
        if nv(3)==1,
            gh{i,j}=gh{i,j}(t(1)-hmax(1):t(1)+hmax(1),t(2)-hmax(2):t(2)+hmax(2));
        else
            gh{i,j}=gh{i,j}(t(1)-hmax(1):t(1)+hmax(1),t(2)-hmax(2):t(2)+hmax(2),t(3)-hmax(3):t(3)+hmax(3));
        end;
    end;
end;

if categ
    nh=fftshift(nh);
    if nv(3)==1,
        nht=nh(t(1)-hmax(1):t(1)+hmax(1),t(2)-hmax(2):t(2)+hmax(2));
    else
        nht=nh(t(1)-hmax(1):t(1)+hmax(1),t(2)-hmax(2):t(2)+hmax(2),t(3)-hmax(3):t(3)+hmax(3));
    end;
end

% Display graph

if display && (length(size(gh{1,1}))==2) % if display and we are in 2D
    figure
    for i=1:nvar,
        for j=1:nvar,
            if isempty(gh{i,j}),
                continue;
            end;
            subplot(nvar,nvar,j+nvar*(i-1));
            pcolor(Xmesh,Ymesh,gh{i,j});
%             imagesc(gh{i,j})
            title(['Var. ',num2str(i),' vs. Var. ',num2str(j)]);
            axis equal
            axis([min(Xmesh) max(Xmesh) min(Ymesh) max(Ymesh)]);
            shading flat;
            if icode==6 || icode==7,
                caxis([0 1]);
%             else
%                 caxis([0 gmax]);
            end;
            colorbar
        end;
    end;
end;

if display && (length(size(gh{1,1}))~=2)
    warning('Plots cannot be drawn if dimension is not 2D');
end;

%space2=memory;
%mem=space2.MemUsedMATLAB-space1.MemUsedMATLAB;
