function G = TwoPointMaster(corrtype,periodicity,masked,memory,varargin)
% Master function for calculating two point statistics. A by product of
% ongoing computational materials science research at MINED@Gatech.
%
% Author: Ahmet Cecen
% Author Web: ahmetcecen.github.io
% Group Web: mined.gatech.edu
%
% This function changes functionality according to the following formula:
%
% corrtype    - Auto('a') or Cross Correlation('c')
% periodicity - Periodic('p') or Non-Periodic('n')
% masked      - Masked('m') or Raw('r')
% memory      - Memory Mapped('m') or Fast('f')
%
% For example:
%
% G = TwoPointMaster('c','n','m','m',...) would calculate the Masked Non
% Periodic Cross Correlation using the memory mapping.
% 
% After the functionality strings, an appropriate number of variables need
% to be entered in the following order:
%
% H1   -  (Mandatory) First Microstructure.
% H2   -  (If Cross Correlation) Second Microstructure.
% Mask -  (If Masked) Mask Matrix same size as H1.
%
% Examples:
%
% G = TwoPointMaster('c','n','m','m',H1,H2,Mask) - > This is a masked cross
% correlation.
%
% G = TwoPointMaster('a','n','m','m',H1,Mask) - > This is a masked auto
% correlation.
%
% G = TwoPointMaster('a','n','r','m',H1) - > This is a raw auto
% correlation.

% Case Division
sw=strcat(corrtype,periodicity,masked,memory);

switch sw
    case 'aprf'
        H1 = varargin{1};
        G = fftshift(ifftn(fftn(H1).*conj(fftn(H1))))./numel(H1);
    case 'aprm'
        H1 = varargin{1};
        save('TwoPointData.mat','H1');
        G = memacp2pt('TwoPointData.mat');
    case 'apmf'
        H1 = varargin{1};
        Mask = varargin{2};
        G = fftshift(ifftn(fftn(H1.*Mask).*conj(fftn(H1.*Mask))))...
            ./fftshift(ifftn(fftn(Mask).*conj(fftn(Mask))));
    case 'apmm'
        H1 = varargin{1};
        Mask = varargin{2};
        save('TwoPointData.mat','H1','Mask');
        G = memacp2ptmask('TwoPointData.mat');       
    case 'cprf'
        H1 = varargin{1};
        H2 = varargin{2};
        G = fftshift(ifftn(fftn(H1).*conj(fftn(H2))))./numel(H1);
    case 'cprm'
        H1 = varargin{1};
        H2 = varargin{2};
        save('TwoPointData.mat','H1','H2');
        G = memccp2pt('TwoPointData.mat');
    case 'cpmf'
        H1 = varargin{1};
        H2 = varargin{2};
        Mask = varargin{3};
        G = fftshift(ifftn(fftn(H1.*Mask).*conj(fftn(H2.*Mask))))...
            ./fftshift(ifftn(fftn(Mask).*conj(fftn(Mask))));
    case 'cpmm'
        H1 = varargin{1};
        H2 = varargin{2};
        Mask = varargin{3};
        save('TwoPointData.mat','H1','H2','Mask');
        G = memccp2ptmask('TwoPointData.mat');       
    case 'anrf'
        H1 = varargin{1};
        Base = padarray(ones(size(H1)),size(H1)-1,0,'post');
        H1 = padarray(H1,size(H1)-1,0,'post');
        G = fftshift(ifftn(fftn(H1).*conj(fftn(H1))))...
            ./fftshift(ifftn(fftn(Base).*conj(fftn(Base))));
    case 'anrm'
        H1 = varargin{1};
        save('TwoPointData.mat','H1');
        G = memacnp2pt('TwoPointData.mat');               
    case 'anmf'
        H1 = varargin{1};
        Mask = varargin{2};
        Base = padarray(Mask,size(H1)-1,0,'post');
        H1 = padarray(H1,size(H1)-1,0,'post');
        G = fftshift(ifftn(fftn(H1.*Base).*conj(fftn(H1.*Base))))...
            ./fftshift(ifftn(fftn(Base).*conj(fftn(Base))));    
    case 'anmm'
        H1 = varargin{1};
        Mask = varargin{2};
        save('TwoPointData.mat','H1','Mask');
        G = memacnp2ptmask('TwoPointData.mat');       
    case 'cnrf'
        H1 = varargin{1};
        H2 = varargin{2};
        Base = padarray(ones(size(H1)),size(H1)-1,0,'post');
        H1 = padarray(H1,size(H1)-1,0,'post');
        H2 = padarray(H2,size(H2)-1,0,'post');
        G = fftshift(ifftn(fftn(H1).*conj(fftn(H2))))...
            ./fftshift(ifftn(fftn(Base).*conj(fftn(Base))));
    case 'cnrm'
        H1 = varargin{1};
        H2 = varargin{2};
        save('TwoPointData.mat','H1','H2');
        G = memccnp2pt('TwoPointData.mat');       
    case 'cnmf'
        H1 = varargin{1};
        H2 = varargin{2};
        Mask = varargin{3};
        Base = padarray(Mask,size(H1)-1,0,'post');
        H1 = padarray(H1,size(H1)-1,0,'post');
        H2 = padarray(H2,size(H2)-1,0,'post');
        G = fftshift(ifftn(fftn(H1.*Base).*conj(fftn(H2.*Base))))...
            ./fftshift(ifftn(fftn(Base).*conj(fftn(Base))));    
    case 'cnmm'
        H1 = varargin{1};
        H2 = varargin{2};
        Mask = varargin{3};
        save('TwoPointData.mat','H1','H2','Mask');
        G = memccnp2ptmask('TwoPointData.mat');  
end
        
