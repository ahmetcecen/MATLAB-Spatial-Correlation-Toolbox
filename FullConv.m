function G = FullConv(corrtype,periodicity,masked,varargin)
% Full convolution wrapper for two point statistics. A by product of
% ongoing computational materials science research at MINED@Gatech.
%
% Author: Ahmet Cecen
% Author Web: ahmetcecen.github.io
% Group Web: mined.gatech.edu
%
% This function changes functionality according to the following formula:
%
% corrtype         - Auto('a') or Cross Correlation('c')
% periodicity      - Periodic('p') or Non-Periodic('n')
% masked          - Masked('m') or Raw('r')
%
% For example:
%
% G = TwoPointMaster('c','n','m',...) would calculate the Masked Non
% Periodic Cross Correlation using the memory mapping.
% 
% After the functionality strings, an appropriate number of variables need
% to be entered in the following order:
%
% H1        -  (Mandatory) First Microstructure.
% H2        -  (If Cross Correlation) Second Microstructure.
% Mask     -  (If Masked) Mask Matrix same size as H1.
%
% Examples:
%
% G = TwoPointMaster('c','n','m',H1,H2,Mask) - > This is a masked cross
% correlation.
%
% G = TwoPointMaster('a','n','m',H1,Mask) - > This is a masked auto
% correlation.
%
% G = TwoPointMaster('a','n','r',H1) - > This is a raw auto
% correlation.

% Case Division
sw=strcat(corrtype,periodicity,masked);

switch sw

    case 'apr'
        H1 = varargin{1};
        G = fftshift(ifftn(fftn(H1).*conj(fftn(H1))));
		
    case 'apm'
        H1 = varargin{1};
        Mask = varargin{2};
        G = fftshift(ifftn(fftn(H1.*Mask).*conj(fftn(H1.*Mask))));

    case 'cpr'
        H1 = varargin{1};
        H2 = varargin{2};
        G = fftshift(ifftn(fftn(H1).*conj(fftn(H2))));

    case 'cpm'
        H1 = varargin{1};
        H2 = varargin{2};
        Mask = varargin{3};
        G = fftshift(ifftn(fftn(H1.*Mask).*conj(fftn(H2.*Mask))));
  
    case 'anr'
        H1 = varargin{1};
        H1 = padarray(H1,size(H1)-1,0,'post');
        G = fftshift(ifftn(fftn(H1).*conj(fftn(H1))));
        
    case 'anm'
        H1 = varargin{1};
        Mask = varargin{2};
        Base = padarray(Mask,size(H1)-1,0,'post');
        H1 = padarray(H1,size(H1)-1,0,'post');
        G = fftshift(ifftn(fftn(H1.*Base).*conj(fftn(H1.*Base))));    
    
    case 'cnr'
        H1 = varargin{1};
        H2 = varargin{2};
        H1 = padarray(H1,size(H1)-1,0,'post');
        H2 = padarray(H2,size(H2)-1,0,'post');
        G = fftshift(ifftn(fftn(H1).*conj(fftn(H2))));
   
    case 'cnm'
        H1 = varargin{1};
        H2 = varargin{2};
        Mask = varargin{3};
        Base = padarray(Mask,size(H1)-1,0,'post');
        H1 = padarray(H1,size(H1)-1,0,'post');
        H2 = padarray(H2,size(H2)-1,0,'post');
        G = fftshift(ifftn(fftn(H1.*Base).*conj(fftn(H2.*Base))));    

end
        
