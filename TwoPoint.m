function GG = TwoPoint(corrtype,cutoff,periodicity,varargin)
% Calculates the two point statistics using the straightforward convolution
% method using FFTs. Use it if you have small datasets or an abundance of
% memory. 
%
% GG = TwoPoint('auto',cutoff,periodicity,H1) will calculate the
% autocorrelation of H1.
%
% GG = TwoPoint('auto',cutoff,periodicity,H1,M1) will calculate the
% autocorrelation of H1 with a mask M1.
%
% GG = TwoPoint('cross',cutoff,periodicity,H1,H2) will calculate the
% crosscorrelation of H1 with H2.
%
% GG = TwoPoint('cross',cutoff,periodicity,H1,H2,M1) will calculate the
% crosscorrelation of H1 with H2, with a uniform mask M1.
%
% GG = TwoPoint('cross',cutoff,periodicity,H1,H2,M1,M2) will calculate the
% crosscorrelation of H1 with mask M1 and H2 with mask M2.


memtype = 'full';

switch periodicity
    
    case 'periodic'
        
        switch corrtype
            
            case 'auto'
                
                if nargin == 4
                    
                    % Periodic 2pt Auto X
                    H1 = varargin{1};
                    GG = CorrMaster(memtype,'auto',cutoff,double(H1));
                    GG = GG./numel(H1);
                    
                elseif nargin == 5
                    
                    % Periodic 2pt Auto Masked 
                    H1 = varargin{1};
                    M1 = varargin{2};
                    GG = CorrMaster(memtype,corrtype,cutoff,double(H1.*M1));
                    BB = CorrMaster(memtype,'auto',cutoff,double(M1));
                    GG = GG./BB;
                    
                end
                
            case 'cross'
                
                if nargin == 5
                    
                    % Periodic 2pt Cross X
                    H1 = varargin{1};
                    H2 = varargin{2};
                    GG = CorrMaster(memtype,corrtype,cutoff,double(H1),double(H2));
                    GG = GG./numel(H1);
                    
                elseif nargin == 6
                    
                    % Periodic 2pt Cross Masked
                    H1 = varargin{1};
                    H2 = varargin{2};
                    M1 = varargin{3};
                    GG = CorrMaster(memtype,corrtype,cutoff,double(H1.*M1),double(H2.*M1));
                    BB = CorrMaster(memtype,'auto',cutoff,double(M1));
                    GG = GG./BB;
                    
                elseif nargin == 7
                    
                    % Periodic 2pt Cross Bi-Masked
                    H1 = varargin{1};
                    H2 = varargin{2};
                    M1 = varargin{3};
                    M2 = varargin{4};
                    GG = CorrMaster(memtype,corrtype,cutoff,double(H1.*M1),double(H2.*M2));
                    BB = CorrMaster(memtype,corrtype,cutoff,double(M1),double(M2));
                    GG = GG./BB;
                    
                end
                
        end
        
    case 'nonperiodic'
        
        switch corrtype
            
            case 'auto'
                
                if nargin == 4
                    
                    % Non-Periodic 2pt Auto X
                    H1 = varargin{1};
                    GG = CorrMaster(memtype,corrtype,cutoff,padarray(double(H1),repmat(cutoff,[1 ndims(H1)]),0,'post'));
                    BB = CorrMaster(memtype,'auto',cutoff,padarray(ones(size(H1)),repmat(cutoff,[1 ndims(H1)]),0,'post'));
                    GG = GG./BB;
                    
                elseif nargin == 5
                    
                    % Non-Periodic 2pt Cross Masked
                    H1 = varargin{1};
                    M1 = varargin{2};
                    GG = CorrMaster(memtype,corrtype,cutoff,padarray(double(H1.*M1),repmat(cutoff,[1 ndims(H1)]),0,'post'));
                    BB = CorrMaster(memtype,'auto',cutoff,padarray(double(M1),repmat(cutoff,[1 ndims(M1)]),0,'post'));
                    GG = GG./BB;
                    
                end
                
            case 'cross'
                
                if nargin == 5
                    
                    % Non-Periodic 2pt Cross X
                    H1 = varargin{1};
                    H2 = varargin{2};                    
                    GG = CorrMaster(memtype,corrtype,cutoff,padarray(double(H1),repmat(cutoff,[1 ndims(H1)]),0,'post'),padarray(double(H2),repmat(cutoff,[1 ndims(H2)]),0,'post'));
                    BB = CorrMaster(memtype,'auto',cutoff,padarray(ones(size(H1)),repmat(cutoff,[1 ndims(H1)]),0,'post'));
                    GG = GG./BB;
                    
                elseif nargin == 6
                    
                    % Non-Periodic 2pt Cross Masked
                    H1 = varargin{1};
                    H2 = varargin{2};
                    M1 = varargin{3};
                    GG = CorrMaster(memtype,corrtype,cutoff,padarray(double(H1.*M1),repmat(cutoff,[1 ndims(H1)]),0,'post'),padarray(double(H2.*M1),repmat(cutoff,[1 ndims(H2)]),0,'post'));
                    BB = CorrMaster(memtype,'auto',cutoff,padarray(double(M1),repmat(cutoff,[1 ndims(M1)]),0,'post'));
                    GG = GG./BB;
                    
                elseif nargin == 7
                    
                    % Non-Periodic 2pt Cross Bi-Masked
                    H1 = varargin{1};
                    H2 = varargin{2};
                    M1 = varargin{3};
                    M2 = varargin{4};
                    GG = CorrMaster(memtype,corrtype,cutoff,padarray(double(H1.*M1),repmat(cutoff,[1 ndims(H1)]),0,'post'),padarray(double(H2.*M2),repmat(cutoff,[1 ndims(H2)]),0,'post'));
                    BB = CorrMaster(memtype,corrtype,cutoff,padarray(double(M1),repmat(cutoff,[1 ndims(M1)]),0,'post'),padarray(double(M2),repmat(cutoff,[1 ndims(M2)]),0,'post'));
                    GG = GG./BB;
                    
                end
                
        end
        
end

end





