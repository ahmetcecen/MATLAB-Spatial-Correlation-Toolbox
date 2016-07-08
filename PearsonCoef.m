function GG = PearsonCoef(corrtype,cutoff,periodicity,varargin)
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
                    H1 = H1 - mean(H1(:)); % Subtract mean
                    GG = CorrMaster(memtype,'auto',cutoff,double(H1));
                    GG = GG./sum(H1(:).^2);
                    
                elseif nargin == 5
                    
                    % Periodic 2pt Auto Masked 
                    H1 = varargin{1};
                    M1 = varargin{2};
                    H1 = H1 - mean(H1(:).*M1(:)); % Subtract mean
                    GG = CorrMaster(memtype,corrtype,cutoff,double(H1.*M1));
                    sigma1 = CorrMaster(memtype,'cross',cutoff,double(H1.^2.*M1),double(M1));
                    GG = GG./sqrt(sigma1.*sigma1(end:-1:1,end:-1:1,end:-1:1));
                    
                end
                
            case 'cross'
                
                if nargin == 5
                    
                    % Periodic 2pt Cross X
                    H1 = varargin{1};
                    H2 = varargin{2};
                    H1 = H1 - mean(H1(:)); % Subtract mean
                    H2 = H2 - mean(H2(:)); % Subtract mean
                    GG = CorrMaster(memtype,corrtype,cutoff,double(H1),double(H2));
                    GG = GG./sqrt(sum(H1(:).^2)*sum(H2(:).^2));
                    
                elseif nargin == 6
                    
                    % Periodic 2pt Cross Masked
                    H1 = varargin{1};
                    H2 = varargin{2};
                    M1 = varargin{3};
                    H1 = H1 - mean(H1(:).*M1(:)); % Subtract mean
                    H2 = H2 - mean(H2(:).*M1(:)); % Subtract mean
                    GG = CorrMaster(memtype,corrtype,cutoff,double(H1.*M1),double(H2.*M1));
                    sigma1 = CorrMaster(memtype,'cross',cutoff,double(H1.^2.*M1),double(M1));
                    sigma2 = CorrMaster(memtype,'cross',cutoff,double(M1),double(H2.^2.*M1));
                    GG = GG./sqrt(sigma1.*sigma2);
                    
                elseif nargin == 7
                    
                    % Periodic 2pt Cross Bi-Masked
                    H1 = varargin{1};
                    H2 = varargin{2};
                    M1 = varargin{3};
                    M2 = varargin{4};
                    H1 = H1 - mean(H1(:).*M1(:)); % Subtract mean
                    H2 = H2 - mean(H2(:).*M2(:)); % Subtract mean
                    GG = CorrMaster(memtype,corrtype,cutoff,double(H1.*M1),double(H2.*M2));
                    sigma1 = CorrMaster(memtype,'cross',cutoff,double(H1.^2.*M1),double(M1));
                    sigma2 = CorrMaster(memtype,'cross',cutoff,double(M2),double(H2.^2.*M2));
                    GG = GG./sqrt(sigma1.*sigma2);
                    
                end
                
        end
        
    case 'nonperiodic'
        
        switch corrtype
            
            case 'auto'
                
                if nargin == 4
                    
                    % Non-Periodic 2pt Auto X
                    H1 = varargin{1};
                    H1 = H1 - mean(H1(:)); % Subtract mean
                    GG = CorrMaster(memtype,corrtype,cutoff,padarray(double(H1),repmat(cutoff,[1 ndims(H1)]),0,'post'));
                    sigma1 = CorrMaster(memtype,'cross',cutoff,padarray(double(H1.^2),repmat(cutoff,[1 ndims(H1)]),0,'post'),...
                                        padarray(ones(size(H1)),repmat(cutoff,[1 ndims(H1)]),0,'post'));
                    GG = GG./sqrt(sigma1.*sigma1(end:-1:1,end:-1:1,end:-1:1));
                    
                elseif nargin == 5
                    
                    % Non-Periodic 2pt Cross Masked
                    H1 = varargin{1};
                    M1 = varargin{2};
                    H1 = H1 - mean(H1(:).*M1(:)); % Subtract mean
                    GG = CorrMaster(memtype,corrtype,cutoff,padarray(double(H1.*M1),repmat(cutoff,[1 ndims(H1)]),0,'post'));
                    sigma1 = CorrMaster(memtype,'cross',cutoff,...
                                        padarray(double(H1.^2.*M1),repmat(cutoff,[1 ndims(H1)]),0,'post'),...
                                        padarray(double(M1),repmat(cutoff,[1 ndims(M1)]),0,'post'));
                    GG = GG./sqrt(sigma1.*sigma1(end:-1:1,end:-1:1,end:-1:1));
                    
                end
                
            case 'cross'
                
                if nargin == 5
                    
                    % Non-Periodic 2pt Cross X
                    H1 = varargin{1};
                    H2 = varargin{2};
                    H1 = H1 - mean(H1(:)); % Subtract mean
                    H2 = H2 - mean(H2(:)); % Subtract mean
                    GG = CorrMaster(memtype,corrtype,cutoff,...
                                    padarray(double(H1),repmat(cutoff,[1 ndims(H1)]),0,'post'),...
                                    padarray(double(H2),repmat(cutoff,[1 ndims(H2)]),0,'post'));
                    sigma1 = CorrMaster(memtype,'cross',cutoff,...
                                        padarray(double(H1.^2),repmat(cutoff,[1 ndims(H1)]),0,'post'),...
                                        padarray(ones(size(H1)),repmat(cutoff,[1 ndims(H1)]),0,'post'));
                    sigma2 = CorrMaster(memtype,'cross',cutoff,...
                                        padarray(ones(size(H1)),repmat(cutoff,[1 ndims(H1)]),0,'post'),...
                                        padarray(double(H2.^2),repmat(cutoff,[1 ndims(H2)]),0,'post'));
                    GG = GG./sqrt(sigma1.*sigma2);
                    
                elseif nargin == 6
                    
                    % Non-Periodic 2pt Cross Masked
                    H1 = varargin{1};
                    H2 = varargin{2};
                    M1 = varargin{3};
                    H1 = H1 - mean(H1(:).*M1(:)); % Subtract mean
                    H2 = H2 - mean(H2(:).*M1(:)); % Subtract mean
                    GG = CorrMaster(memtype,corrtype,cutoff,padarray(double(H1.*M1),repmat(cutoff,[1 ndims(H1)]),0,'post'),padarray(double(H2.*M1),repmat(cutoff,[1 ndims(H2)]),0,'post'));
                    sigma1 = CorrMaster(memtype,'cross',cutoff,...
                                        padarray(double(H1.^2.*M1),repmat(cutoff,[1 ndims(H1)]),0,'post'),...
                                        padarray(double(M1),repmat(cutoff,[1 ndims(M1)]),0,'post'));
                    sigma2 = CorrMaster(memtype,'cross',cutoff,...
                                        padarray(double(M1),repmat(cutoff,[1 ndims(M1)]),0,'post'),...
                                        padarray(double(H2.^2.*M1),repmat(cutoff,[1 ndims(H2)]),0,'post'));
                    GG = GG./sqrt(sigma1.*sigma2);
                    
                elseif nargin == 7
                    
                    % Non-Periodic 2pt Cross Bi-Masked
                    H1 = varargin{1};
                    H2 = varargin{2};
                    M1 = varargin{3};
                    M2 = varargin{4};
                    H1 = H1 - mean(H1(:).*M1(:)); % Subtract mean
                    H2 = H2 - mean(H2(:).*M2(:)); % Subtract mean
                    GG = CorrMaster(memtype,corrtype,cutoff,padarray(double(H1.*M1),repmat(cutoff,[1 ndims(H1)]),0,'post'),padarray(double(H2.*M2),repmat(cutoff,[1 ndims(H2)]),0,'post'));
                    sigma1 = CorrMaster(memtype,'cross',cutoff,...
                                        padarray(double(H1.^2.*M1),repmat(cutoff,[1 ndims(H1)]),0,'post'),...
                                        padarray(double(M1),repmat(cutoff,[1 ndims(M1)]),0,'post'));
                    sigma2 = CorrMaster(memtype,'cross',cutoff,...
                                        padarray(double(M2),repmat(cutoff,[1 ndims(M2)]),0,'post'),...
                                        padarray(double(H2.^2.*M2),repmat(cutoff,[1 ndims(H2)]),0,'post'));
                    GG = GG./sqrt(sigma1.*sigma2);
                    
                end
                
        end
        
end

end





