function G = memccp2pt(filename)
% Periodic two point statistics with memory mapping. The filename must 
% link to a mat file with H1, H2 variables. 

% Map the Data
Data = matfile(filename,'Writable',true);

% Find and Save P(H1H2R)
H1 = Data.H1;
H1 = fftn(H1);

H2 = Data.H2;
H2 = fftn(H2);

H1 = H1.*conj(H2);
clearvars H2;
H1 = fftshift(ifftn(H1));

Data.G = H1./numel(H1);
G = H1./numel(H1);
clearvars H1 H2
