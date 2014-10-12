function G = memacnp2pt(filename)
% Nonperiodic two point statistics with memory mapping. The filename must 
% link to a mat file with H1 variables. 

% Map the Data
Data = matfile(filename,'Writable',true);

% Find P(R)
sizeH = size(Data.H1);

Base = padarray(ones(sizeH),sizeH-1,0,'post');
Base = fftn(Base);
Base = Base.*conj(Base);
Base = fftshift(ifftn(Base));

Data.PR = Base;
clearvars Base

% Find and Save P(H1H2R)
H1 = Data.H1;
H1 = padarray(H1,sizeH-1,0,'post');
H1 = fftn(H1);

H1 = H1.*conj(H1);
H1 = fftshift(ifftn(H1));

Data.PH1H2R = H1;
clearvars H1

% Find G
G = Data.PH1H2R./Data.PR;
Data.G=G;