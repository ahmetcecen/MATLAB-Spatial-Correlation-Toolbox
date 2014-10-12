function G = memacp2ptmask(filename)
% Periodic two point statistics with memory mapping. The filename must 
% link to a mat file with H1 and Mask variables. 

% Map the Data
Data = matfile(filename,'Writable',true);

% Find and Save P(R)

Mask = Data.Mask;
Mask = fftn(Mask);
Mask = Mask.*conj(Mask);
Mask = fftshift(ifftn(Mask));

Data.PR = Mask;
clearvars Mask

% Find and Save P(H1H2R)
H1 = Data.H1.*Data.Mask;
H1 = fftn(H1);
H1 = H1.*conj(H1);
H1 = fftshift(ifftn(H1));

Data.PH1H2R = H1;
clearvars H1

% Find G
G = Data.PH1H2R./Data.PR;
Data.G=G;