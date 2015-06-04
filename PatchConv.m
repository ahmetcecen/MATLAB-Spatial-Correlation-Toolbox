function GG=PatchConv(cutoff,winmulti,corrtype,DataFile,varargin)
% Experimental function for patched convolution. A by product of
% ongoing computational materials science research at MINED@Gatech.
%
% Author: Ahmet Cecen
% Author Web: ahmetcecen.github.io
% Group Web: mined.gatech.edu
%
% Inputs:
%    -cutoff: Maximum vector size as measured from center pixel to center
%             pixel.
%    -winmulti: The patched windowsize as an integer multiplier of cutoff.
%    -corrtype: 'auto' or 'cross'
%    -DataFile: String refering to a filename of a v7.3 mat file. The file
%    needs to have 2 variables:
%           + H1 is the input structure, ALREADY PADDED by cutoff-1  
%             on all 6 sides.
%           + Hsize is size(H1), this saves memory later.
%    -DataFile2: Optional for cross correlation, same rules as DataFile.
%
% Entering Signals Have to All be the Same Size
%
% The current gridding scheme may fail at obscenely high number of windows.
% 

% Memory Map
m=matfile(DataFile);

% Window Indexing Variables
DataSize=m.Hsize-[2*(cutoff-1),2*(cutoff-1),2*(cutoff-1)];
nx=DataSize(1);
ny=DataSize(2);
nz=DataSize(3);

xtra=mod(nx,(winmulti*cutoff));
ytra=mod(ny,(winmulti*cutoff));
ztra=mod(nz,(winmulti*cutoff));

xbc=floor(nx/(winmulti*cutoff));
ybc=floor(ny/(winmulti*cutoff));
zbc=floor(nz/(winmulti*cutoff));

xbuf=floor(xtra/xbc);
ybuf=floor(ytra/ybc);
zbuf=floor(ztra/zbc);

xwinsize=(winmulti*cutoff)+xbuf;
ywinsize=(winmulti*cutoff)+ybuf;
zwinsize=(winmulti*cutoff)+zbuf;

xbuftra=mod(xtra,xbuf);
ybuftra=mod(ytra,ybuf);
zbuftra=mod(ztra,zbuf);

xwinlist=[xwinsize*ones(xbc-1,1);(xwinsize+xbuftra)];
ywinlist=[ywinsize*ones(ybc-1,1);(ywinsize+ybuftra)];
zwinlist=[zwinsize*ones(zbc-1,1);(zwinsize+zbuftra)];

% Main Loop
GG = zeros((1+2*(cutoff-1)),(1+2*(cutoff-1)),(1+2*(cutoff-1))); 


switch corrtype 
	case 'auto'
	
		for xx=1:xbc
			for yy=1:ybc
				for zz=1:zbc
				
					% Grab H1 Window
					HH1 = m.H1(((sum(xwinlist(1:xx)))-xwinlist(xx)+1):(sum(xwinlist(1:xx))+2*(cutoff-1)),...
						((sum(ywinlist(1:yy)))-ywinlist(yy)+1):(sum(ywinlist(1:yy))+2*(cutoff-1)),...
						((sum(zwinlist(1:zz)))-zwinlist(zz)+1):(sum(zwinlist(1:zz))+2*(cutoff-1)));
						
					% Create H2 Mask
					MM1 = padarray(ones(xwinlist(xx),ywinlist(yy),zwinlist(zz)),[cutoff-1,cutoff-1,cutoff-1],0,'both');        
					
					% Compute Convolution
					G = fftshift(ifftn(fftn(HH1).*conj(fftn(HH1.*MM1))));
					
					G = G((floor((size(G,1)/2)+1)-(cutoff-1)):(floor((size(G,1)/2)+1)+(cutoff-1)),...
						(floor((size(G,2)/2)+1)-(cutoff-1)):(floor((size(G,2)/2)+1)+(cutoff-1)),...
						(floor((size(G,3)/2)+1)-(cutoff-1)):(floor((size(G,3)/2)+1)+(cutoff-1)));
				 
					GG = GG + G;           
	
				end
			end
		end
		
	case 'cross'

		n=matfile(varargin{1});
	
		for xx=1:xbc
			for yy=1:ybc
				for zz=1:zbc
				
					% Grab H1 Window
					HH1 = m.H1(((sum(xwinlist(1:xx)))-xwinlist(xx)+1):(sum(xwinlist(1:xx))+2*(cutoff-1)),...
						((sum(ywinlist(1:yy)))-ywinlist(yy)+1):(sum(ywinlist(1:yy))+2*(cutoff-1)),...
						((sum(zwinlist(1:zz)))-zwinlist(zz)+1):(sum(zwinlist(1:zz))+2*(cutoff-1)));
						
					% Grab H2 Window
					HH2 = n.H1(((sum(xwinlist(1:xx)))-xwinlist(xx)+1):(sum(xwinlist(1:xx))+2*(cutoff-1)),...
						((sum(ywinlist(1:yy)))-ywinlist(yy)+1):(sum(ywinlist(1:yy))+2*(cutoff-1)),...
						((sum(zwinlist(1:zz)))-zwinlist(zz)+1):(sum(zwinlist(1:zz))+2*(cutoff-1)));
						
					% Create H2 Mask
					MM1 = padarray(ones(xwinlist(xx),ywinlist(yy),zwinlist(zz)),[cutoff-1,cutoff-1,cutoff-1],0,'both');

					% Compute Convolution
					G = fftshift(ifftn(fftn(HH1).*conj(fftn(HH2.*MM1))));
					
					G = G((floor((size(G,1)/2)+1)-(cutoff-1)):(floor((size(G,1)/2)+1)+(cutoff-1)),...
						(floor((size(G,2)/2)+1)-(cutoff-1)):(floor((size(G,2)/2)+1)+(cutoff-1)),...
						(floor((size(G,3)/2)+1)-(cutoff-1)):(floor((size(G,3)/2)+1)+(cutoff-1)));
				 
					GG = GG + G;           
					
				end
			end
		end
		
end

