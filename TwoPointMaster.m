function GG = TwoPointMaster(memtype,corrtype,varargin)

%G = TwoPointMaster('full','auto',H1)

%G = TwoPointMaster('full','cross',H1,H2)
 
%G = TwoPointMaster('patched','auto',DataFile,cutoff,winmulti)

%G = TwoPointMaster('patched','cross',DataFile,DataFile2,cutoff,winmulti)
 
%G = TwoPointMaster(memtype,corrtype,varargin)
  
%(cutoff,winmulti,corrtype,DataFile,varargin)

switch memtype
	
	case 'full'
	
		switch corrtype 

			case 'auto'
				H1 = varargin{1};
				GG = fftshift(ifftn(fftn(H1).*conj(fftn(H1))));
				
				if ismatrix(H1)
					GG = GG((floor((size(GG,1)/2)+1)-(cutoff-1)):(floor((size(GG,1)/2)+1)+(cutoff-1)),...
								(floor((size(GG,2)/2)+1)-(cutoff-1)):(floor((size(GG,2)/2)+1)+(cutoff-1)));
				elseif ndims(H1)==3
					GG = GG((floor((size(GG,1)/2)+1)-(cutoff-1)):(floor((size(GG,1)/2)+1)+(cutoff-1)),...
								(floor((size(GG,2)/2)+1)-(cutoff-1)):(floor((size(GG,2)/2)+1)+(cutoff-1)),...
								(floor((size(GG,3)/2)+1)-(cutoff-1)):(floor((size(GG,3)/2)+1)+(cutoff-1)));
				else
					error('Incorrect Number of Dimensions!');
				end
				
			case 'cross'
				H1 = varargin{1};
				H2 = varargin{2};
				GG = fftshift(ifftn(fftn(H1).*conj(fftn(H2))));			
							
				if ismatrix(H1)
					GG = GG((floor((size(GG,1)/2)+1)-(cutoff-1)):(floor((size(GG,1)/2)+1)+(cutoff-1)),...
								(floor((size(GG,2)/2)+1)-(cutoff-1)):(floor((size(GG,2)/2)+1)+(cutoff-1)));
				elseif ndims(H1)==3
					GG = GG((floor((size(GG,1)/2)+1)-(cutoff-1)):(floor((size(GG,1)/2)+1)+(cutoff-1)),...
								(floor((size(GG,2)/2)+1)-(cutoff-1)):(floor((size(GG,2)/2)+1)+(cutoff-1)),...
								(floor((size(GG,3)/2)+1)-(cutoff-1)):(floor((size(GG,3)/2)+1)+(cutoff-1)));
				else
					error('Incorrect Number of Dimensions!');
				end							
							
							
		end
		
	case 'patched'
	
		% Window Indexing Variables
		
        m=matfile(DataFile);
		
        % If 3D
        if length(Hsize)==3
			DataSize=m.Hsize-[2*(cutoff-1),2*(cutoff-1),2*(cutoff-1)];
			nz=DataSize(3);		
			ztra=mod(nz,(winmulti*cutoff));		
			zbc=floor(nz/(winmulti*cutoff));
			zbuf=floor(ztra/zbc);
			zwinsize=(winmulti*cutoff)+zbuf;
			zbuftra=mod(ztra,zbuf);
			zwinlist=[zwinsize*ones(zbc-1,1);(zwinsize+zbuftra)];
			GG = zeros((1+2*(cutoff-1)),(1+2*(cutoff-1)),(1+2*(cutoff-1))); 
        end

        if length(Hsize)==2
            DataSize=m.Hsize-[2*(cutoff-1),2*(cutoff-1)];
			GG = zeros((1+2*(cutoff-1)),(1+2*(cutoff-1))); 
        end
	
		nx=DataSize(1);
		ny=DataSize(2);

		xtra=mod(nx,(winmulti*cutoff));
		ytra=mod(ny,(winmulti*cutoff));

		xbc=floor(nx/(winmulti*cutoff));
		ybc=floor(ny/(winmulti*cutoff));

		xbuf=floor(xtra/xbc);
		ybuf=floor(ytra/ybc);

		xwinsize=(winmulti*cutoff)+xbuf;
		ywinsize=(winmulti*cutoff)+ybuf;

		xbuftra=mod(xtra,xbuf);
		ybuftra=mod(ytra,ybuf);

		xwinlist=[xwinsize*ones(xbc-1,1);(xwinsize+xbuftra)];
		ywinlist=[ywinsize*ones(ybc-1,1);(ywinsize+ybuftra)];
		

		
		% Main Loop
		
		switch corrtype 
		
			case 'auto'
				% Memory Map
                
				if length(Hsize)==2
				
					for xx=1:xbc
						for yy=1:ybc
							
							% Grab H1 Window
							HH1 = m.H1(((sum(xwinlist(1:xx)))-xwinlist(xx)+1):(sum(xwinlist(1:xx))+2*(cutoff-1)),...
								((sum(ywinlist(1:yy)))-ywinlist(yy)+1):(sum(ywinlist(1:yy))+2*(cutoff-1)));
								
							% Create H2 Mask
							MM1 = padarray(ones(xwinlist(xx),ywinlist(yy)),[cutoff-1,cutoff-1],0,'both');        
							
							% Compute Convolution
							G = fftshift(ifftn(fftn(HH1).*conj(fftn(HH1.*MM1))));
							
							G = G((floor((size(G,1)/2)+1)-(cutoff-1)):(floor((size(G,1)/2)+1)+(cutoff-1)),...
								(floor((size(G,2)/2)+1)-(cutoff-1)):(floor((size(G,2)/2)+1)+(cutoff-1)));
						 
							GG = GG + G;           
			
						end
					end
								
				elseif length(Hsize)==3
				
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
				
				else
					error('Incorrect Number of Dimensions!');

				end
				
			case 'cross'
				% Memory Map
				n=matfile(DataFile2);
				
				if length(Hsize)==2
				
					for xx=1:xbc
						for yy=1:ybc
							
							% Grab H1 Window
							HH1 = m.H1(((sum(xwinlist(1:xx)))-xwinlist(xx)+1):(sum(xwinlist(1:xx))+2*(cutoff-1)),...
								((sum(ywinlist(1:yy)))-ywinlist(yy)+1):(sum(ywinlist(1:yy))+2*(cutoff-1)));
								
							% Grab H2 Window
							HH2 = n.H1(((sum(xwinlist(1:xx)))-xwinlist(xx)+1):(sum(xwinlist(1:xx))+2*(cutoff-1)),...
								((sum(ywinlist(1:yy)))-ywinlist(yy)+1):(sum(ywinlist(1:yy))+2*(cutoff-1)));
								
							% Create H2 Mask
							MM1 = padarray(ones(xwinlist(xx),ywinlist(yy)),[cutoff-1,cutoff-1],0,'both');

							% Compute Convolution
							G = fftshift(ifftn(fftn(HH1).*conj(fftn(HH2.*MM1))));
							
							G = G((floor((size(G,1)/2)+1)-(cutoff-1)):(floor((size(G,1)/2)+1)+(cutoff-1)),...
								(floor((size(G,2)/2)+1)-(cutoff-1)):(floor((size(G,2)/2)+1)+(cutoff-1)));
						 
							GG = GG + G;           
								
						end
					end					
				
				elseif length(Hsize)==3
				
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
				
				else
					error('Incorrect Number of Dimensions!');				
				end
				
		end
end









