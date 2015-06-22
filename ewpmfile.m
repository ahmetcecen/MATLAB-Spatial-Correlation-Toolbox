function ewpmfile(datafile1,varname1,datafile2,varname2,savefile)
% Takes the Element-Wise Product of two very large memory mapped files 
% slice by slice.
% A by-product of ongoing computational materials science research 
% at MINED@Gatech.(http://mined.gatech.edu/)
%
% Author: Ahmet Cecen
% Contact: ahmetcecen@gmail.com

%Map Data1
m=matfile(datafile1);
datasize=size(m,varname1);

%Map Data2
k=matfile(datafile2);

if length(datasize)==3
    
    %Map Savefile
    H1 = zeros(2,2,2); save(savefile,'H1','-v7.3');
    n = matfile(savefile,'Writable',true);

    nx = datasize(1);
    ny = datasize(2);
    nz = datasize(3);

    % Progress Bar Initialize
    progress=0;
    reverseStr = '';

    for ii=1:nz
        eval(['current=m.',varname1,'(:,:,ii);'])
        eval(['current2=k.',varname2,'(:,:,ii);'])
        n.H1(1:nx,1:ny,ii)=current.*current2;

        % Progress Bar
        progress=progress+1;
        fprintf(reverseStr);
        msg=fprintf('Progress = %.2f %%',100*progress/(nz+2*cutoff-2));                                
        reverseStr = repmat(sprintf('\b'), 1, msg);
    end

elseif length(datasize)==2

    %Map Savefile
    H1 = zeros(2,2); save(savefile,'H1','-v7.3');
    n = matfile(savefile,'Writable',true);

    nx = datasize(1);
    ny = datasize(2);

    % Progress Bar Initialize
    progress=0;
    reverseStr = '';

    for ii=1:ny
        eval(['current=m.',varname1,'(:,ii);'])
        eval(['current2=k.',varname2,'(:,ii);'])
        n.H1(1:nx,ii)=current.*current2;

        % Progress Bar
        progress=progress+1;
        fprintf(reverseStr);
        msg=fprintf('Progress = %.2f %%',100*progress/(ny+2*cutoff-2));                                
        reverseStr = repmat(sprintf('\b'), 1, msg);
    end
    
end
