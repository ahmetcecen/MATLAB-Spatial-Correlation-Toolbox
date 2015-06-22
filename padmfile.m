function padmfile(datafile,varname,cutoff,savefile)
% Pads a very large memory mapped file with zeros slice by slice.
% A by-product of ongoing computational materials science research 
% at MINED@Gatech.(http://mined.gatech.edu/)
%
% Author: Ahmet Cecen
% Contact: ahmetcecen@gmail.com

%Map Data
m=matfile(datafile);
datasize=size(m,varname);

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

    for ii=1:cutoff-1
        n.H1(1:(nx+2*(cutoff-1)),1:(ny+2*(cutoff-1)),ii)=zeros(nx+2*(cutoff-1),ny+2*(cutoff-1),1);

        % Progress Bar
        progress=progress+1;
        fprintf(reverseStr);
        msg=fprintf('Progress = %.2f %%',100*progress/(nz+2*cutoff-2));                                
        reverseStr = repmat(sprintf('\b'), 1, msg);
    end

    for ii=1:nz
        eval(['current=m.',varname,'(:,:,ii);'])
        n.H1(1:(nx+2*(cutoff-1)),1:(ny+2*(cutoff-1)),cutoff+ii-1)=padarray(current,[cutoff-1 cutoff-1],0,'both');

        % Progress Bar
        progress=progress+1;
        fprintf(reverseStr);
        msg=fprintf('Progress = %.2f %%',100*progress/(nz+2*cutoff-2));                                
        reverseStr = repmat(sprintf('\b'), 1, msg);
    end

    for ii=cutoff-1+nz+1:cutoff-1+nz+cutoff-1
        n.H1(1:(nx+2*(cutoff-1)),1:(ny+2*(cutoff-1)),ii)=zeros(nx+2*(cutoff-1),ny+2*(cutoff-1),1);

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

    for ii=1:cutoff-1
        n.H1(1:(nx+2*(cutoff-1)),ii)=zeros(nx+2*(cutoff-1),1);

        % Progress Bar
        progress=progress+1;
        fprintf(reverseStr);
        msg=fprintf('Progress = %.2f %%',100*progress/(ny+2*cutoff-2));                                
        reverseStr = repmat(sprintf('\b'), 1, msg);
    end

    for ii=1:ny
        eval(['current=m.',varname,'(:,ii);'])
        n.H1(1:(nx+2*(cutoff-1)),cutoff+ii-1)=padarray(current,[cutoff-1 0],0,'both');

        % Progress Bar
        progress=progress+1;
        fprintf(reverseStr);
        msg=fprintf('Progress = %.2f %%',100*progress/(ny+2*cutoff-2));                                
        reverseStr = repmat(sprintf('\b'), 1, msg);
    end

    for ii=cutoff-1+ny+1:cutoff-1+ny+cutoff-1
        n.H1(1:(nx+2*(cutoff-1)),ii)=zeros(nx+2*(cutoff-1),1);

        % Progress Bar
        progress=progress+1;
        fprintf(reverseStr);
        msg=fprintf('Progress = %.2f %%',100*progress/(ny+2*cutoff-2));                                
        reverseStr = repmat(sprintf('\b'), 1, msg);
    end 
    
end