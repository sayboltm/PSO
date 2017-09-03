function [ null ] = matlabpoolStart( num_cores )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if verLessThan('matlab','8.2.0')
    if matlabpool('size')==0
        set(findResource(),'ClusterSize',num_cores);
        %matlabpool open int(num_cores)
        matlabpool('open', 'local', num_cores)
    end
else
    if isempty(gcp('nocreate'))
        parpool;
    end
end

end

