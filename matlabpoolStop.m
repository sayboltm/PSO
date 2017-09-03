function [ null ] = matlabpoolStop(  )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if verLessThan('matlab','8.2.0')
    if matlabpool('size') > 0
        %set(findResource(),'ClusterSize',8);
        matlabpool close
    end
else
    if isempty(gcp('nocreate'))
        delete(gcp('nocreate'));
    end
end

end

