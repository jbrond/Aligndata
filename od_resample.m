function [aligneddata] = od_resample(data,SF,offset,drift,interpolationMethod)
%od_resample Is a function that resamples the data based on an estimated
%offset and drift
%   

%Do we need to remove data from the start?
if offset>-1
    if offset==0
        offset = 1;
    end
    
    %Copy the data from raw    
    offsetdata = data(offset:end,:);
else
    %Need to append
    offset = abs(offset);
    offsetdata = [zeros(offset,3) data(1:end,:)];
end

%The original time
orgTime = lstime(0,offsetdata,SF);
%Number of samples in the data series
Nsamples = length(offsetdata);
%Determining the total drift at the last sample from the drift
totalDrift = orgTime(end)*drift/86400+orgTime(end);

%Generating the new resampleded time points
adjustedTime = linspace(0,totalDrift,Nsamples);

%Applying the resampling
aligneddata = interp1(orgTime,offsetdata,adjustedTime,interpolationMethod);

end

