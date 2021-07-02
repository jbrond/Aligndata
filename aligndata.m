function [ table, aligneddata, aux, stats ] = aligndata( data1, data2 )
% aligndata is function to temporal align acceleration 
%   Input:
%	data1 : OmGui wave struct (use loadomwav)
%	data2 : OmGui wave struct
% interpolationMethod : Interpolation method for the alignment. Default is
% 'previous'
%
%   Output:
%	table - table with the raw temporal offset values used in the estimating the drift and offset
%	aligneddata - The aligned data
%	aux - The aligned aux data (temperature, light, battery)
%	stats - Descriptive statistics of the model performance for fitting bias and offset


% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

if nargin<3
    disp('No interpolation method specified. Using ''previous'''); 
    interpolationMethod = 'previous';
end

aligneddata = [];
aux = [];
table = [];
stats = [];
   
%Making sure we do have data2
if nargin<2
    disp('No data2 provided');
    return;
end

%Making sure that both time series are based on the same sample rate
if data1.INFO.SampleRate~=data2.INFO.SampleRate
    disp('Both time series must have the same sample rate');
    return;
end

%Sampling frequency
SF = data1.INFO.SampleRate;

%Default nosie filter
[B,A] = butter(4,[0.1 7]./(SF/2));

%Window length (hours)
wndlength = 1;

%Maximum residual distance
maxResidualDistance = 1.0;

%data block used with cross correlation
ccWndSize = round(33.33/(1/SF));
%
dayt = array2table([data2.INFO.orgstart]','VariableNames',{'Day'});

daysamplesstart = 1;
validwnd = 1;

%Data block in number of samples
wndsamples = SF*wndlength*60*60;

%Restrict the data to the shortest data
NallData = min([length(data1.raw),length(data2.raw)]);

%Number of windows in the data
Nwnds = floor((NallData-wndsamples)/wndsamples);

unsyncedsamples = NaN(Nwnds,3);

%
%Stage I
%

for n=1:Nwnds
    
    %Data1 is hip
    h24 = data1.raw(daysamplesstart:daysamplesstart+wndsamples,:);
    w24 = data2.raw(daysamplesstart:daysamplesstart+wndsamples,:);
    
    mactivityw = mean(abs(vmag(w24)-1.0));
    mactivityh = mean(abs(vmag(h24)-1.0));
    
    %Only do cross correlation if we have enough movement
    if (mactivityw>0.01 && mactivityh>0.01)
        
        %Reduce noise in the data
        h = abs(filtfilt(B,A,trunc(vmag(h24)-1,0.068)));
        w = abs(filtfilt(B,A,trunc(vmag(w24)-1,0.068)));

        %Cross correlation coefficient        
        [acora,lag] = xcorr(w,h,ccWndSize,'coeff');
        
        %Find the maximum correlation coeffient
        [~,I] = max(acora);
    
        %Retrieve the coeffient for later regression analysis
        unsyncedsamples(validwnd,1) = data1.INFO.orgstart+datenum(0,0,0,0,0,1)/SF*daysamplesstart;
        unsyncedsamples(validwnd,2) = lag(I);
        unsyncedsamples(validwnd,3) = mactivityw;

        validwnd = validwnd + 1;
    end

    daysamplesstart = daysamplesstart + wndsamples;
    
end

%Remove the last entry
unsyncedsamples = unsyncedsamples(1:validwnd-1,:);

%There is not enough activity in the file to do the alignment
if (validwnd<4)
    aligneddata = data2.raw;
    aux = data2.aux;
    table = [];
    stats = [];
    return;
end

%
%Stage II
%

%Linear fit of the unsynced samples

%Genrate the table
table = array2table(unsyncedsamples,'VariableNames',{'Day','lag','MAD'});

%Do the regression analysis
mdl = fitlm(table,'lag~Day');

%Find the standardized residual that are to far from the mean
I = find(abs(mdl.Residuals.Standardized)>maxResidualDistance);
%Remove them from the data set
table.lag(I) = NaN;

str = sprintf('Excluding %d data points',round(length(I)));
disp(str);

%Excluding more data?
while isempty(I)==0
    mdl = fitlm(table,'lag~Day');
    str = sprintf('Excluding %d number of data points',round(length(I)));
    disp(str);
    maxResidualDistance = exp(maxResidualDistance); %1, 2.7, 14...
    I = find(abs(mdl.Residuals.Standardized)>maxResidualDistance);
    table.lag(I) = NaN;
end

%Estimate the offset
psamples = predict(mdl,dayt);

%Store the results from the regression analysis in stats
stats.mdl = mdl;
stats.time.offset = psamples(1)*1/SF;
stats.time.slope = 1/SF*mdl.Coefficients.Estimate(2);
stats.offset = psamples(1);
stats.slope = mdl.Coefficients.Estimate(2);

str = sprintf('Offset: %2.3f slope %2.3f',psamples(1)*1/SF,1/SF*mdl.Coefficients.Estimate(2));
disp(str);

%
%Stage III
%
%Aligning the acceleration
aligneddata = od_resample(data2.raw,SF,round(stats.offset),stats.time.slope,interpolationMethod);

%Aligning the temperature, light and battery
%This data is sampled at 1 Hz so we need to scale down the offset and the
%drift
offset = floor(stats.offset/SF);
drift = stats.time.slope/SF;

alignedtemp = od_resample(data2.temp,1,offset,drift,interpolationMethod);
alignedlight = od_resample(data2.light,1,offset,drift,interpolationMethod);
alignedbat = od_resample(data2.bat,1,offset,drift,interpolationMethod);

Nsamples = length(alignedtemp);

aux = zeros(Nsamples*SF,1);
aux(1:SF:end) = alignedtemp;
aux(2:SF:end) = alignedlight;
aux(3:SF:end) = alignedbat;

end

