function [time, x, y, z, etaEns] = readSurfaceElevationNonEquidistantEnsembleAverage( tString, period, dt, startTime, extrap )
% [time, x, y, z, etaEns] = readSurfaceElevationNonEquidistantEnsembleAverage( tString, period, dt, startTime, extrap )
% 
% This function returns the ensemble averaging of the surface elevation
% computed using the waves2Foam utility "surfaceElevation". This method
% assumes non-equidistant time steps in the surfaceElevation.dat file.
%
% The function only works on linux/unix machines.
%
% The input variables are the following:
%
% tString:   A string denoting the folder name, where surfaceElevation.dat
%            is found. Specifically: <rootCase>/surfaceElevation/<tString>
%
% period:    The period over which, the ensemble averaging is performed.
%
% dt:        The delta time. E.g. period = 6.0 s, dt = 1.0 s, gives six
%            instances per period, where the ensemble average is computed.
%
% startTime: Start the ensemble averaging for t >= startTime.
%
% extrap:    As time is non-equistant a interpolation is utilised. If
%            extrap is set to "true", extrapolation is allowed. 
%            [Default = false].
%
% Niels Gjoel Jacobsen
% Technical University of Denmark, 8th of June 2012.
%

if nargin < 5
    extrap = false;
end

% Remove header from surfaceElevation.dat
unix(sprintf('grep -v Time surfaceElevation/%s/surfaceElevation.dat > temporarySurfaceFile.dat',tString));

% Load the data and remove the temporary file
X = load('temporarySurfaceFile.dat');
!rm temporarySurfaceFile.dat

% Extract spatial coordinates and time
time = X(4:end,1);

x    = X(1,2:end);
y    = X(2,2:end);
z    = X(3,2:end);

% Extract instantaneous surface elevations
eta      = X(4:end,2:end);

% Eliminate times, which are close to each other. This is a potential
% problem, if the writing accuracy is so small the times in
% surfaceElevation.dat are not distinct.
I = find( diff(time) < 1e-5 );

time(I + 1)  = [];
eta(I + 1,:) = [];

% Find the number of instances per period over which to perform the
% ensemble averaging
N  = int32(period / dt);

% Initialise return field
etaEns = zeros(N+1,length(x));

% Make interpolation time field
tInterp = (startTime:dt:time(end))';

% Ensemble averaging
for j=1:length(x)
    if extrap
        disp('using extrap')
        etaInterp = interp1(time, eta(:,j), tInterp, 'linear', 'extrap');
    else
        etaInterp = interp1(time, eta(:,j), tInterp);
    end
    
    for i=1:N
        etaEns(i,j) = mean( etaInterp(i:N:end) );
    end
    
end

% Make the first and last field identical (notice, above looping is only 
% performed over N, whereas the size of etaEns is N+1)
etaEns(end,:) = etaEns(1,:);

% Make time for the ensemble averaged field
time = (0:double(N)) * dt;









