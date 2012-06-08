function [time, x, y, z, eta] = readSurfaceElevation( tString )
% [time, x, y, z, eta] = readSurfaceElevation( tString )
% 
% This function returns the surface elevation computed using the waves2Foam 
% utility "surfaceElevation". 
%
% The function only works on linux/unix machines.
%
% The input variable is the following:
%
% tString:   A string denoting the folder name, where surfaceElevation.dat
%            is found. Specifically: <rootCase>/surfaceElevation/<tString>
%            [Default = '0']
%
% Niels Gjoel Jacobsen
% Technical University of Denmark, 8th of June 2012.
%

if nargin == 0
    time = '0';
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
eta  = X(4:end,2:end);

% Eliminate times, which are close to each other. This is a potential
% problem, if the writing accuracy is so small the times in
% surfaceElevation.dat are not distinct, and is thus treated here before
% subsequent data processing (especially interpolation on non-distinct data
% sets make the program fail).
I = find( diff(time) < 1e-5 );

time(I + 1)  = [];
eta(I + 1,:) = [];


