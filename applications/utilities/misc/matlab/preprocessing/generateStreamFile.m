function [k,etaNew] = generateStreamFile(H,h,T,EorS,uEorS,N)
% GENERATESTREAMFILE(H,h,T,EorS,uEorS, N)
% 
% This function computes the stream function coefficients and output them
% to a file, which can be pasted into waveProperties
%
% H:     Wave height [m]
% h:     Water depth [m]
% T:     Wave period [s]
% EorS:  String. Either 'Stokes' or 'Eulerian' drift.
% uEorS: Magnitude of the Stokes or Eulerian drifts. E.g. 'Stokes' and 0
%        gives a zero Stokes drift, e.g. a closed wave flume.
% N:     The number of stream function coefficients. N=1 yields Airy
%        theory.
%
% Modified from a stream function script by Harry Bingham, Technical
% University of Denmark.
%
% Niels G. Jacobsen, Technical University of Denmark
%

close all

nsteps = 10;
g = 9.81;
if nargin == 5
    N = 15;
end

% Solving for the stream function coefficients
[eta, B, Q, c, k, R, uBar] = StreamFunctionCoefficientsPeriod(N,H,h,T,uEorS,EorS,nsteps,g);

% Solving for the surface elevation coefficients
matrix = ones(N+1);

for i=0:N
    matrix(i+1,1:end-1) = cos((1:N) * i * pi / N);
end

A = matrix\eta';

% Writing the stream function coefficients in the format needed for
% waves2Foam
fid = fopen('streamFunctionCoefficient','wt');
fprintf(fid,'    waveType\tstreamFunction;\n');
fprintf(fid,'    N\t\t%.0f;\n',N);
fprintf(fid,'    depth\t%f;\n',h);
fprintf(fid,'    omega\t%f;\n',2 * pi / T);
fprintf(fid,'    phi\t\t%f;\n',0);
fprintf(fid,'    waveNumber\t (%f 0.0 0.0);\n',k);
fprintf(fid,'    uBar\t%g;\n',uBar);
fprintf(fid,'    A\t\tnonuniform List<scalar>\t%.0f\n    (\n',N);
for i=1:N
    fprintf(fid,'        %g\n',A(i));
end
fprintf(fid,'    );\n');
fprintf(fid,'    B\t\tnonuniform List<scalar>\t%.0f\n    (\n',N);
for i=1:N
    fprintf(fid,'        %g\n',B(i));
end
fprintf(fid,'    );\n');
fclose(fid);

% Stretch the figure
figure

hpp = get(gcf,'PaperPosition');
hpp(4) = 2 * hpp(4);
set(gcf,'PaperPosition',hpp);

hp = get(gcf,'Position');
hp(4) = 2 * hp(4);
set(gcf,'Position',hp);

% Plotting A and B coefficients and the surface elevation
subplot(3,1,1), hold on, grid on
plot(abs(A(1:end-1)),'linewidth',2), xlim([1,N])
set(gca,'YScale','Log','Fontname','Times','Fontsize',15)
ylabel('|A|'), xlabel('N')

subplot(3,1,2), hold on, grid on
plot(abs(B),'linewidth',2), xlim([1,N])
set(gca,'YScale','Log','Fontname','Times','Fontsize',15)
ylabel('|B|'), xlabel('N')

t = linspace(0,T,1000);
etaT = zeros(size(t)) + A(end);
for j=1:N
    etaT = etaT + A(j) * cos(- j * 2 * pi / T * t);
end

subplot(3,1,3), hold on, grid on
plot(t / T, etaT,'linewidth',2), xlim([0 1])
set(gca,'Fontname','Times','Fontsize',15)
xlabel('t/T, [-]'), ylabel('\eta, [m]')


