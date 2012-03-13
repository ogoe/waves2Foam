function [k,etaNew] = generateStreamFile(H,h,T,EorS,uEorS,N)
% GENERATESTREAMFILE(H,h,T,EorS,uEorS)
close all

nsteps = 10;
g = 9.81;
if nargin == 5
    N = 15;
end

[eta, B, Q, c, k, R, uBar] = StreamFunctionCoefficientsPeriod(N,H,h,T,uEorS,EorS,nsteps,g);

% size(eta)
% size(B)

j=1:N;
for i=0:N
   U(i+1) = -uBar + sum(B .* ...
          cosh(j * k * (eta(i+1) + h)) ./ cosh(j * k*h) ...
          .* cos(j * i * pi / N));
    V(i+1) = sum(B .* ...
          sinh(j * k * (eta(i+1) + h)) ./ cosh(j * k*h) ...
          .* sin(j * i * pi / N));
end

matrix = ones(N+1);

for i=0:N
    matrix(i+1,1:end-1) = cos((1:N) * i * pi / N);
end

L = 2 * pi / k

A = matrix\eta';

sum(A);

x = linspace(0,2*pi / k,1000);
t = linspace(0,T,1000);
etaNew = zeros(size(x)) + A(end);
etaT = etaNew;
for j=1:N
    etaNew = etaNew + A(j) * cos(j * k * x);
    etaT = etaT + A(j) * cos(- j * 2 * pi / T * t);
end

fid = fopen('streamCoeffs','wt');
fprintf(fid,'%f\n',k);
fprintf(fid,'%f\n',2 * pi / T);
fprintf(fid,'%f\n',uBar);
for i=1:N
    fprintf(fid,'%g\n%g\n',A(i),B(i));
end
fclose(fid);

fid = fopen('streamNewApproach','wt');
fprintf(fid,'    waveType\tstreamFunction;\n');
fprintf(fid,'    N\t\t%.0f;\n',N);
fprintf(fid,'    depth\t%f;\n',h);
fprintf(fid,'    omega\t%f;\n',2 * pi / T);
fprintf(fid,'    phi\t\t%f;\n',0);
fprintf(fid,'    waveNumber\twaveNumber [0 -1 0 0 0 0 0] (%f 0.0 0.0);\n',k);
fprintf(fid,'    uBar\t%g;\n',uBar);
fprintf(fid,'    A\t\tnonuniform List<scalar>\t%.0f\n    (\n',N);
for i=1:N
    fprintf(fid,'    %g\n',A(i));
end
fprintf(fid,'    );\n');
fprintf(fid,'    B\t\tnonuniform List<scalar>\t%.0f\n    (\n',N);
for i=1:N
    fprintf(fid,'    %g\n',B(i));
end
fprintf(fid,'    );\n');
fclose(fid);



figure, plot(A),title('A')
