function cnoidalFirst(h, H, T, stream, g)
% cnoidal(h, H, T, g)
%
% Solves for the parameter "m" used in the elliptic integrals and the
% jacobian elliptic functions. Bi-section is used above Matlab's solvers,
% as the elliptic functions and integrals are undefined for m >= 1 and the
% solvers are not readily available for constraint optimization.
%
% Writes a file, STREAM, with the necessary parameters to be used together
% with the OpenFOAM implementation of this wave theory.
%
% h = depth
% H = wave height
% T = wave period
% stream = file name (default = 'cnoidal.dat')
% g = acceleration due to gravity
%
% For the solution to m, equation 5.2(20), 5.2(22) and 5.2(23) are solved
% in combination. These equations are from
% 
% Hydrodynamics of Coastal Regions
% Ib. A. Svendsen and Ivar G. Jonsson
% Den Private Ingenioerfond
% Technical University of Denmark
% 1982
%
% Author: Niels G. Jacobsen
% 2010-05-26
%

close all

if nargin < 5
    g = 9.81;
end
if nargin < 4
    stream = 'cnoidal.dat';
end

m0 = 1e-15;
m2 = 1 - 1e-15;

m   = linspace(m0,m2,2001);
err = zeros(size(m));

for i=1:size(m,2)
    err(i) = evalExpression(m(i), h, H, T, g);
end

I      = find(abs(imag(err)) > 1e-14);
err(I) = [];
m(I)   = [];

min(m)
I      = find(err == max(err));
m0     = m(I);
m0
m1     = 0.5 * (m0 + m2);
plot(m, real(err))
figure

i = 0;
maxIter = 5000;

while true
    err0 = evalExpression(m0, h, H, T, g);
    err1 = evalExpression(m1, h, H, T, g);
    err2 = evalExpression(m2, h, H, T, g);

    if i == 0
        if (sign(err0) == sign(err2))
            disp('Error! The error function has only one sign!')
            i = -1;
            break
        end
    end

    if (sign(err0) == sign(err1))
        m0 = m1;
    else
        m2 = m1;
    end



    if abs(err1) < 1e-10 || i > maxIter
        break
    end

    m1     = 0.5 * (m0 + m2);
    i = i + 1;
end

if i~=-1
    m      = m1;
    [K,E]  = ellipke(m1);
    A      = 2 / m - 1 - 3 / m * E / K;
    etaMin = (1 / m * (1 - E / K) - 1) * H;
    L      = sqrt( 16 * m * K^2 * h^3 / (3 * H));
    c      = sqrt(g * h * (1 + A * H / h));

    fid = fopen(stream, 'wt');

    fprintf(fid,'    waveType\tcnoidalFirst;\n');
    fprintf(fid,'    depth\t%.4f;\n',h);
    fprintf(fid,'    height\t%.6f;\n',H);
    fprintf(fid,'    omega\t%.8f;\n', 2 * pi / T);
    fprintf(fid,'    length\t%.8f;\n',L);
    fprintf(fid,'    celerity\t%.8f;\n',c);
    fprintf(fid,'    direction\t(     );\n',L);
    fprintf(fid,'    m\t\t%.14f;\n',m);

    fclose(fid);
    !cat cnoidal.dat

    if true
        x            = linspace(0,3, 1001);
        [sn, cn, dn] = ellipj(- 2 * K * x, m);
        eta          = etaMin + H * cn.^2;

        subplot(2,1,1)
        plot(x * L, etaMin + H * cn.^2);
        ylabel('\eta, [m]');

        detadx = (eta(3:end) - eta(1:end-2)) / ((x(3) - x(1)) * L);
        subplot(2,1,2)
        plot(x(2:end-1) * L, detadx);
        ylabel('d\eta / dx [-]'); xlabel('x, [m]')

    end
end





function err = evalExpression(m, h, H, T, g)

[K,E] = ellipke(m);

A     = 2 / m - 1 - 3 / m * E / K;

err   = T * sqrt(g / h) * sqrt(1 + H / h * A) - sqrt( 16 * h / (3 * H) * m * K^2);


