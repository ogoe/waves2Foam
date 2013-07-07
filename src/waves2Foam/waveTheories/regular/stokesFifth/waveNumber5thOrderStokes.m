function [k,eta] = waveNumber5thOrderStokes(H,T,h,EulerOrStokes,valOfEuOrStok)

format long g
format compact

g = 9.81;
kInit = fsolve(@(k) (2 * pi / T)^2 - g * k * tanh(k * h),1);

EulerOrStokes

if strcmp(EulerOrStokes,'Euler')
    
elseif strcmp(EulerOrStokes,'Stokes')
   k = fsolve(@(k) stokes(k,H,T,h,valOfEuOrStok,g),kInit);
else
    return;
end

L = 2 * pi / k;
Ursell = H * L^2 / h^3;


testCoeffs = false;
if testCoeffs == true
    kd = 0.753982;
    k = kd / h;
end

S = sech(2 * k * h);

% A-coefficients
A11 = 1 / sinh(k * h);
A22 = 3 * S^2 / (2 * (1 - S)^2);
A31 = (-4 - 20 * S + 10 * S^2 - 13 * S^3) / (8 * sinh(k * h) * (1 - S)^3);
A33 = (-2 * S^2 + 11 * S^3) / (8 * sinh(k * h) * (1 - S)^3);
A42 = (12 * S - 14 * S^2 - 264 * S^3 - 45 * S^4 - 13 * S^5) / (24 * (1 - S)^5);
A44 = (10 * S^3 - 174 * S^4 + 291 * S^5 + 278 * S^6) / (48 * (3 + 2 * S) * (1 - S)^5);
A51 = (-1184 + 32 * S + 13232 * S^2 + 21712 * S^3 + 20940 * S^4 + 12554 * S^5 - 500 * S^6 - 3341 * S^7 - 670 * S^8) / (64 * sinh(k * h) * (3 + 2 * S) * (4 + S) * (1 - S)^6);
A53 = (4 * S + 105 * S^2 + 198 * S^3 - 1376 * S^4 - 1302 * S^5 - 117 * S^6 + 58 * S^7) / (32 * sinh(k * h) * (3 + 2 * S) * (1 - S)^6);
A55 = (-6 * S^3 + 272 * S^4 - 1552 * S^5 + 852 * S^6 + 2029 * S^7 + 430 * S^8) / (64 * sinh(k * h) * (3 + 2 * S) * (4 + S) * (1 - S)^6);

% B-coefficients
B22 = coth(k * h) * (1 + 2 * S) / (2 * (1 - S));
B31 = -3 * (1 + 3 * S + 3 * S^2 + 2 * S^3) / (8 * (1 - S)^3);
B42 = coth(k * h) * (6 - 26 * S - 182 * S^2 - 204 * S^3 - 25 * S^4 + 26 * S^5) / (6 * (3 + 2 * S) * (1 - S)^4);
B44 = coth(k * h) * (24 + 92 * S + 122 * S^2 + 66 * S^3 + 67 * S^4 + 34 * S^5) / (24 * (3 + 2 * S) * (1 - S)^4);
B53 = 9 * (132 + 17 * S - 2216 * S^2 - 5897 * S^3 - 6292 * S^4 - 2687 * S^5 + 194 * S^6 + 467 * S^7 + 82 * S^8) / (128 * (3 + 2 * S) * (4 + S) * (1 - S)^6);
B55 = 5 * (300 + 1579 * S + 3176 * S^2 + 2949 * S^3 + 1188 * S^4 + 675 * S^5 + 1326 * S^6 + 827 * S^7 + 130 * S^8) / (384 * (3 + 2 * S) * (4 + S) * (1 - S)^6);


% C-coefficients
C0 = sqrt(tanh(k * h));
C2 = sqrt(tanh(k * h)) * (2 + 7 * S^2) / (4 * (1 - S)^2);
C4 = sqrt(tanh(k * h)) * (4 + 32 * S - 116 * S^2 - 400 * S^3 - 71 * S^4 + 146 * S^5) / (32 * (1 - S)^5);

% D-coefficients
D2 = -sqrt(coth(k * h)) / 2;
D4 = sqrt(coth(k * h)) * (2 + 4 * S + S^2 + 2 * S^3) / (8 * (1 - S)^3);

% E-coefficients
E2 = tanh(k * h) * (2 + 2 * S + 5 * S^2) / (4 * (1 - S)^2);
E4 = tanh(k * h) * (8 + 12 * S - 152 * S^2 - 308 * S^3 - 42 * S^4 + 77 * S^5) / (32 * (1 - S)^5);

epsilon = k * H / 2;
kx = k * linspace(0,L,1000);

eta = (k * h + epsilon * cos(kx) + epsilon^2 * B22 * cos(2 * kx) + epsilon^3 * B31 * (cos(kx) - cos(3 * kx)) + epsilon^4 * (B42 * cos(2 * kx) + B44 * cos(4 * kx)) + epsilon^5 * (-(B53 + B55) * cos(kx) + B53 * cos(3*kx) + B55 * cos(5 * kx))) / k - h;



function f = stokes(k,H,T,h,cS,g)

S = sech(2 * k * h);
C0 = sqrt(tanh(k * h));
C2 = sqrt(tanh(k * h)) * (2 + 7 * S^2) / (4 * (1 - S)^2);
C4 = sqrt(tanh(k * h)) * (4 + 32 * S - 116 * S^2 - 400 * S^3 - 71 * S^4 + 146 * S^5) / (32 * (1 - S)^5);
D2 = -sqrt(coth(k * h)) / 2;
D4 = sqrt(coth(k * h)) * (2 + 4 * S + S^2 + 2 * S^3) / (8 * (1 - S)^3);

f = sqrt(k / g) * cS - 2 * pi / (T * sqrt(g * k)) + C0 + (k * H / 2)^2 * (C2 + D2 / (k * h)) ...
    + (k * H / 2)^4 * (C4 + D4 / (k * h));




function euler(k,H,T,h,cE,g)