function [eta, B, Q, c, k, R, uBar] = StreamFunctionCoefficientsPeriod(N,H,h,T,uEorS,EorS,nsteps,g)
% This implementation of the streamfunction follows David R. Fuhrman's 
% lecture notes from the course 'Linear and non-linear wave theory' given
% at the Technical University of Denmark.

% Initial quess on wave number (2 * pi / 1.56 T^2 is the deep water wave number)
k0 = abs(fzero(@(x) (2*pi/T)^2 - g*x*tanh(x*h), 2 * pi / (1.56 * T^2)));

% Initial wave height
H0 = H/nsteps;

% Initial linear quess on the surface elevation
% x0(1:N+1) = H0/2 * cos(k0 * (0:N) * pi / N)
x0(1:N+1) = H0/2 * cos((0:N) * pi / N);

% Initial quess on B_j
x0(N+2:2*N+1) = [pi* H0 / T,zeros(1,N-1)];

% Initial quess on Q
x0(2*N+2) = - 2 * pi / T / k0 * h;

% Initial guess on celerity
x0(2*N+3) = 2 * pi / T / k0;

% Initial quess on k
x0(2*N+4) = k0;

% initial quess on the bernoulli constant
x0(2*N+5) = g/(2*k0)*tanh(k0*h);

% Initial quess on uBar 
x0(2*N+6) = 2 * pi / T / k0;

[x0', solveStream(x0, N, h, H0, T, g, uEorS, EorS)'];

for i=1:nsteps
    x = fsolve(@(x) solveStream(x,N,h,H0 * i,T,g,uEorS,EorS),x0);
    x0 = x;
    [x(1) x(2 * N + 4)];
end

plot(x(1:N+1))
figure
plot(x(N+2:2*N+1))
    
eta = x(1:N+1);
B = x(N+2:2*N+1);
Q = x(2*N+2);
c = x(2*N+3);
k = x(2*N+4);
R = x(2*N+5);
uBar = x(2*N+6);


function F = solveStream(sol,N,h,H,T,g,uEorS,EorS)
% eta  = sol(1:N+1)
% B    = sol(N+2:2N+1)
% Q    = sol(2N+2)
% c    = sol(2N+3)
% k    = sol(2N+4)
% R    = sol(2N+5)
% uBar = sol(2N+6)
sol = (sol(:))';
j = 1:N;

if strcmp(EorS,'Stokes')
    % Q + c h - uS h = 0
    F(1) = sol(2*N+2) + sol(2*N+3) * h - uEorS * h;
else
    F(1) = -sol(2*N+3) + uEorS + sol(2*N+6);
end
% eta_0 - eta_N - H = 0
F(2) = sol(1) - sol(N+1) - H;

% 0.5 (eta_0 + eta_N) + sum_(j=2)^(N-1) eta_j = 0
F(3) = 0.5 * (sol(1)+sol(N+1)) + sum(sol(2:N));

% 2 * pi / (T * c) - k = 0
F(4) = 2 * pi / (T * sol(2*N+3)) - sol(2*N+4);

for i=0:N
    psi = -sol(2*N+6) * (sol(i+1) + h) + sum(sol(N+2:2*N+1) ./ (j*sol(2*N+4)).* ...
          sinh(j * sol(2*N+4) * (sol(i+1) + h)) ./ cosh(j * sol(2*N+4)*h) ...
          .* cos(j * i * pi / N));
    % psi - Q = 0
    F(5+i) = psi - sol(2*N+2);
 
    % 0.5(psi_z^2 + psi_x^2) + g eta - R = 0    
    U = -sol(2*N+6) + sum(sol(N+2:2*N+1) .* ...
          cosh(j * sol(2*N+4) * (sol(i+1) + h)) ./ cosh(j * sol(2*N+4)*h) ...
          .* cos(j * i * pi / N));
    V = sum(sol(N+2:2*N+1) .* ...
          sinh(j * sol(2*N+4) * (sol(i+1) + h)) ./ cosh(j * sol(2*N+4)*h) ...
          .* sin(j * i * pi / N));
    F(6+N+i) = g * sol(i+1) + 0.5 * (U^2 + V^2) - sol(2*N+5);
end
    


