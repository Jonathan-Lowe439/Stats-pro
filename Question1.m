% EM strong test (Euler - Maruyama)

rng(67)  % pseudorandom seed

% setting up stochastic equation as in section 5 of paper 
% dX = f(X)dt + g(X)dW over range t = 0:1
% for simplicity I am using dX = a*X dt + b*x dW

a=2;
b=1;


T = 1;    % range
N = 2^14;  % number of nodes
dt = T/N;  % step size 
dW = sqrt(dt)*randn(1, N);
W = cumsum(dW); % discretized Brownian path 

Initial_condition = 1;


t_list = dt:dt:T; % array of t values separated by steps of dt


% array of the true values of X
Xtrue = Initial_condition*exp((a-0.5*b^2)*t_list+b*W);


% EM strong test
% use a bigger dt as we are approximating


R = 8;
Dt = R*dt;
L = N/R;
xaprox = zeros(1,L);
xem_current = Initial_condition;

for i = 1:L
    step_em = sum(dW(R*(i-1)+1:R*i));
    xem_current = xem_current+Dt*a*xem_current+b*xem_current*step_em;
    xaprox(i) = xem_current;
end

t_aprox_list = 0:Dt:T;



figure;
hold on
plot([0,t_list], [Initial_condition,Xtrue], 'LineWidth', 2);
plot(t_aprox_list, [Initial_condition,xaprox], 'LineWidth', 1.5);
grid on;
xlabel('t');
ylabel('X(t)');
title('Euler-Maruyama Sample Path on [0, 1]');