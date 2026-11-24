% Simple ODE f(x) = X - X^
% d(X(t))/dt = X - X^3


X0 = [0.2, 0.5, 0.9, -0.2, -0.5, -0.9, 0, 2, 5, -2, -5]; 
N = 1000;
t = linspace(0,3,N);
X = zeros(size(T));

figure; close all;
for i = 1:length(X0)
    if ismember(X0(i), [0,1,-1])
        disp("Undefined")
    elseif abs(X0(i)) < 1
        C = X0(i)^2/(1-X0(i)^2);
        X = sign(X0(i)) * sqrt( C*exp(2*t) ./ (1+C*exp(2*t)) );
    else
        C = X0(i)^2/(X0(i)^2-1);
        X = sign(X0(i)) * sqrt( C*exp(2*t) ./ (C*exp(2*t)-1) );
    end

plot(t, X, 'LineWidth', 1.5); hold on;
end
plot(t,ones(N),'--','Color', [0 0 0], 'LineWidth', 2); hold on;
plot(t,-1.*ones(N), '--','Color', [0 0 0], 'LineWidth', 2); hold on;
xlabel('t');
ylabel('X(t)');
grid on;

