% SDE of f(X) = X - X^3, g(X) = c, for some constant c


randn('state',100)
X0 = [0.2, 0.5, 0.9, -0.2, -0.5, -0.9, 0, 2, 5, -2, -5];
M = length(X0);
c = 2;
T = 3; N = 2^12; dt = T/N;
dW = sqrt(dt)*randn(1,N);
W = cumsum(dW);

R = 2; Dt = R*dt; L = N/R;        % L EM steps of size Dt = R*dt
Xem = zeros(M,L);                 % preallocate for efficiency
Xtemp = X0;
for j = 1:L
    Winc = sum(dW(R*(j-1)+1:R*j));
    Xtemp = Xtemp + Dt*(Xtemp-Xtemp.^3) + c*Winc;
    Xem(:,j) = Xtemp;
end



close all;
plot([0:Dt:T],[X0(:),Xem]), hold off
xlabel('t','FontSize',12)
ylabel('X','FontSize',16,'Rotation',0,'HorizontalAlignment','right')