

%A continuous model of the Keller-Segel system used within the hybrid 
%Model.

Dn = 1/10;
Ds = 0;
k2 = 10^-3;

dt = 0.0001;
tend = 0.01;
Xlen = 3;
J = 60;

dx = Xlen/J;
steps = tend/dt;

N = zeros(steps,J);
S = ones(steps,J);
preDif = zeros(1,J);

constK = 2;
idx = 2:1:(J-1);

N(1,J/2) = 10;

for t = 1:1:2
    preDif(1,idx) = constK./S(t,idx) .* N(t,idx) .* (S(t,idx+1)-S(t,idx-1))./(2*dx);
    N(t+1,idx) = N(t,idx) +  ...
                 dt .* (Dn .* ( N(t,idx) + N(t,idx+1) - 2*N(t,idx-1)) ...
                    - (preDif(1,idx+1) - preDif(1,idx-1))./(2*dx));
    S(t+1,idx) = Ds .* (S(t,idx+1) - S(t,idx-1))./(2*dx) - 10^(-3).*N(t,idx); 
end

figure(1);clf;surf(N);
figure(2);clf;surf(S);