

%A continuous model of the Keller-Segel system used within the hybrid 
%Model.

Dn = 1/10;
Ds = 1.5;
k2 = 10^-3;

dt = 0.0001;
tend = 1;
Xlen = 1;
J = 20;

dx = Xlen/J;
steps = tend/dt
mu = dt/dx^2

N = zeros(int32(steps),J);
S = zeros(int32(steps),J);
preDif = zeros(1,J);

constK = 2;
ks = 0.9;
idx = 2:1:(J-1);


x = meshgrid(linspace(0,Xlen,J));
%x = 15*(x - x.^2).*exp(-50 .*((x - 0.5).^2));
%N(1,:) = x(1,:) + abs(0.01*randn(1,length(x(1,:))));

N(1,(J/2-2):(J/2+2)) = 2;

S(1,:) = 1;


for t = 1:1:steps
    
    betaP = (N(t,idx+1) .* (constK./S(t,idx+1)) + (N(t,idx) .* (constK ./ S(t,idx)))) ./ 2;
    betaN = (N(t,idx-1) .* (constK./S(t,idx-1)) + (N(t,idx) .* (constK ./ S(t,idx)))) ./ 2;
        
    N(t+1,idx) = N(t,idx) +  ...
                 dt/dx^2 .* (Dn .* ( N(t,idx-1) + N(t,idx+1) - 2.*N(t,idx)) ...
                     - (betaP .* (S(t,idx+1) - S(t,idx)) ...
                     -  betaN .* (S(t,idx)- S(t,idx-1))));
    S(t+1,idx) = S(t,idx) + dt*( Ds .* (S(t,idx+1) + S(t,idx-1) - 2*S(t,idx))./(dx^2) - ks.*N(t,idx)); 
    
    %Treat the boundarys:
    N(t+1,1)   = N(t,1) + dt/dx^2 * (Dn*( 2*N(t,2) - 2*N(t,1))) - 2*betaP(1)*(S(t,2)-S(t,1));
    N(t+1,end) = N(t,end) + dt/dx^2 * (Dn*(2*N(t,end-1) - 2*N(t,end))) - 2*betaP(end)*(S(t,end-1)-S(t,end));
    
    S(t+1,1)   = S(t,1) + dt * ( (2*S(t,2) - 2*S(t,1))./dx^2 - ks.*N(t,1));
    S(t+1,end) = S(t,end) + dt * ( (2*S(t,end-1) - 2*S(t,end))/dx^2 - ks.*N(t,end));
    
   

end

figure(1);clf;surf(N);shading('flat');
figure(2);clf;surf(S);shading('flat');