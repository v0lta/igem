
% compute the model from "Spatio-Temporal Patterns Generated by 
% Salmonella typhimurium "

%% load parameters.
alpha = 3.0;
beta  = 2.0;
gamma = 0.2;
%gamma = 0;
rho   = 0.03;
Dn    = 0.1;
Dc    = 0.3;
% Bifurcataion parameter.
s = 1;
%pattern propagation.
r = 0.9;
%r = 0;

length = 3;
tend = 0.1;

dt = 0.005;
J = 40;
idx = 2:(J-1);

dx = length/J

steps = floor(tend/dt)

%there are only cells at the center.
[x,y] = meshgrid(linspace(0,1,J));
%N = 15*(x - x.^2).*(y-y.^2).*exp(-50 .*((x - 0.5).^2 + (y - 0.5).^2 ));

N = 0.*x + 0.*y;
N((J/2-10):(J/2+10),(J/2-10):(J/2+10)) = 1;
Nnew = N;


%transmitters are zero initially.
C = 0.*x+0.*y;
Cnew = C;
iChem = C;

%% start computation!
for t = 1:steps
    
      %compute the inner chemotaxis term.
      iChem(idx,idx) = N(idx,idx) ./ ((1 + beta.*C(idx,idx)).^2) ...
              .* ((C(idx+1,idx) + C(idx-1,idx) + C(idx,idx+1) + ...
              C(idx,idx-1))./(2*dx)); 
            
                  %diffusion term.
      Nnew(idx,idx) = N(idx,idx) + dt.*( Dn./(dx^2).*(N(idx+1,idx) + N(idx-1,idx) + ...
                  N(idx,idx+1) + N(idx,idx-1) - 4.*N(idx,idx))  ...
                  ... %chemotaxis term.
                  - (alpha./(2*dx)).*(iChem(idx+1,idx) + iChem(idx-1,idx) ...
                  + iChem(idx,idx+1) + iChem(idx,idx-1)) ...
                  ... %logistic growth, proliferation term.
                  + rho.*N(idx,idx).*( 1  - N(idx,idx)./s));
        
      %Only the growth model.
      %Nnew(i,j) = N(i,j) + dt*( rho*N(i,j)*( 1  - N(i,j)/s));
  
                    %diffusion term
      Cnew(idx,idx) = C(idx,idx) + dt*( Dc/(dx^2) .* (C(idx+1,idx) + ...
                  C(idx-1,idx) + C(idx,idx+1) + ...
                  C(idx,idx-1) - 4.*C(idx,idx))  ...
                  ... %production term
                  + s.*N(idx,idx)./(1 + gamma.*N(idx,idx))  ...
                  ... %consumption term
                  - (C(idx,idx).*(N(idx,idx).^r))./(1 + gamma.*C(idx,idx)) );
    
   N = Nnew;
   C = Cnew;
end
figure(1)
surf(real(N));shading('flat');
figure(2)
surf(real(C));shading('flat');

