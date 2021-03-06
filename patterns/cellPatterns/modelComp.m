
%% Constants
DRhoMin = 10;
DRhoMax = 450;
Dh = 400;
Dn = 800;
n0 = 15*10^8;
RhoS = 15*10^8;
gamma = 0.7;
alpha = 1; %missing?????!?
beta = 1.04;
kn = 1;
Kn = 10^9;
Kh = 4*10^8;

% sumlation parameters
dt = 0.01;
dx = 0.1;
tEnd = 0.1;
gridLength = 1;

%% test the muFunction.

h = 4*10^7:10^5:4*10^9;
for i = 1:length(h)
    sols(i) = muFunction(h(i),DRhoMin,DRhoMax);
end

plot(h,sols);


%% numerical simulation.
%Initial cell density.
Rnew = ones(floor(gridLength/dx),floor(gridLength/dx));
R = Rnew;
%AHL-concentrations
Hnew = zeros(floor(gridLength/dx),floor(gridLength/dx));
H = Hnew;
%Nutrient-levels.
Nnew = ones(floor(gridLength/dx),floor(gridLength/dx));
N = Nnew;

for t = 1:1:floor(tEnd/dt)
    
    muR = muFunction(H,DRhoMin,DRhoMax) * R;
          
    for i = 2:1:(floor(gridLength/dx)-1)
        for j = 2:1:(floor(gridLength/dx)-1)
                
            deltaR = (muR(i+1,j) + muR(i-1,j) + muR(i,j+1) + muR(i,j-1) - 4*muR(i,j))/dx^2; 
            deltaH = (H(i+1,j) + H(i-1,j) + H(i,j+1) + H(i,j-1) - 4*H(i,j))/dx^2;
            deltaN = (N(i+1,j) + N(i-1,j) + N(i,j+1) + N(i,j-1) - 4*N(i,j))/dx^2;
                    
            Rnew(i,j) = R(i,j) + dt*(deltaR + gamma*N(i,j)^2*R(i,j) / (N(i,j)^2 + Kn^2));
            Hnew(i,j) = H(i,j) + dt*(Dh*deltaH + alpha*R(i,j)  - beta*H(i,j));
            Nnew(i,j) = N(i,j) + dt*(Dn*deltaN - (kn * gamma * N(i,j)^2 * R(i,j))/(N(i,j)^2 + Kn^2));           
        end
    end
    R = Rnew;
    H = Hnew;
    N = Nnew;
            
end 

figure(1); clf; surf(H); title('H')
figure(2); clf; surf(R); title('R')
figure(3); clf; surf(N); title('N')

