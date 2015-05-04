
% dn = rN(1 - N/K)

% growth rate 
r = 2.1;

% carrying capacity.
K = (10^8);

%cell number:
N = 0:10^1:2*K;
dN = zeros(1,length(N));

for i = 1:length(N)
   dN(i) = r*N(i)*(1 - N(i)/K);
end

figure(1)
plot(N,dN)
xlabel('# cells (N)')
ylabel('population growth rate (dN)')
grid on

%solve with forward euler and plot different trajectories:
for e = 0:1:9;
    step = 0.001;
    tend = 15;

    N = ones(1,floor(tend/step));
    N(1) = 10^e;

    for i = 1:floor(tend/step)
       
        N(i+1) = N(i) + step*(r*N(i)*(1 - N(i)/K));
    
    end

    figure(2)
    plot(0:step:tend,N)
    xlabel('time (t)')
    ylabel('# cells (N)')
    axis([0 10 0 10^8.25])
    hold on
    grid on
end