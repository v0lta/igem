% The matlab source code used to compute the 1d simulation at the top of
% the page. The 2015 Ku Leuven Igem team.

%% Load constants.
Da = 0.072*10^(-3);     %1*10^(-6)*3600; %(cm)^2/h;
Db = 2.376*10^(-3);     %5*10^(-6)*3600; %(cm)^2/h;
Dr = 26.46*10^(-3);     %(cm)^2/h; 0.02628 = 0.2628*10^(-1) = 2.628*10^(-2)
Dh = 50*10^(-3);        %(cm)^2/h; 

dt = 0.15; %h
tend = 80; %h
Xlen = 8.0; %cm
J = 40; %#

dx = Xlen/J; %cm
steps = tend/dt %#
mu = dt/dx^2 %s/cm^2

A = zeros(int32(steps),J);
B = zeros(int32(steps),J);
R = zeros(int32(steps),J);
H = zeros(int32(steps),J);

preDif = zeros(1,J);

Kc = 0.00095;
gamma = 10^(-5);
As = 1.0 * 10^2; % (cl)^(-1);
Bs = 1.0 * 10^2; % (cl)^(-1);
kr = 5.4199*10^(-4); %fmol/h
kh = 17.9*10^(-4);   %fnmol/h
KAbsorbR = 1/80;
KAbsorbH = 1/48;
idx = 2:1:(J-1);

%% Set intial condition.
x = meshgrid(linspace(0,Xlen,J));
x = -100*(x - x.^2).*exp(-1 .*((x - 4).^2));
%Cell density of celly type A.
A(1,:) = x(1,:) + abs(0.01*randn(1,length(x(1,:))));
%A(1,(J/2-2):(J/2+2)) = 2;
%Cell density of cell type B.
B(1,:) = x(1,:) + abs(0.01*randn(1,length(x(1,:))));
%B(1,(J/2-2):(J/2+2)) = 2;

%Concentration of the replellant.
R(1,:) = 0.01;
H(1,:) = 0.01;

%% run simulation.
for t = 1:1:steps
    
    %% Compute the next time step of the A-Type cells:
    %A(t+1,:) = A(t,:);
    A(t+1,idx) = A(t,idx) + (dt.* (Da/dx^2  .* ( A(t,idx-1) + A(t,idx+1) - 2.*A(t,idx))) ...
                 +gamma.* A(t,idx) .* (1 - A(t,idx)./As));
    
    %% Compute the next time step of the B-Type cells:
    
    %Xn =  -B(t,:) .* (Db.*Kc1./(Kc2 + Kc3.*H(t,:).*R(t,:)));
    Xn = -B(t,:) .* Kc .* H(t,:)./(R(t,:));
   
    betaP = (Xn(idx+1) + Xn(idx))./2; 
    betaN = (Xn(idx-1) + Xn(idx)) ./ 2;
        
    B(t+1,idx) = B(t,idx) +  ...
                 dt .* (1/dx^2  .* (Db.*( B(t,idx-1) + B(t,idx+1) - 2.*B(t,idx)) ...
                     - (betaP .* (R(t,idx+1) - R(t,idx)) ...
                     -  betaN .* (R(t,idx)- R(t,idx-1))))...
                     + gamma .* B(t,idx) .* (1 - B(t,idx)./Bs));
    
    %% Compute the repellant:
    R(t+1,idx) = R(t,idx) + dt*( Dr .* (R(t,idx+1) + R(t,idx-1) -  ...
                 2*R(t,idx))./(dx^2) + kr.*A(t,idx) - KAbsorbR * R(t,idx)); 
             
    %% Compute the AHL    
    H(t+1,idx) = H(t,idx) + dt*( Dh .* (H(t,idx+1) + H(t,idx-1) -  ...
                 2*H(t,idx))./(dx^2) + kh.*A(t,idx) - KAbsorbH * H(t,idx)); 
              
             
       
    %% Symmetric boundary conditions:
    %A
    %A(t+1,1)   = A(t,1) + dt/dx^2 * (Da*( A(t,2)+A(t,end) - 2*A(t,1)));
    %A(t+1,end) = B(t,end) + dt/dx^2 * (Da*(A(t,end-1)+A(t,1) - 2*A(t,end)));
    A(t+1,1) = 0;
    A(t+1,end) = 0;
    
    %B
    betaP    = (Xn(2) + Xn(1))./2;
    betaBeam = (Xn(end) + Xn(1))./ 2;
    B(t+1,1)   = B(t,1) + ...
                     dt .* (Db/dx^2 .* ( B(t,end) + B(t,2) - 2.*B(t,1)) ...
                     - (betaP .* (R(t,2) - R(t,1)) ...
                     -  betaBeam .* (R(t,1)- R(t,end)))) ...
                    + gamma .* B(t,1) .* (1 - B(t,1)./Bs);
    betaN   = (Xn(end-1) + Xn(end))./2;
    B(t+1,end) = B(t,end) + ...
                     dt .* (Db/dx^2 .* ( B(t,end-1) + B(t,1) - 2.*B(t,end)) ...
                     - (betaBeam .* (R(t,1) - R(t,end)) ...
                     -  betaN .* (R(t,end)- R(t,end-1))) ...
                     + gamma .* B(t,end) .* (1 - B(t,end)./Bs));
                 
    %Repellant             
    R(t+1,1) = R(t,1) + dt*( Dr .* (R(t,2) + R(t,end)...
                 - 2*R(t,1))./(dx^2) - kr.*A(t,1));
    R(t+1,end) = R(t,end) + dt*( Dr .* (R(t,1) ...
                 + R(t,end-1) - 2*R(t,end))./(dx^2) ...
                 + kr.*A(t,end) - KAbsorbR * R(t,end));
             
    %AHL
    H(t+1,1) = H(t,1) + dt*( Dh .* (H(t,2) + H(t,end)...
                 - 2*H(t,1))./(dx^2) - kh.*A(t,1));
    H(t+1,end) = H(t,end) + dt*( Dh .* (H(t,1) ...
                 + H(t,end-1) - 2*H(t,end))./(dx^2) ...
                 + kh.*A(t,end) - KAbsorbH * H(t,end));
    

end

figure(1);clf;surf(A);shading('flat');xlabel('x');ylabel('time');title('cell A');zlabel('cell density');
figure(2);clf;surf(B);shading('flat');xlabel('x');ylabel('time');title('cell B');zlabel('cell density');
figure(3);clf;surf(R);shading('flat');xlabel('x');ylabel('time');title('repellant');zlabel('concentration');
figure(4);clf;surf(H);shading('flat');xlabel('x');ylabel('time');title('AHL');zlabel('concentration');
%figure(5);clf;surf(R-H);shading('flat');xlabel('x');ylabel('time');title('repellant-AHL');

%% make movie:
 vidObj=VideoWriter('simulation_Hours2.avi');
 set(vidObj,'FrameRate',10);
 open(vidObj);
 
 for t = 1:2:steps
     figure(1);clf;
     plot(A(t,:));xlabel('x');ylabel('#cells and concentration*10');
     hold on;
     plot(B(t,:))
     plot(R(t,:).*10);
     plot(H(t,:).*10);
     ylim([0 1300]);
     legend('A','B','rep','AHL')
     
    % Display the time.
    hours = dt*t;
    mFigure = gcf;
    % Create a uicontrol of type "text"
    mTextBox = uicontrol('style','text')
    string = [num2str(hours) ' h'];
    set(mTextBox,'string',string);
    mTextBox.Position = [85 360 80 20];
    mTextBox.BackgroundColor = [1 1 1];
    mTextBox.FontWeight = 'bold';
    writeVideo(vidObj,getframe(gcf)); 
    hours
 end
 vidObj.close();


