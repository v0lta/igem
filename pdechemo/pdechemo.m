function pdechemo

m = 0;
x = linspace(0,1,100);
tend = 2;
steps = 400;
t = linspace(0,tend,steps);
% dt = (2-0)/(100-1);
% tend = 2;
% steps = tend/dt;

sol = pdepe(m,@pdechemopde,@pdechemoic,@pdechemobc,x,t);
u1 = sol(:,:,1);
u2 = sol(:,:,2);

figure
surf(x,t,u1)
shading('flat')
title('u1(x,t)')
xlabel('Distance x')
ylabel('Time t')

figure
surf(x,t,u2)
shading('flat')
title('u2(x,t)')
xlabel('Distance x')
ylabel('Time t')

% make movie:

vidObj=VideoWriter('simulation_chemo.avi');
set(vidObj,'FrameRate',1);
open(vidObj);

for t = 1:1:steps
    figure(3);clf;
    plot(u1(t,:));
    xlabel('x');
    ylabel('cells and concentration');
    hold on;
    plot(u2(t,:));
    ylim([0 4]);
    writeVideo(vidObj,getframe(gcf)); 
    t
end
vidObj.close();