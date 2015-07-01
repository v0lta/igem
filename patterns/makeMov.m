
vidObj = VideoWriter('Movie.avi');
vidObj.FrameRate=30;
open(vidObj);

figure(1);
title(['Concentration of u']);
for i = 1:1:499
    str = ['u' num2str(i)];
    surface(eval(str));shading('flat');axis off;
    A(i) = getframe;
    disp(i)
end
writeVideo(vidObj,A);
close(vidObj);
%movie(A);