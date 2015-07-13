vidObj=VideoWriter('test5.avi');
set(vidObj,'FrameRate',5);
open(vidObj);

fig=figure(1);
for i=1:N+1
	model.plot(i,fig);
	ylim([0 150]);
	%ylim([0 1.2]);
	writeVideo(vidObj,getframe());
	clf;
end

close(fig);
vidObj.close();
