%% Make videos

save(filename);
disp('Saving workspace and videos');
vidObj3D=VideoWriter([filename '_3D.avi']);
set(vidObj3D,'FrameRate',framerate);
open(vidObj3D);

vidObj2D=VideoWriter([filename '_2D.avi']);
set(vidObj2D,'FrameRate',framerate);
open(vidObj2D);

nFrames=model.getlength();
fig=figure();
set(fig,'units','normalized','outerposition',[0 0 1 1]);
%fig=figure(1);
for i=1:nFrames
	%3D
	model.plot3D(i,fig);
	%ylim([-5 75]);
	writeVideo(vidObj3D,getframe(fig));
	clf;
end

for i=1:nFrames
	%2D
	model.plot2D(i,fig);
	%ylim([-5 75]);
	writeVideo(vidObj2D,getframe(fig));
	clf;
end

close(fig);
vidObj3D.close();
vidObj2D.close();

beep on;
beep;
beep off;
