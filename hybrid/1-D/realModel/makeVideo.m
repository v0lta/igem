vidObj=VideoWriter(filename);
set(vidObj,'FrameRate',framerate);
open(vidObj);

nFrames=model.getlength();
fig=figure('Position',[100 100 1000 1000]);
for i=1:nFrames
	model.plot(i,fig);
	ylim([-5 75]);
	writeVideo(vidObj,getframe());
	clf;
end

close(fig);
vidObj.close();
