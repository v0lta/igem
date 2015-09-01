function runanalysis(previewOrSave,filename,framerate,scaling)

	%filename='functiontest1';
	eval([filename '_arguments']);
	eval([filename '_AHLArray']);
	eval([filename '_leucineArray']);
	eval([filename '_rhoAArray']);
	eval([filename '_rhoBArray']);
	eval([filename '_coordinateAMatrix']);
	eval([filename '_coordinateBMatrix']);

	domain.x=linspace(0,XLength,Jx);
	domain.y=linspace(0,YLength,Jy);

	[X,Y]=meshgrid(domain.x,domain.y);
	domainGrid.X=X;
	domainGrid.Y=Y;

	paramAnal.scaling=scaling;			%scaling for plotting concentrations
	paramAnal.framerate=framerate;	%framerate

	analObject=analyzer(paramAnal,AHLArray,leucineArray,rhoAArray,rhoBArray,coordinateAMatrix,coordinateBMatrix,domain,domainGrid);

	switch previewOrSave
	case 'preview'
		disp('Preview of simulation');
		analObject.preview();
	case 'save'
		disp('Saving workspace and videos');
		t2=tic;
		save(filename);
		analObject.makevideoparallel(filename);
		disp('Workspace and videos saved');
		toc(t2);
	otherwise
		warning('Unknown setting, defaulting to preview');
		disp('Preview of simulation');
		analObject.preview();
	end
end
