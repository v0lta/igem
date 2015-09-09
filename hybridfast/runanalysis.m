function runanalysis(previewOrSave,filename,framerate,scaling)

	%disp(filename);
	%eval([filename '_arguments;']);
	%eval([filename '_AHLArray;']);
	%eval([filename '_leucineArray;']);
	%eval([filename '_rhoAArray;']);
	%eval([filename '_rhoBArray;']);
	%eval([filename '_coordinateAMatrix;']);
	%eval([filename '_coordinateBMatrix;']);

	%load([filename '_data.mat']);

	%XLength=200;
	%YLength=200;
	%Jx=1001;
	%Jy=1001;

	%domain.x=linspace(0,XLength,Jx);
	%domain.y=linspace(0,YLength,Jy);

	%[X,Y]=meshgrid(domain.x,domain.y);
	%domainGrid.X=X;
	%domainGrid.Y=Y;

	paramAnal.scaling=scaling;			%scaling for plotting concentrations
	paramAnal.framerate=framerate;	%framerate
	%paramAnal.modulo=modulo;	%framerate

	%analObject=analyzer(paramAnal,AHLArray,leucineArray,rhoAArray,rhoBArray,coordinateAMatrix,coordinateBMatrix,domain,domainGrid);
	analObject=fileanalyzer(filename,paramAnal);

	switch previewOrSave
	case 'preview'
		disp('Preview of simulation');
		analObject.preview();
	case 'save'
		disp('Saving videos');
		t2=tic;
		analObject.makevideo(filename);
		disp('Videos saved');
		toc(t2);
	case '2D'
		disp('Saving 2D video');
		t2=tic;
		analObject.make2Dvideo(filename);
		disp('Video saved');
		toc(t2);
	case '3D'
		disp('Saving 3D video');
		t2=tic;
		analObject.make3Dvideo(filename);
		disp('Video saved');
		toc(t2);
	otherwise
		warning('Unknown setting, defaulting to preview');
		disp('Preview of simulation');
		analObject.preview();
	end
end
