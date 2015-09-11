function outputstr=displaytime(seconds)

%hours
hours=floor(seconds/3600);
seconds=seconds-hours*3600;

%minutes
minutes=floor(seconds/60);
seconds=seconds-minutes*60;

outputstr=[num2str(hours) ' hours, ' num2str(minutes) ' minutes and ' num2str(seconds) ' seconds'];
