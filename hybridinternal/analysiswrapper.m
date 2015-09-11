%% Written bij KU Leuven iGEM team
%% Model I
%close all;clear all;clc;
close all;clear all;

%previewOrSave='preview';
previewOrSave='save';
%previewOrSave='2D';
%previewOrSave='3D';
%previewOrSave='interaction';

%filename='real_large_domain';
%filename='toytest1';
%filename='speedupshort1';
%filename='toyfull';
%filename='revert_test';
%filename='periodic_field_test';
%filename='decoupled_timesteps_test';
filename='toy_finished_test';

list=ls([filename '*_data.mat']);
[k,~]=size(list);
filename=[filename num2str(k)];

%filename='revert_test';
%filename='10_sep_real_test';
filename='11_sep_real_test';

%framerate=25;
framerate=5;
scaling=1;

runanalysis(previewOrSave,filename,framerate,scaling);

beep on;beep;pause(1);beep;beep off;
