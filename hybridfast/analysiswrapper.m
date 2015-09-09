%% Written bij KU Leuven iGEM team
%% Model I
%close all;clear all;clc;
close all;clear all;

%previewOrSave='preview';
previewOrSave='save';
%previewOrSave='2D';
%filename='real_large_domain';
%filename='toytest1';
%filename='speedupshort1';
%filename='toyfull';
%filename='revert_test';
filename='periodic_test';

%list=ls([filename '*_data.mat']);
%[k,~]=size(list);
%filename=[filename num2str(k)];

%filename='revert_test';

%framerate=25;
framerate=15;
scaling=20;

runanalysis(previewOrSave,filename,framerate,scaling);
