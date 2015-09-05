%% Written bij KU Leuven iGEM team
%% Model I
%close all;clear all;clc;
close all;clear all;

%previewOrSave='preview';
%previewOrSave='save';
previewOrSave='2D';
%filename='real_large_domain';
%filename='toytest1';
%filename='speedupshort1';
%filename='toyfull';
filename='toy_no_attr';

list=ls([filename '*_data.mat']);
[k,~]=size(list);
filename=[filename num2str(k)];

filename='toyfull5';

framerate=25;
%framerate=1;
scaling=20;

runanalysis(previewOrSave,filename,framerate,scaling);
