%% Written bij KU Leuven iGEM team
%% Model I
%close all;clear all;clc;
close all;clear;

%previewOrSave='preview';
previewOrSave='save';
%filename='real_large_domain';
%filename='toytest1';
%filename='speedupshort1';
%filename='toyfull';
filename='real_large_test';

list=ls([filename '*_data.mat']);
[k,~]=size(list);
filename=[filename num2str(k)];

filename='15_sep_real_gaussian_test';

%framerate=25;
framerate=1;
%scaling=20;
scaling=1;

runanalysis(previewOrSave,filename,framerate,scaling);
