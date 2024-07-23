clear all
close all
clc
%%
% Author: Yun Zhao (Monash University, Australia)
% Email: yun.zhao@monash.edu
% -------------------------------------------------------------------------
% Here we demonstrate dynamic stability analysis, and dynamic chaos
% analysis of the epileptogenic zones during seizures.
%
% We provide two 16-channel iEEG recordings in which seizures started at 60
% second and ended 10 seconds before the end of recordings. A multi-region
% model consisting of 16 neural mass models is fitted to the data.
% Parameters of the model are estimated and are used in stability and chaos
% analyzes.
%
% Showcase 1: Dynamic cortical stability
%
% This showcase calculates and shows the time-evolving Jacobi's eigenvalue 
% spectrum of the model which indicates dynamic stability of the 
% epileptogenic regions.
% 
% Showcase 2: Dynamic cortical chaos
%
% This showcase calculates and shows the time-evolving Lyapunov spectrum of
% the model which indicates the dynamic chaos of the epileptogenic regions.
% -------------------------------------------------------------------------
% 
% Please run 'parameter estimation' before running 'stability analysis' or
% 'chaos analysis'.

%% Run showcase to show dynamic cortical stability
data_file = 3; % choose from two iEEG files 1, 2 or 3

%% Need to run 'parameter estimation' first then choose from 'stability analysis' and 'chaos analysis'
tasks = ['stability analysis'];
main(data_file, tasks)
