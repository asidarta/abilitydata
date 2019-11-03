
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot different tasks as a function of joint-spaces. In other words, 
% how different tasks involve different variation in each of the joint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

% Define path and filename.
mypath  = 'D:\extracted angle\UL angle\';
head    = readtable(strcat(mypath,'Head.txt'), 'HeaderLines',5);
head.Color = repmat('r',[size(head,1),1]);;
back    = readtable(strcat(mypath,'back.txt'), 'HeaderLines',5);
back.Color = repmat('b',[size(back,1),1]);
grasp   = readtable(strcat(mypath,'Grasp.txt'), 'HeaderLines',5);
grasp.Color = repmat('k',[size(grasp,1),1]);
key     = readtable(strcat(mypath,'key.txt'), 'HeaderLines',5);
key.Color = repmat('y',[size(key,1),1]);
lateral = readtable(strcat(mypath,'lateral.txt'), 'HeaderLines',5);
lateral.Color = repmat('m',[size(lateral,1),1]);
mouth   = readtable(strcat(mypath,'Mouth.txt'), 'HeaderLines',5);
mouth.Color = repmat('c',[size(mouth,1),1]);
towel   = readtable(strcat(mypath,'Towel.txt'), 'HeaderLines',5);
towel.Color = repmat('g',[size(towel,1),1]);

T = [head; back; grasp; key; lateral; mouth; towel];

key.Properties.VariableNames = {'time','Cx1','Cx2','Cx3','elbow1','elbow2','elbow3',...
     'shoulder1','shoulder2','shoulder3','Tx1','Tx2','Tx3','wrist1','wrist2','wrist3','Color'};
grasp.Properties.VariableNames = {'time','Cx1','Cx2','Cx3','elbow1','elbow2','elbow3',...
     'shoulder1','shoulder2','shoulder3','Tx1','Tx2','Tx3','wrist1','wrist2','wrist3','Color'};
mouth.Properties.VariableNames = {'time','Cx1','Cx2','Cx3','elbow1','elbow2','elbow3',...
     'shoulder1','shoulder2','shoulder3','Tx1','Tx2','Tx3','wrist1','wrist2','wrist3','Color'};
head.Properties.VariableNames = {'time','Cx1','Cx2','Cx3','elbow1','elbow2','elbow3',...
     'shoulder1','shoulder2','shoulder3','Tx1','Tx2','Tx3','wrist1','wrist2','wrist3','Color'};

plot3(key.shoulder1, key.shoulder2, key.shoulder3, '.', 'color', 'r');
hold on;
plot3(grasp.shoulder1, grasp.shoulder2, grasp.shoulder3, '.', 'color', 'k');
plot3(head.shoulder1, head.shoulder2, head.shoulder3, '.', 'color', 'g');
plot3(mouth.shoulder1, mouth.shoulder2, mouth.shoulder3, '.', 'color', 'b');


