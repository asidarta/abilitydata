
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following code examines grasping pattern of LATERAL and BOX GRASP tasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all

% Define path and filename.
mypath = 'D:\extracted angle\';
%myfolder = {'SN011\', 'SN021\', 'SN032\', 'SN044\', 'SN050\'};
%myfile = {'grasp_R', 'lateral_R', 'keysit_R', 'head_R'};
myfolder = {'SN064\', 'SN065\', 'SN067\'};
myfile = {'lateral', 'lateral_cube',};

% Declare empty arrays. Note the shape later on, it will be 2D array.
grasp=[]; lateral=[]; keysit=[]; head=[];

for i = 1:length(myfolder)
    h = strcat(mypath,myfolder{i},myfile{1},'.txt');
    temp = readtable(h, 'HeaderLines',0);
    aperture = sqrt( (temp.Thumb_pos_X-temp.Index_pos_X).^2 + ...
                     (temp.Thumb_pos_Y-temp.Index_pos_Y).^2 + ...
                     (temp.Thumb_pos_Z-temp.Index_pos_Z).^2 );
    aperture = reshape(aperture, [200 length(aperture)/200]);
    grasp = [grasp mean(aperture,2)];
end

for i = 1:length(myfolder)
    h = strcat(mypath,myfolder{i},myfile{2},'.txt');
    temp = readtable(h, 'HeaderLines',0);
    aperture = sqrt( (temp.Thumb_pos_X-temp.Index_pos_X).^2 + ...
                     (temp.Thumb_pos_Y-temp.Index_pos_Y).^2 + ...
                     (temp.Thumb_pos_Z-temp.Index_pos_Z).^2 );
    aperture = reshape(aperture, [200 length(aperture)/200]);
    lateral = [lateral mean(aperture,2)];
end

if(0)
for i = 1:length(myfolder)
    h = strcat(mypath,myfolder{i},myfile{3},'.txt');
    temp = readtable(h, 'HeaderLines',0);
    aperture = sqrt( (temp.Thumb_pos_X-temp.Index_pos_X).^2 + ...
                     (temp.Thumb_pos_Y-temp.Index_pos_Y).^2 + ...
                     (temp.Thumb_pos_Z-temp.Index_pos_Z).^2 );
    aperture = reshape(aperture, [200 length(aperture)/200]);
    keysit = [keysit mean(aperture,2)];
end

for i = 1:length(myfolder)
    h = strcat(mypath,myfolder{i},myfile{4},'.txt');
    temp = readtable(h, 'HeaderLines',0);
    aperture = sqrt( (temp.Thumb_pos_X-temp.Index_pos_X).^2 + ...
                     (temp.Thumb_pos_Y-temp.Index_pos_Y).^2 + ...
                     (temp.Thumb_pos_Z-temp.Index_pos_Z).^2 );
    aperture = reshape(aperture, [200 length(aperture)/200]);
    head = [head mean(aperture,2)];
end
end


% Use 1-way ANOVA spm1d to compare time-series to check the main effect of
% different task. Note that we should change the format of input data!!
mymean = [grasp lateral keysit head]';
mytask = [repmat(1,[length(myfolder) 1]);
    repmat(2,[length(myfolder) 1]);
    repmat(3,[length(myfolder) 1]);
    repmat(4,[length(myfolder) 1])];
F = spm1d.stats.anova1(mymean, mytask);
Fi = F.inference(0.001)
%Fi.plot()




%--------------------------------------------------------------------------
if(0)
grasp=[]; lateral=[];
for i = 1:length(myfolder)
    h = strcat(mypath,myfolder{i},'grasp_R.txt');
    temp = readtable(h, 'HeaderLines',0);
    aperture = sqrt( (temp.Thumb_pos_X-temp.Index_pos_X).^2 + ...
                     (temp.Thumb_pos_Y-temp.Index_pos_Y).^2 + ...
                     (temp.Thumb_pos_Z-temp.Index_pos_Z).^2 );
    aperture = reshape(aperture, [200 length(aperture)/200]);
    grasp = [grasp aperture(:,5)];
end

for i = 1:length(myfolder)
    h = strcat(mypath,myfolder{i},'lateral_R.txt');
    temp = readtable(h, 'HeaderLines',0);
    aperture = sqrt( (temp.Thumb_pos_X-temp.Index_pos_X).^2 + ...
                     (temp.Thumb_pos_Y-temp.Index_pos_Y).^2 + ...
                     (temp.Thumb_pos_Z-temp.Index_pos_Z).^2 );
    aperture = reshape(aperture, [200 length(aperture)/200]);
    lateral = [lateral aperture(:,5)];
end
end

plot((lateral),'*');hold on;
plot((grasp));
plot(mean(keysit,2));
plot(mean(head,2)-5);

% Use spm1d to perform POST-HOC with Bonferroni correction.
t12        = spm1d.stats.ttest2(lateral', head');
t13        = spm1d.stats.ttest2(lateral', keysit');
t14        = spm1d.stats.ttest2(lateral', grasp');
t23        = spm1d.stats.ttest2(head', keysit');
t24        = spm1d.stats.ttest2(head', grasp');
t34        = spm1d.stats.ttest2(keysit', grasp');
% inference:
alpha      = 0.05;
nTests     = 4;
p_critical = spm1d.util.p_critical_bonf(alpha, nTests);
t12i       = t12.inference(p_critical, 'two_tailed',true);
t13i       = t13.inference(p_critical, 'two_tailed',true);
t14i       = t14.inference(p_critical, 'two_tailed',true);
t23i       = t23.inference(p_critical, 'two_tailed',true);
t24i       = t24.inference(p_critical, 'two_tailed',true);
t34i       = t34.inference(p_critical, 'two_tailed',true);

%(2) Plot:
close all
subplot(231);  t12i.plot();  ylim([-40 40]);  title('Cylinder - Head')
subplot(232);  t13i.plot();  ylim([-40 40]);  title('Cylinder - Key')
subplot(233);  t14i.plot();  ylim([-40 40]);  title('Cylinder - Cube')
subplot(234);  t23i.plot();  ylim([-40 40]);  title('Head - Key')
subplot(235);  t24i.plot();  ylim([-40 40]);  title('Head - Cube')
subplot(236);  t34i.plot();  ylim([-40 40]);  title('Key - Cube')
