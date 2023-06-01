BWD='/p01-hdd/dsb/asever/For_Aleksandar/Both_TimePoints';
addpath('/p01-hdd/dsb/asever/freesurfer/matlab/lme');
addpath('/p01-hdd/dsb/asever/freesurfer/matlab');

%% read in the relevant tables and data
% read spatially smoothed data
[Y_lhthick,mri_lhthick] = fs_read_Y('lh.thickness_sm10.mgh');
[Y_rhthick,mri_rhthick] = fs_read_Y('rh.thickness_sm10.mgh');
[Y_lharea,mri_lharea] = fs_read_Y('lh.area_sm10.mgh');
[Y_rharea,mri_rharea] = fs_read_Y('rh.area_sm10.mgh');

% read fsaverage's spherical surface and cortex label
lhsphere = fs_read_surf('fsaverage/surf/lh.sphere');
rhsphere = fs_read_surf('fsaverage/surf/rh.sphere');
lhcortex = fs_read_label('fsaverage/label/lh.cortex.label');
rhcortex = fs_read_label('fsaverage/label/rh.cortex.label');

%% construct the Design matrix based on the Qdec table

% load the Qdec table
Qdec = fReadQdec('long.qdec.2tps.manually.dat');
Qdec = rmQdecCol(Qdec,1);  % (removes first column)

sID = Qdec(2:end,1); %(grabs the subjects' IDs)
for i = 1:numel(sID);
    sID{i} = strrep(sID{i}, '_Template', '');
end;
sID = str2double(sID);

Qdec = rmQdecCol(Qdec,1);  %(removes the subjects'ID column)

% converts groups & sex to numerical 
Qdec(2:end,1) = cellstr(num2str(grp2idx(Qdec(2:end, 1))));
Qdec(2:end,4) = cellstr(num2str(grp2idx(Qdec(2:end, 4))));

% convert remaining Qdec table to matrix
M = str2double(Qdec);
M(1, :) = [];

% sorts the data and creates vector ni which represents number of 
% observations per each subject
[M,Y_lhthick,ni_lhthick] = sortData(M,1,Y_lhthick,sID);  
[M,Y_rhthick,ni_rhthick] = sortData(M,1,Y_rhthick,sID);  
[M,Y_lharea,ni_lharea] = sortData(M,1,Y_lharea,sID);  
[M,Y_rharea,ni_rharea] = sortData(M,1,Y_rharea,sID);  

% content of the design matrix                             
Intercept = ones(length(M),1);
Time = M(:,2);
Age = M(:,3);
Sex = zeros(length(M),1);
Sex(M(:,4)==2,1)=1;
ICV = M(:,5);
Group = zeros(length(M),2);
Group(M(:,1)==3,1)=1;
Group(M(:,1)==1,2)=1;
Interaction = [Time.*Group(:,1),Time.*Group(:,2)];

% final design matrix X. Reference (=0) for Sex is female. Reference for
% Group is NoTrain. Column 5 is 1 where subject is trained. Column 6 is 1
% where subject is literate. lastly the interaction with time for those two
% group columns.
X_thick = [Intercept, Time, Age, Sex, Group, Interaction];
X_area = [Intercept, Time, Age, Sex, ICV, Group, Interaction];

% the following analyses can also be done with area measures. For the
% master's thesis only thickness was analysed

%% Estimate model with one random effect (Intercept)

% lme_mass_fit_init computes initial temporal covariance component
% estimates. [1 2] indicates a random intercept & slope(years) model.
[Theta0_lhthick_1, Re_lhthick_1] = lme_mass_fit_EMinit(X_thick,[1],Y_lhthick,ni_lhthick,lhcortex);
[Theta0_rhthick_1, Re_rhthick_1] = lme_mass_fit_EMinit(X_thick,[1],Y_rhthick,ni_rhthick,rhcortex);

% lme_mass_Rgw uses the above estimates to segment the brain into
% homogeneous regions of vertices with similar covariance parameters.
[Rgs_lhthick_1,RgMeans_lhthick_1] = lme_mass_RgGrow(lhsphere,Re_lhthick_1,Theta0_lhthick_1,lhcortex);
[Rgs_rhthick_1,RgMeans_rhthick_1] = lme_mass_RgGrow(rhsphere,Re_rhthick_1,Theta0_rhthick_1,rhcortex);

% lme_mass_fit_Rgw fits the spatiotemporal model using the segmentation and
% covariance estimates.
stats_lhthick_1 = lme_mass_fit_Rgw(X_thick,[1],Y_lhthick,ni_lhthick,Theta0_lhthick_1,Rgs_lhthick_1,lhsphere);
stats_rhthick_1 = lme_mass_fit_Rgw(X_thick,[1],Y_rhthick,ni_rhthick,Theta0_rhthick_1,Rgs_rhthick_1,rhsphere);


%% Estimate model with two random effects (Intercept, Time)

% lme_mass_fit_init computes initial temporal covariance component
% estimates. [1 2] indicates a random intercept & slope(years) model.
[Theta0_lhthick, Re_lhthick] = lme_mass_fit_EMinit(X_thick,[1 2],Y_lhthick,ni_lhthick,lhcortex);
[Theta0_rhthick, Re_rhthick] = lme_mass_fit_EMinit(X_thick,[1 2],Y_rhthick,ni_rhthick,rhcortex);

% lme_mass_Rgw uses the above estimates to segment the brain into
% homogeneous regions of vertices with similar covariance parameters.
[Rgs_lhthick,RgMeans_lhthick] = lme_mass_RgGrow(lhsphere,Re_lhthick,Theta0_lhthick,lhcortex);
[Rgs_rhthick,RgMeans_rhthick] = lme_mass_RgGrow(rhsphere,Re_rhthick,Theta0_rhthick,rhcortex);


% lme_mass_fit_Rgw fits the spatiotemporal model using the segmentation and
% covariance estimates.
stats_lhthick = lme_mass_fit_Rgw(X_thick,[1 2],Y_lhthick,ni_lhthick,Theta0_lhthick,Rgs_lhthick,lhsphere);
stats_rhthick = lme_mass_fit_Rgw(X_thick,[1 2],Y_rhthick,ni_rhthick,Theta0_rhthick,Rgs_rhthick,rhsphere);

%% compare the models with one and two random effects
LR_pval_lh = lme_mass_LR(stats_lhthick,stats_lhthick_1,1);
dvtx_lh = lme_mass_FDR2(LR_pval_lh,ones(1,length(LR_pval_lh)),lhcortex,0.05,0);
LR_pval_rh = lme_mass_LR(stats_rhthick,stats_rhthick_1,1);
dvtx_rh = lme_mass_FDR2(LR_pval_rh,ones(1,length(LR_pval_rh)),rhcortex,0.05,0);
% length of dvtx indicates how many vertices survive the above correction
% only 30.76% (lh) and 29.63% (rh) survived the above correction.
% Therefore, only some vertices survived this comparison and we assumed
% that the model with one random effect on the Intercept is significantly
% better than the model with two random effects (Intercept and Slope(time).
length(lhcortex); %149955
length(rhcortex); %149926
length(dvtx_lh); %46126
length(dvtx_rh); %44430

%% the null hypothesis of no group differences in the rate of change over time
% among the three groups was tested with the following contrast
CM_thick.C = [zeros(2,6) [1 0; -1 1]]
CM_thick.C = [0 0 0 1 0 0 0 0]


%% Statistics with one random effect
F_lhstats_thick_1 = lme_mass_F(stats_lhthick_1,CM_thick)
F_rhstats_thick_1 = lme_mass_F(stats_rhthick_1,CM_thick)

%% write significance map to disc
fs_write_fstats(F_lhstats_thick_1,mri_lhthick,'sig_lhthick_sex.mgh','sig');
fs_write_fstats(F_rhstats_thick_1,mri_rhthick,'sig_rhthick_sex.mgh','sig');

%% multiple comparison correction
P = [ F_lhstats_thick_1.pval(lhcortex) F_rhstats_thick_1.pval(rhcortex) ];
G = [ F_lhstats_thick_1.sgn(lhcortex) F_rhstats_thick_1.sgn(rhcortex) ];
[detvtx_both,sided_pval_both,pth_both] = lme_mass_FDR2(P,G,[],0.05,0);
pcor = -log10(pth_both)
%% extract the vertices that have a p-value p =< 0.0001
pvalueslh = F_lhstats_thick_1.pval;
verticeslh = find(pvalueslh < 0.0001);
pvalues_threshlh = pvalueslh(verticeslh);

sign_verticeslh = [verticeslh; pvalues_threshlh];

pvaluesrh = F_rhstats_thick_1.pval;
verticesrh = find(pvaluesrh < 0.0001);
pvalues_threshrh = pvaluesrh(verticesrh);

sign_verticesrh = [verticesrh; pvalues_threshrh];

%rh: 148868 (pval: 0.0016) lh: 79752 (pval: 2.6 e⁻5), 38425 (pval: .00024), 119434 (pval
%.00013),  44237 (plval: .0002)
% first we have to work with the matrix M which has the groups saved.
% Groups: 1= literate, 2=notrain, 3=train. 1 for baseline and 2 for second
% time point
M2 = [sID, M];
M2(M2(:, 3) == 0, 3) = 1;
M2(M2(:, 3) == 0.5, 3) = 2;
M2(:,4:6) = [];

% write down the vertices for plotting in R
M3 = [M2, Y_rhthick(:,148868), Y_lhthick(:,79752), Y_lhthick(:,38425), Y_lhthick(:,119434), Y_lhthick(:,44237)]; 
csvwrite('sign_roi_vertices.csv', M3);





