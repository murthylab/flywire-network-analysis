% Scatterplots, distributions using codex data format (Albert Lin)
directory = 'C:\Users\alber\Dropbox\AL Murthy Lab\FlyWire Data';
% mac
% directory = '/Users/Albert/Dropbox/AL Murthy Lab/FlyWire Data';
cd(directory)

%% load data from csvs
% neurons in neurons.csv and neuropil synapse table are in the same order,
% thresholded
data_directory = 'C:\Users\alber\Dropbox\AL Murthy Lab\FlyWire Data\Arie Data exports\630';

opts = detectImportOptions(fullfile(data_directory,'neurons.csv'));
disp([opts.VariableNames' opts.VariableTypes'])

opts = setvartype(opts,{'root_id'},'int64');

T_n = readtable(fullfile(data_directory,'neurons.csv'),opts); 

opts2 = detectImportOptions(fullfile(data_directory,'neuropil_synapse_table.csv'));
disp([opts2.VariableNames' opts2.VariableTypes'])

opts2 = setvartype(opts2,{'root_id'},'int64');

T_np = readtable(fullfile(data_directory,'neuropil_synapse_table.csv'),opts2); 

%% load sven's data sheet
data_directory_2 = 'C:\Users\alber\Dropbox\AL Murthy Lab\FlyWire Data\python_571_noKCcorr';
opts3 = detectImportOptions(fullfile(data_directory_2,'cell_stats_571.csv'));
disp([opts3.VariableNames' opts3.VariableTypes'])

opts3 = setvartype(opts3,{'root_id'},'int64');

T_stats = readtable(fullfile(data_directory_2,'cell_stats_571.csv'),opts3); 

%% load colors from spreadsheet
np_color_import = readcell('Amy neuropil color key.xlsx','Sheet','nps','Range','A2:C76');
np_ordered = np_color_import(:,1);
np_hex_colors = np_color_import(:,2);
np_names = np_color_import(:,3);

nt_hex_colors = readcell('Amy neuropil color key.xlsx','Sheet','nts','Range','B1:B6');

superclass_import = readcell('Amy neuropil color key.xlsx','Sheet','superclasses','Range','A1:B11');
superclasses = superclass_import(:,1);
superclass_hex_colors = superclass_import(:,2);

%% nps alphabetical
all_regions = {'AL_L', 'AL_R', 'AME_L', 'AME_R', 'AMMC_L', 'AMMC_R', 'AOTU_L',...
           'AOTU_R', 'ATL_L', 'ATL_R', 'AVLP_L', 'AVLP_R', 'BU_L', 'BU_R',...
           'CAN_L', 'CAN_R', 'CRE_L', 'CRE_R', 'EB', 'EPA_L', 'EPA_R', 'FB',...
           'FLA_L', 'FLA_R', 'GA_L', 'GA_R', 'GNG', 'GOR_L', 'GOR_R', 'IB_L',...
           'IB_R', 'ICL_L', 'ICL_R', 'IPS_L', 'IPS_R', 'LAL_L', 'LAL_R',...
           'LH_L', 'LH_R', 'LOP_L', 'LOP_R', 'LO_L', 'LO_R', 'MB_CA_L',...
           'MB_CA_R', 'MB_ML_L', 'MB_ML_R', 'MB_PED_L', 'MB_PED_R', 'MB_VL_L',...
           'MB_VL_R', 'ME_L', 'ME_R', 'NO', 'PB', 'PLP_L', 'PLP_R',...
           'PRW', 'PVLP_L', 'PVLP_R', 'SAD', 'SCL_L', 'SCL_R', 'SIP_L',...
           'SIP_R', 'SLP_L', 'SLP_R', 'SMP_L', 'SMP_R', 'SPS_L', 'SPS_R',...
           'VES_L', 'VES_R', 'WED_L', 'WED_R'};

all_nts = {'ach', 'gaba','glut', 'da', 'oct', 'ser'};

%%
nt_colors = hex2rgb(nt_hex_colors);
np_colors = hex2rgb(np_hex_colors);
superclass_colors = hex2rgb(superclass_hex_colors);

%% to do: extract neuron info as required
neuron_IDs = table2array(T_n(:,1));
out_nts = table2array(T_n(:,4));
n_lengths = table2array(T_n(:,16))./1e3; % units nm to um
n_surfarea = table2array(T_n(:,17))./1e6; % units nm^2 to um
n_vol = table2array(T_n(:,18))./1e9; % units nm^3 to um

%% load KC list for overwrite
% load KC file
data_directory_3 = 'C:\Users\alber\Dropbox\AL Murthy Lab\FlyWire Data\python_571';
opts4 = detectImportOptions(fullfile(data_directory_3,'KC_list_571.csv'));
disp([opts4.VariableNames' opts4.VariableTypes'])

opts4 = setvartype(opts4,{'Var2'},'int64');

T_KC = readtable(fullfile(data_directory_3,'KC_list_571.csv'),opts4); 
KC_list = table2array(T_KC(2:end,2));

%% Overwrite KCs to ach
[vals,idxs]=intersect(neuron_IDs,KC_list);

out_nts(idxs) = {'ACH'};

%% find root IDs from stats sheet all in stats order
stats_neuron_IDs = table2array(T_stats(:,2));
stats_strahler = table2array(T_stats(:,6));
stats_branchpoints = table2array(T_stats(:,7));
stats_extent = table2array(T_stats(:,10));

%% reorder
reordered_idx = zeros(size(neuron_IDs));

for i = 1:length(neuron_IDs)
    reordered_idx(i) = find(stats_neuron_IDs == neuron_IDs(i));
end

n_strahler = zeros(size(reordered_idx));
n_branchpoints = zeros(size(reordered_idx));
n_extent = zeros(size(reordered_idx));


for i = 1:length(neuron_IDs)
         n_strahler(i) = stats_strahler(reordered_idx(i));
         n_branchpoints(i) = stats_branchpoints(reordered_idx(i));
         n_extent(i) = stats_extent(reordered_idx(i))/1e3;
end


%% histogram of length by nt
hist_output_dir = 'C:\Users\alber\Dropbox\AL Murthy Lab\FlyWire Data\Neuron Feature Scatterplots';
mkdir(hist_output_dir)

for i = 1:6
histogram(log10(n_lengths(strcmpi(out_nts,all_nts{i}))),'FaceColor',nt_colors(i,:),...
    'BinWidth',0.04,'EdgeColor','w');
 xlim([1 5])
 set(gca,'Xtick',[1 2 3 4 5])
 set(gca,'XtickLabel',[10 10^2 10^3 10^4 10^5])
ylabel('number of neurons')
xlabel('neuron path length (\mum)')
title(all_nts{i})
set(gcf,'position',[1000,100,400,300])
set(gcf,'renderer','Painters')
saveas(gcf,strcat(hist_output_dir,'\length_',all_nts{i}),'epsc')
saveas(gcf,strcat(hist_output_dir,'\length_',all_nts{i}),'png')
close all
end

%% cdf of length stacked

for i = 1:6
h = cdfplot(n_lengths(strcmpi(out_nts,all_nts{i})));
set(h,'color',nt_colors(i,:),'linewidth',2);
hold on
set(gca,'xscale','log')  
xlim([10 10^5])
%set(gca,'yscale','log')
xlabel('neuron path length (\mum)')
ylabel('cdf')
end

legend(all_nts,'location','northwest')
title('')
set(gcf,'position',[1000,100,350,300])
set(gcf,'renderer','Painters')
saveas(gcf,strcat(hist_output_dir,'\length_cdf'),'epsc')
saveas(gcf,strcat(hist_output_dir,'\length_cdf'),'png')
%saveas(gcf,strcat(output_dir,'\in_deg_cdf'),'png')
% loglog logx logy

%% histogram of surfarea by nt

for i = 1:6
histogram(log10(n_surfarea(strcmpi(out_nts,all_nts{i}))),'FaceColor',nt_colors(i,:),...
    'BinWidth',0.05,'EdgeColor','w');
 xlim([1 6])
 set(gca,'Xtick',[1 2 3 4 5 6])
 set(gca,'XtickLabel',[10 10^2 10^3 10^4 10^5 10^6])
ylabel('number of neurons')
xlabel('surface area (\mum^2)')
title(all_nts{i})
set(gcf,'position',[1000,100,400,300])
set(gcf,'renderer','Painters')
saveas(gcf,strcat(hist_output_dir,'\surfarea_',all_nts{i}),'epsc')
saveas(gcf,strcat(hist_output_dir,'\surfarea_',all_nts{i}),'png')
close all
end

%% cdf of surfarea stacked
% or cdf comparison

% degree cdf
for i = 1:6
h = cdfplot(n_surfarea(strcmpi(out_nts,all_nts{i})));
set(h,'color',nt_colors(i,:),'linewidth',2);
hold on
set(gca,'xscale','log')  
xlim([10 10^6])
%set(gca,'yscale','log')
xlabel('surface area (\mum^2)')
ylabel('cdf')
end

%legend(all_nts,'location','northwest')
title('')
set(gcf,'position',[1000,100,350,300])
set(gcf,'renderer','Painters')
saveas(gcf,strcat(hist_output_dir,'\surfarea_cdf'),'epsc')
saveas(gcf,strcat(hist_output_dir,'\surfarea_cdf'),'png')


%% histogram of volume by nt

for i = 1:6
histogram(log10(n_vol(strcmpi(out_nts,all_nts{i}))),'FaceColor',nt_colors(i,:),...
    'BinWidth',0.03,'EdgeColor','w');
 xlim([1 4])
 set(gca,'Xtick',[1 2 3 4])
 set(gca,'XtickLabel',[10 10^2 10^3 10^4])
ylabel('number of neurons')
xlabel('volume (\mum^3)')
title(all_nts{i})

set(gcf,'position',[1000,100,400,300])
set(gcf,'renderer','Painters')
saveas(gcf,strcat(hist_output_dir,'\volume_',all_nts{i}),'epsc')
saveas(gcf,strcat(hist_output_dir,'\volume_',all_nts{i}),'png')
close all
end

%% cdf of volume stacked
% or cdf comparison

% degree cdf
for i = 1:6
h = cdfplot(n_vol(strcmpi(out_nts,all_nts{i})));
set(h,'color',nt_colors(i,:),'linewidth',2);
hold on
set(gca,'xscale','log')  
xlim([10 10^4])
%set(gca,'yscale','log')
xlabel('volume (\mum^3)')
ylabel('cdf')
end

%legend(all_nts,'location','northwest')
title('')
set(gcf,'position',[1000,100,350,300])
set(gcf,'renderer','Painters')
saveas(gcf,strcat(hist_output_dir,'\volume_cdf'),'epsc')
saveas(gcf,strcat(hist_output_dir,'\volume_cdf'),'png')

%% histogram of extent by nt

for i = 1:6
histogram(log10(n_extent(strcmpi(out_nts,all_nts{i}))),'FaceColor',nt_colors(i,:),...
    'BinWidth',0.02,'EdgeColor','w');
 xlim([1 3])
 set(gca,'Xtick',[1 2 3])
 set(gca,'XtickLabel',[10 10^2 10^3])
ylabel('number of neurons')
xlabel('extent (\mum)')
title(all_nts{i})

set(gcf,'position',[1000,100,400,300])
set(gcf,'renderer','Painters')
saveas(gcf,strcat(hist_output_dir,'\extent_',all_nts{i}),'epsc')
saveas(gcf,strcat(hist_output_dir,'\extent_',all_nts{i}),'png')
close all
end

%% cdf of extent stacked
% or cdf comparison

% degree cdf
for i = 1:6
h = cdfplot(n_extent(strcmpi(out_nts,all_nts{i})));
set(h,'color',nt_colors(i,:),'linewidth',2);
hold on
set(gca,'xscale','log')  
xlim([10 10^3])
%set(gca,'yscale','log')
xlabel('extent (\mum)')
ylabel('cdf')
end

%legend(all_nts,'location','northwest')
title('')
set(gcf,'position',[1000,100,350,300])
set(gcf,'renderer','Painters')
saveas(gcf,strcat(hist_output_dir,'\extent_cdf'),'epsc')
saveas(gcf,strcat(hist_output_dir,'\extent_cdf'),'png')

%% histogram of branchpoint by nt

for i = 1:6
histogram(log10(n_branchpoints(strcmpi(out_nts,all_nts{i}))),'FaceColor',nt_colors(i,:),...
    'BinWidth',0.037,'EdgeColor','w'); %about 100 bins
 xlim([0 3.4])
 set(gca,'Xtick',[0 1 2 3])
 set(gca,'XtickLabel',[1 10 10^2 10^3])
ylabel('number of neurons')
xlabel('branch point count')
title(all_nts{i})

set(gcf,'position',[1000,100,400,300])
set(gcf,'renderer','Painters')
saveas(gcf,strcat(hist_output_dir,'\branchpoint_',all_nts{i}),'epsc')
saveas(gcf,strcat(hist_output_dir,'\branchpoint_',all_nts{i}),'png')
close all
end

%% cdf of branch point stacked
% or cdf comparison

% degree cdf
for i = 1:6
h = cdfplot(n_branchpoints(strcmpi(out_nts,all_nts{i})));
set(h,'color',nt_colors(i,:),'linewidth',2);
hold on
set(gca,'xscale','log')  
xlim([1 2*10^3])
%set(gca,'yscale','log')
xlabel('branch point count')
ylabel('cdf')
end

%legend(all_nts,'location','northwest')
title('')
set(gcf,'position',[1000,100,350,300])
set(gcf,'renderer','Painters')
saveas(gcf,strcat(hist_output_dir,'\branchpoint_cdf'),'epsc')
saveas(gcf,strcat(hist_output_dir,'\branchpoint_cdf'),'png')

%% histogram of strahler by nt
% max is 9
for i = 1:6
histogram(n_strahler(strcmpi(out_nts,all_nts{i})),'FaceColor',nt_colors(i,:),'BinEdges',linspace(-0.5,9.5,11),'EdgeColor','w');
xlim([0 10])
set(gca,'Xtick',[0 1 2 3 4 5 6 7 8 9])
set(gca,'XtickLabel',[0 1 2 3 4 5 6 7 8 9])
ylabel('number of neurons')
xlabel('Strahler number')
title(all_nts{i})

set(gcf,'position',[1000,100,400,300])
set(gcf,'renderer','Painters')
saveas(gcf,strcat(hist_output_dir,'\strahler_',all_nts{i}),'epsc')
saveas(gcf,strcat(hist_output_dir,'\strahler_',all_nts{i}),'png')
close all
end

%% cdf of strahler stacked
% or cdf comparison

% degree cdf
for i = 1:6
h = cdfplot(n_strahler(strcmpi(out_nts,all_nts{i})));
set(h,'color',nt_colors(i,:),'linewidth',2);
hold on
%set(gca,'xscale','log')  
%xlim([1 2*10^3])
%set(gca,'yscale','log')
xlabel('branch point count')
ylabel('cdf')
end

%legend(all_nts,'location','northwest')
title('')
set(gcf,'position',[1000,100,350,300])
set(gcf,'renderer','Painters')
saveas(gcf,strcat(hist_output_dir,'\strahler_cdf'),'epsc')
saveas(gcf,strcat(hist_output_dir,'\strahler_cdf'),'png')
%% scatterplot by nt lengths vs area log fixed scale
for i = 1%:6

%scatter(in_deg_by_in_deg(strcmp(out_nt_in_deg,all_nts{i})),out_deg_by_in_deg(strcmp(out_nt_in_deg,all_nts{i})),'.','markeredgecolor',nt_colors(i,:))
scatterhist(log10(n_lengths(strcmpi(out_nts,all_nts{i}))),log10(n_surfarea(strcmpi(out_nts,all_nts{i}))),...
    'color',nt_colors(i,:),'Marker','.','Direction','out','NBins',100);

xlabel('in degree')
ylabel('out degree')

% xlim([-0.2 4])
% ylim([-0.2 4])
% 
% 
% set(gca,'Xtick',[0 1 2 3 4])
% set(gca,'XtickLabel',[1 10 100 1000 10000])
% set(gca,'Ytick',[0 1 2 3 4])
% set(gca,'YtickLabel',[1 10 100 1000 10000])

title(all_nts{i})
%axis square
%hold on

set(gcf,'position',[2000,100,800,800])

%set(gcf,'renderer','Painters')

%saveas(gcf,strcat(output_dir,'\all_neurons_scatter_log_',all_nts{i}),'epsc')
%saveas(gcf,strcat(output_dir,'\all_neurons_scatter_log_',all_nts{i}),'png')
%close all
end

