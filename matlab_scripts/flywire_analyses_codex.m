% Using codex data format (Albert Lin)
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

%% load colors from spreadsheet
np_color_import = readcell('Amy neuropil color key.xlsx','Sheet','nps','Range','A2:C79');
np_ordered = np_color_import(:,1);
np_hex_colors = np_color_import(:,2);
np_names = np_color_import(:,3);

nt_hex_colors = readcell('Amy neuropil color key.xlsx','Sheet','nts','Range','B1:B6');

superclass_import = readcell('Amy neuropil color key.xlsx','Sheet','superclasses','Range','A2:B12');
superclasses = superclass_import(:,1);
superclass_hex_colors = superclass_import(:,2);

%% nps alphabetical
all_regions = {'AL_L', 'AL_R', 'AME_L', 'AME_R', 'AMMC_L', 'AMMC_R', 'AOTU_L',...
           'AOTU_R', 'ATL_L', 'ATL_R', 'AVLP_L', 'AVLP_R', 'BU_L', 'BU_R',...
           'CAN_L', 'CAN_R', 'CRE_L', 'CRE_R', 'EB', 'EPA_L', 'EPA_R', 'FB',...
           'FLA_L', 'FLA_R', 'GA_L', 'GA_R', 'GNG', 'GOR_L', 'GOR_R', 'IB_L',...
           'IB_R', 'ICL_L', 'ICL_R', 'IPS_L', 'IPS_R', 'LAL_L', 'LAL_R','LA_L','LA_R',...
           'LH_L', 'LH_R', 'LOP_L', 'LOP_R', 'LO_L', 'LO_R', 'MB_CA_L',...
           'MB_CA_R', 'MB_ML_L', 'MB_ML_R', 'MB_PED_L', 'MB_PED_R', 'MB_VL_L',...
           'MB_VL_R', 'ME_L', 'ME_R', 'NO','OCG', 'PB', 'PLP_L', 'PLP_R',...
           'PRW', 'PVLP_L', 'PVLP_R', 'SAD', 'SCL_L', 'SCL_R', 'SIP_L',...
           'SIP_R', 'SLP_L', 'SLP_R', 'SMP_L', 'SMP_R', 'SPS_L', 'SPS_R',...
           'VES_L', 'VES_R', 'WED_L', 'WED_R'};

all_nts = {'ach', 'gaba','glut', 'da', 'oct', 'ser'};

%% directory for general plots
alt_output_dir = fullfile(directory,'matlab_630_R');
mkdir(alt_output_dir)
%%
nt_colors = hex2rgb(nt_hex_colors);
np_colors = hex2rgb(np_hex_colors);
superclass_colors = hex2rgb(superclass_hex_colors);

%% to do: extract neuron info as required

% not all neurons in T_n are also in T_np: find the np ones in the T_n list

neuron_IDs_T_n = table2array(T_n(:,1));
out_nts_T_n = table2array(T_n(:,4));

neuron_IDs = table2array(T_np(:,1));
% extract nts for T_np neurons
clear out_nts
for i = 1:length(neuron_IDs)
    idx = find(neuron_IDs_T_n == neuron_IDs(i));
out_nts{i} = out_nts_T_n{idx};
end

%% load KC list for overwrite NOT for 630
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
%% fractional projectome
% compute fraction neurons incoming/outgoing for each neuron
% make matrix: neuron x incoming fraction, neuron x outgoing fraction
% for each neuron, multiply incoming fractions by outgoing fractions (75 x
% 75 matrix)
% add that together to generate overall matrix (add differently for
% different nts)

%% determine neuron fractions by degree
tic

in_degree = table2array(T_np(:,3));
out_degree = table2array(T_np(:,5));
in_neuron_frac = zeros(length(neuron_IDs),length(all_regions));
out_neuron_frac = zeros(length(neuron_IDs),length(all_regions));


% synapses
in_syn = table2array(T_np(:,2));
out_syn = table2array(T_np(:,4));
in_syn_frac = zeros(length(neuron_IDs),length(all_regions));
out_syn_frac = zeros(length(neuron_IDs),length(all_regions));


% extracting in alphabetical order first since that follows the file format
% later, swap matrices into np_ordered order

for i = 1:length(neuron_IDs)
        input_neuron_num = table2array(T_np(i,84:161)); 
        output_neuron_num = table2array(T_np(i,240:317)); 
        input_syn_num = table2array(T_np(i,6:83)); 
        output_syn_num = table2array(T_np(i,162:239));

        in_neuron_frac(i,:) = input_neuron_num./in_degree(i);
        out_neuron_frac(i,:) = output_neuron_num./out_degree(i);
        in_syn_frac(i,:) = input_syn_num./in_syn(i);
        out_syn_frac(i,:) = output_syn_num./out_syn(i);
end

toc
%% fraction input and output for each neuron by nt (sum up the nts to get full projectome)

output_dir = fullfile(directory,strcat('Fractional_Projectome_t5'));
mkdir(output_dir)

%%
clear project_to_matrix_nt project_to_matrix_synweights_nt project_to_matrix_synnum_nt
project_to_matrix_nt = zeros(length(all_regions),length(all_regions),6);

project_to_matrix_synweights_nt = zeros(length(all_regions),length(all_regions),6);

project_to_matrix_synnum_nt = zeros(length(all_regions),length(all_regions),6);

debug_AL_counter = 0;

tic

for i = 1:length(neuron_IDs) % neurons
    clear temp_neuron_matrix temp_syn_matrix temp_synnum_matrix k
temp_neuron_matrix = in_neuron_frac(i,:)' *out_neuron_frac(i,:); % input on x axis, output on y
temp_neuron_matrix(isnan(temp_neuron_matrix)) = 0;

temp_syn_matrix = in_syn_frac(i,:)' *out_syn_frac(i,:); % input on x axis, output on y
temp_syn_matrix(isnan(temp_syn_matrix)) = 0;

temp_synnum_matrix = (in_syn_frac(i,:)' *out_syn_frac(i,:)).*(in_syn(i)+out_syn(i)); % input on x axis, output on y, multiplied by total synapse number
temp_synnum_matrix(isnan(temp_synnum_matrix)) = 0;

if temp_syn_matrix(1,1)>0
debug_AL_counter = debug_AL_counter+1;
end

k = find(strcmpi(all_nts,out_nts{i})); % which nt is this neuron

project_to_matrix_nt(:,:,k) = project_to_matrix_nt(:,:,k) + temp_neuron_matrix;
project_to_matrix_synweights_nt(:,:,k) = project_to_matrix_synweights_nt(:,:,k) + temp_syn_matrix;

project_to_matrix_synnum_nt(:,:,k) = project_to_matrix_synnum_nt(:,:,k) + temp_synnum_matrix;

end

toc
%% reorder
reordered_idx = zeros(size(all_regions));

for i = 1:length(all_regions)
    reordered_idx(i) = find(strcmp(np_ordered,all_regions{i}));
end

ordered_project_to_matrix_nt = zeros(size(project_to_matrix_nt));
ordered_project_to_matrix_syn_nt = zeros(size(project_to_matrix_nt));
ordered_project_to_matrix_synnum_nt = zeros(size(project_to_matrix_nt));

for i = 1:length(np_ordered) % row
    for j = 1:length(np_ordered) % column
        ordered_project_to_matrix_nt(reordered_idx(i),reordered_idx(j),:) = project_to_matrix_nt(i,j,:);
        ordered_project_to_matrix_syn_nt(reordered_idx(i),reordered_idx(j),:) = project_to_matrix_synweights_nt(i,j,:);
        ordered_project_to_matrix_synnum_nt(reordered_idx(i),reordered_idx(j),:) = project_to_matrix_synnum_nt(i,j,:);
    end
end

%% print all neurons projecting to each neuropil even in part SKIP THESE
% for winner-take-all use flywire_analyses_nt code
% projecting to means outputs point to this neuropil

np_csv_folder = fullfile(directory,'neuropil_fractional_csvs_571\');
mkdir(np_csv_folder)

project_to_neuron_num = zeros(size(all_regions));

% neuron fracs are unordered here
for i = 1:length(all_regions)
clear out_neurons
neuropil = all_regions{i}; 

idx = find(out_syn_frac(:,i)>0); % >0 to exclude NaN

out_neurons = neuron_IDs(idx);

project_to_neuron_num(i) = length(idx);

if ~isempty(out_neurons)
writematrix(out_neurons',strcat(np_csv_folder,neuropil,'_project_to_neurons.csv'))
end

end

%% print all neurons projecting to each neuropil even in part
% for winner-take-all use flywire_analyses_nt code
% projecting to means outputs point to this neuropil

project_from_neuron_num = zeros(size(all_regions));

% neuron fracs are unordered here
for i = 1:length(all_regions)
clear out_neurons
neuropil = all_regions{i}; 

idx = find(in_syn_frac(:,i)>0);

in_neurons = neuron_IDs(idx);

project_from_neuron_num(i) = length(idx);

if ~isempty(in_neurons)
writematrix(in_neurons',strcat(np_csv_folder,neuropil,'_project_from_neurons.csv'))
end

end


%% stacked bar plots input

output_dir_np_proj = fullfile(directory,strcat('Fractional_NP_Projections_t5'));
mkdir(output_dir_np_proj)

for i = 1:length(np_ordered)
for k = 1:length(all_nts)
    
    neuropil = np_ordered{i};
    neuropil_print = np_names{i};
nt = all_nts{k};

bar(ordered_project_to_matrix_syn_nt(:,i,k),'facecolor',nt_colors(k,:),'edgecolor','none')
set(gca,'xtick',1:78,'xticklabel',np_names)
 xtickangle(90)
set(gcf,'position',[100,100,1800,300])
set(gca,'FontSize',15)
ylabel('total frac. weights')
xlabel('input neuropils of neurons')

title(strcat(nt,{' '}, 'projections to',{' '}, neuropil_print))

saveas(gcf,fullfile(output_dir_np_proj,strcat(neuropil,'_',nt,'_input_linear.png')));
%close all

close all
%

%[neuropils,~,idx] = unique(out_nps);
% hist(idx,unique(idx))
% set(gca,'xtick',1:length(neuropils),'xticklabel',neuropils)
% xtickangle(90)
 
end

% make combo plot stacked bar plot
h = bar(reshape(ordered_project_to_matrix_syn_nt(:,i,:),78,6),'stacked','edgecolor','none');%,'facecolor',color_matrix)
for k = 1:6
h(k).FaceColor = nt_colors(k,:);
end

ax = gca;
set(ax,'xtick',1:78,'xticklabel',np_names)
 xtickangle(90)
for m = 1:length(np_names)
    ax.XTickLabel{m} = [strcat('\color[rgb]{',num2str(np_colors(m,:)),'} '),ax.XTickLabel{m}];
end
set(gcf,'position',[100,100,1800,300])
set(gca,'FontSize',15)
ylabel('total frac. weights')
xlabel('input neuropils of neurons')

title(strcat('projections to',{' '}, neuropil_print))
legend(all_nts,'Location','eastoutside')

saveas(gcf,strcat(output_dir_np_proj,'\',neuropil,'_all_nts_input_linear'),'epsc')
saveas(gcf,strcat(output_dir_np_proj,'\',neuropil,'_all_nts_input_linear'),'png')

%saveas(gcf,fullfile(output_dir,strcat(neuropil,'_all_nts_max_in_linear.png')));
close all

end

%% stacked bar plots output

output_dir_np_proj = fullfile(directory,strcat('Fractional_NP_Projections_t5'));
mkdir(output_dir_np_proj)

for i = 1:length(np_ordered)
for k = 1:length(all_nts)
    
    neuropil = np_ordered{i};
    neuropil_print = np_names{i};
nt = all_nts{k};

bar(ordered_project_to_matrix_syn_nt(i,:,k),'facecolor',nt_colors(k,:),'edgecolor','none')
set(gca,'xtick',1:78,'xticklabel',np_names)
 xtickangle(90)
set(gcf,'position',[100,100,1800,300])
set(gca,'FontSize',15)
ylabel('total frac. weights')
xlabel('output neuropils of neurons')

title(strcat(nt,{' '}, 'projections from',{' '}, neuropil_print))

saveas(gcf,fullfile(output_dir_np_proj,strcat(neuropil,'_',nt,'_output_linear.png')));
%close all

close all
%

%[neuropils,~,idx] = unique(out_nps);
% hist(idx,unique(idx))
% set(gca,'xtick',1:length(neuropils),'xticklabel',neuropils)
% xtickangle(90)
 
end

% make combo plot stacked bar plot
h = bar(reshape(ordered_project_to_matrix_syn_nt(i,:,:),78,6),'stacked','edgecolor','none');%,'facecolor',color_matrix)
for k = 1:6
h(k).FaceColor = nt_colors(k,:);
end

ax = gca;
set(ax,'xtick',1:78,'xticklabel',np_names)
 xtickangle(90)
for m = 1:length(np_names)
    ax.XTickLabel{m} = [strcat('\color[rgb]{',num2str(np_colors(m,:)),'} '),ax.XTickLabel{m}];
end

set(gcf,'position',[100,100,1800,300])
set(gca,'FontSize',15)
ylabel('total frac. weights')
xlabel('output neuropils of neurons')

title(strcat('projections from',{' '}, neuropil_print))
legend(all_nts,'Location','eastoutside')

saveas(gcf,strcat(output_dir_np_proj,'\',neuropil,'_all_nts_output_linear'),'epsc')
saveas(gcf,strcat(output_dir_np_proj,'\',neuropil,'_all_nts_output_linear'),'png')

%saveas(gcf,fullfile(output_dir,strcat(neuropil,'_all_nts_max_in_linear.png')));
close all

end


%% compute number of internal, external incoming, external outgoing fractions

ordered_project_to_matrix_syn = sum(ordered_project_to_matrix_syn_nt,3);

% outgoing neurons are neurons for which maxin is this neuropil
internal_num = zeros(size(np_ordered));
external_outgoing_num = zeros(size(np_ordered));
external_incoming_num = zeros(size(np_ordered));

for i=1:length(internal_num)
    internal_num(i) = ordered_project_to_matrix_syn(i,i);
    external_outgoing_num(i) = sum(ordered_project_to_matrix_syn(i,:))-internal_num(i);
    external_incoming_num(i) = sum(ordered_project_to_matrix_syn(:,i))-internal_num(i);
end

num_for_hist = [internal_num external_outgoing_num external_incoming_num];


%totals = sum(num_for_hist,2);
%frac_for_hist = num_for_hist./totals;

fig = figure(1);
b = bar(num_for_hist,'edgecolor','none'); %,'stacked'
b(1).FaceColor = [0.9 0.8 0.5];
b(2).FaceColor = [0.5 0.6 0.8];
b(3).FaceColor = [0.75 0.4 0.5];

set(gca,'xtick',1:78,'xticklabel',np_names)
 xtickangle(90)
 
ylabel('total fractional weights')
legend('internal','external outgoing','external incoming','Location','southoutside')
set(gcf,'position',[100,100,1500,700])
set(gca,'FontSize',15)

saveas(fig,strcat(alt_output_dir,'\internal_num'),'epsc')
saveas(fig,strcat(alt_output_dir,'\internal_num'),'png')

%% compute fraction of internal, external incoming, external outgoing fractions

totals = sum(num_for_hist,2);
frac_for_hist = num_for_hist./totals;

fig = figure(1);
b = bar(frac_for_hist,'stacked','edgecolor','none'); %,'stacked'
b(1).FaceColor = [0.9 0.8 0.5];
b(2).FaceColor = [0.5 0.6 0.8];
b(3).FaceColor = [0.75 0.4 0.5];

ax = gca;
set(ax,'xtick',1:78,'xticklabel',np_names)
 xtickangle(90)
for m = 1:length(np_names)
    ax.XTickLabel{m} = [strcat('\color[rgb]{',num2str(np_colors(m,:)),'} '),ax.XTickLabel{m}];
end

ylim([0 1])
 
ylabel('fraction of total weights')
legend('internal','external outgoing','external incoming','Location','southoutside')
set(gcf,'position',[100,100,1500,700])
set(gca,'FontSize',15)

saveas(fig,strcat(alt_output_dir,'\internal_frac'),'epsc')
saveas(fig,strcat(alt_output_dir,'\internal_frac'),'png')
%% split up internal by nt

% for each 
internal_num_nt = zeros(78,6);
external_outgoing_num_nt = zeros(78,6);
external_incoming_num_nt = zeros(78,6);

for i=1:length(internal_num)
    internal_num_nt(i,:) = ordered_project_to_matrix_syn_nt(i,i,:);
    external_outgoing_num_nt(i,:) = reshape(sum(ordered_project_to_matrix_syn_nt(i,:,:)),[1 6])-internal_num_nt(i,:);
    external_incoming_num_nt(i,:) = reshape(sum(ordered_project_to_matrix_syn_nt(:,i,:)),[1 6])-internal_num_nt(i,:);
end

% internal hist
internal_totals = sum(internal_num_nt,2);
internal_frac_for_hist = internal_num_nt./internal_totals;

% make combo plot stacked bar plot
h = bar(internal_frac_for_hist,'stacked','edgecolor','none');%,'facecolor',color_matrix)
for k = 1:6
h(k).FaceColor = nt_colors(k,:);
end

ax = gca;
set(ax,'xtick',1:78,'xticklabel',np_names)
 xtickangle(90)
for m = 1:length(np_names)
    ax.XTickLabel{m} = [strcat('\color[rgb]{',num2str(np_colors(m,:)),'} '),ax.XTickLabel{m}];
end

set(gcf,'position',[100,100,1500,300])
set(gca,'FontSize',15)
ylabel('frac. of weights')
xlabel('neuropils')

title('fraction of internal weights by nt')
legend(all_nts,'Location','eastoutside')

saveas(gcf,strcat(alt_output_dir,'\','nt_frac_internal'),'epsc')
saveas(gcf,strcat(alt_output_dir,'\','nt_frac_internal'),'png')


%% split up ext incoming by nt

% ext hist
ext_inc_totals = sum(external_incoming_num_nt,2);
ext_inc_frac_for_hist = external_incoming_num_nt./ext_inc_totals;

% make combo plot stacked bar plot
h = bar(ext_inc_frac_for_hist,'stacked','edgecolor','none');%,'facecolor',color_matrix)
for k = 1:6
h(k).FaceColor = nt_colors(k,:);
end

ax = gca;
set(ax,'xtick',1:78,'xticklabel',np_names)
 xtickangle(90)
for m = 1:length(np_names)
    ax.XTickLabel{m} = [strcat('\color[rgb]{',num2str(np_colors(m,:)),'} '),ax.XTickLabel{m}];
end

set(gcf,'position',[100,100,1500,300])
set(gca,'FontSize',15)
ylabel('frac. of weights')
xlabel('neuropils')

title('fraction of external incoming weights by nt')
legend(all_nts,'Location','eastoutside')

saveas(gcf,strcat(alt_output_dir,'\','nt_frac_ext_inc'),'epsc')
saveas(gcf,strcat(alt_output_dir,'\','nt_frac_ext_inc'),'png')

%% split up ext outgoing by nt


% internal hist
ext_out_totals = sum(external_outgoing_num_nt,2);
ext_out_frac_for_hist = external_outgoing_num_nt./ext_out_totals;

% make combo plot stacked bar plot
h = bar(ext_out_frac_for_hist,'stacked','edgecolor','none');%,'facecolor',color_matrix)
for k = 1:6
h(k).FaceColor = nt_colors(k,:);
end

ax = gca;
set(ax,'xtick',1:78,'xticklabel',np_names)
 xtickangle(90)
for m = 1:length(np_names)
    ax.XTickLabel{m} = [strcat('\color[rgb]{',num2str(np_colors(m,:)),'} '),ax.XTickLabel{m}];
end

set(gcf,'position',[100,100,1500,300])
set(gca,'FontSize',15)
ylabel('frac. of weights')
xlabel('neuropils')

title('fraction of external outgoing weights by nt')
legend(all_nts,'Location','eastoutside')

saveas(gcf,strcat(alt_output_dir,'\','nt_frac_ext_out'),'epsc')
saveas(gcf,strcat(alt_output_dir,'\','nt_frac_ext_out'),'png')

%% sum over whole brain: numbers of int/ext outgoing/ext incoming

wholebrain_int_nt = sum(internal_num_nt,1);
wholebrain_ext_in_nt = sum(external_incoming_num_nt,1);
wholebrain_ext_out_nt = sum(external_outgoing_num_nt,1);

wholebrain_num_for_hist = [wholebrain_int_nt; wholebrain_ext_in_nt];

h = bar(wholebrain_num_for_hist,'stacked','edgecolor','none');%,'facecolor',color_matrix)
for k = 1:6
h(k).FaceColor = nt_colors(k,:);
end
ax = gca;
set(ax,'xtick',1:2,'xticklabel',{'internal','external'})
 xtickangle(90)


set(gcf,'position',[100,100,400,600])
set(gca,'FontSize',20)
ylabel('total fractional weights')

%title('fraction of external outgoing weights by nt')
legend(all_nts,'Location','eastoutside')

saveas(gcf,strcat(alt_output_dir,'\','nt_frac_int_ext_wholebrain'),'epsc')
saveas(gcf,strcat(alt_output_dir,'\','nt_frac_int_ext_wholebrain'),'png')



%% project matrices by nt

for k = 1:length(all_nts)
nt = all_nts{k};
figure(1)
imagesc(ordered_project_to_matrix_nt(:,:,k))
axis square
set(gca,'ytick',1:78,'yticklabel',np_names)
set(gca,'xtick',1:78,'xticklabel',np_names)
 xtickangle(90)
 ylabel('input neuropils','fontsize',15)
 xlabel('output neuropils','fontsize',15)
 title(nt,'fontsize',15)

 % construct custom colormap
 max_color = nt_colors(k,:);
 map = ones(256,3);
 map(1,:)=[0.95 0.95 0.95];
 for i=2:256
    map(i,:)= map(1,:) - (map(1,:) - max_color).*i/256;
 end

colormap(map);
%colormap(colorcet('L12')); % set the colorscheme, diverging
c = colorbar;
c.Label.String='total fractional weights';
c.Label.FontSize = 15;
set(gcf, 'Position',  [100, 0, 1000, 1000])
%set(gca,'ColorScale','log')

saveas(gcf,strcat(output_dir,'\',nt,'_lin_proj_unclustered'),'epsc')
saveas(gcf,strcat(output_dir,'\',nt,'_lin_proj_unclustered'),'png')

close all

end

%% project matrices by nt log
%log_ordered_project_to_matrix_nt = log10(ordered_project_to_matrix_nt);

for k = 1:length(all_nts)
nt = all_nts{k};
figure(1)
imagesc(ordered_project_to_matrix_nt(:,:,k))
axis square
set(gca,'ytick',1:78,'yticklabel',np_names)
set(gca,'xtick',1:78,'xticklabel',np_names)
 xtickangle(90)
 ylabel('input neuropils','fontsize',15)
 xlabel('output neuropils','fontsize',15)
 title(nt,'fontsize',15)

 % construct custom colormap
 max_color = nt_colors(k,:);
 map = ones(256,3);
 map(1,:)=[0.95 0.95 0.95];
 for i=2:256
    map(i,:)= map(1,:) - (map(1,:) - max_color).*i/256;
 end

colormap(map);
%colormap(colorcet('L12')); % set the colorscheme, diverging
c = colorbar;
c.Label.String='total fractional weights';
c.Label.FontSize = 15;
set(gca,'ColorScale','log')
set(gcf, 'Position',  [100, 0, 1000, 1000])
%set(gca,'ColorScale','log')

saveas(gcf,strcat(output_dir,'\',nt,'_log_proj_unclustered'),'epsc')
saveas(gcf,strcat(output_dir,'\',nt,'_log_proj_unclustered'),'png')

close all

end

%% total projectome neuron count
ordered_project_to_matrix = sum(ordered_project_to_matrix_nt,3);

figure(2)
imagesc(ordered_project_to_matrix)
axis square
set(gca,'ytick',1:78,'yticklabel',np_names)
set(gca,'xtick',1:78,'xticklabel',np_names)
 xtickangle(90)
 ylabel('input neuropils','fontsize',15)
 xlabel('output neuropils','fontsize',15)
colormap(colorcet('L1','reverse',true)); % set the colorscheme, diverging
c = colorbar;
c.Label.String='total fractional weights';
c.Label.FontSize = 15;
set(gcf, 'Position',  [100, 0, 1000, 1000])
set(gca,'ColorScale','log')

saveas(gcf,strcat(output_dir,'\','proj_unclustered'),'epsc')
saveas(gcf,strcat(output_dir,'\','proj_unclustered'),'png')

%% project matrices by nt

for k = 1:length(all_nts)
nt = all_nts{k};
figure(1)
imagesc(ordered_project_to_matrix_nt(:,:,k))
axis square
set(gca,'ytick',1:78,'yticklabel',np_names)
set(gca,'xtick',1:78,'xticklabel',np_names)
 xtickangle(90)
 ylabel('input neuropils','fontsize',15)
 xlabel('output neuropils','fontsize',15)
 title(nt,'fontsize',15)

 % construct custom colormap
 max_color = nt_colors(k,:);
 map = ones(256,3);
 map(1,:)=[0.95 0.95 0.95];
 for i=2:256
    map(i,:)= map(1,:) - (map(1,:) - max_color).*i/256;
 end

colormap(map);
%colormap(colorcet('L12')); % set the colorscheme, diverging
c = colorbar;
c.Label.String='total fractional weights';
c.Label.FontSize = 15;
set(gcf, 'Position',  [100, 0, 1000, 1000])
%set(gca,'ColorScale','log')

saveas(gcf,strcat(output_dir,'\',nt,'_lin_proj_unclustered'),'epsc')
saveas(gcf,strcat(output_dir,'\',nt,'_lin_proj_unclustered'),'png')

close all

end

%% project matrices by nt log synapses
%log_ordered_project_to_matrix_nt = log10(ordered_project_to_matrix_nt);

for k = 1:length(all_nts)
nt = all_nts{k};
figure(1)
imagesc(ordered_project_to_matrix_syn_nt(:,:,k))
axis square
set(gca,'ytick',1:75,'yticklabel',np_names)
set(gca,'xtick',1:75,'xticklabel',np_names)
 xtickangle(90)

ax = gca;
for m = 1:length(np_names)
    ax.XTickLabel{m} = [strcat('\color[rgb]{',num2str(np_colors(m,:)),'} '),ax.XTickLabel{m}];
    ax.YTickLabel{m} = [strcat('\color[rgb]{',num2str(np_colors(m,:)),'} '),ax.YTickLabel{m}];
end

 ylabel('input neuropils','fontsize',15)
 xlabel('output neuropils','fontsize',15)
 title(nt,'fontsize',15)

 % construct custom colormap
 max_color = nt_colors(k,:);
 map = ones(256,3);
 map(1,:)=[0.95 0.95 0.95];
 for i=2:256
    map(i,:)= map(1,:) - (map(1,:) - max_color).*i/256;
 end

colormap(map);
%colormap(colorcet('L12')); % set the colorscheme, diverging
c = colorbar;
c.Label.String='total fractional weights';
c.Label.FontSize = 15;
set(gca,'ColorScale','log')
set(gcf, 'Position',  [100, 0, 1000, 1000])
%set(gca,'ColorScale','log')

saveas(gcf,strcat(output_dir,'\',nt,'_synweights_log_proj_unclustered'),'epsc')
saveas(gcf,strcat(output_dir,'\',nt,'_synweights_log_proj_unclustered'),'png')

close all

end

%% ach gaba glut combining colors WIP

rgb_matrix = zeros(size(ordered_project_to_matrix_syn_nt(:,:,1:3)));
% rgb
combo_colors = [0 0 1; 1 0 0; 0 1 0];

%combo_colors = nt_colors(1:3,:);

 maps = ones(256,3,3);
% save color maps
for k=1:3
 max_color = combo_colors(k,:);
 maps(1,:,k)=[0.95 0.95 0.95];
 for i=2:256
    maps(i,:,k)= maps(1,:,k) - (maps(1,:,k) - max_color).*i/256;
 end

end

for i=1:length(np_ordered)
    for j=1:length(np_ordered)
        colors_temp = zeros(3,3); % each row is a nt's rgb values
        for k = 1:3
            % normalize by the element value, then get its color, log scale
            el_val = log10(ordered_project_to_matrix_syn_nt(i,j,k)/max(max(ordered_project_to_matrix_syn_nt(:,:,k))));
            color_range_min = -5;
            color_range_max = 1; % logscale color bar limits manually set
            idx_temp = round((el_val-color_range_min)*256/(color_range_max-color_range_min));
            color_idx = max(idx_temp,1);
            colors_temp(k,:) = maps(color_idx,:,k);
        end
        % blend colors avg
        %rgb_matrix(i,j,:) = 1/3*(colors_temp(1,:)+colors_temp(2,:)+colors_temp(3,:));
        % blend colors quadratic
        rgb_matrix(i,j,:) = 1/3*(colors_temp(1,:).^2+colors_temp(2,:).^2+colors_temp(3,:).^2);
    end
end

imagesc(rgb_matrix)
axis square
set(gca,'ytick',1:75,'yticklabel',np_names)
set(gca,'xtick',1:75,'xticklabel',np_names)
 xtickangle(90)
title('RGB composite of ach (blue), gaba (red), glut (green)','fontsize',15)

ax = gca;
for m = 1:length(np_names)
    ax.XTickLabel{m} = [strcat('\color[rgb]{',num2str(np_colors(m,:)),'} '),ax.XTickLabel{m}];
    ax.YTickLabel{m} = [strcat('\color[rgb]{',num2str(np_colors(m,:)),'} '),ax.YTickLabel{m}];
end

set(gcf, 'Position',  [100, 0, 1000, 1000])

% could generate legend by grabbing log-scale bars off of individual
% matrices plotted with combo colors

saveas(gcf,strcat(output_dir,'\','rgb_white_synweights_log_proj_unclustered'),'epsc')
saveas(gcf,strcat(output_dir,'\','rgb_white_synweights_log_proj_unclustered'),'png')

%% ach gaba glut "rgb"
imagesc(ordered_project_to_matrix_syn_nt(:,:,1:3))
axis square
set(gca,'ytick',1:75,'yticklabel',np_names)
set(gca,'xtick',1:75,'xticklabel',np_names)
 xtickangle(90)
 ylabel('input neuropils','fontsize',15)
 xlabel('output neuropils','fontsize',15)
 title('R=ach, G=gaba, B=glut','fontsize',15)

set(gca,'ColorScale','log')
set(gcf, 'Position',  [100, 0, 1000, 1000])
%set(gca,'ColorScale','log')

ax = gca;
for m = 1:length(np_names)
    ax.XTickLabel{m} = [strcat('\color[rgb]{',num2str(np_colors(m,:)),'} '),ax.XTickLabel{m}];
    ax.YTickLabel{m} = [strcat('\color[rgb]{',num2str(np_colors(m,:)),'} '),ax.YTickLabel{m}];
end

saveas(gcf,strcat(output_dir,'\','rgb_synweights_log_proj_unclustered'),'epsc')
saveas(gcf,strcat(output_dir,'\','rgb_synweights_log_proj_unclustered'),'png')



%% total projectome neuron count synapses
ordered_project_to_matrix_syn = sum(ordered_project_to_matrix_syn_nt,3);

figure(1)
imagesc(ordered_project_to_matrix_syn)
axis square
set(gca,'ytick',1:75,'yticklabel',np_names)
set(gca,'xtick',1:75,'xticklabel',np_names)
 xtickangle(90)
 ylabel('input neuropils','fontsize',15)
 xlabel('output neuropils','fontsize',15)
colormap(colorcet('L1','reverse',true)); % set the colorscheme, diverging
c = colorbar;
c.Label.String='total fractional weights';
c.Label.FontSize = 15;
set(gcf, 'Position',  [100, 0, 1000, 1000])
set(gca,'ColorScale','log')

ax = gca;
for m = 1:length(np_names)
    ax.XTickLabel{m} = [strcat('\color[rgb]{',num2str(np_colors(m,:)),'} '),ax.XTickLabel{m}];
    ax.YTickLabel{m} = [strcat('\color[rgb]{',num2str(np_colors(m,:)),'} '),ax.YTickLabel{m}];
end

saveas(gcf,strcat(output_dir,'\','synweights_proj_unclustered'),'epsc')
saveas(gcf,strcat(output_dir,'\','synweights_proj_unclustered'),'png')

%% project matrices by nt log synapse numbers
%log_ordered_project_to_matrix_nt = log10(ordered_project_to_matrix_nt);

for k = 1:length(all_nts)
nt = all_nts{k};
figure(1)
imagesc(ordered_project_to_matrix_synnum_nt(:,:,k))
axis square
set(gca,'ytick',1:75,'yticklabel',np_ordered)
set(gca,'xtick',1:75,'xticklabel',np_ordered)
 xtickangle(90)
 ylabel('input neuropils','fontsize',15)
 xlabel('output neuropils','fontsize',15)
 title(nt,'fontsize',15)

 % construct custom colormap
 max_color = nt_colors(k,:);
 map = ones(256,3);
 map(1,:)=[0.95 0.95 0.95];
 for i=2:256
    map(i,:)= map(1,:) - (map(1,:) - max_color).*i/256;
 end

colormap(map);
%colormap(colorcet('L12')); % set the colorscheme, diverging
c = colorbar;
c.Label.String='total fractional synapses';
c.Label.FontSize = 15;
set(gca,'ColorScale','log')
set(gcf, 'Position',  [100, 0, 1000, 1000])
%set(gca,'ColorScale','log')

saveas(gcf,strcat(output_dir,'\',nt,'_synnum_log_proj_unclustered'),'epsc')
saveas(gcf,strcat(output_dir,'\',nt,'_synnum_log_proj_unclustered'),'png')

close all

end

%% total projectome neuron count synapse numbers
ordered_project_to_matrix_synnum = sum(ordered_project_to_matrix_synnum_nt,3);

figure(1)
imagesc(ordered_project_to_matrix_synnum)
axis square
set(gca,'ytick',1:75,'yticklabel',np_ordered)
set(gca,'xtick',1:75,'xticklabel',np_ordered)
 xtickangle(90)
 ylabel('input neuropils','fontsize',15)
 xlabel('output neuropils','fontsize',15)
colormap(colorcet('L1','reverse',true)); % set the colorscheme, diverging
c = colorbar;
c.Label.String='total fractional synapses';
c.Label.FontSize = 15;
set(gcf, 'Position',  [100, 0, 1000, 1000])
set(gca,'ColorScale','log')

saveas(gcf,strcat(output_dir,'\','synnum_proj_unclustered'),'epsc')
saveas(gcf,strcat(output_dir,'\','synnum_proj_unclustered'),'png')


%% ach - gaba normalized

achgaba_diff_norm = ordered_project_to_matrix_syn_nt(:,:,1)./sum(sum(ordered_project_to_matrix_syn_nt(:,:,1)))...
    -ordered_project_to_matrix_syn_nt(:,:,2)./sum(sum(ordered_project_to_matrix_syn_nt(:,:,2)));

achgaba_diff = ordered_project_to_matrix_syn_nt(:,:,1)-ordered_project_to_matrix_syn_nt(:,:,2);

log_achgaba_diff_norm = zeros(size(achgaba_diff));

%log pulls negatives
for i = 1:length(np_ordered)
    for j = 1:length(np_ordered)
        if achgaba_diff_norm(i,j) ~= 0 && ~isnan(achgaba_diff_norm(i,j))
log_achgaba_diff_norm(i,j) = log(abs(achgaba_diff_norm(i,j)))*sign(achgaba_diff_norm(i,j));
        end
        end
end

nt = all_nts{k};
figure(1)
imagesc(log_achgaba_diff_norm)
axis square
set(gca,'ytick',1:75,'yticklabel',np_names)
set(gca,'xtick',1:75,'xticklabel',np_names)
 xtickangle(90)
 ylabel('input neuropils','fontsize',15)
 xlabel('output neuropils','fontsize',15)
 title('ach - gaba normalized','fontsize',15)

 ax = gca;
for m = 1:length(np_names)
    ax.XTickLabel{m} = [strcat('\color[rgb]{',num2str(np_colors(m,:)),'} '),ax.XTickLabel{m}];
    ax.YTickLabel{m} = [strcat('\color[rgb]{',num2str(np_colors(m,:)),'} '),ax.YTickLabel{m}];
end

 % construct custom colormap
 max_color = nt_colors(1,:);
 map = ones(256,3);
map(128,:)=[0.95 0.95 0.95];
min_color =nt_colors(2,:);
 for i=1:127
    map(i,:)= min_color - (min_color - map(128,:)).*i/128;
 end
 for i=129:256
    map(i,:)= map(128,:) - (map(128,:) - max_color).*(i-128)/128;
 end
colormap(map);

%colormap(colorcet('D1','reverse',1))
caxis([-30 30])

c = colorbar;
c.Label.String='signed log(norm. difference)';
c.Label.FontSize = 15;
%set(gca,'ColorScale','log')
set(gcf, 'Position',  [100, 0, 1000, 1000])
%set(gca,'ColorScale','log')

saveas(gcf,strcat(output_dir,'\','synweights_log_achgabadiff'),'epsc')
saveas(gcf,strcat(output_dir,'\','synweights_log_achgabadiff'),'png')



%% reorder projectome by L/R
% establish new order
reordered_idx_LR = zeros(1,300);
% center, place left ones at 100, right ones at 200, then remove all zero
% entries
for i = 1:75
    neuropil = np_ordered{i};
    if contains(neuropil,'_L')
        reordered_idx_LR(i + 100) = i;
    elseif contains(neuropil,'_R')
        reordered_idx_LR(i + 200) = i;
    else
        reordered_idx_LR(i) = i;
    end
end

reordered_idx_LR(reordered_idx_LR==0) = [];

np_colors_LR = zeros(size(np_colors));

% reorder correlation matrix and neuron labels by L/R
LR_synweight_matrix = zeros(size(ordered_project_to_matrix_syn_nt));
for i = 1:length(np_ordered) % row
    LR_np{i} = np_ordered{reordered_idx_LR(i)};
    LR_names{i} = np_names{reordered_idx_LR(i)};
    np_colors_LR(i,:) = np_colors(reordered_idx_LR(i),:);
    for j = 1:length(np_ordered) % column
        LR_synweight_matrix(i,j,:) = ordered_project_to_matrix_syn_nt(reordered_idx_LR(i),reordered_idx_LR(j),:);
    end
end

%% project matrices by nt log synapses LR

for k = 1:length(all_nts)
nt = all_nts{k};
figure(1)
imagesc(LR_synweight_matrix(:,:,k))
hold on
yline(7.5,'r:','linewidth',3)
hold on
xline(7.5,'r:','linewidth',3)
hold on
yline(41.5,'b:','linewidth',3)
hold on
xline(41.5,'b:','linewidth',3)

axis square
set(gca,'ytick',1:75,'yticklabel',LR_names)
set(gca,'xtick',1:75,'xticklabel',LR_names)
 xtickangle(90)
 ylabel('input neuropils','fontsize',15)
 xlabel('output neuropils','fontsize',15)
 title(nt,'fontsize',15)

 ax = gca;
for m = 1:length(np_names)
    ax.XTickLabel{m} = [strcat('\color[rgb]{',num2str(np_colors_LR(m,:)),'} '),ax.XTickLabel{m}];
    ax.YTickLabel{m} = [strcat('\color[rgb]{',num2str(np_colors_LR(m,:)),'} '),ax.YTickLabel{m}];
end

 % construct custom colormap
 max_color = nt_colors(k,:);
 map = ones(256,3);
 map(1,:)=[0.95 0.95 0.95];
 for i=2:256
    map(i,:)= map(1,:) - (map(1,:) - max_color).*i/256;
 end

colormap(map);
%colormap(colorcet('L12')); % set the colorscheme, diverging
c = colorbar;
c.Label.String='total fractional weights';
c.Label.FontSize = 15;
set(gca,'ColorScale','log')
set(gcf, 'Position',  [100, 0, 1000, 1000])
%set(gca,'ColorScale','log')

saveas(gcf,strcat(output_dir,'\',nt,'_synweights_log_proj_LR'),'epsc')
saveas(gcf,strcat(output_dir,'\',nt,'_synweights_log_proj_LR'),'png')

close all

end

%% total projectome neuron count synapses LR
LR_synweight_matrix_sum = sum(LR_synweight_matrix,3);

figure(1)
imagesc(LR_synweight_matrix_sum)
hold on
yline(7.5,'r:','linewidth',3)
hold on
xline(7.5,'r:','linewidth',3)
hold on
yline(41.5,'b:','linewidth',3)
hold on
xline(41.5,'b:','linewidth',3)

axis square
set(gca,'ytick',1:75,'yticklabel',LR_names)
set(gca,'xtick',1:75,'xticklabel',LR_names)
 xtickangle(90)

  ax = gca;
for m = 1:length(np_names)
    ax.XTickLabel{m} = [strcat('\color[rgb]{',num2str(np_colors_LR(m,:)),'} '),ax.XTickLabel{m}];
    ax.YTickLabel{m} = [strcat('\color[rgb]{',num2str(np_colors_LR(m,:)),'} '),ax.YTickLabel{m}];
end

 ylabel('input neuropils','fontsize',15)
 xlabel('output neuropils','fontsize',15)
colormap(colorcet('L1','reverse',true)); % set the colorscheme, diverging
c = colorbar;
c.Label.String='total fractional weights';
c.Label.FontSize = 15;
set(gcf, 'Position',  [100, 0, 1000, 1000])
set(gca,'ColorScale','log')

saveas(gcf,strcat(output_dir,'\','synweights_proj_LR'),'epsc')
saveas(gcf,strcat(output_dir,'\','synweights_proj_LR'),'png')

%% total projectome neuron count synapses LR without parentheses
for i = 1:length(LR_names)
    if strcmp(LR_names{i}(end),')')
    LR_names_noparens{i} = LR_names{i}(1:end-3);
    else
    LR_names_noparens{i} = LR_names{i};   
    end
end

figure(1)
imagesc(LR_synweight_matrix_sum)
hold on
yline(7.5,'r:','linewidth',3)
hold on
xline(7.5,'r:','linewidth',3)
hold on
yline(41.5,'b:','linewidth',3)
hold on
xline(41.5,'b:','linewidth',3)

axis square
set(gca,'ytick',1:75,'yticklabel',LR_names_noparens)
set(gca,'xtick',1:75,'xticklabel',LR_names_noparens)
 xtickangle(90)

  ax = gca;
for m = 1:length(np_names)
    ax.XTickLabel{m} = [strcat('\color[rgb]{',num2str(np_colors_LR(m,:)),'} '),ax.XTickLabel{m}];
    ax.YTickLabel{m} = [strcat('\color[rgb]{',num2str(np_colors_LR(m,:)),'} '),ax.YTickLabel{m}];
end

 ylabel('input neuropils','fontsize',15)
 xlabel('output neuropils','fontsize',15)
colormap(colorcet('L1','reverse',true)); % set the colorscheme, diverging
c = colorbar;
c.Label.String='total fractional weights';
c.Label.FontSize = 15;
set(gcf, 'Position',  [100, 0, 1000, 1000])
set(gca,'ColorScale','log')

saveas(gcf,strcat(output_dir,'\','synweights_proj_LR_noparens'),'epsc')
saveas(gcf,strcat(output_dir,'\','synweights_proj_LR_noparens'),'png')

%% test plot matrix before ordering
figure(1)
imagesc(project_to_matrix_nt(:,:,1))
axis square
set(gca,'ytick',1:75,'yticklabel',all_regions)
set(gca,'xtick',1:75,'xticklabel',all_regions)
 xtickangle(90)
 ylabel('input neuropils of neurons','fontsize',15)
 xlabel('output neuropils of neurons','fontsize',15)
colormap(colorcet('L1','reverse',true)); % set the colorscheme, diverging
c = colorbar;
c.Label.String='synapse number';
c.Label.FontSize = 15;
set(gcf, 'Position',  [100, 0, 1000, 1000])
set(gca,'ColorScale','log')

%% cartoon matrices (ABCD)

all_inputs = [ 10 10 20 0; 0 5 0 15; 20 0 0 5];
all_outputs = [0 20 0 30; 10 0 0 0; 0 15 5 10];

for i = 1:3

input = all_inputs(i,:);
output = all_outputs(i,:);

inputfrac = input./sum(input);
outputfrac = output./sum(output);

project_matrix_cartoon(:,:,i) = inputfrac'*outputfrac;
imagesc(project_matrix_cartoon(:,:,i))
axis square
set(gca,'ytick',1:4,'yticklabel',[{'A'} {'B'} {'C'} {'D'}])
set(gca,'xtick',1:4,'xticklabel',[{'A'} {'B'} {'C'} {'D'}])
%xtickangle(90)
%ylabel('input neuropils of neurons','fontsize',15)
% xlabel('output neuropils of neurons','fontsize',15)
colormap(colorcet('L1','reverse',true)); % set the colorscheme, diverging
set(gcf, 'Position',  [300, 300, 200, 200])
%set(gca,'ColorScale','log')

saveas(gcf,strcat(output_dir,'\','cartoon',num2str(i)),'epsc')
saveas(gcf,strcat(output_dir,'\','cartoon',num2str(i)),'png')
close all
end

%%

imagesc(sum(project_matrix_cartoon,3))
axis square
set(gca,'ytick',1:4,'yticklabel',[{'A'} {'B'} {'C'} {'D'}])
set(gca,'xtick',1:4,'xticklabel',[{'A'} {'B'} {'C'} {'D'}])
%xtickangle(90)
%ylabel('input neuropils of neurons','fontsize',15)
% xlabel('output neuropils of neurons','fontsize',15)
colormap(colorcet('L1','reverse',true)); % set the colorscheme, diverging
set(gcf, 'Position',  [300, 300, 200, 200])
%set(gca,'ColorScale','log')

saveas(gcf,strcat(output_dir,'\','cartoon_all'),'epsc')
saveas(gcf,strcat(output_dir,'\','cartoon_all'),'png')