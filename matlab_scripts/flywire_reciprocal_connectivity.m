% Using codex data format (Albert Lin)
directory = 'DIRECTORY';

cd(directory)

%% load data from csvs
% neurons in neurons.csv and neuropil synapse table are in the same order,
% thresholded
data_directory = 'DATA DIRECTORY';

opts = detectImportOptions(fullfile(data_directory,'neurons.csv'));
disp([opts.VariableNames' opts.VariableTypes'])

opts = setvartype(opts,{'root_id'},'int64');

T_n = readtable(fullfile(data_directory,'neurons.csv'),opts); 

opts2 = detectImportOptions(fullfile(data_directory,'neuropil_synapse_table.csv'));
disp([opts2.VariableNames' opts2.VariableTypes'])

opts2 = setvartype(opts2,{'root_id'},'int64');

T_np = readtable(fullfile(data_directory,'neuropil_synapse_table.csv'),opts2); 

%% load colors from spreadsheet
np_color_import = readcell('neuropil color key.xlsx','Sheet','nps','Range','A2:C79');
np_ordered = np_color_import(:,1);
np_hex_colors = np_color_import(:,2);
np_names = np_color_import(:,3);

nt_hex_colors = readcell('neuropil color key.xlsx','Sheet','nts','Range','B1:B6');

superclass_import = readcell('neuropil color key.xlsx','Sheet','superclasses','Range','A2:B12');
superclasses = superclass_import(:,1);
superclass_hex_colors = superclass_import(:,2);

%% nps alphabetical with 3 new ones for 630
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

%% load array of classification for flow
opts5 = detectImportOptions(fullfile(data_directory,'classification.csv'));
disp([opts5.VariableNames' opts5.VariableTypes'])

opts5 = setvartype(opts5,{'root_id'},'int64');

T_class = readtable(fullfile(data_directory,'classification.csv'),opts5); 

neuron_IDs_T_class = table2array(T_class(:,1));
flow_T_class = table2array(T_class(:,2));
class_T_class = table2array(T_class(:,3));

clear flow
for i = 1:length(neuron_IDs)
    idx = find(neuron_IDs_T_class == neuron_IDs(i));
flow{i} = flow_T_class{idx};
classes{i} = class_T_class{idx};
end

% these variables got moved to another sheet
%flow = table2array(T_n(:,12));
%classes = table2array(T_n(:,13));


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

%%
[max_in_deg, max_in_idx] = maxk(in_degree,10);
[max_out_deg, max_out_idx] = maxk(out_degree,10);

top10_in_deg_ns = neuron_IDs(max_in_idx)

top10_out_deg_ns = neuron_IDs(max_out_idx)


% avg synapses per incoming connection:
%sum(in_syn)/sum(in_degree)

%sum(out_syn)/sum(out_degree)
%% correlation of in and out degree
filter = (in_degree >= 5) & (out_degree >= 5);
[R,P,RL,RU] = corrcoef(in_degree(filter),out_degree(filter))

%% scatterplot of in degree vs out degree across neurons

output_dir = fullfile(directory,'Degree Scatterplots 630');
mkdir(output_dir)

%h = scatterhist(in_deg_by_in_deg,out_deg_by_in_deg,'Direction','out',...
%    'Color','k','LineStyle',{'-','-'},'Marker','.');


s = scatterhist(log10(in_degree),log10(out_degree),'Direction','out',...
    'Color','k','LineStyle','none','Marker','.');

s(2).Children.EdgeColor = 'w';
s(2).Children.BinWidth = 0.07;
s(3).Children.EdgeColor = 'w';
s(3).Children.BinWidth = 0.07;
%s(2).GridLineStyle = 'none';
xlabel('in degree')
ylabel('out degree')


xlim([-0.2 4.2])
ylim([-0.2 4.2])

axis square

set(gca,'Xtick',[0 1 2 3 4])
set(gca,'XtickLabel',[1 10 100 1000 10000])
set(gca,'Ytick',[0 1 2 3 4])
set(gca,'YtickLabel',[1 10 100 1000 10000])

%h = axes;
%get(h,'XTickLabel')
%h(2).Children(1).YScale = 'log';

%xlim([0 3000])
%ylim([0 3000])

%set(gca,'xscale','log')
%set(gca,'yscale','log')


set(gcf,'position',[2000,100,800,800])
set(gca,'fontsize',20)

set(gcf,'renderer','Painters')

%saveas(gcf,strcat(output_dir,'\all_neurons_scatter'),'epsc')
%saveas(gcf,strcat(output_dir,'\all_neurons_scatter'),'png')

%% load reciprocal lists
% load KC file
data_directory_rc = 'DATA DIRECTORY';
opts_allpairs = detectImportOptions(fullfile(data_directory_rc,'v630-all-reciprocal-pairs-s1.csv'));
disp([opts_allpairs.VariableNames' opts_allpairs.VariableTypes'])

opts_allpairs = setvartype(opts_allpairs,{'n1'},'int64');
opts_allpairs = setvartype(opts_allpairs,{'n2'},'int64');

T_rc = readtable(fullfile(data_directory_rc,'v630-all-reciprocal-pairs-s1.csv'),opts_allpairs); 
rc_ns = table2array(T_rc(1:end,1:2));

%% distributions of numbers of neurons participating in recurrent connections

output_dir = 'OUTPUT DIRECTORY';

rc_ns_reshape = reshape(rc_ns,[],1);
[unique_rns,ia,ic] = unique(rc_ns_reshape);
rconn_counts = accumarray(ic,1);

rconn_ns = [unique_rns, rconn_counts];
rconn_ns_ordered = sortrows(rconn_ns,2,'descend'); % row 2 is recurrent degree

rconn_ns_ids_ordered = rconn_ns_ordered(:,1);
rconns_ordered = cast(rconn_ns_ordered(:,2),'double');

% order by rconn_counts

edges = 0:0.1:4;

%histogram(in_degree,10.^edges,'FaceColor','r','EdgeColor','r')
%hold on
%histogram(out_degree,10.^edges,'FaceColor','b','EdgeColor','b')
%hold on
fig = histogram(rconns_ordered,10.^edges,'FaceColor','k','EdgeColor','k');
set(gca,'xscale','log')
set(gca,'yscale','log')
%xlim([0 5000])
ylim([0 10^5])
xlabel('number of reciprocal connections')
ylabel('neuron count')
%title(strcat('neurons for which',{' '},neuropil,' is the max-in region'))

set(gcf,'position',[100,100,300,300])
set(gca,'FontSize',14)

set(gcf,'renderer','Painters')

saveas(fig,strcat(output_dir,'\recurrent_degree_dist'),'epsc')
saveas(fig,strcat(output_dir,'\recurrent_degree_dist'),'png')

%%
histogram(in_degree,10.^edges,'FaceColor','r','EdgeColor','r')
hold on
histogram(out_degree,10.^edges,'FaceColor','b','EdgeColor','b')
hold on
fig = histogram(rconns_ordered,10.^edges,'FaceColor','k','EdgeColor','k');
set(gca,'xscale','log')
%set(gca,'yscale','log')
%xlim([0 5000])
xlabel('no. of reciprocal connections')
ylabel('neuron count')
%title(strcat('neurons for which',{' '},neuropil,' is the max-in region'))

set(gcf,'position',[100,100,500,500])
set(gca,'FontSize',14)

legend("in degree","out degree","reciprocal degree")

set(gcf,'renderer','Painters')

saveas(fig,strcat(output_dir,'\recurrent_degree_dist_in_out_deg_semilog'),'epsc')
saveas(fig,strcat(output_dir,'\recurrent_degree_dist_in_out_deg_semilog'),'png')

%% build individual vectors that contains other relevant info (nts, normal degree, highest neuropil in and out fractions)

clear rc_ns_nt rc_ns_in_np rc_ns_out_np

rc_ns_in_degree = zeros(size(unique_rns));
rc_ns_out_degree = zeros(size(unique_rns));
rc_ns_nt = {};
rc_ns_flow = {};

rc_ns_max_in_frac = zeros(size(unique_rns));
rc_ns_max_out_frac = zeros(size(unique_rns));
rc_ns_in_np = {};
rc_ns_out_np = {};

n = 0;

for i = 1:length(unique_rns)
    root_id = rconn_ns_ordered(i,1);
    idx = find(neuron_IDs == root_id);
    if ~isempty(idx)
    rc_ns_in_degree(i) = in_degree(idx);
    rc_ns_out_degree(i) = out_degree(idx);
    rc_ns_nt{i} = out_nts{idx};
    rc_ns_flow{i} = flow{idx};

    [rc_ns_max_in_frac(i),np_in_idx] = max(in_neuron_frac(idx,:)); %nps not reordered
    rc_ns_in_np{i} = all_regions{np_in_idx};

    [rc_ns_max_out_frac(i),np_out_idx] = max(out_neuron_frac(idx,:)); %nps not reordered
    rc_ns_out_np{i} = all_regions{np_out_idx};

    else
    %disp(root_id);
    rc_ns_in_degree(i) = NaN;
    rc_ns_out_degree(i) = NaN;
    n = n+1;
    end
end

n

%% identify the neurotransmitters of those neurons, distribution by nt

for i = 1:6

nt_idx = strcmpi(rc_ns_nt,all_nts{i});
edges = 0:0.1:4; 

%histogram(in_degree,10.^edges,'FaceColor','r','EdgeColor','r')
%hold on
%histogram(out_degree,10.^edges,'FaceColor','b','EdgeColor','b')
%hold on
histogram(rconns_ordered(nt_idx),10.^edges,'FaceColor',nt_colors(i,:),'EdgeColor','k','FaceAlpha',1);
set(gca,'xscale','log')
set(gca,'yscale','log')
%xlim([0 5000])
ylim([0 10^5])
xlabel('no. of reciprocal connections')
ylabel('neuron count')
title(strcat(all_nts{i},{' '},'neurons'))

set(gcf,'position',[100,100,350,220])
set(gca,'FontSize',14)

%legend("in degree","out degree","reciprocal degree")

set(gcf,'renderer','Painters')

saveas(gcf,strcat(output_dir,strcat('\recurrent_degree_dist_',all_nts{i})),'epsc')
saveas(gcf,strcat(output_dir,strcat('\recurrent_degree_dist_',all_nts{i})),'png')
close all

end

%% identify the neurotransmitters of those neurons, distribution by nt semilog

for i = 1:6

nt_idx = strcmpi(rc_ns_nt,all_nts{i});
edges = 0:0.1:4;

%histogram(in_degree,10.^edges,'FaceColor','r','EdgeColor','r')
%hold on
%histogram(out_degree,10.^edges,'FaceColor','b','EdgeColor','b')
%hold on
histogram(rconns_ordered(nt_idx),10.^edges,'FaceColor',nt_colors(i,:),'EdgeColor','k','FaceAlpha',1);
set(gca,'xscale','log')
%xlim([0 5000])
ylim([0 3*10^4])
xlabel('no. of reciprocal connections')
ylabel('neuron count')
title(strcat(all_nts{i},{' '},'neurons'))

set(gcf,'position',[100,100,350,220])
set(gca,'FontSize',14)

%legend("in degree","out degree","reciprocal degree")

set(gcf,'renderer','Painters')

saveas(gcf,strcat(output_dir,strcat('\recurrent_degree_dist_semilog_',all_nts{i})),'epsc')
saveas(gcf,strcat(output_dir,strcat('\recurrent_degree_dist_semilog_',all_nts{i})),'png')
close all

end


%% fraction of in degrees that are recurrent

% some are over, why?

rc_ns_in_recfrac = rconns_ordered./rc_ns_in_degree;

edges = 0:0.02:1;

fig = histogram(rc_ns_in_recfrac,edges,'FaceColor',[0.8 0.5 0.4],'EdgeColor','k','FaceAlpha',1);
%set(gca,'xscale','log')
%set(gca,'yscale','log')
xlim([0 1])
xlabel('reciprocal connections/incoming connections')
ylabel('neuron count')
%title(strcat('neurons for which',{' '},neuropil,' is the max-in region'))

set(gcf,'position',[100,100,500,400])
set(gca,'FontSize',12)

set(gcf,'renderer','Painters')

saveas(fig,strcat(output_dir,'\recurrent_in_frac'),'epsc')
saveas(fig,strcat(output_dir,'\recurrent_in_frac'),'png')

%% fraction of out degrees that are recurrent

rc_ns_out_recfrac = rconns_ordered./rc_ns_out_degree;

edges = 0:0.02:1;

fig = histogram(rc_ns_out_recfrac,edges,'FaceColor',[0.4 0.6 1],'EdgeColor','k','FaceAlpha',1);
%set(gca,'xscale','log')
%set(gca,'yscale','log')
xlim([0 1])
xlabel('reciprocal connections/outgoing connections')
ylabel('neuron count')
%title(strcat('neurons for which',{' '},neuropil,' is the max-in region'))

set(gcf,'position',[100,100,500,400])
set(gca,'FontSize',12)

set(gcf,'renderer','Painters')

saveas(fig,strcat(output_dir,'\recurrent_out_frac'),'epsc')
saveas(fig,strcat(output_dir,'\recurrent_out_frac'),'png')

%% compare recurrent degree to in and out degree fraction (scatterplots)
% remember that a recurrent connection counts as both an in degree and and
% out degree
% we could ask what fraction of degree is accounted for by recurrent
% connections

s = scatterhistogram(rc_ns_in_recfrac,rc_ns_out_recfrac,'color',[0 0 0],'MarkerStyle','.');
s.ScatterPlotLocation = 'NorthEast';
s.ScatterPlotProportion = 0.8;
s.XHistogramDirection = 'down';
s.YHistogramDirection = 'left';
s.HistogramDisplayStyle = 'bar';
s.LineStyle = 'none';
%scatterhandle = findobj(s(1))

xlabel('reciprocal connections/incoming connections')
ylabel('reciprocal connections/outgoing connections')
%xlabel('in degree')
%ylabel('out degree')
%s(2).Children(1).YScale = 'log';

xlim([0 1])
ylim([0 1])

set(gca,'FontSize',12)

%set(gca,'xscale','log')
%set(gca,'yscale','log')




%% compare recurrent degree to in and out degree fraction heatmap
% remember that a recurrent connection counts as both an in degree and and
% out degree
% we could ask what fraction of degree is accounted for by recurrent
% connections

Xedges = 0:0.05:1;
Yedges = 0:0.05:1;

h = histogram2(rc_ns_in_recfrac,rc_ns_out_recfrac,Xedges,Yedges,'DisplayStyle','tile','ShowEmptyBins','on');
h.EdgeColor = 'none';
set(gca,'ColorScale','log')
cb1 = colorbar;
cb1.Label.String = 'neuron counts';

%scatterhandle = findobj(s(1))

xlabel('reciprocal conn./incoming conn.')
ylabel('reciprocal conn./outgoing conn.')
%xlabel('in degree')
%ylabel('out degree')
%s(2).Children(1).YScale = 'log';

xlim([0 1])
ylim([0 1])
axis square

set(gca,'FontSize',20)

set(gcf,'position',[100,100,600,600])

set(gcf,'renderer','Painters')

saveas(gcf,strcat(output_dir,'\infrac_vs_outfrac_heatmap'),'epsc')
saveas(gcf,strcat(output_dir,'\infrac_vs_outfrac_heatmap'),'png')

%% try scatterplot of in vs out, colored by recurrence
f = figure;

s = scatter(rc_ns_in_degree,rc_ns_out_degree,3,rconns_ordered,'filled');

colormap(colorcet('L12','reverse',false)); % set the colorscheme, diverging
c = colorbar;
set(gca,'ColorScale','log')


xlabel('in degree')
ylabel('out degree')
c.Label.String='no. reciprocal connections';
c.Label.FontSize = 15;

xlim([0 6500])
ylim([0 6500])

set(gca,'FontSize',15)

%set(gca,'xscale','log')
%set(gca,'yscale','log')


set(gcf,'position',[100,100,600,600])

set(gcf,'renderer','Painters')

saveas(gcf,strcat(output_dir,'\in_vs_out_deg_rec_lin'),'epsc')
saveas(gcf,strcat(output_dir,'\in_vs_out_deg_rec_lin'),'png')

%% scatterplot of recurrence x 2 vs total degree
f = figure;

error_filter = (rc_ns_in_degree+rc_ns_out_degree) >= rconns_ordered*2;

s = scatter(rc_ns_in_degree(error_filter)+rc_ns_out_degree(error_filter),rconns_ordered(error_filter)*2,10,'k','filled',...
    'MarkerFaceAlpha',0.2);
hold on
xline = 0:5;
plot(10.^xline,10.^xline,'k','linewidth',2)
hold on
plot(10.^xline,2*10.^xline,'--k')
hold on
plot(2*10.^xline,10.^xline,'--k')

xlabel('total degree')
ylabel('2x reciprocal degree')


xlim([0 12000])
ylim([0 12000])

set(gca,'FontSize',13)

set(gca,'xscale','log')
set(gca,'yscale','log')


set(gcf,'position',[100,100,600,600])

set(gcf,'renderer','Painters')

% works but takes forever and kind of crashes matlab
%saveas(gcf,strcat(output_dir,'\deg_vs_rec'),'pdf')

%works okay but is slow
saveas(gcf,strcat(output_dir,'\deg_vs_rec'),'svg')

% doesn't work for scatter
%plot2svg(strcat(output_dir,'\deg_vs_rec.svg'))

saveas(gcf,strcat(output_dir,'\deg_vs_rec'),'png')

%% count of neurons
rc_ns_intrinsic_filter = strcmp(rc_ns_flow,'intrinsic');

sum((rc_ns_in_degree(rc_ns_intrinsic_filter)+rc_ns_out_degree(rc_ns_intrinsic_filter))/2 < rconns_ordered(rc_ns_intrinsic_filter)*2)
%% scatterplot of recurrence x 2 vs total degree by neurotransmitter
f = figure;

error_filter = (rc_ns_in_degree+rc_ns_out_degree) >= rconns_ordered*2;


for i = 1:6
nt_idx = strcmpi(rc_ns_nt,all_nts{i});

s = scatter(rc_ns_in_degree(error_filter & nt_idx')+rc_ns_out_degree(error_filter & nt_idx'),...
    rconns_ordered(error_filter& nt_idx')*2,15,nt_colors(i,:),'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3);
hold on

end
xline = 0:5;
plot(10.^xline,10.^xline,'k','linewidth',2)
hold on
plot(10.^xline,2*10.^xline,'--k')
hold on
plot(2*10.^xline,10.^xline,'--k')

xlabel('total degree')
ylabel('2x reciprocal degree')

xlim([0 13000])
ylim([0 13000])

%legend(all_nts,'location','northwest')

set(gca,'FontSize',13)

set(gca,'xscale','log')
set(gca,'yscale','log')


set(gcf,'position',[100,100,600,600])

set(gcf,'renderer','Painters')

%works okay but is slow (eps doesn't work for transparency)
saveas(gcf,strcat(output_dir,'\deg_vs_rec_allnt'),'svg')

saveas(gcf,strcat(output_dir,'\deg_vs_rec_allnt'),'png')

%% scatterplot of recurrence x 2 vs total degree by individual neurotransmitter
f = figure;

error_filter = (rc_ns_in_degree+rc_ns_out_degree) >= rconns_ordered*2;


for i = 1:6
nt_idx = strcmpi(rc_ns_nt,all_nts{i});

s = scatter(rc_ns_in_degree(error_filter & nt_idx')+rc_ns_out_degree(error_filter & nt_idx'),...
    rconns_ordered(error_filter& nt_idx')*2,15,nt_colors(i,:),'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3);
hold on 

xline = 0:5;
plot(10.^xline,10.^xline,'k','linewidth',2)
hold on
plot(10.^xline,2*10.^xline,'--k')
hold on
plot(2*10.^xline,10.^xline,'--k')

xlabel('total degree')
ylabel('2x reciprocal degree')


title(all_nts{i})


xlim([0 13000])
ylim([0 13000])


set(gca,'FontSize',13)

set(gca,'xscale','log')
set(gca,'yscale','log')


set(gcf,'position',[100,100,600,600])

set(gcf,'renderer','Painters')

%works okay but is slow (eps doesn't work for transparency)
saveas(gcf,strcat(output_dir,strcat('\deg_vs_rec_',all_nts{i})),'svg')

saveas(gcf,strcat(output_dir,strcat('\deg_vs_rec_',all_nts{i})),'png')

close all
end


%% distribution of max in single-neuropil fraction of recurrent neurons (break this by nt)

edges = 0:0.02:1;

% remove 0s
fig = histogram(nonzeros(rc_ns_max_in_frac),edges,'FaceColor',[0.8 0.5 0.4],'EdgeColor','k','FaceAlpha',1);
%set(gca,'xscale','log')
%set(gca,'yscale','log')
xlim([0 1])
xlabel('max frac. of incoming connections in a single neuropil')
ylabel('neuron count')
title(strcat('neurons participating in reciprocal connections'))

set(gcf,'position',[100,100,500,400])
set(gca,'FontSize',12)

set(gcf,'renderer','Painters')

saveas(gcf,strcat(output_dir,'\max_in_frac'),'epsc')
saveas(gcf,strcat(output_dir,'\max_in_frac'),'png')

%% distribution of max out single-neuropil fraction of recurrent neurons (break this by nt)
edges = 0:0.02:1;

% remove 0s
fig = histogram(nonzeros(rc_ns_max_out_frac),edges,'FaceColor',[0.4 0.6 1],'EdgeColor','k','FaceAlpha',1);
%set(gca,'xscale','log')
%set(gca,'yscale','log')
xlim([0 1])
xlabel('max frac. of outgoing connections in a single neuropil')
ylabel('neuron count')
title(strcat('neurons participating in reciprocal connections'))

set(gcf,'position',[100,100,500,400])
set(gca,'FontSize',12)

set(gcf,'renderer','Painters')

saveas(gcf,strcat(output_dir,'\max_out_frac'),'epsc')
saveas(gcf,strcat(output_dir,'\max_out_frac'),'png')

%% set a cutoff for "highly reciprocal" neurons, histogram of neuropil frac outgoing
% more than 50 reciprocal connections
thresh = 98/2;
high_filter = rconns_ordered >= thresh;
sum(high_filter)

edges = 0:0.02:1;
fig = histogram(rc_ns_max_out_frac(high_filter),edges,'FaceColor',[0.4 0.6 1],'EdgeColor','k','FaceAlpha',1);
xlim([0 1])
xlabel('max frac. of outgoing connections in a single neuropil')
ylabel('neuron count')
title(strcat('neurons with >=',num2str(thresh), 'reciprocal connections'))

set(gcf,'position',[100,100,500,400])
set(gca,'FontSize',12)

set(gcf,'renderer','Painters')

saveas(gcf,strcat(output_dir,'\max_out_frac_',num2str(thresh)),'epsc')
saveas(gcf,strcat(output_dir,'\max_out_frac_',num2str(thresh)),'png')

%% histogram of neuropil frac incoming
edges = 0:0.02:1;
fig = histogram(rc_ns_max_in_frac(high_filter),edges,'FaceColor',[0.8 0.5 0.4],'EdgeColor','k','FaceAlpha',1);
xlim([0 1])
xlabel('max frac. of incoming connections in a single neuropil')
ylabel('neuron count')
title(strcat('neurons with >=',num2str(thresh), 'reciprocal connections'))

set(gcf,'position',[100,100,500,400])
set(gca,'FontSize',12)

set(gcf,'renderer','Painters')

saveas(gcf,strcat(output_dir,'\max_in_frac_',num2str(thresh)),'epsc')
saveas(gcf,strcat(output_dir,'\max_in_frac_',num2str(thresh)),'png')


%% do neuropils (based on fraction of connections)
% if both are high in the same neuropil, then it is tiling
thresh = 37;
high_filter = (rc_ns_in_degree+rc_ns_out_degree) >= thresh;
highly_rec_filter = ((rc_ns_in_degree+rc_ns_out_degree)/2 < rconns_ordered*2);

% TO DO: set a threshold for fraction of connections before doing this
frac_thresh = 0.5;

high_fracin_filter = rc_ns_max_in_frac >= frac_thresh;
high_fracout_filter = rc_ns_max_out_frac >= frac_thresh;

high_in_np = rc_ns_in_np(high_filter & high_fracin_filter & high_fracout_filter & highly_rec_filter & rc_ns_intrinsic_filter');
high_out_np = rc_ns_out_np(high_filter & high_fracin_filter & high_fracout_filter & highly_rec_filter & rc_ns_intrinsic_filter');
high_ns_ids = rconn_ns_ids_ordered(high_filter & high_fracin_filter & high_fracout_filter & highly_rec_filter & rc_ns_intrinsic_filter');
high_ns_nt = rc_ns_nt(high_filter & high_fracin_filter & high_fracout_filter & highly_rec_filter & rc_ns_intrinsic_filter');

internal_high_filter = strcmp(high_in_np,high_out_np);

% number of "internal" reciprocal neurons
sum(internal_high_filter)

internal_high_nps = high_in_np(internal_high_filter);
internal_high_ns = high_ns_ids(internal_high_filter);
internal_high_nts = high_ns_nt(internal_high_filter);

[unique_high_internal_nps,ia,ic] = unique(internal_high_nps);
unique_high_internal_np_counts = accumarray(ic,1);

% bar plot tiling neurons
% bar(unique_high_internal_np_counts)
% xlabel('neuropils')
% set(gca,'xtick',1:length(unique_high_internal_nps),'xticklabel',unique_high_internal_nps)
% xtickangle(90)
% ylabel('neuron count')
% title(strcat('tiling neurons: frac. of neuropil conn. >= ',num2str(frac_thresh),'& reciprocal deg. >= ',num2str(thresh)))


% sort for the purposes of plotting histogram with names and colors
unique_high_internal_np_counts_forplotting = zeros(size(unique_high_internal_np_counts));
unique_high_internal_nps_forplotting = {};
unique_high_internal_npnames_forplotting = {};

idx = 1;
for i = 1:length(np_ordered)
if sum(strcmp(unique_high_internal_nps,np_ordered{i})) > 0
    data_idx = find(strcmp(unique_high_internal_nps,np_ordered{i}));
    unique_high_internal_np_counts_forplotting(idx) = unique_high_internal_np_counts(data_idx);

    unique_high_internal_nps_forplotting{idx} = np_ordered{i};
    unique_high_internal_npnames_forplotting{idx} = np_names{i};

    unique_high_np_colors(idx,:) = np_colors(i,:);

    idx = idx+1;
end
end

% find the proper printing version of these neuropils
bar(unique_high_internal_np_counts_forplotting,'facecolor','k','edgecolor','none')
%xlabel('neuropils')
set(gca,'xtick',1:length(unique_high_internal_npnames_forplotting),'xticklabel',unique_high_internal_npnames_forplotting)
xtickangle(90)

ax = gca;
for m = 1:length(unique_high_internal_npnames_forplotting)
     ax.XTickLabel{m} = [strcat('\color[rgb]{',num2str(unique_high_np_colors(m,:)),'} '),ax.XTickLabel{m}];
end


set(gca,'FontSize',11)

ylabel('neuron count','FontSize',15)
%title(strcat('frac. of np conn. >= ',num2str(frac_thresh),'& rec. deg. >= ',num2str(thresh)))

set(gcf,'position',[100,100,800,400])


set(gcf,'renderer','Painters')

saveas(gcf,strcat(output_dir,strcat('\tiling_reciprocal_neurons_bar')),'epsc')
saveas(gcf,strcat(output_dir,strcat('\tiling_reciprocal_neurons_bar')),'png')

% print this plot
%% PRINT THESE AS CSV

%%
neuropil = "LO_R";
% find neurons associated with this neuropil
np_idx = find(strcmp(unique_high_internal_nps,neuropil));
internal_high_ns(find(ic == np_idx))

% scatterplot in frac out frac

%% nts of internal high neurons
%internal_high_ns = high_ns_ids(internal_high_filter);
internal_high_nts = high_ns_nt(internal_high_filter);

% pie chart of nt composition of 4 categories

nt_num_internal_high = zeros(1,6);

nt_frac_internal_high = zeros(1,6);


for i = 1:6

nt_num_internal_high(i) = sum(strcmpi(internal_high_nts,all_nts{i}));
nt_frac_internal_high(i) = nt_num_internal_high(i)/length(internal_high_nts);

end



explode = [0 0 0 1 1 1];
% Create pie charts
p1 = pie([nt_frac_internal_high,1-sum(nt_frac_internal_high)]); %,explode
for i=1:2:13
p1(i).EdgeColor = 'none';
p1(i+1).FontSize = 30;
end
colormap([nt_colors;0.8 0.8 0.8])
%title('all neurons','FontSize',20)


% Create legend
nt_legend = all_nts;
nt_legend{7} = 'uncertain';

lgd = legend(nt_legend,'FontSize',30,'location','eastoutside');
%lgd.Layout.Tile = 'east';


% nt fractions for all neurons
set(gcf,'position',[100,100,1000,600])

set(gcf,'renderer','Painters')

saveas(gcf,strcat(output_dir,'\internal_nt'),'epsc')
saveas(gcf,strcat(output_dir,'\internal_nt'),'png')


%% apply filters using flow and classes

n_intrinsic_filter = strcmp(flow,'intrinsic');

factor = 5;

% look at degree scatterplot for integrators and broadcasters
% to do: shade each region
% scatterplot

scatter(in_degree(n_intrinsic_filter),out_degree(n_intrinsic_filter),10,'filled','markerfacecolor',[0.5 0.5 0.5]);
hold on
% plot circle
r = 37; % double the reciprocal degree
circ_x = 1:0.001:r;
circ_y = (r^2 - circ_x.^2).^0.5;
plot(circ_x,circ_y,'k','linewidth',2)
hold on

% number of neurons in the circle
small_neuron_number = sum((in_degree(n_intrinsic_filter)+out_degree(n_intrinsic_filter)) < r)

% plot boundary for "broadcast" (few to many) neurons
% out_degree > 40, out_degree 5x in_degree
broadcast_number = sum(((in_degree(n_intrinsic_filter)+out_degree(n_intrinsic_filter)) >= r) & (out_degree(n_intrinsic_filter) >= in_degree(n_intrinsic_filter)*factor))
y_broadcast = r:7000;
x_broadcast = y_broadcast./factor;
plot(x_broadcast,y_broadcast,'k','linewidth',2)
hold on

% plot boundary for "integrator" (many to few) neurons
% in_degree > 40, out_degree 5x in_degree
integrate_number = sum(((in_degree(n_intrinsic_filter)+out_degree(n_intrinsic_filter)) >= r) & (in_degree(n_intrinsic_filter) >= out_degree(n_intrinsic_filter)*factor))
x_integrate = r:7000;
y_integrate = x_integrate./factor;
plot(x_integrate,y_integrate,'k','linewidth',2)

prolific_number = length(in_degree(n_intrinsic_filter))-small_neuron_number-broadcast_number-integrate_number

xlabel('in degree')
ylabel('out degree')


xlim([1 7000])
ylim([1 7000])

set(gca,'FontSize',16)

set(gca,'xscale','log')
set(gca,'yscale','log')


set(gcf,'position',[100,100,600,600])

set(gcf,'renderer','Painters')

saveas(gcf,strcat(output_dir,'\in_vs_out_deg_phase'),'epsc')
saveas(gcf,strcat(output_dir,'\in_vs_out_deg_phase'),'png')

%% plot scatterplot by neurotransmitter
%% look at degree scatterplot for integrators and broadcasters
% to do: shade each region

for i = 1:6
nt_idx = strcmpi(out_nts,all_nts{i});

s = scatter(in_degree(n_intrinsic_filter & nt_idx),out_degree(n_intrinsic_filter & nt_idx),...
    10,nt_colors(i,:),'filled');

hold on

end

% plot circle
r = 37; % double the reciprocal degree
circ_x = 1:0.001:r;
circ_y = (r^2 - circ_x.^2).^0.5;
plot(circ_x,circ_y,'k','linewidth',2)
hold on

% number of neurons in the circle
%small_neuron_number = sum((in_degree(n_intrinsic_filter)+out_degree(n_intrinsic_filter)) < r)

% plot boundary for "broadcast" (few to many) neurons
% out_degree > 40, out_degree 5x in_degree
%broadcast_number = sum(((in_degree(n_intrinsic_filter)+out_degree(n_intrinsic_filter)) > r) & (out_degree(n_intrinsic_filter) >= in_degree(n_intrinsic_filter)*5))
y_broadcast = r:7000;
x_broadcast = y_broadcast./factor;
plot(x_broadcast,y_broadcast,'k','linewidth',2)
hold on

% plot boundary for "integrator" (many to few) neurons
% in_degree > 40, out_degree 5x in_degree
%integrate_number = sum(((in_degree(n_intrinsic_filter)+out_degree(n_intrinsic_filter)) > r) & (in_degree(n_intrinsic_filter) >= out_degree(n_intrinsic_filter)*5))
x_integrate = r:7000;
y_integrate = x_integrate./factor;
plot(x_integrate,y_integrate,'k','linewidth',2)

%prolific_number = length(in_degree(n_intrinsic_filter))-small_neuron_number-broadcast_number-integrate_number

xlabel('in degree')
ylabel('out degree')

legend(all_nts,'location','northwest')

xlim([1 7000])
ylim([1 7000])

set(gca,'FontSize',16)

set(gca,'xscale','log')
set(gca,'yscale','log')


set(gcf,'position',[100,100,600,600])

set(gcf,'renderer','Painters')

saveas(gcf,strcat(output_dir,'\in_vs_out_deg_phase_nt'),'epsc')
saveas(gcf,strcat(output_dir,'\in_vs_out_deg_phase_nt'),'png')

%% print each nt individually for supplement

%% plot by superclass
% look at degree scatterplot for integrators and broadcasters
% to do: shade each region
% intrinsic classes only
for i = 6:9
class_idx = contains(classes,superclasses{i});

s = scatter(in_degree(n_intrinsic_filter & class_idx),out_degree(n_intrinsic_filter & class_idx),...
    10,superclass_colors(i,:),'filled');

hold on

end

% plot circle
r = 37; % double the reciprocal degree
circ_x = 1:0.001:r;
circ_y = (r^2 - circ_x.^2).^0.5;
plot(circ_x,circ_y,'k','linewidth',2)
hold on

% number of neurons in the circle
%small_neuron_number = sum((in_degree(n_intrinsic_filter)+out_degree(n_intrinsic_filter)) < r)

% plot boundary for "broadcast" (few to many) neurons
% out_degree > 40, out_degree 5x in_degree
%broadcast_number = sum(((in_degree(n_intrinsic_filter)+out_degree(n_intrinsic_filter)) > r) & (out_degree(n_intrinsic_filter) >= in_degree(n_intrinsic_filter)*5))
y_broadcast = r:7000;
x_broadcast = y_broadcast./factor;
plot(x_broadcast,y_broadcast,'k','linewidth',2)
hold on

% plot boundary for "integrator" (many to few) neurons
% in_degree > 40, out_degree 5x in_degree
%integrate_number = sum(((in_degree(n_intrinsic_filter)+out_degree(n_intrinsic_filter)) > r) & (in_degree(n_intrinsic_filter) >= out_degree(n_intrinsic_filter)*5))
x_integrate = r:7000;
y_integrate = x_integrate./factor;
plot(x_integrate,y_integrate,'k','linewidth',2)

%prolific_number = length(in_degree(n_intrinsic_filter))-small_neuron_number-broadcast_number-integrate_number

xlabel('in degree')
ylabel('out degree')

% superclasses(6:9)
legend({'optic', 'visual projection', 'visual centrifugal','central'},'location','northwest')

xlim([1 7000])
ylim([1 7000])

set(gca,'FontSize',16)

set(gca,'xscale','log')
set(gca,'yscale','log')


set(gcf,'position',[100,100,600,600])

set(gcf,'renderer','Painters')

saveas(gcf,strcat(output_dir,'\in_vs_out_deg_phase_class'),'epsc')
saveas(gcf,strcat(output_dir,'\in_vs_out_deg_phase_class'),'png')

%% pull neuron lists in many-to-few, few-to-many, with in-degree, out-degree, nt, class

broadcast_filter = ((in_degree+out_degree) >= r)...
    & (out_degree >= in_degree*factor)  & n_intrinsic_filter';

broadcast_neurons = neuron_IDs(broadcast_filter);


integrate_filter = ((in_degree+out_degree) >= r)...
    & (in_degree >= out_degree*factor)  & n_intrinsic_filter';

integrate_neurons = neuron_IDs(integrate_filter);

rich_club_filter = ((in_degree+out_degree) >= r);

rich_club_neurons = neuron_IDs(rich_club_filter);

intrinsic_balanced_filter = rich_club_filter & ~broadcast_filter & ~integrate_filter & n_intrinsic_filter';

intrinsic_balanced_neurons = neuron_IDs(intrinsic_balanced_filter);


%% print integrate and broadcast neurons
% for winner-take-all use flywire_analyses_nt code
% projecting to means outputs point to this neuropil

largescale_folder = fullfile(directory,'largescale_csvs_630\');
mkdir(largescale_folder)

project_to_neuron_num = zeros(size(all_regions));


writematrix(integrate_neurons,strcat(largescale_folder,'integrate_neurons.csv'))

writematrix(broadcast_neurons,strcat(largescale_folder,'broadcast_neurons.csv'))

writematrix(rich_club_neurons,strcat(largescale_folder,'rich_club_neurons.csv'))

writematrix(intrinsic_balanced_neurons,strcat(largescale_folder,'intrinsic_balanced_neurons.csv'))

%% comparisons between "all neurons", "all rich club neurons", and broadcasters and integrators
% pie chart of nt composition of 4 categories

nt_num_all = zeros(1,6);
nt_num_rich = zeros(1,6);
nt_num_integrate = zeros(1,6);
nt_num_broadcast = zeros(1,6);
nt_num_balanced = zeros(1,6);

nt_frac_all = zeros(1,6);
nt_frac_rich = zeros(1,6);
nt_frac_integrate = zeros(1,6);
nt_frac_broadcast = zeros(1,6);
nt_frac_balanced = zeros(1,6);

for i = 1:6

nt_num_all(i) = sum(strcmpi(out_nts,all_nts{i}));
nt_frac_all(i) = nt_num_all(i)/length(neuron_IDs);

nt_num_rich(i) = sum(strcmpi(out_nts(rich_club_filter),all_nts{i}));
nt_frac_rich(i) = nt_num_rich(i)/sum(rich_club_filter);

nt_num_integrate(i) = sum(strcmpi(out_nts(integrate_filter),all_nts{i}));
nt_frac_integrate(i) = nt_num_integrate(i)/sum(integrate_filter);

nt_num_broadcast(i) = sum(strcmpi(out_nts(broadcast_filter),all_nts{i}));
nt_frac_broadcast(i) = nt_num_broadcast(i)/sum(broadcast_filter);

nt_num_balanced(i) = sum(strcmpi(out_nts(intrinsic_balanced_filter),all_nts{i}));
nt_frac_balanced(i) = nt_num_balanced(i)/sum(intrinsic_balanced_filter);

end

t = tiledlayout(1,4,'TileSpacing','compact');

explode = [0 0 0 1 1 1];
% Create pie charts
ax1 = nexttile;
p1 = pie(ax1,[nt_frac_all,1-sum(nt_frac_all)]); %,explode
for i=1:2:13
p1(i).EdgeColor = 'none';
p1(i+1).FontSize = 30;
end
colormap([nt_colors;0.8 0.8 0.8])
%title('all neurons','FontSize',20)



ax2 = nexttile;
p2=pie(ax2, [nt_frac_rich,1-sum(nt_frac_rich)]);
for i=1:2:13
p2(i).EdgeColor = 'none';
p2(i+1).FontSize = 30;
end
%title('rich club neurons','FontSize',20)

ax3 = nexttile;
p3=pie(ax3, [nt_frac_integrate,1-sum(nt_frac_integrate)]);
for i=1:2:13
p3(i).EdgeColor = 'none';
p3(i+1).FontSize = 30;
end
%title('integrator neurons','FontSize',20)

ax4 = nexttile;
p4=pie(ax4, [nt_frac_broadcast,1-sum(nt_frac_broadcast)]);
for i=1:2:13
p4(i).EdgeColor = 'none';
p4(i+1).FontSize = 30;
end
%title('broadcaster neurons','FontSize',20)

% ax5 = nexttile;
% p5=pie(ax5, [nt_frac_balanced,1-sum(nt_frac_balanced)]);
% for i=1:2:13
% p5(i).EdgeColor = 'none';
% end
% title('large balanced neurons')

% Create legend
nt_legend = all_nts;
nt_legend{7} = 'uncertain';

lgd = legend(nt_legend,'FontSize',30);
lgd.Layout.Tile = 'east';


% nt fractions for all neurons
set(gcf,'position',[100,100,2500,600])

%set(gcf,)

set(gcf,'renderer','Painters')

saveas(gcf,strcat(output_dir,'\richclub_nt'),'epsc')
saveas(gcf,strcat(output_dir,'\richclub_nt'),'png')

%% pie charts by the four intrinsic classes
class_num_all = zeros(1,4);
class_num_rich = zeros(1,4);
class_num_integrate = zeros(1,4);
class_num_broadcast = zeros(1,4);
class_num_balanced = zeros(1,4);

class_frac_all = zeros(1,4);
class_frac_rich = zeros(1,4);
class_frac_integrate = zeros(1,4);
class_frac_broadcast = zeros(1,4);
class_frac_balanced = zeros(1,4);

for i = 1:4

class_num_all(i) = sum(strcmpi(classes,superclasses{i+5}));
class_frac_all(i) = class_num_all(i)/sum(n_intrinsic_filter);

class_num_rich(i) = sum(strcmpi(classes(rich_club_filter & n_intrinsic_filter'),superclasses{i+5}));
class_frac_rich(i) = class_num_rich(i)/sum(rich_club_filter & n_intrinsic_filter');

class_num_integrate(i) = sum(strcmpi(classes(integrate_filter),superclasses{i+5}));
class_frac_integrate(i) = class_num_integrate(i)/sum(integrate_filter);

class_num_broadcast(i) = sum(strcmpi(classes(broadcast_filter),superclasses{i+5}));
class_frac_broadcast(i) = class_num_broadcast(i)/sum(broadcast_filter);

class_num_balanced(i) = sum(strcmpi(classes(intrinsic_balanced_filter),superclasses{i+5}));
class_frac_balanced(i) = class_num_balanced(i)/sum(intrinsic_balanced_filter);

end

t = tiledlayout(1,4,'TileSpacing','compact');

explode = [0 0 0 1 1 1];
% Create pie charts
ax1 = nexttile;
p1 = pie(ax1,class_frac_all); %,explode
for i=1:2:7
p1(i).EdgeColor = 'none';
p1(i+1).FontSize = 30;
end
colormap(superclass_colors(6:9,:))
%title('all intrinsic neurons','FontSize',20)



ax2 = nexttile;
p2=pie(ax2, class_frac_rich);
for i=1:2:7
p2(i).EdgeColor = 'none';
p2(i+1).FontSize = 30;
end
%title('intrinsic rich club neurons','FontSize',20)

ax3 = nexttile;
p3=pie(ax3, class_frac_integrate);
for i=1:2:7
p3(i).EdgeColor = 'none';
p3(i+1).FontSize = 30;
end
%title('integrator neurons','FontSize',20)

ax4 = nexttile;
p4=pie(ax4, class_frac_broadcast);
for i=1:2:7
p4(i).EdgeColor = 'none';
p4(i+1).FontSize = 30;
end
%title('broadcaster neurons','FontSize',20)
% 
% % ax5 = nexttile;
% % p5=pie(ax5, [nt_frac_balanced,1-sum(nt_frac_balanced)]);
% % for i=1:2:13
% % p5(i).EdgeColor = 'none';
% % end
% % title('large balanced neurons')
% 
% Create legend

lgd = legend({'optic', 'vis. proj.', 'vis. cent. ','central'},'FontSize',30);
lgd.Layout.Tile = 'east';

set(gcf,'position',[100,100,2500,600])

%set(gcf,)

set(gcf,'renderer','Painters')

saveas(gcf,strcat(output_dir,'\richclub_class'),'epsc')
saveas(gcf,strcat(output_dir,'\richclub_class'),'png')




%% track input and output sides, neuropil number
% find which entries in all neurons are L and R
% currently manual overwrite L/R!!!!
left_regions = strfind(all_regions,'_R');
left_idcs = find(~cellfun(@isempty,left_regions));
right_regions = strfind(all_regions,'_L');
right_idcs = find(~cellfun(@isempty,right_regions));

% track input L, input both, input R for all neurons
% same for outputs
% track total neuropil number inputs and outputs
%in_neuron_frac(:,left_idcs)

input_np_num = zeros(size(in_degree));
output_np_num = zeros(size(in_degree));

for i = 1:length(neuron_IDs)
    input_np_num(i) = nnz(in_neuron_frac(i,:));
    output_np_num(i) = nnz(out_neuron_frac(i,:));
    
    left_in_temp = nnz(in_neuron_frac(i,left_idcs));
    right_in_temp = nnz(in_neuron_frac(i,right_idcs));

    if left_in_temp ~= 0 && right_in_temp == 0
          input_side{i} = 'left';
    elseif left_in_temp == 0 && right_in_temp ~= 0
          input_side{i} = 'right';
    elseif left_in_temp ~= 0 && right_in_temp ~= 0
           input_side{i} = 'both';
    else
            input_side{i} = 'central';
    end

    left_out_temp = nnz(out_neuron_frac(i,left_idcs));
    right_out_temp = nnz(out_neuron_frac(i,right_idcs));

    if left_out_temp ~= 0 && right_out_temp == 0
          output_side{i} = 'left';
    elseif left_out_temp == 0 && right_out_temp ~= 0
          output_side{i} = 'right';
    elseif left_out_temp ~= 0 && right_out_temp ~= 0
           output_side{i} = 'both';
    else
        output_side{i} = 'central';
    end
end


%% do rich club neurons cross hemispheres more than others? input sides

hemi_in_num_all = zeros(1,4);
hemi_in_num_rich = zeros(1,4);
hemi_in_num_integrate = zeros(1,4);
hemi_in_num_broadcast = zeros(1,4);
hemi_in_num_balanced = zeros(1,4);

hemi_in_frac_all = zeros(1,4);
hemi_in_frac_rich = zeros(1,4);
hemi_in_frac_integrate = zeros(1,4);
hemi_in_frac_broadcast = zeros(1,4);
hemi_in_frac_balanced = zeros(1,4);

sides = {'left','right','both','central'};

for i = 1:4

hemi_in_num_all(i) = sum(strcmpi(input_side,sides{i}));
hemi_in_frac_all(i) = hemi_in_num_all(i)/length(neuron_IDs);

hemi_in_num_rich(i) = sum(strcmpi(input_side(rich_club_filter),sides{i}));
hemi_in_frac_rich(i) = hemi_in_num_rich(i)/sum(rich_club_filter);

hemi_in_num_integrate(i) = sum(strcmpi(input_side(integrate_filter),sides{i}));
hemi_in_frac_integrate(i) = hemi_in_num_integrate(i)/sum(integrate_filter);
 
hemi_in_num_broadcast(i) = sum(strcmpi(input_side(broadcast_filter),sides{i}));
hemi_in_frac_broadcast(i) = hemi_in_num_broadcast(i)/sum(broadcast_filter);
% 
hemi_in_num_balanced(i) = sum(strcmpi(input_side(intrinsic_balanced_filter),sides{i}));
hemi_in_frac_balanced(i) = hemi_in_num_balanced(i)/sum(intrinsic_balanced_filter);

end

t = tiledlayout(1,4,'TileSpacing','compact');

explode = [0 0 0 1 1 1];
% Create pie charts
ax1 = nexttile;
p1 = pie(ax1,hemi_in_frac_all); %,explode
for i=1:2:7
p1(i).EdgeColor = 'none';
p1(i+1).FontSize = 30;
end
colormap([0.8 0.3 0.3;...
    0.2 0.7 0.7;...
    0.4 0.8 0.3;...
    0.9 0.9 0.3])
%title('inputs: all neurons','FontSize',20)

ax2 = nexttile;
p2=pie(ax2, hemi_in_frac_rich);
for i=1:2:7
p2(i).EdgeColor = 'none';
p2(i+1).FontSize = 30;
end
%title('rich club neurons','FontSize',20)

ax3 = nexttile;
p3=pie(ax3, hemi_in_frac_integrate);
for i=1:2:7
p3(i).EdgeColor = 'none';
p3(i+1).FontSize = 30;
end
%title('integrator neurons','FontSize',20)

ax4 = nexttile;
p4=pie(ax4, hemi_in_frac_broadcast);
for i=1:2:7
p4(i).EdgeColor = 'none';
p4(i+1).FontSize = 30;
end
%title('broadcaster neurons','FontSize',20)

% ax5 = nexttile;
% p5=pie(ax5, [nt_frac_balanced,1-sum(nt_frac_balanced)]);
% for i=1:2:13
% p5(i).EdgeColor = 'none';
% end
% title('large balanced neurons')

% Create legend
lgd = legend({'left','right','both','central only'},'FontSize',30);
lgd.Layout.Tile = 'east';


% nt fractions for all neurons
set(gcf,'position',[100,100,2500,600])

%set(gcf,)

set(gcf,'renderer','Painters')

saveas(gcf,strcat(output_dir,'\richclub_inputs'),'epsc')
saveas(gcf,strcat(output_dir,'\richclub_inputs'),'png')


%% do rich club neurons cross hemispheres more than others? output sides


hemi_out_num_all = zeros(1,4);
hemi_out_num_rich = zeros(1,4);
hemi_out_num_integrate = zeros(1,4);
hemi_out_num_broadcast = zeros(1,4);
hemi_out_num_balanced = zeros(1,4);

hemi_out_frac_all = zeros(1,4);
hemi_out_frac_rich = zeros(1,4);
hemi_out_frac_integrate = zeros(1,4);
hemi_out_frac_broadcast = zeros(1,4);
hemi_out_frac_balanced = zeros(1,4);

sides = {'left','right','both','central'};

for i = 1:4

hemi_out_num_all(i) = sum(strcmpi(output_side,sides{i}));
hemi_out_frac_all(i) = hemi_out_num_all(i)/length(neuron_IDs);

hemi_out_num_rich(i) = sum(strcmpi(output_side(rich_club_filter),sides{i}));
hemi_out_frac_rich(i) = hemi_out_num_rich(i)/sum(rich_club_filter);

hemi_out_num_integrate(i) = sum(strcmpi(output_side(integrate_filter),sides{i}));
hemi_out_frac_integrate(i) = hemi_out_num_integrate(i)/sum(integrate_filter);
 
hemi_out_num_broadcast(i) = sum(strcmpi(output_side(broadcast_filter),sides{i}));
hemi_out_frac_broadcast(i) = hemi_out_num_broadcast(i)/sum(broadcast_filter);
% 
hemi_out_num_balanced(i) = sum(strcmpi(output_side(intrinsic_balanced_filter),sides{i}));
hemi_out_frac_balanced(i) = hemi_out_num_balanced(i)/sum(intrinsic_balanced_filter);

end

t = tiledlayout(1,4,'TileSpacing','compact');

explode = [0 0 0 1 1 1];
% Create pie charts
ax1 = nexttile;
p1 = pie(ax1,hemi_out_frac_all); %,explode
for i=1:2:7
p1(i).EdgeColor = 'none';
p1(i+1).FontSize = 30;
end
colormap([0.8 0.3 0.3;...
    0.2 0.7 0.7;...
    0.4 0.8 0.3;...
    0.9 0.9 0.3])
%title('outputs: all neurons','FontSize',11)

ax2 = nexttile;
p2=pie(ax2, hemi_out_frac_rich);
for i=1:2:7
p2(i).EdgeColor = 'none';
p2(i+1).FontSize = 30;
end
%title('rich club neurons','FontSize',11)

ax3 = nexttile;
p3=pie(ax3, hemi_out_frac_integrate);
for i=1:2:7
p3(i).EdgeColor = 'none';
p3(i+1).FontSize = 30;
end
%title('integrator neurons','FontSize',11)

ax4 = nexttile;
p4=pie(ax4, hemi_out_frac_broadcast);
for i=1:2:7
p4(i).EdgeColor = 'none';
p4(i+1).FontSize = 30;
end
%title('broadcaster neurons','FontSize',11)

% ax5 = nexttile;
% p5=pie(ax5, [nt_frac_balanced,1-sum(nt_frac_balanced)]);
% for i=1:2:13
% p5(i).EdgeColor = 'none';
% end
% title('large balanced neurons')

% Create legend
lgd = legend({'left','right','both','central only'},'FontSize',30);
lgd.Layout.Tile = 'east';


% nt fractions for all neurons
set(gcf,'position',[100,100,2500,600])

%set(gcf,)

set(gcf,'renderer','Painters')

saveas(gcf,strcat(output_dir,'\richclub_outputs'),'epsc')
saveas(gcf,strcat(output_dir,'\richclub_outputs'),'png')

% print examples of rich club neurons that are balanced, integrators, or
% broadcasters (6, plus more in supplement?)
%try visualizing all broadcasters and all integrators
