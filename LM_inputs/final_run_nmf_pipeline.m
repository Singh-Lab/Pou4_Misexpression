% Pipeline to calculate Cell-type Specific Solved Coefficients.
% Refer to Methods for detailed description.
% It also visualizes Coefficients at Single-Cell and Population levels.

% Inputs: Average WT Cell Type (av_mat) and OE Expression Vectors (oe_vec)
%         generated from "final_pou4_analysis_pipeline.Rmd".
% Outputs: Various statistics and visualizations. 

% Change oe_csv Input.
av_vec_csv = "t_av_all_wt.csv";
oe_csv = "oe4_data.csv";

% Run pipeline.
[p_val_list, sp_corr_list, a_list, b_list, c_list, p_list, wt_list, qual_cell_list, unqual_cell_list, weird_cells, p_av_mat_init, oe_vec] = final_nmf_pipeline(av_vec_csv, oe_csv);

% Get the the min and max coefficient value among all the lists.
all_list_mat = [a_list, b_list, c_list, p_list, wt_list];
x_min = min(all_list_mat)
x_max = max(all_list_mat)

% weird_cells is only relevant to oe 4.
a_list(weird_cells) = [];
b_list(weird_cells) = [];
c_list(weird_cells) = [];
p_list(weird_cells) = [];
wt_list(weird_cells) = [];

% Get number of qualified cells and total percentage.
num_qual_cells = size(qual_cell_list);
total_cells = size(qual_cell_list)+size(unqual_cell_list);
perc_qual_cells = size(qual_cell_list)/total_cells;

figure
histogram(a_list, 'Normalization','probability')

hold on 
histogram(b_list,  'Normalization','probability')

hold on 
histogram(c_list,  'Normalization','probability')

hold on 
histogram(p_list, 'Normalization','probability')

% hold on 
% histogram(wt_list, 'EdgeAlpha', 0.4, 'Normalization','probability')
% ylim([0,1]);

xlabel('Coefficient Value/Percent Contribution');
ylabel('Percentage');
legend({'y = aATEN','y = BTN', 'y = CESN','y = PSC'});
% legend({'y = aATEN','y = BTN', 'y = CESN','y = PSC', 'y = WT'})
title('Distribution of Sensory Coefficients');
% title('Distribution of Sensory Neuron and WT Coefficients')

% Stacked Barplots.
% Sort based on non-decreasing BTN coefficients/
% "increasing character" of other Sensory Cell Types.  
[sorted_b_list, new_indices] = sort(b_list, 'descend');
sorted_p_list = p_list(new_indices); 
sorted_a_list = a_list(new_indices);
sorted_c_list = c_list(new_indices);

% All.
figure
y = [sorted_b_list; sorted_p_list; sorted_a_list; sorted_c_list]';
bar(y, 'stacked')
xlabel('Number of Cells');
ylabel('Coefficient Value');
legend('BTN', 'PSC', 'aATEN', 'CESN')

% BTN/PSC.
figure
y = [sorted_b_list; sorted_p_list]';
bar(y, 'stacked')
xlabel('Number of Cells');
ylabel('Coefficient Value');
legend('BTN', 'PSC')

% BTN/aATEN.
figure
y = [sorted_b_list; sorted_a_list]';
bar(y, 'stacked')
xlabel('Number of Cells');
ylabel('Coefficient Value');
legend('BTN', 'aATEN')

% BTN/CESN.
figure
y = [sorted_b_list; sorted_c_list]';
bar(y, 'stacked')
xlabel('Number of Cells');
ylabel('Coefficient Value');
legend('BTN', 'CESN')

% Get average delta between BTN coefficients and those of other Sensory
% Cell Types.
avg_btn_aaten_delta = sum(b_list-a_list)/length(b_list)
avg_btn_cesn_delta = sum(b_list-c_list)/length(b_list)
avg_btn_psc_delta = sum(b_list-p_list)/length(b_list)

% Summary Statistic- average of coefficients. 
aaten_mean = mean(a_list)
btn_mean = mean(b_list)
cesn_mean = mean(c_list)
psc_mean= mean(p_list)
