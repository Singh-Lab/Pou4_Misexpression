function [p_val_list, sp_corr_list, a_list, b_list, c_list, p_list, wt_list, qual_cell_list, unqual_cell_list, weird_cells, av_mat, oe_vec] = final_nmf_pipeline(av_vec_csv, oe_csv)

% Inputs:
% 1. av_vec_csv ("t_av_all_wt.csv") 
%    This is the csv file that has the matrix of Average WT Cell Type expression vectors. 
%    It has dims 248 (DEGs) x 5 (4 Sensory Cell Types and WT Epidermis)
%    
%    Note: 248 corresponds to 50 DEGs of all WT Cell Types, except for
%          CESNs, which only have 48 DEGs.

% 2. oe_csv (e.g. "oe0_data.csv")
%    It has the expression vectors of cells in an OE Subcluster. 
%    e.g. oe_vec from "oe0_data.csv" has dims 248 x 69, because oe0
%    subcluster has 69 cells.

% Outputs:
% 1. p_val_list = list of p_values < 0.05 from spearman correlation between
%    expression vector of a single oe cell and its solved expression vector. 
%    (Methods) 

% 2. sp_corr_list = list of spearman correlation values. 

% 3. *_list = list of solved WT Cell Type Coefficients.
%             ("Cell-type Specific Solved Coefficients" in Methods)
%    a corresponds to aATEN, b to BTN, c to CESN, p to PSC, and wt
%    to WT Epidermis.

% 4. qual_cell_list = list of cells that had significant spearman
%    correlations.
% 5. unqual_cell_list = complement of above. 
% 6. weird_cells = cells that had a spearman correlation > 1.
%    These were present only in subcluster 4 (untransformed WT Epi cells) 
    
av_mat = (csvread(av_vec_csv));

oe_vec = csvread(oe_csv);

    p_val_list = [];
    sp_corr_list = [];
    a_list = [];
    b_list = [];
    c_list = [];
    p_list = [];
    wt_list = [];
    qual_cell_list = [];
    unqual_cell_list = [];
    % for oe 4 
    weird_cells = [];

% % Do for all cells. 
num_cells = size(oe_vec);
num_cells = num_cells(2);

for i = 1:num_cells
    oe_vec_2 = oe_vec(1:248, i);
    
    oe_fit = lsqnonneg(av_mat, oe_vec_2);
    [oe_corr, oe_pval] = corr(oe_vec_2, av_mat*oe_fit, 'Type', 'Spearman');
    
    if oe_pval < 0.05
        p_val_list = [p_val_list oe_pval];
        sp_corr_list = [sp_corr_list oe_corr];
        a_list = [a_list oe_fit(1)];
        b_list = [b_list oe_fit(2)];
        c_list = [c_list oe_fit(3)];
        p_list = [p_list oe_fit(4)];
        wt_list = [wt_list oe_fit(5)];
        qual_cell_list = [qual_cell_list i];
    else
        unqual_cell_list = [unqual_cell_list i];
    end
    % Relevant only for oe 4.
    if oe_fit(5) > 1
        weird_cells = [weird_cells i];
    end
end
