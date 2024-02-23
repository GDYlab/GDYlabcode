% Code used is in-built and custom MATLAB routines (tested using MATLAB versions R2015a or R2020a) 
% No further dependency, installation or special hardware is required
% Details of inputs and outputs provided in each script

% List of codes:

% Track end similarity, Shape of population neural trajectory and deviation from linearity of the trajectory
% For Main Figures 3a-d and 7c-d; and Supplementary Figures 4a-c and 7d,f-h

track_end_similarity_computation.m
shape_of_popneural_trajectory.m
deviation_from_linearity_assessment.m

% Cell-pair correlations across all sleep frames
% For Figures 4a-c and 8b left

cellpaircorrelations_allframes.m

% Ensemble frame-pair correlations and clustering of sleep frame pairs
% For Main Figures 4d-g, 6b, 6e-f , 8b (right), 8c and Supplementary Figures 5, 6e

ensembleframepaircorrelations.m
framepairsimilarity_normalization_by_shuffles.m
clustering_of_frames.m
silhouette_calc_bestcliteration.m
plotting_detecting_clustermats.m
withinacrossclustercorrs.m

% Similarity of population rates across tracks
% For Main Figure 5c

population_dynamics_across_tracks.m

% Specificity of frames in individual clusters versus tracks
% For Main Figures 6b and 8d

specificity_of_clusters_for_tracks.m

% Applying two simultaneous criteria for sequence assessment
% For Supplementary Figures 2c-d, 6d,6f and 8a,8d

wcorrs_and_maxjump_comparisons.m

% Remaining codes used in the manuscript are commonly-used routines in the field based on prior work cited in the paper
