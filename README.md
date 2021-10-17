# READ ME FOR CODE TO ANALYZE CELL VELOCITIES IN COLLECTIVE CELL MIGRATION EXPERIMENTS

*Written by Notbohm Research Group, University of Wisconsin-Madison.* https://notbohm.ep.wisc.edu 

This document gives information about common code to plot data from experiments
in collective cell migration. The scripts in this repository are customized for analyzing
a scratch wound assay, though they will work for other assays as well (perhaps with minor modifications).

**IMPORTANT: These scripts give examples of how to plot outputs. You may have to modify them
or write your own scripts to analyze the quantity of interest in your experiment.** 

## Starting Point

Begin by computing displacements, e.g., using FIDIC ( https://github.com/jknotbohm/FIDIC ).

## Text Files Required

In addition to running FIDIC and having the the output mat file in your directory of interest,
the scripts in this repository require three txt files:
- **ExperimentalSettings.txt**: Only the first and 7th lines are required for analyzing/plotting
cell velocities. First line: pixel size (in m). 7th line: 1 or 2, where 1 identifies that the 
geometry of the region containing cells is a strip and 2 identifies that it's an island. For a 
scratch wound assay, use 1.
- **TimeIncrement.txt**: Time between images (min)
- **time_points_start_end.txt**: The time points at which to start and end the analysis. 
Enter an empty array [] to begin at first time point and/or end at last time point. (no units)

Examples are included in this repository with "_example" appended at the
end. See the notes in these files for more information.

## List of Files to Run

After computing cell displacements over time using FIDIC (link to repository above), run the 
following script

``find_boundary_scratch_single``: Identifies the boundary of the cell island and saves an 
image identifying regions with and without cells. This script identifies the cell boundaries
based on a threshold and other settings set by the user. You may have to iterate to get good
results. At a minimum, run for the first and last time points.

After running ``find_boundary_scratch_single``, it is possible to run ``wound_areas``.

Run ``compute_cell_trajectories``, which interpolates the displacement data to stitch 
together approximate trajectories traveled by different cells. See the documentation at the start
of the code for more information.

Now it is possible to run the remaining m files, 
``plot_cell_trajectories``,
``plot_MSD``,
``plot_PathRatio``,
``plot_D2min``,
``vel_autocorr_nogrid``, and
``plot_directionality``.

## Additional Subfunctions

``make_fig.m``: Subfunction used to make figure window with desired size and settings.

``autocorr_nogrid.m`` and ``autocorr_angle_nogrid.m``: functions called by
vel_autocorr_nogrid.m.

Some scripts may require other m files available from the Matlab
file repository. See comments of each script for more information.

## Running as a Batch

``find_boundary_scratch_single`` typically has to be run manually, but for all other m files,
I run as a batch. Each function name is called in a for loop, which loops through
different directories, with each directory containing one data set, e.g.,

>`for 1:length(dirlist)`   
>`cd(dirlist{k});`   
>`wound_areas;`   
>`compute_cell_trajectories;`   
>`plot_cell_trajectories;`   
>`plot_MSD;`   
>`plot_PathRatio;`   
>`plot_D2min;`   
>`vel_autocorr_nogrid;`   
>`plot_directionality;`   
>`end`



## Matlab Version

This software was tested on Matlab version R2019b, though other recent versions should work. 
Note that versions earlier than 2014b likely will **not** work, because graphics functions changed
significantly in that version.