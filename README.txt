These directories contain the raw data and behavioral analyses for the paper:

            Michael Jigo, Marisa Carrasco; Attention alters spatial resolution by modulating second-order processing. 
            Journal of Vision 2018;18(7):2. doi: https://doi.org/10.1167/18.7.2.

-----------------------------------------
---------- Recreate figures -------------
-----------------------------------------

All necessary functions are located in the "analyses" directory

To recreate a figure run the command within the square brackets:
   - Figure 2    [figure2]
   - Figure 3    [figure3]
   - Figure 4    [figure4]



-----------------------------------------
------------- Data analyses -------------
-----------------------------------------

All necessary functions are located in the "analyses" directory for the respective experiment.

To re-run behavioral data analyses, run the command "do_analyses".



-----------------------------------------
---------- Directory structure ----------
-----------------------------------------

Directories are structured as follows:
Key: |--  = directory branch
     |--> = description of files contained in directory

|-- analyses
|   |-> functions that perform analyses and recreate figures in the paper.
|   |-- helperfun
|   |   |--> functions (.m) that facilitate the main analyses
|   |    
|-- data
|   |-- subdirectories for each subject (e.g., AJF)
|   |   |-- endo
|   |   |   |-- cue
|   |   |   |   |--> raw stimfiles from main cueing experiment
|   |   |   |-- prac
|   |   |   |   |--> raw stimfiles from practice session
|   |   |   |-- thresh
|   |   |   |   |--> raw stimfiles from thresholding session
|   |   |-- exo (organization identical to endo)
|   |    
|   |-- bootstrap
|   |   |--> directories containing bootstrapped samples for each subject in exo and endo conditions
|   |   |--> outputted from bootstrap_observer_data
|   |    
|   |-- parsed_data
|   |   |--> analyzed data (i.e., the output from parse_observer_data)
|   |    
|-- stim
|   |-> functions that create and display experimental protocol
|-- manuscript
|   |-> final manuscript files (and drafts)
