# LiLS-PSD-analysis
The codes that analyze the LiLS PSD data for PROSPECT Experiment

1. ProcessOneRun_CsData_20180213.py

processes the runs with both Cs137 and AmBe sources. Input file name is needed directly in the script.

2. ProcessOneRun_AmBeData_20180213_batch.py 

processes the runs with only the AmBe source. Input file name is set as parameter, but need to check the folder name inside the code.

3. drawCal.C

processes the Cs137 runs and prepares the energy scales for later analysis.

4. fitPSD.C

does the fit to the PSD distributions and extract the quantities for the QA/QC.

5. the rest of the C scripts

plot histograms and results. Will need to write codes separately if I need to plot the data in different ways. 
