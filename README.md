# Electrolysis modelling toolbox
10/2021


VERSION

Beta 0.7


COPYRIGHT AND LICENSING

All rights reserved!
Redistribution of the Beta version not allowed!



GENERAL DESCRIPTION

This toolbox is a Matlab tool for electrolysis modelling, for both PEM and 
alkaline systems. Main functionality is in parametrization of the UI curve
based on measured data. This toolbox aims to simplify the process and 
provide easy-to-use commands for quickly determining the UI curve parameter 
values and their margins of uncertainty.



INSTALLATION

There are three options for preparing this toolbox ready for usage in any 
Matlab project:

1. If you have admin rights and the file electrolyzerModellingToolbox.mltbx
    - Double click on the .mltbx file to initiate installation.
2. If you are allowed to modify the Matlab search path defined in pathdef.m
	- Open the directory of the toolbox in Matlab
    - Run script installEModel.m to add all the necessary folders to the
        search path and save it in pathdef.m.
3. If you are not allowed to modify permanently the Matlab search path
    - Copy the script startup.m from the toolbox to any project you want to
        use it.
    - Replace the value of the variable toolpath to include the full path
        to the toolbox.
    - When matlab is opened in the project directory, the startup.m script
        is run automatically. Alternatively, if the project directory is 
        switched from another one while Matlab is running, you can manually
        run startup.m.



AUTHORS AND CONTACTS

Pietari Puranen
Junior researcher at the Lappeenranta-Lahti University of Technology LUT
Email: Pietari.Puranen@lut.fi

Lauri Järvinen
Junior researcher at the Lappeenranta-Lahti University of Technology LUT
Email: Lauri.Jarvinen@lut.fi



AKNOWLEDGEMENTS

Marchov Chain Monte Carlo (MCMC) tools provided by Marko Laine have been 
used for the determination of uncertainty in curve fitting when using the 
particleswarm optimization method. The licensing for the tools can be found
under ./Utils/mcmcstat/LICENSE.txt

