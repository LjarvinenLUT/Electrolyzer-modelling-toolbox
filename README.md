### GENERAL DESCRIPTION

This toolbox is a Matlab tool for electrolysis modelling, for both PEM and 
alkaline systems. Main functionality is in parametrization of the UI curve
based on measured data. This toolbox aims to simplify the process and 
provide easy-to-use commands for quickly determining the UI curve parameter 
values and their margins of uncertainty.

---

### INSTALLATION

MATLAB needs to be installed on the computer for the user to be able to
install the electrolyzer modelling toolbox. Following steps install the
toolbox for that matlab installation:

1. Download the provided electrolyzerModellingToolbox.mltbx file from the
   releases page.
    - Double click on the .mltbx file to initiate installation.
2. Once the installation completes the toolbox functions are usable in all
   matlab projects

---

### AUTHORS AND CONTACTS

Pietari Puranen
Junior researcher at the Lappeenranta-Lahti University of Technology LUT
Email: Pietari.Puranen@lut.fi

Lauri Järvinen
Junior researcher at the Lappeenranta-Lahti University of Technology LUT
Email: Lauri.Jarvinen@lut.fi

---

### AKNOWLEDGEMENTS

Marchov Chain Monte Carlo (MCMC) tools provided by Marko Laine have been 
used for the determination of uncertainty in curve fitting when using the 
particleswarm optimization method. The licensing for the tools can be found
under ./Utils/mcmcstat/LICENSE.txt

The Academy of Finland is acknowledged for the main financial support of
the Research of Power Quality Effect on Water Electrolyzer Operation (POQELYZER)
project. The project partners are also thanked for their financial contributions.
