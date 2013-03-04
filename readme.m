%% README
% A tutorial on how to use p53Cinema and a discussion of its pitfalls.
%% *Setting up the microscope*
% I have found that the complexity of code rapidly increases when a program
% is made to accommodate many different options, settings, or preferences.
% Perhaps, the most frustrating part of this complexity is that almost all
% of the extra code has nothing to do with the core code that was written
% to solve a particular problem in the first place. This code, hereto
% referred to as bloat code, can mask the true purpose of a code and make
% its interpretation difficult to a programmer trying to read it.
% Additionally, if the bloat code was not written well it can be very
% difficult to modify or extract the core code to have it repurposed.
%
% On the other hand, if an operator's preferences tend to change, perhaps
% due to conditions that vary from experiment to experiment, then having a
% versatile program is useful. Also, if the program is to be distributed to
% a large number of people that together have a wide range of preferences
% then having a program that offers flexibility is essential. When is it
% worth the investment in bloat code? For p53Cinema it was decided that the
% bloat code should be kept to a minimum due to the fact that the desired
% tasks a software should perform are apt to evolve relatively quickly and
% members of the Lahav Lab are all proficient-enough programmers that any
% trivial adjustment can be directly made in the code.
%
% One way to minimize the neccessity of bloat code is to minimize the
% variety of inputs that p53Cinema will handle. Performing experiments with
% a set of standard microscope settings is one way this minimization can be
% accomplished. Standardization can also improve the repeatability of
% experiments and edify new lab members who are learning microscopy for the
% first time. The overall goal of the following settings is to reduce the
% amount of exposure necessary to fill the dynamic range of the camera.
% Therefore, please use the following microscope settings when performing
% time-lapse fluorescence microscopy:
%
%% Recommended Microscope Settings
% Please use the following settings for time-lapse fluorescence
% microscopy... 
%
% # *Use a binning of _2_* The binning can help in several ways. It can
% reduce the file size (albeit this is not always an improvement, e.g. this
% is a terrible idea if doing smFISH where detection of the PSF is vital).
% It can improve signal to noise by reducing the read noise (which is
% typically insignificant if not doing single fluorophore work). It reduces
% the length of exposure necessary to fill the dynamic range of the camera;
% this last point is probably the most significant advantage to binning.
% # *Use maximum gain* When a signal is low the gain in a camera can be
% used to increase the dynamic range. In otherwords, the gain can be used
% to amplify a signal so that it stretches over all of the gray levels when
% undergoing analog to digital conversion. The downside to using gain is
% that it introduces shot noise into the signal. The upside is that the
% full dynamic range of the camera is still accessible even when using a
% very low exposure.
% # *Use a 20x objective* A 20x objective provides a large field of view.
% The human cell lines used in the lab are large enough that the nucleus is
% readily identifiable at 20x, including with a binning of 2. A large field
% of view is essential when imaging for more than 24 hours as cells tend to
% move around. Furthermore, the field of view can become very crowed after
% a division or two, so have a larger field view lends itself to plating
% fewer cells; this also means there are more cells to measure per image.
% The disadvantage to 20x is that more out of focus light and
% autofluorescence is measured in the focal volume, which can reduce the
% signal to noise ratio. Also, 20x is fine for measuring a diffuse nuclear
% signal, but not sub-nuclear structures. For instance, when Ketki was
% measuring the dynamics of p53bp51 foci she had to use 60x to get enough
% spatial resolution to identify individual foci.
% # *A note on exposures* The length of exposure for one microscope does
% not translate to the length of exposure on another microscope.
% Differences in absorbance of the light path and intensity of the light
% source mean that the flux of photons varies between the scopes for any
% given wavelength. Also, the exposure for different colors of light is not
% comparable. This is because the intensity of the light source varies
% across the spectrum and the fluorescent molecules being excited can
% differ in number per target and quantum efficiency. Therefore, the
% exposure must be tuned for each channel on each microscope individually.
% # *Another note about exposure*
%% Flat-field and dark-field correction
% Flatfield and darkfield images can be measured directly from the
% microscope. It is highly recommended that these images are collected for
% each experiment to achieve the highest quality of data for
% quantification.