%% README
% A tutorial on how to use p53Cinema and a discussion of its pitfalls.
%% Setting up the microscope
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
% first time. Therefore, please use the following microscope settings when
% performing time-lapse fluorescence microscopy:
%
%% *Microscope Settings*
% Please use the following settings...
%
% # *Binning* Use a binning of _2_. The binning