July 24, 2024

After much ado, Robyn is able to compile iscam again (thanks Chris Grandin).

The 2018 exe, which is the basis for all the P cod updates since 2018, did not report all the objective function components (priors, qvec or pvec).

Robyn added these to the rep file so they can be compared, along with nlvec_dd, with the obj function components from iscam.

Unfortunately, for the moment, all the projection_model components have been commented out as they were causing errors. Will fix later.

So there is a new rep file (for which all mpd values match the old exe), which also reports priors, qvec and pvec.

This was built in C:\GitHub\gfiscam-wsl2-copy. Open a command window in the root directory then 'make clean' then 'make dist'. 
The new exe appears in C:\GitHub\gfiscam-wsl2-copy/build/dist/bin
This can then be copied into C:\GitHub\DDRTMB\devs\0_1a_5ABCD_BASE_2020 and run.
Then put the new rep file into data-raw (with a new file name) and make sure to update the path in devs/plots.r
