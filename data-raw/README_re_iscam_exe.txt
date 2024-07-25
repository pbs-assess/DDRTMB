July 24, 2024

NOTES FOR TESTING RTMB OUTPUTS AGAINST ISCAM MPD OUTPUTS

After much ado, Robyn is able to compile iscam again (thanks Chris Grandin).

The 2018 exe, which is the basis for all the P cod updates since 2018, did not report all the objective function components (priors, qvec or pvec).

In iscam.tpl, Robyn added these to the rep file so they could be compared, along with nlvec_dd, with the obj function components from iscam.

Unfortunately, for the moment, all the projection_model components had to be commented out of iscam.tpl as they were causing errors. 
Doesn't matter because we already have the original projections outputted to file.

So there is a new rep and par file in this folder (for which all mpd values exactly match those from the old exe).
However, the new rep file also reports the objective function components: priors, qvec and pvec.

Notes for Robyn to remember what to do if we need to make changes to iscam and compile again:
The new exe was built in C:\GitHub\gfiscam-wsl2-copy on Robyn's local computer. 
1. Unless the updates are specifically for projections, edit the iscam.tpl file that is in C:\GitHub\gfiscam-wsl2-copy\src\main
2. If not, drop the updated iscam.tpl into C:\GitHub\gfiscam-wsl2-copy\src\main, noting you will have to solve compile errors related to the projection_dd model
3. Open a command window in the root directory of C:\GitHub\gfiscam-wsl2-copy then 'make clean' then 'make dist'. 
4. The new exe appears in C:\GitHub\gfiscam-wsl2-copy\build\dist\bin
5. This can then be copied into C:\GitHub\DDRTMB\devs\0_1a_5ABCD_BASE_2020 and run using the batch file.
6. Then put the new rep and par file into data-raw.
NOTE there is a hardwired path to statsLib.h in GLOBALS_SECTION

Note also that Robyn had to build ADMB from source. ADMB_HOME now points to C:GitHub\admb\build\admb.
To build ADMB from source if ADMB updates: 
1. Update the GitHub admb repo on local computer. 
2. Open command window in GitHub\admb root directory. 
3. Type make. Wait ages (about 40 mins).
4. In Environment Variables, make sure ADMB_HOME points to C:GitHub\admb\build\admb
5. In Environment Variables, make sure C:\GitHub\admb\build\admb\bin is on the Path
6. [To update Environment Variables, go to C:\Windows\System32. Right click on SystemPropertiesAdvanced.exe. Open with iprivilege management.
[probably need to delete the existing build folder]
BUT hopefully won't need to do it again

