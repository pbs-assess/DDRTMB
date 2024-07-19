July 17, 2024

Robyn is now unable to compile iscam.tpl on Windows and is relying on Chris Grandin to compile it.

The make files depend on a Linux setup and/or ADMB compiled from source, which I can't do without VisualStudio.

Any executable that Chris builds on his Linux system will not run on Windows.

Therefore, the report file in this folder is from a version of devs/iscam.tpl compiled by Chris (with gratitude).

The only difference is this report file includes reporting of all the objective function components: nlvec_dd, priors, qvec and pvec (the original version only reports nlvec_dd).

NOTE: if we do want to do any more test runs (which include priors, qvec and pvec), we will have to run iscam on Chris's server. 

Otherwise, the original exe will still run on Windows, it just won't have those objects reported out.