# Note on HTPB Formula and Heat of Formation
The HTPB heat of formation in thermo_spg.inp is likely wrong, but it shouldn't affect calculations too much.

### Taken from an email chain from Dr. Stephen Heister at Purdue in 2020

Rocketfolk:

The heat of formation for HTPB is in error and you should use a value of -250 cal/mol. ...  See exchange with Dave McGrath – a very senior alum who has worked in the solid propellant industry for nearly 4 decades:

Steve

Great question and one that has had lots of discussion over the years.  I did a quick look at several thermochemical program ingredient files (I can’t access the current one for CEA600 but did look at some output files).  We have at least 10 slightly different ways we did it for different grades of HTPB, both with and without incorporation of various typical cure agent ratios. It encompasses different manufacturers as well and generations of data.

All of them have some small amount of nitrogen in the formula and heats of formation vary from slightly positive to significantly negative.

Dave