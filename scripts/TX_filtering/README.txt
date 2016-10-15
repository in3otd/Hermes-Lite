Tested with GNU Octave Version 4.0.3

The script need need the signal package, which needs the control package.
To install these packages:
  pkg install -forge control
  pkg install -forge signal

If autoload was not selected during the installation, you may need to load the signal package explicitly:
  pkg load signal

To read S-parameters files, the SXPParse.m function  from the sbox Toolbox by Tudor Dima is needed. Ask the author for a copy or dowload it from this unofficial repository https://github.com/taoyilee/HFSS_API/blob/master/sbox/SXPParse.m

