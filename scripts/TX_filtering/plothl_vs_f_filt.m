#!/usr/bin/octave -q

pkg load signal

# filter bank : f_max, filter_type, order, cutoff_freq [, [stopband_attenuation] [passband ripple]]
#   f_max : the LPF will be used up to this frequency
#   filter_type : 'Butterworth', 'ChebyshevII', 'file' (can be .s2p) 
#   order : filter order
#   cutoff_freq : frequency where the filter response is 3 dB down
#   stopband_attenuation : ChebyshevII only, minimun stopband attenuation
#   passband_ripple : ChebyshevI only, passband ripple
# filters need to be in ascending f_max order
#filters = {
#  {2.5e6, 'Butterworth', 5, 3e6} # 160 m
#  {4.5e6, 'Butterworth', 5, 5.5e6} # 80 m
#  {10e6, 'ChebyshevI', 5, 11e6, 0.5} # 40 m / 60 m
#  {17e6, 'Butterworth', 5, 20e6} # 20 m / 30 m
#  {24e6, 'ChebyshevI', 5, 30e6, 0.5} # 15 m / 17 m
#  {31e6, 'ChebyshevII', 5, 37e6, 20} # 10 m / 12 m
#};

# bank similar to G11 TX LPFs
#filters = {
#  {2.5e6, 'ChebyshevII', 5, 2.65e6, 45}
#  {4.7e6, 'ChebyshevII', 5, 5.1e6, 45}
#  {8e6, 'ChebyshevII', 5, 8.8e6, 45}
#  {12e6, 'ChebyshevII', 5, 13.5e6, 45}
#  {19e6, 'ChebyshevII', 5, 22.3e6, 45}
#  {34e6, 'ChebyshevII', 5, 34.5e6, 45}
#};

# G11 TX LPFs bank
filters = {
  {2.5e6, 'file', 'G11_bank_1.txt'}
  {4.7e6, 'file', 'G11_bank_2.txt'}
  {8e6, 'file', 'G11_bank_3.txt'}
  {12e6, 'file', 'G11_bank_4.txt'}
  {19e6, 'file', 'G11_bank_5.txt'}
  {34e6, 'file', 'G11_bank_6.s2p'}
};

# sampling frequency (used to compute spurs location)
global fs = 79.872e6;

# define the HF bands limits - used for annotating the graphs
function plot_HF_bands(Ypos)
  HF_bands = {
    {'160 m', 1.8, 2.0};
    {'80 m', 3.5, 4.0};
    {'60 m', 5.3305, 5.4035};
    {'40 m', 7.0, 7.3};
    {'30 m', 10.1, 10.15};
    {'20 m', 14.0, 14.35};
    {'17 m', 18.068, 18.168};
    {'15 m', 21.0, 21.45};
    {'12 m', 24.89, 24.99};
    {'10 m', 28.0, 29.7}
  };
  hold on;
  for bidx = 1:length(HF_bands)
    bdata = HF_bands{bidx};
    bname = bdata{1};
    bstart = bdata{2}*1e6; # in Hz
    bstop = bdata{3}*1e6; # in Hz
    plot([bstart, bstop], [Ypos, Ypos], 'g', 'LineWidth',4);
    text(0.5*(bstart+bstop), Ypos, [' ' bname], 'horizontalalignment', 'left', 'verticalalignment', 'middle', 'fontsize', 8, 'rotation', 90);
  endfor
  hold off;
endfunction

# draw a limit line at 'level' dB on a graph
function draw_limit_line(level)
  hold on;
  plot(xlim(), [level, level], 'r', 'LineWidth',3);
  hold off;
endfunction

# apply a Butterworth filter
function out = filt_Butterworth(in, frqs, n, f0)
  # ideal Butterworth filter
  # n : order
  # f0 : cutoff frequency
  Amax = 50; # maximum final attenuation
  [b, a] = butter (n, f0 ,'s');
  G_dB = 20*log10(abs(freqs(b, a, frqs)));
  G_dB = merge(G_dB < -Amax, -Amax, G_dB);
  out = in + G_dB;
endfunction

# apply a type I (equiripple in the passband) Chebyshev filter
function out = filt_ChebyshevI(in, frqs, n, f0, rp)
  # ideal (type I) Chebyshev filter
  # n : order
  # f0 : cutoff frequency
  # rp : passband ripple
  
  Amax = 50; # maximum final attenuation  
  
  # compute epsilon from the passband ripple
  eps = sqrt(10^(rp/10)-1);
  # compute the the ripple cutoff frequency from the 3 dB cutoff frequency
  f0 = f0 / cosh(acosh(1/eps)/n);    
  [b, a] = cheby1(n, rp, f0, 's');
  
  G_dB = 20*log10(abs(freqs(b, a, frqs)));
  G_dB = merge(G_dB < -Amax, -Amax, G_dB);
  out = in + G_dB;
endfunction

# apply a type II (equiripple in the stopband) Chebyshev filter
function out = filt_ChebyshevII(in, frqs, n, f0, Rs)
  # ideal inverse (type II) Chebyshev filter
  # n : order
  # f0 : cutoff frequency
  # Rs : stopband attenuation
  
  Amax = 50; # maximum final attenuation  
  # compute epsilon from the stopband attenuation
  eps = 1 / sqrt(10^(Rs/10)-1);
  # compute the the stopband cutoff frequency from the 3 dB cutoff frequency
  f0 = f0 * cosh(acosh(1/eps)/n);
  [b, a] = cheby2(n, Rs, f0, 's');
  
  G_dB = 20*log10(abs(freqs(b, a, frqs)));
  G_dB = merge(G_dB < -Amax, -Amax, G_dB);
  out = in + G_dB;
endfunction

# apply a filter function from file
#   if file ends in .s2p it's read as a Touchstone file
#   (only S21 will be used)
#   otherwise it is assumed to contain two data columns
#   containg the frequency and the filter response in dB
function out = filt_file(in, freqs, fname)
  Amax = 50; # maximum final attenuation  
  [dir, name, fext] = fileparts(fname);
  if strcmp(fext, '.s2p')
   [ff, fr] = SXPParse(fname, '/dev/null'); # suppress output text
   ff = ff';
   fr = squeeze(fr(2,1,:)); # extract S21
   fr = 20 *log10(abs(fr)); # in dB
  else
    fdata = load(fname);
    ff = fdata(:,1);
    fr = fdata(:,2);
  endif
  
  G_dB = interp1 (ff, fr, freqs, 'spline');
  G_dB = merge(G_dB < -Amax, -Amax, G_dB);
  out = in + G_dB;
endfunction

# filter the TX spectrum using the provided filter bank
function ftxdata = LPF(txdata, filters)
  global fs;
  idxmin = 1;
  for fidx = 1:length(filters)
    idxmax = lookup(txdata(:,1), filters{fidx}{1});
  
    switch (filters{fidx}{2})
      case "Butterworth"
        BW = @(in, f) filt_Butterworth(in, f, filters{fidx}{3}, filters{fidx}{4});
      case "ChebyshevI"
        BW = @(in, f) filt_ChebyshevI(in, f, filters{fidx}{3}, filters{fidx}{4}, filters{fidx}{5});        
      case "ChebyshevII"
        BW = @(in, f) filt_ChebyshevII(in, f, filters{fidx}{3}, filters{fidx}{4}, filters{fidx}{5});
      case "file"
        BW = @(in, f) filt_file(in, f, filters{fidx}{3});
    endswitch
  
    ftxdata(idxmin:idxmax,1) = txdata(idxmin:idxmax,1); # TX frequency
    frqs = txdata(idxmin:idxmax,1); # TX frequency
    ftxdata(idxmin:idxmax,2) = BW(txdata(idxmin:idxmax,2), frqs);
    ftxdata(idxmin:idxmax,3) = BW(txdata(idxmin:idxmax,3), 2*frqs);
    ftxdata(idxmin:idxmax,4) = BW(txdata(idxmin:idxmax,4), 3*frqs);
    ftxdata(idxmin:idxmax,5) = BW(txdata(idxmin:idxmax,5), 4*frqs);
    ftxdata(idxmin:idxmax,6) = BW(txdata(idxmin:idxmax,6), 5*frqs);

    ftxdata(idxmin:idxmax,7) = BW(txdata(idxmin:idxmax,7), fs-frqs);
    ftxdata(idxmin:idxmax,8) = BW(txdata(idxmin:idxmax,8), fs+frqs);

    ftxdata(idxmin:idxmax,9) = BW(txdata(idxmin:idxmax,9), 2*fs-5*frqs);
    ftxdata(idxmin:idxmax,10) = BW(txdata(idxmin:idxmax,10), 2*fs-4*frqs);
    ftxdata(idxmin:idxmax,11) = BW(txdata(idxmin:idxmax,11), 2*fs-3*frqs);
    ftxdata(idxmin:idxmax,12) = BW(txdata(idxmin:idxmax,12), 2*fs-2*frqs);
    ftxdata(idxmin:idxmax,13) = BW(txdata(idxmin:idxmax,13), 2*fs-frqs);

    ftxdata(idxmin:idxmax,14) = BW(txdata(idxmin:idxmax,14), 2*fs);

    ftxdata(idxmin:idxmax,15) = BW(txdata(idxmin:idxmax,9), 2*fs+frqs);
    ftxdata(idxmin:idxmax,16) = BW(txdata(idxmin:idxmax,10), 2*fs+2*frqs);
    ftxdata(idxmin:idxmax,17) = BW(txdata(idxmin:idxmax,11), 2*fs+3*frqs);
    ftxdata(idxmin:idxmax,18) = BW(txdata(idxmin:idxmax,12), 2*fs+4*frqs);
    ftxdata(idxmin:idxmax,19) = BW(txdata(idxmin:idxmax,13), 2*fs+5*frqs);
    idxmin = idxmax + 1;
  endfor
endfunction

########################################

graphics_toolkit("gnuplot")

# output graphs basename
datafile = argv(){1};
[pathstr, fbasename, ext] = fileparts(datafile);

# load curve with Pout around 37 dBm (5W)
txdata = load(datafile);
# data format is:
#   spot_level
#   f_{tx}, 2*f_{tx}, 3*f_{tx}, 4*f_{tx}, 5*f_{tx},
#   f_s-f_{tx}, f_s+f_{tx},
#   2*f_s-5*f_{tx}, 2*f_s-4*f_{tx}, 2*f_s-3*f_{tx}, 2*f_s-2*f_{tx}, 2*f_s-f_{tx},
#   2*f_s,
#   2*f_s+f_{tx}, 2*f_s+2*f_{tx}, 2*f_s+3*f_{tx}, 2*f_s+4*f_{tx}, 2*f_s+5*f_{tx}]

# apply LPF bank to the TX signal
ftxdata = LPF(txdata, filters);
frqs = txdata(:,1); # TX frequency

fidx = 1; # figures index

# plot relative output power values (dBc) with limit line
figure(fidx++)
plot(frqs, ftxdata(:,3:6)-ftxdata(:,2));
draw_limit_line(-50);
yr = ylim();
plot_HF_bands(yr(1)+0.02*(yr(2)-yr(1))); # place band markers slightly above x axis
grid on;
#ylim([-70 -30]);
xlim([0 31e6]);
#set(gca, 'ytick', [-90:10:10]);
xlabel("TX frequency, Hz");
ylabel('Pout, dBc');
legend('2f_{tx}', '3f_{tx}', '4f_{tx}', '5f_{tx}', 'location', 'southwest');
legend("boxoff");
legend ("location", "northeast")
title('Harmonics');
print([fbasename '_harmonics_vs_frequency_LPF_dBc.png'], '-dpng', '-F:12');

# plot fundamental output power (dBm)
figure(fidx++)
plot(frqs, ftxdata(:,2), 'b;after LPF bank;'); # filtered output
hold on;
plot(frqs, txdata(:,2), 'r--;PA output;'); # unfiltered output
yr = ylim();
plot_HF_bands(yr(1)+0.02*(yr(2)-yr(1))); # place band markers slightly above x axis
grid on;
#ylim([-90 15]);
xlim([0 31e6]);
#set(gca, 'ytick', [-90:10:10]); 
xlabel("TX frequency, Hz");
ylabel('Pout, dBm');
legend('location', 'northwest');
legend("boxoff");
title('Fundamental output power');
print([fbasename '_fundamental_vs_frequency_LPF.png'], '-dpng', '-F:12');

##########
# plot absolute output power values (dBm)
figure(fidx++)
plot(frqs, ftxdata(:,2:6));
yr = ylim();
plot_HF_bands(yr(1)+0.02*(yr(2)-yr(1))); # place band markers slightly above x axis
grid on;
#ylim([-90 15]);
xlim([0 31e6]);
#set(gca, 'ytick', [-90:10:10]); 
xlabel("TX frequency, Hz");
ylabel('Pout, dBm');
legend('f_{tx}', '2f_{tx}', '3f_{tx}', '4f_{tx}', '5f_{tx}', 'location', 'northwest');
legend("boxoff");
title('Harmonics after LPF');
#print([fbasename '_harmonics_vs_frequency.png'], '-dpng', '-F:12');

figure(fidx++)
plot(frqs, ftxdata(:,7:8));
yr = ylim();
plot_HF_bands(yr(1)+0.02*(yr(2)-yr(1))); # place band markers slightly above x axis
grid on;
#ylim([-90 15]);
xlim([0 31e6]);
#set(gca, 'ytick', [-90:10:10]); 
xlabel("TX frequency, Hz");
ylabel('Pout, dBm');
legend('f_s-f_{tx}', 'f_s+f_{tx}', 'location', 'northwest');
legend("boxoff");
title('Images around f_s');
#print([fbasename '_images_fs_vs_frequency.png'], '-dpng', '-F:12');

figure(fidx++)
plot(frqs, ftxdata(:,9:14));
yr = ylim();
plot_HF_bands(yr(1)+0.02*(yr(2)-yr(1))); # place band markers slightly above x axis
grid on;
#ylim([-90 15]);
xlim([0 31e6]);
#set(gca, 'ytick', [-90:10:10]); 
xlabel("TX frequency, Hz");
ylabel('Pout, dBm');
legend('2f_s-5f_{tx}', '2f_s-4f_{tx}', '2f_s-3f_{tx}', '2f_s-2f_{tx}', '2f_s-f_{tx}', '2f_s', 'location', 'northwest');
legend("boxoff");
title('Images below 2f_s');
#print([fbasename '_images_2fs_low_vs_frequency.png'], '-dpng', '-F:12');

figure(fidx++)
plot(frqs, ftxdata(:,14:19));
yr = ylim();
plot_HF_bands(yr(1)+0.02*(yr(2)-yr(1))); # place band markers slightly above x axis
grid on;
#ylim([-90 15]);
xlim([0 31e6]);
#set(gca, 'ytick', [-90:10:10]); 
xlabel("TX frequency, Hz");
ylabel('Pout, dBm');
legend('2f_s', '2f_s+f_{tx}', '2f_s+2f_{tx}', '2f_s+3f_{tx}', '2f_s+4f_{tx}', '2f_s+5f_{tx}', 'location', 'northwest')
legend("boxoff");
title('Images above 2f_s');
#print([fbasename '_images_2fs_high_vs_frequency.png'], '-dpng', '-F:12');

##########
# plot relative output power values (dBc)
figure(fidx++)
plot(frqs, ftxdata(:,3:6)-ftxdata(:,2));
yr = ylim();
plot_HF_bands(yr(1)+0.02*(yr(2)-yr(1))); # place band markers slightly above x axis
grid on;
#ylim([-70 -30]);
xlim([0 31e6]);
#set(gca, 'ytick', [-90:10:10]);
xlabel("TX frequency, Hz");
ylabel('Pout, dBc');
legend('2f_{tx}', '3f_{tx}', '4f_{tx}', '5f_{tx}', 'location', 'southwest');
legend("boxoff");
legend ("location", "northeast")
title('Harmonics');
#print([fbasename '_harmonics_vs_frequency_dBc.png'], '-dpng', '-F:12');

figure(fidx++)
plot(frqs, ftxdata(:,7:8)-ftxdata(:,2));
yr = ylim();
plot_HF_bands(yr(1)+0.02*(yr(2)-yr(1))); # place band markers slightly above x axis
grid on;
#ylim([-60 -40]);
xlim([0 31e6]);
#set(gca, 'ytick', [-90:10:10]); 
xlabel("TX frequency, Hz");
ylabel('Pout, dBc');
legend('f_s-f_{tx}', 'f_s+f_{tx}', 'location', 'northwest');
legend("boxoff");
title('Images around f_s');
#print([fbasename '_images_fs_vs_frequency_dBc.png'], '-dpng', '-F:12');

figure(fidx++)
plot(frqs, ftxdata(:,9:14)-ftxdata(:,2));
yr = ylim();
plot_HF_bands(yr(1)+0.02*(yr(2)-yr(1))); # place band markers slightly above x axis
grid on;
#ylim([-90 0]);
xlim([0 31e6]);
#set(gca, 'ytick', [-90:10:10]); 
xlabel("TX frequency, Hz");
ylabel('Pout, dBc');
legend('2f_s-5f_{tx}', '2f_s-4f_{tx}', '2f_s-3f_{tx}', '2f_s-2f_{tx}', '2f_s-f_{tx}', '2f_s', 'location', 'southwest');
legend("boxoff");
title('Images below 2f_s');
#print([fbasename '_images_2fs_low_vs_frequency_dBc.png'], '-dpng', '-F:12');

figure(fidx++)
plot(frqs, ftxdata(:,14:19)-ftxdata(:,2));
yr = ylim();
plot_HF_bands(yr(1));
grid on;
#ylim([-90 0]);
xlim([0 31e6]);
#set(gca, 'ytick', [-90:10:10]); 
xlabel("TX frequency, Hz");
ylabel('Pout, dBc');
legend('2f_s', '2f_s+f_{tx}', '2f_s+2f_{tx}', '2f_s+3f_{tx}', '2f_s+4f_{tx}', '2f_s+5f_{tx}', 'location', 'southwest')
legend("boxoff");
title('Images above 2f_s');
#print([fbasename '_images_2fs_high_vs_frequency_dBc.png'], '-dpng', '-F:12');
