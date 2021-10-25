%%PLOT: compare downsampled DVARS/FD and physioparams.
figure()
subplot(4,1,1)
plot(dvars_std, 'y','LineWidth',2)
title("DVARS")
subplot(4,1,2)
plot(FD_power, 'g','LineWidth',2)
title("FD displacement")
subplot(4,1,3)
%plot(phys_data.physio.ons_secs.t,phys_data.physio.ons_secs.fr, 'g','LineWidth',2)
%hold on
plot(resp, 'blue','LineWidth',2)
title("respiratory trace (a.u.), downsampled")
subplot(4,1,4)
%plot(phys_data.physio.ons_secs.t,phys_data.physio.ons_secs.c, 'blue','LineWidth',2)
%hold on
plot(oxy, 'b','LineWidth',2)
title("PPU trace (a.u.), downsampled")

%%PLOT original data compared to downsampled data
figure()
subplot(4,1,1)
plot(phys_data.physio.ons_secs.fr, 'g','LineWidth',2)
title("resp data, original time resolution, filtered, cropped to acquisition window")
subplot(4,1,2)
plot(resp, 'blue','LineWidth',2)
title("resp data, resampled to TR, filtered")
subplot(4,1,3)
plot(phys_data.physio.ons_secs.c, 'blue','LineWidth',2)
title("cardiac data, original time resolution, filtered, cropped to acquisition window")
subplot(4,1,4)
plot(oxy, 'b','LineWidth',2)
title("cardiac data, resampled to TR, filtered")

%%PLOT HRV and RVT and mean FD and DVARS
figure()
subplot(4,1,1)
plot(phys_data.physio.ons_secs.hr(6:end), 'r','LineWidth',2)
title("HR")
subplot(4,1,2)
plot(phys_data.physio.ons_secs.rvt(6:end), 'r','LineWidth',2)
title("RVR")
subplot(4,1,3)
plot(FD_power, 'g','LineWidth',2)
title("FD displacement")
subplot(4,1,4)
plot(dvars_std, 'y','LineWidth',2)
title("DVARS")

