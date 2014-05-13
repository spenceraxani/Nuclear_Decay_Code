import numpy as np
from scipy.signal import lombscargle
import matplotlib.pyplot as plt


timepoints = np.loadtxt('timesequence', usecols=(0,), delimiter=",")

# Coalesce the repeated times into the `times` and `counts` arrays.
times, inv = np.unique(timepoints, return_inverse=True)
counts = np.bincount(inv)
print(counts)
# Check the sample spacing--is the sampling uniform?
dt = np.diff(times)
dts, inv = np.unique(dt, return_inverse=True)
dt_counts = np.bincount(inv)
#print dts
# Inspection of `dts` shows that the smallest dt is about 1.7, and there
# are many near multiples of 1.74,  but the sampling is not uniform,
# so we'll analyze the spectrum using lombscargle.


# First remove the mean from the data.  This is not essential; it just
# removes the large value at the 0 frequency that we don't care about.
counts0 = counts - counts.mean()

# Minimum interarrival time.
dt_min = dt.min()

# Total time span.
trange = times[-1] - times[0]

# --- Lomb-Scargle calculation ---
num_ls_freqs = 16000
ls_min_freq = 1.0 / trange
ls_max_freq = 1.0 / dt_min
freqs = np.linspace(ls_min_freq, ls_max_freq, num_ls_freqs)
ls_pgram = lombscargle(times, counts0, 2*np.pi*freqs)

ls_peak_k = ls_pgram.argmax()
ls_peak_freq = freqs[ls_peak_k]
print "ls_peak_freq  =", ls_peak_freq


# --- FFT calculation of the binned data ---
# Assume the Lomb-Scargle calculation gave a good estimate
# of the fundamental frequency.  Use a bin size for the histogram
# of timepoints that oversamples that period by m.
m = 8
nbins = int(m * ls_peak_freq * trange + 0.5)
hist, bin_edges = np.histogram(timepoints, bins=nbins, density=True)
delta = bin_edges[1] - bin_edges[0]

fft_coeffs = np.fft.rfft(hist - hist.mean())
fft_freqs = np.fft.fftfreq(hist.size, d=delta)[:fft_coeffs.size]
# Hack to handle the case when hist.size is even.  `fftfreq` puts
# -nyquist where we want +nyquist.
fft_freqs[-1] = abs(fft_freqs[-1])

fft_peak_k = np.abs(fft_coeffs).argmax()
fft_peak_freq = fft_freqs[fft_peak_k]
print "fft_peak_freq =", fft_peak_freq


# --- Lomb-Scargle plot ---
plt.figure(1)
plt.clf()
plt.plot(freqs, ls_pgram)
plt.title('Spectrum computed by Lomb-Scargle')
plt.annotate("%6.4f Hz" % ls_peak_freq, 
             xy=(ls_peak_freq, ls_pgram[ls_peak_k]),
             xytext=(10, -10), textcoords='offset points')
plt.xlabel('Frequency (Hz)')
plt.grid(True)


# --- FFT plot ---
plt.figure(2)
plt.clf()
plt.plot(fft_freqs, np.abs(fft_coeffs)**2)
plt.annotate("%6.4f Hz" % fft_peak_freq,
             xy=(fft_peak_freq, np.abs(fft_coeffs[fft_peak_k])**2),
             xytext=(10, -10), textcoords='offset points')
plt.title("Spectrum computed by FFT")
plt.grid(True)

plt.show()