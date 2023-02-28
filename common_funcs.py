import numpy as np
import ROOT

noisey_channels = np.array([ i for i in range(128) if i%4 in [0, 3]])
clean_channels = np.array([ i for i in range(128) if i%4 in [1, 2]])

attenuation_factor = 0.455462
adc_to_volts_factor = 1.9/16384.0

kernels = np.array([eval(x) for x in open("kernels.txt").readlines()])
kernels = [k/np.sum(np.abs(k)) for k in kernels]

kernel2 = [-0.11557788944719505, 0.025125628140813205, -0.4045226130656374, 0.22613065326549986, 0.15326633165750536, -0.04648241205995873, -2.1469849246232116, -8.247487437185555, -12.341708542713604, -7.153266331658415, -3.4874055415621115, -3.1368221941993397, -2.892676767676676, -1.4266750948163462, -0.9816687737038592, -1.2249683143218135, -0.840966921119616, -0.3613231552162688, -0.5318471337577648, -0.5650510204077364]
#kernel2 = np.array([-3.18, -13.18, -40.18, -33.18, -16.18,  -9.18, -10.18, -10.18, -3.18], dtype=np.float64)
kernel2 -= np.sum(kernel2)/float(len(kernel2))
kernel2 = kernel2/np.sum(np.abs(kernel2))

def calculate_charge(wf, baseline):
    corr_factor = -1*adc_to_volts_factor*(1000.0/50.0)/attenuation_factor
    charge = corr_factor*np.trapz(wf-baseline, dx=2);
    return charge

def rolling_window(a, window):
    pad = np.ones(len(a.shape), dtype=np.int32)
    pad[-1] = window-1
    pad = list(zip(pad, np.zeros(len(a.shape), dtype=np.int32)))
    a = np.pad(a, pad,mode='reflect')
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)

# Values determined for VAR_WINDOW = 25
#def cut_value(var):
#    return 8.425/44.314*var - 0.901

# Values determined for VAR_WINDOW = 6
def cut_value(var):
    return 0.040559*var +0.516458

def find_pulses(wf, threshold=23, reset=6, window=6):
    var_wf = np.var(rolling_window(wf, window), axis=-1)
    ret = []
    start = 0
    end = 0
    bline, conv, conv2 = None, None, None
    #locs = np.where(var_wf > threshold)[0]
    #start_locs = np.concatenate((locs[0:1], locs[1:][np.diff(locs)>1]))

    while np.count_nonzero(var_wf[end:] > threshold):
        if(bline is None):
            bline, _ = find_baseline(wf)
            conv = np.convolve(kernels[0], wf-bline, mode="same") 
            conv2 = np.convolve(kernel2, wf-bline, mode="same")

        start = end + np.where(var_wf[end:] > threshold)[0][0]
        end = len(var_wf)
        try:
            end = start + np.where(var_wf[start:] < reset)[0][0]
            #end = min(end+2, len(var_wf))
        except IndexError:
            pass # end = len(var_wf)

        # Apply cuts to remove noise bursts
        c = np.max(conv[start:end])
        c2 = np.max(conv2[start:end])
        v = np.max(var_wf[start:end])
        cut = cut_value(v)
        if(v > 150 or (c<cut and c2 > 3)):
            # TODO! the "actual_start" should be replaced by something dynamic instead of a fixed number
            actual_start = max(start-3, 0)
            ret.append((actual_start, bline, wf[actual_start:end]))
    return ret

def find_baseline(wf, window=50):
    data_windows = rolling_window(wf, window)
    var_wf = np.var(data_windows, axis=-1)
    mean_wf = np.mean(data_windows, axis=-1)
    min_loc = np.argmin(var_wf)
    return mean_wf[min_loc], np.sqrt(var_wf[min_loc])

def plot_it(wf, dx=1):
    g = ROOT.TGraph()
    for i, v in enumerate(wf):
        g.SetPoint(i, i*dx, v)
    return g

def plot_var(wf, dx=1):
    var_wf = np.var(rolling_window(wf, 25), axis=-1)
    return plot_it(var_wf, dx)

def read_waveforms(fn, include_times=False):
        f = ROOT.TFile(fn, "READ")
        t= f.Get("tree")
        for ev in t:
            wfs = np.reshape(list(ev.FADC), (128, ev.nsample))
            if(include_times):
                times = np.array(list(ev.TimeTag))
                yield wfs, times
            else:
                yield wfs

def correlation_lags(in1_len, in2_len, mode="full"):
   if mode == "full":
       # the output is the full discrete linear convolution
       # of the inputs. (Default)
       lags = np.arange(-in2_len + 1, in1_len)
   elif mode == "same":
       # the output is the same size as `in1`, centered
       # with respect to the 'full' output.
       # calculate the full output
       lags = np.arange(-in2_len + 1, in1_len)
       # determine the midpoint in the full output
       mid = lags.size // 2
       # determine lag_bound to be used with respect
       # to the midpoint
       lag_bound = in1_len // 2
       # calculate lag ranges for even and odd scenarios
       if in1_len % 2 == 0:
           lags = lags[(mid-lag_bound):(mid+lag_bound)]
       else:
           lags = lags[(mid-lag_bound):(mid+lag_bound)+1]
   elif mode == "valid":
       # the output consists only of those elements that do not
       # rely on the zero-padding. In 'valid' mode, either `in1` or `in2`
       # must be at least as large as the other in every dimension.

       # the lag_bound will be either negative or positive
       # this let's us infer how to present the lag range
       lag_bound = in1_len - in2_len
       if lag_bound >= 0:
           lags = np.arange(lag_bound + 1)
       else:
           lags = np.arange(lag_bound, 1)
   return lags
