import yt
import matplotlib.pyplot as plt
import numpy as np
from sys import argv

def generate_frbs(ds, width):
    frbs = yt.SlicePlot(ds, 'z', [('boxlib', 'conv_actual'), ('boxlib', 'conv_adiabatic'), ('boxlib', 'conv_ledoux')], width=width).frb

    act = frbs.data[('boxlib', 'conv_actual')]
    ad = frbs.data[('boxlib', 'conv_adiabatic')]
    led = frbs.data[('boxlib', 'conv_ledoux')]
    
    schwarz_cri = act - ad
    ledoux_cri = act - led
    
    return act, ad, led, schwarz_cri, ledoux_cri

def plot_criteria(ds, s, l, width):
    plt.figure(2)

    center = [int(act.shape[0]/2), int(act.shape[1]/2)]
    # take a 4 average at least
    # rightward, downward, leftward, upward
    s_1d = np.mean([s[center[0]:, center[1]], s[center[0], center[1]:], np.flip(s[center[0], :center[1]]), np.flip(s[:center[0], center[1]])], axis=0)
    l_1d = np.mean([l[center[0]:, center[1]], l[center[0], center[1]:], np.flip(l[center[0], :center[1]]), np.flip(l[:center[0], center[1]])], axis=0)
    radius = width[0]/2. * np.linspace(0., 1., num=len(s_1d))

    plt.plot(radius, s_1d, 'k',label='schwarzchild criterion')
    plt.plot(radius, l_1d, 'r--',label='Ledoux criterion')
    plt.ylim(-0.4, 0.2)
    plt.legend()
    left, right = plt.xlim(300, 600)
    plt.hlines(0., left, right, color='k', linestyle='--', alpha=0.2)
    
    plt.savefig(f"plots_conv_profile/{ds.basename}_criteria.png")
    
def plot_gradients(ds, act, ad, led, width):
    plt.figure(1, figsize= (10, 6))
    
    center = [int(act.shape[0]/2), int(act.shape[1]/2)]
    # take a 4 average at least
    # rightward, downward, leftward, upward
    act_1d = np.mean([act[center[0]:, center[1]], act[center[0], center[1]:], np.flip(act[center[0], :center[1]]), np.flip(act[:center[0], center[1]])], axis=0)
    ad_1d = np.mean([ad[center[0]:, center[1]], ad[center[0], center[1]:], np.flip(ad[center[0], :center[1]]), np.flip(ad[:center[0], center[1]])], axis=0)
    led_1d = np.mean([led[center[0]:, center[1]], led[center[0], center[1]:], np.flip(led[center[0], :center[1]]), np.flip(led[:center[0], center[1]])], axis=0)
    radius = width[0]/2. * np.linspace(0., 1., num=len(ad_1d))

    plt.plot(radius, act_1d, label='actual grad')
    plt.plot(radius, ad_1d, label='adiabt grad')
    plt.plot(radius, led_1d, label='adiabat+comp grad')
    plt.ylim(0., 0.7)
    plt.legend()
    left, right = plt.xlim()
    #plt.hlines(0., left, right, color='k', linestyle='--', alpha=0.2)
    
    plt.savefig(f"plots_conv_profile/{ds.basename}_gradients.png")
    
def _ad_excess(field, data):
    return  data[('boxlib', 'conv_actual')] - data[('boxlib', 'conv_adiabatic')]

def _ad_excess_led(field, data):
    return  data[('boxlib', 'conv_actual')] - data[('boxlib', 'conv_ledoux')]



if __name__ == '__main__':
    #load and add mass field
    fname = argv[1]
    
    if len(argv) == 3:
        width = (argv[2], 'km')
    else:
        width = (1.5e3, 'km')
    
    ds = yt.load(fname, hint='amrex')
    
    ds.add_field(name=("boxlib", "ad_excess"),
                function=_ad_excess,
                take_log=False,
                units = 'dimensionless',
                display_name="$\\Delta \\nabla $",
                sampling_type="local")

    ds.add_field(name=("boxlib", "ad_excess_led"),
                function=_ad_excess_led,
                take_log=False,
                units = 'dimensionless',
                display_name="Ledoux $\\Delta \\nabla$",
                sampling_type="local")
    
    act, ad, led, s_cri, l_cri = generate_frbs(ds, width)
    
    plot_gradients(ds, act, ad, led, width)
    plot_criteria(ds, s_cri, l_cri, width)

