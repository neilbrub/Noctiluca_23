import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import MinuteLocator, DateFormatter
from sklearn.linear_model import LinearRegression

import re
import csv
import datetime as dt


def Aanderaa_O2_compensation(meas_o2, temp, pres, sal, ref_sal=0):
    # Equations & coefficients from Aanderaa TD 269 Operating Manual for 4831, June 2017
    # https://www.aanderaa.com/media/pdfs/oxygen-optode-4330-4835-and-4831.pdf
    
    # Salinity
    B0 = -6.24097e-3
    B1 = -6.93498e-3
    B2 = -6.90358e-3
    B3 = -4.29155e-3
    C0 = -3.11680e-7

    ts = np.log((298.15-temp)/(273.15+temp))

    if ref_sal == 0:
        sal_corr_o2 = meas_o2 * np.exp(sal*(B0 + B1*ts + B2*ts**2 + B3*ts**3) + C0*sal**2)
    else:
        sal_corr_o2 = meas_o2 * np.exp((sal-ref_sal)*(B0 + B1*ts + B2*ts**2 + B3*ts**3) + C0*(sal**2 - ref_sal**2))

    # Apply pressure compensation to salinity-compensated values
    corr_o2 = sal_corr_o2 * (1 + (0.032*pres)/1000)

    return corr_o2


def PSU_to_ref_sal(psu): return psu * (35.16504/35)


def read_seabird(fname, pres_thresh=0.5):
    pres = []  # 0
    temps = []  # 1
    conds = []  # 2
    oxys = []  # 3
    sals = []  # 7 (absolute salinity)
    dens = []  # 8

    with open(fname, 'r') as f:
        lines = f.readlines()
        for l_i, line in enumerate(lines):
            if line[0] in ['*', '#']: continue

            line = re.sub('\s+', ',', line)[1:-1]
            vals = line.split(',')

            # Clip to min pressures
            if float(vals[0]) < pres_thresh: continue

            p = float(vals[0])
            temp = float(vals[1])
            cond = float(vals[2])
            oxy = float(vals[3])
            sal = float(vals[7])
            den = float(vals[8])

            pres.append(p)
            temps.append(temp)
            conds.append(cond)
            oxys.append(oxy)
            sals.append(sal)
            dens.append(den)

    return {
        'pres': pres,
        'temps': temps,
        'sals': sals,
        'oxys': oxys,
        'dens': dens,
        'conds': conds
    }


def parse_ctd_csv(cast_fpath, pres_thresh=0.5):
    datetimes = []
    pres = []
    temps = []
    sals = []
    oxys = []
    corr_oxys = []
    turbs = []
    dens = []

    cast_datestr = ''
    cast_timestr = ''
    cast_datetime = None
    date_idx = None
    time_idx = None
    pres_idx = None
    temp_idx = None
    sal_idx = None
    oxy_idx = None

    with open(cast_fpath, newline='') as f:
        reader = csv.reader(f)

        read_headings = False
        data_start = False
        for row in reader:
            if not data_start and not read_headings:
                if row[0] == '[data]':
                    # Next row will contain header
                    read_headings = True
                elif row[0].split('=')[0] == 'date':
                    cast_datestr = row[0].split('=')[1]  # yyyy-mm-dd
                elif row[0].split('=')[0] == 'time':
                    cast_timestr = row[0].split('=')[1]  # HH:MM:SS.ms
                    
            elif read_headings:
                cast_datetime_str = cast_datestr + ':' + cast_timestr[:-3]  # ignore deci-second
                cast_datetime = dt.datetime.strptime(cast_datetime_str, '%Y-%m-%d:%H:%M:%S')

                date_idx = row.index('Date (yyyy-mm-dd)')
                turb_idx = row.index('Turbidity (NTU)')
                time_idx = row.index('Time')
                pres_idx = row.index('Pressure (dBar)')
                temp_idx = row.index('Temperature (C)')
                sal_idx = row.index('Salinity (PSU)')
                oxy_idx = row.index('Aanderaa 4831')
                dens_idx = row.index('Density (kg m-3)')
                read_headings = False
                # Now expect only data rows 
                data_start = True
            
            elif data_start:
                # Merge date & time into datetime object
                sample_time = row[time_idx]
                sample_min = sample_time.split(':')[0]
                sample_sec = sample_time.split(':')[1].split('.')[0]
                sample_dsec = sample_time.split(':')[1].split('.')[1]
                sample_dt = cast_datetime.replace(
                    minute=int(sample_min),
                    second=int(sample_sec),
                    microsecond=int(sample_dsec) * int(1e5))


                # if sample_dt < sample_dt.replace(minute=52, second=45): continue
                if float(row[pres_idx]) < pres_thresh: continue

                datetimes.append(sample_dt)
                p = float(row[pres_idx])
                temp = float(row[temp_idx])
                sal = PSU_to_ref_sal(float(row[sal_idx]))
                oxy = float(row[oxy_idx])
                turb = float(row[turb_idx])
                den = float(row[dens_idx])

                pres.append(p)
                temps.append(temp)
                sals.append(sal)
                oxys.append(oxy)
                turbs.append(turb)
                dens.append(den)

                corr_o2 = Aanderaa_O2_compensation(oxy, temp, p, sal, ref_sal=0)
                corr_oxys.append(corr_o2)
                
                
    return {
        'dts': datetimes,
        'pres': pres,
        'temps': temps,
        'sals': sals,
        'oxys': oxys,
        'corr_oxys': corr_oxys,
        'dens': dens,
        'turbs': turbs
    }

def separate_casts_seabird(down_stop_idx, up_start_idx, sb_data):
    """
    down_stop_idx & up_start_idx - sample indices
    """
    down_pres = []
    down_temps = []
    down_oxys = []
    down_corr_oxys = []
    down_sals = []
    down_dens = []

    up_pres = []
    up_temps = []
    up_oxys = []
    up_corr_oxys = []
    up_sals = []
    up_dens = []

    for i in range(len(sb_data['pres'])):
        if i < down_stop_idx:
            down_pres.append(sb_data['pres'][i])
            down_temps.append(sb_data['temps'][i])
            down_oxys.append(sb_data['oxys'][i])
            down_sals.append(sb_data['sals'][i])
            down_dens.append(sb_data['dens'][i])
        elif i > up_start_idx:
            up_pres.append(sb_data['pres'][i])
            up_temps.append(sb_data['temps'][i])
            up_oxys.append(sb_data['oxys'][i])
            up_sals.append(sb_data['sals'][i])
            up_dens.append(sb_data['dens'][i])
    
    return {
        'down': {
            'pres': down_pres,
            'temps': down_temps,
            'sals': down_sals,
            'oxys': down_oxys,
            'dens': down_dens
        },
        'up': {
            'pres': up_pres,
            'temps': up_temps,
            'sals': up_sals,
            'oxys': up_oxys,
            'dens': up_dens
        }
    }


def separate_casts_aml(down_stop_time, up_start_time, aml_data): 
    """
    down_stop_time & up_start_time - lists of [hour, minute, second]
    dts - datetime objects for each measurement
    """
    down_pres = []
    down_temps = []
    down_sals = []
    down_oxys = []
    down_corr_oxys = []
    down_dens = []

    up_pres = []
    up_temps = []
    up_sals = []
    up_oxys = []
    up_corr_oxys = []
    up_dens = []
    
    for i in range(len(aml_data['dts'])):
        if aml_data['dts'][i] < aml_data['dts'][i].replace(
                hour=down_stop_time[0],
                minute=down_stop_time[1],
                second=down_stop_time[2]):
            down_pres.append(aml_data['pres'][i])
            down_temps.append(aml_data['temps'][i])
            down_oxys.append(aml_data['oxys'][i])
            down_corr_oxys.append(aml_data['corr_oxys'][i])
            down_sals.append(aml_data['sals'][i])
            down_dens.append(aml_data['dens'][i])
    
        elif aml_data['dts'][i] > aml_data['dts'][i].replace(
                hour=up_start_time[0],
                minute=up_start_time[1],
                second=up_start_time[2]):
            up_pres.append(aml_data['pres'][i])
            up_temps.append(aml_data['temps'][i])
            up_oxys.append(aml_data['oxys'][i])
            up_corr_oxys.append(aml_data['corr_oxys'][i])
            up_sals.append(aml_data['sals'][i])
            up_dens.append(aml_data['dens'][i])

    return {
        'down': {
            'pres': down_pres,
            'temps': down_temps,
            'sals': down_sals,
            'oxys': down_oxys,
            'corr_oxys': down_corr_oxys,
            'dens': down_dens
        },
        'up': {
            'pres': up_pres,
            'temps': up_temps,
            'sals': up_sals,
            'oxys': up_oxys,
            'corr_oxys': up_corr_oxys,
            'dens': up_dens
        }
    }


def bin_by_pres(pres, vals=None):
    min_pres = np.floor(min(pres))
    max_pres = np.ceil(max(pres))
    pres_bins = np.arange(min_pres, max_pres, 0.5)

    # Option to just return binned pressured values:
    if not vals: return pres_bins

    val_sums = np.zeros(pres_bins.size)
    val_counts = np.zeros(pres_bins.size)

    for i, val in enumerate(vals):
        p = pres[i]
        target_diff = p - pres_bins
        bin_idx = np.argmin(target_diff[target_diff >= 0])
        val_sums[bin_idx] = val_sums[bin_idx] + val
        val_counts[bin_idx] = val_counts[bin_idx] + 1

    binned_vals = val_sums / np.maximum(val_counts, 1)
    binned_vals[binned_vals==0] = np.nan

    return binned_vals


def main():
    # Read seabird data
    seabird_cast1_fp = 'CTD/Seabird/2023-01-14T185616 SBE0251244_filter_align_ctm_loopteos_10.cnv'
    sb_cast1 = read_seabird(seabird_cast1_fp, pres_thresh=1)

    # Read AML data
    aml_cast1_fp = 'CTD/AML/noctiluca_saturday_cast1.csv'
    aml_cast1 = parse_ctd_csv(aml_cast1_fp, pres_thresh=1)


    # fig, axs = plt.subplots(2, figsize=(12, 8))

    # axs[0].plot(aml_cast1['dts'], aml_cast1['pres'])
    # axs[0].invert_yaxis()
    # axs[0].grid(alpha=0.3)
    # mins = MinuteLocator(interval=1)
    # dateFmt = DateFormatter('%H:%Mf%S')
    # axs[0].xaxis.set_major_locator(mins)
    # axs[0].xaxis.set_major_formatter(dateFmt)
    # axs[0].set_title("AML")
    # axs[0].set_ylim([210, -5])

    # axs[1].plot(np.arange(len(sb_cast1['pres'])), sb_cast1['pres'])
    # axs[1].invert_yaxis()
    # axs[1].grid(alpha=0.3)
    # axs[1].set_title("Seabird")
    # axs[1].set_ylim([210, -5])
    # # axs[1].set_yticklabels([])

    # fig.suptitle("Comparison Cast 1 - Pressure")
    # fig.supxlabel('Time / Sample Idx.')
    # fig.supylabel("Pressure (dBar)")
    # # axs[0].set_ylabel("Pressure (dBar)")

    # # plt.savefig('figs/calib1_pres_vertical.png')
    # fig.tight_layout()
    # plt.show()

    sb_down1_stop_idx = 8000
    sb_up1_start_idx = 8750

    aml_down1_stop_time = [16,12,0]  # H, M, S
    aml_up1_start_time = [16,12,30]

    sb_casts1 = separate_casts_seabird(sb_down1_stop_idx, sb_up1_start_idx, sb_cast1)
    aml_casts1 = separate_casts_aml(aml_down1_stop_time, aml_up1_start_time, aml_cast1)

    # Seabird down
    sb_down1_pres_binned = bin_by_pres(sb_casts1['down']['pres'])
    sb_down1_temps_binned = bin_by_pres(sb_casts1['down']['pres'], sb_casts1['down']['temps'])
    sb_down1_sals_binned = bin_by_pres(sb_casts1['down']['pres'], sb_casts1['down']['sals'])
    sb_down1_oxys_binned = bin_by_pres(sb_casts1['down']['pres'], sb_casts1['down']['oxys'])
    sb_down1_dens_binned = bin_by_pres(sb_casts1['down']['pres'], sb_casts1['down']['dens'])

    # AML down
    aml_down1_pres_binned = bin_by_pres(aml_casts1['down']['pres'])
    aml_down1_temps_binned = bin_by_pres(aml_casts1['down']['pres'], aml_casts1['down']['temps'])
    aml_down1_sals_binned = bin_by_pres(aml_casts1['down']['pres'], aml_casts1['down']['sals'])
    aml_down1_oxys_binned = bin_by_pres(aml_casts1['down']['pres'], aml_casts1['down']['oxys'])
    aml_down1_corr_oxys_binned = bin_by_pres(aml_casts1['down']['pres'], aml_casts1['down']['corr_oxys'])
    aml_down1_dens_binned = bin_by_pres(aml_casts1['down']['pres'], aml_casts1['down']['dens'])

    fig, axs = plt.subplots(2, 2, figsize=(10,10))

    axs[0, 0].plot(sb_down1_temps_binned, sb_down1_pres_binned, label='Seabird', alpha=0.8)
    axs[0, 0].plot(aml_down1_temps_binned, aml_down1_pres_binned, label='AML', alpha=0.7)
    axs[0, 0].legend(loc=4)
    axs[0, 0].invert_yaxis()
    axs[0, 0].grid(alpha=0.5)
    axs[0, 0].set_title("In-Situ Temperature ($\degree$C)")

    axs[0, 1].plot(sb_down1_oxys_binned, sb_down1_pres_binned, label='Seabird', alpha=0.8)
    axs[0, 1].plot(aml_down1_corr_oxys_binned, aml_down1_pres_binned, label='AML (corrected)', alpha=0.7)
    axs[0, 1].plot(aml_down1_oxys_binned, aml_down1_pres_binned, 'k--', label='AML (uncorrected)', alpha=0.6, linewidth=0.5)
    axs[0, 1].legend(loc=4)
    axs[0, 1].invert_yaxis()
    axs[0, 1].grid(alpha=0.5)
    axs[0, 1].set_title("O$_2$ Concentration ($\mu$mol/L)")

    axs[1, 0].plot(sb_down1_sals_binned, sb_down1_pres_binned, label='Seabird (absolute)', alpha=0.8)
    axs[1, 0].plot(aml_down1_sals_binned, aml_down1_pres_binned, label='AML (ref)', alpha=0.7)
    axs[1, 0].legend(loc=3)
    axs[1, 0].invert_yaxis()
    axs[1, 0].grid(alpha=0.5)
    axs[1, 0].set_title("Salinity (g/kg)")

    axs[1, 1].plot(sb_down1_dens_binned, sb_down1_pres_binned, label='Seabird', alpha=0.8)
    axs[1, 1].plot(aml_down1_dens_binned, aml_down1_pres_binned, label='AML', alpha=0.7)
    axs[1, 1].legend(loc=3)
    axs[1, 1].invert_yaxis()
    axs[1, 1].grid(alpha=0.5)
    axs[1, 1].set_title("Density (kg/m$^3$)")

    fig.suptitle("Comparison Cast 1 - Downcast")
    fig.supylabel("Pressure (dBar)")
    fig.tight_layout()

    # plt.savefig('figs/calib1_downcast.png')
    plt.show()



if __name__ == "__main__":
    main()