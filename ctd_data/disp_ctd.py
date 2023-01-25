import csv
import datetime as dt
import matplotlib.pyplot as plt
from matplotlib.dates import MinuteLocator, DateFormatter


def parse_ctd_csv(cast_fpath, pres_thresh=0.5):
    datetimes = []
    pres = []
    temps = []
    sals = []
    oxys = []
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
                pres.append(float(row[pres_idx]))
                temps.append(float(row[temp_idx]))
                sals.append(float(row[sal_idx]))
                oxys.append(float(row[oxy_idx]))
                turbs.append(float(row[turb_idx]))
                dens.append(float(row[dens_idx]))
                
    return {
        'dts': datetimes,
        'pres': pres,
        'temps': temps,
        'sals': sals,
        'oxys': oxys,
        'dens': dens,
        'turbs': turbs
    }


def main():
    cast_file = 'CTD/mytilus_saturday_cast1.csv'

    # Read CTD data
    data = parse_ctd_csv(cast_file)

    # Make the plots
    fig, axs = plt.subplots(1, 2)

    # axs[0].set_title("AML / Seabird Calibration Cast (COM2) - 14/01/23")
    axs[0].plot(data['dts'], data['pres'], label="AML")
    # axs[0].legend()
    axs[0].invert_yaxis()
    axs[0].set_xlabel("Time")
    axs[0].set_ylabel("Pressure (dBar)")
    axs[0].grid(alpha=0.5)
    mins = MinuteLocator(interval=2)
    dateFmt = DateFormatter('%H:%M:%S')
    axs[0].xaxis.set_major_locator(mins)
    axs[0].xaxis.set_major_formatter(dateFmt)

    axs[1].plot(data['temps'], data['pres'])
    axs[1].set_xlabel("Temp")
    axs[1].set_ylabel("Pressure (dBar)")
    axs[1].invert_yaxis()

    plt.show()


if __name__ == "__main__":
    main()