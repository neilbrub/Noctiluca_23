import csv
import datetime as dt
import matplotlib.pyplot as plt
from matplotlib.dates import MinuteLocator, DateFormatter


def main():
    cast_fpath = "CTD/noctiluca_saturday_cast1.csv"
    down_start = '16:42:29'
    down_end = '16:44:40'
    up_start = '16:44:47'
    up_end = '16:46:31'

    down_fpath = 'noctiluca_saturday_cast1_down.csv'
    up_fpath = 'noctiluca_saturday_cast1_up.csv'

    down_datetimes = []
    down_pres = []
    down_temps = []
    down_sals = []
    down_oxys = []
    down_turbs = []

    up_datetimes = []
    up_pres = []
    up_temps = []
    up_sals = []
    up_oxys = []
    up_turbs = []

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


                down_start_min = int(down_start.split(':')[1])
                down_start_sec = int(down_start.split(':')[2])

                down_end_min = int(down_end.split(':')[1])
                down_end_sec = int(down_end.split(':')[2])

                up_start_min = int(up_start.split(':')[1])
                up_start_sec = int(up_start.split(':')[2])

                up_end_min = int(up_end.split(':')[1])
                up_end_sec = int(up_end.split(':')[2])

                if (sample_dt > sample_dt.replace(minute=down_start_min, second=down_start_sec)
                        and sample_dt < sample_dt.replace(minute=down_end_min, second=down_end_sec)):
                    down_datetimes.append(sample_dt)
                    down_pres.append(float(row[pres_idx]))
                    down_temps.append(float(row[temp_idx]))
                    down_sals.append(float(row[sal_idx]))
                    down_oxys.append(float(row[oxy_idx]))
                    down_turbs.append(float(row[turb_idx]))

                elif (sample_dt > sample_dt.replace(minute=up_start_min, second=up_start_sec)
                        and sample_dt < sample_dt.replace(minute=up_end_min, second=up_end_sec)):
                    up_datetimes.append(sample_dt)
                    up_pres.append(float(row[pres_idx]))
                    up_temps.append(float(row[temp_idx]))
                    up_sals.append(float(row[sal_idx]))
                    up_oxys.append(float(row[oxy_idx]))
                    up_turbs.append(float(row[turb_idx]))
            
        fig,ax = plt.subplots(1)
        ax.plot(down_temps, down_pres, label='down')
        ax.plot(up_temps, up_pres, label='up')
        ax.legend()
        ax.invert_yaxis()
        plt.show()

        with open('CTD/isolated_data/' + str(down_fpath), 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Datetime', 'Pressure', 'Temperature', 'Salinity', 'Oxygen', 'Turbidity'])
            for i, time in enumerate(down_datetimes):
                writer.writerow([time, down_pres[i], down_temps[i], down_sals[i], down_oxys[i], down_turbs[i]])

        with open('CTD/isolated_data/' + str(up_fpath), 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Datetime', 'Pressure', 'Temperature', 'Salinity', 'Oxygen', 'Turbidity'])
            for i, time in enumerate(up_datetimes):
                writer.writerow([time, up_pres[i], up_temps[i], up_sals[i], up_oxys[i], up_turbs[i]])



if __name__ == "__main__":
    main()