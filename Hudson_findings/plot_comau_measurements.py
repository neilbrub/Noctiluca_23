import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def main():
    fname = './Hudson_findings/comau_measurements.csv'
    df = pd.read_csv(fname)

    fig, axs = plt.subplots(1, 3)

    fig.suptitle("Mean Hudson measurements in Comau Fjord (1970)")

    depths = np.array(df['depth'])

    temps = np.array(df['temp'])
    temp_ranges = np.array(df['temp_range']) / 2
    axs[0].plot(temps, depths,
        marker='o', markersize=3, linewidth=1, c='k')
    axs[0].plot(temps - temp_ranges, depths, alpha=0.3, linewidth=0.5, c='k')
    axs[0].plot(temps + temp_ranges, depths, alpha=0.3, linewidth=0.5, c='k')
    axs[0].fill_betweenx(depths, temps - temp_ranges, temps + temp_ranges, alpha=0.05, color='k')

    sals = np.array(df['sal'])
    sal_ranges = np.array(df['sal_range']) / 2
    axs[1].plot(sals, depths,
        marker='o', markersize=3, linewidth=1, c='k')
    axs[1].plot(sals - sal_ranges, depths, alpha=0.3, linewidth=0.5, c='k')
    axs[1].plot(sals + sal_ranges, depths, alpha=0.3, linewidth=0.5, c='k')
    axs[1].fill_betweenx(depths, sals - sal_ranges, sals + sal_ranges, alpha=0.05, color='k')

    o2s = np.array(df['o2'])
    o2_ranges = np.array(df['o2_range']) / 2
    axs[2].plot(o2s, depths,
        marker='o', markersize=3, linewidth=1, c='k')
    axs[2].plot(o2s - o2_ranges, depths, alpha=0.3, linewidth=0.5, c='k')
    axs[2].plot(o2s + o2_ranges, depths, alpha=0.3, linewidth=0.5, c='k')
    axs[2].fill_betweenx(depths, o2s - o2_ranges, o2s + o2_ranges, alpha=0.05, color='k')


    axs[0].set_xlim([9, 16])
    axs[1].set_xlim([25, 34])
    axs[2].set_xlim([2, 9])

    axs[0].set_title("Temperature (°C)")
    axs[1].set_title("Salinity (%)")
    axs[2].set_title("Dissolved Oxygen (mL/L)")

    axs[0].set_ylabel("Depth (m)")

    for idx, ax in enumerate(axs):
        ax.invert_yaxis()
        ax.spines.top.set_visible(False)
        ax.spines.right.set_visible(False)
        if idx > 0:
            ax.spines.left.set_visible(False)
        ax.grid(alpha=0.4)

    plt.show()


if __name__ == "__main__":
    main()