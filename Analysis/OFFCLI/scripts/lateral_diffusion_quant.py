import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from sklearn.utils import resample
import argparse
import os
import json

def analyze_diffusion(xvg_file, window_size=10000, window_step=200, min_r2=0.99, title=None, out_dir=None):
    data = np.loadtxt(xvg_file, comments=["@", "#"])
    time = data[:, 0]  # in ps or ns
    msd = data[:, 1]   # in nm^2 (ps/nm^2 from xvg)

    # best_slope = None
    best_r2 = -np.inf
    best_window = None

    # find  most stable/linear region
    for i in range(0, len(msd) - window_size, window_step):
        time_window = time[i:i + window_size]
        msd_window = msd[i:i + window_size]
        slope, intercept, r_value, p_value, std_err = linregress(time_window, msd_window)
        
        if r_value**2 >= min_r2:
            if r_value**2 > best_r2:
                best_r2 = r_value**2
                # best_slope = slope
                best_window = (i, i + window_size)
        
    best_time_window = time[best_window[0]:best_window[1]]
    best_msd_window = msd[best_window[0]:best_window[1]]

    slope, intercept, r_value, x, std_err = linregress(best_time_window, best_msd_window)
    DL = 0.25 * slope # 2D einstein
    DL_error = 0.25 * std_err  # SE of DL

    results = { # result dictionary for json
        "diffusion_coefficient": float(DL),
        "std_error_regression": float(DL_error),
        "linear_region_start_time": float(time[best_window[0]]),
        "linear_region_end_time": float(time[best_window[1] - 1]),
        "r_squared": float(best_r2),
    }

    plt.figure(figsize=(10, 6))
    plt.plot(time, msd, label='MSD', color='#0333b0')
    fit_time = best_time_window
    fit_msd = slope * fit_time + intercept
    plt.plot(fit_time, fit_msd, 'r--', label=f'Linear Fit (R² = {best_r2:.4f})')
    
    # time range annotations. can just comment these out if preferred
    start_time = fit_time[0]
    end_time = fit_time[-1]
    plt.vlines(start_time, 0, slope * start_time + intercept, colors='r', linestyles='dashed', alpha=0.7)
    plt.vlines(end_time, 0, slope * end_time + intercept, colors='r', linestyles='dashed', alpha=0.7)
    plt.annotate(f'{start_time:.1f} ps', xy=(start_time, slope * start_time + intercept/2), xytext=(-15, -30), textcoords='offset points', arrowprops=dict(arrowstyle='->', color='r', alpha=0.7))
    plt.annotate(f'{end_time:.1f} ps', xy=(end_time, slope * end_time + intercept/2), xytext=(15, -30), textcoords='offset points', arrowprops=dict(arrowstyle='->', color='r', alpha=0.7))
    
    plt.xlabel('Time (ps)')
    plt.ylabel('MSD (nm²)')
    plt.legend()
    
    if title:
        plt.title(title)
    else:
        plt.title('MSD with Linear Fit')

    plt.tight_layout()
    
    if out_dir:
        plot_path = os.path.join(out_dir, 'fit_diffusion.png')
        json_path = os.path.join(out_dir, 'diffusion_results.json')

        plt.savefig(plot_path, dpi=300) #
        
        with open(json_path, 'w') as f:
            json.dump(results, f, indent=4)
    else:
        plt.show()
        
    return results

def main():
    parser = argparse.ArgumentParser(description="Analyze lateral diffusion from MSD data")
    parser.add_argument("--xvg", help="Path to XVG file")
    parser.add_argument("--window-size", type=int, default=10000, help="Default: 10000 = 100 ns")
    parser.add_argument("--window-step", type=int, default=200, help="Default: 200 = 2 ns")
    parser.add_argument("--min-r2", type=float, default=0.99, help="Minimum linear fit value (default: 0.99)")
    parser.add_argument("--title")
    parser.add_argument("--out_dir")
    
    args = parser.parse_args()
    
    xvg_file = args.xvg
    if not xvg_file:
        default_dir = os.path.join("outputs", "msd_output", "xvg_files") # checks msd_plot.py output for default (if scripts are ran in order)
        if os.path.exists(default_dir):
            xvg_files = [f for f in os.listdir(default_dir) if f.endswith('.xvg')]
            if xvg_files:
                xvg_file = os.path.join(default_dir, xvg_files[0])
                print(f"Using default XVG file: {xvg_file}")
            else:
                print(f"No XVG files found in {default_dir}")
                return
        else:
            print(f"No XVG file provided and default directory {default_dir} not found")
            return
    
    if args.out_dir and not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir, exist_ok=True)
    
    results = analyze_diffusion(
        xvg_file=xvg_file,
        window_size=args.window_size,
        window_step=args.window_step,
        min_r2=args.min_r2,
        title=args.title,
        out_dir=args.out_dir
    )
###### Result check ######
    # print(f"Diffusion coefficient (D_L): {results['diffusion_coefficient']:.4e} nm²/ps")
    # print(f"Standard error from linear regression: {results['std_error_regression']:.4e} nm²/ps")
    # print(f"Best linear region: Time window {results['linear_region_start_time']:.1f} to {results['linear_region_end_time']:.1f} ps")
    # print(f"R² of linear fit: {results['r_squared']:.6f}")

if __name__ == "__main__":
    main()