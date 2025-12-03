import matplotlib.pyplot as plt
import matplotlib.ticker as tkr

import numpy as np
import pandas as pd 

pd.set_option('display.max_rows', None)
#pd.set_option('display.max_columns', None)

# maths:
def get_coord(x0, v, t):
    return [x0[0] + v[0]*t, x0[1] + v[1]*t, x0[2] + v[2]*t]

def dist(p, q):
    return np.sqrt((p[0]-q[0])**2 + (p[1]-q[1])**2 + (p[2]-q[2])**2)

def hyperbola_minimum(star_pos, star_velocity):
    P, Q, R = star_pos
    A, B, C = star_velocity
    
    Vn = (A*A + B*B + C*C)
    t_min = 0.0 # not moving?

    if Vn != 0.0: 
        t_min = -(A*P + B*Q + C*R)/Vn

    d_min = np.sqrt((P+A*t_min)**2 + (Q+B*t_min)**2 + (R+C*t_min)**2)
    return t_min, d_min

# data utils:
def prepare_data(df):
    sorted_df = df.sort_values(by='dist')[['id', 'hip', 'gl', 'proper', 'dist', 'x', 'y', 'z', 'vx', 'vy', 'vz']]

    light_years_per_parsec = 3.261563 # parsecs to light years everywhere
    for k in ['dist', 'x', 'y', 'z', 'vx', 'vy', 'vz']:
        sorted_df[k] = light_years_per_parsec * sorted_df[k]

    for k in ['vx', 'vy', 'vz']:
        sorted_df[k] = 100 * sorted_df[k] # velocity to ly / century

    sorted_df = sorted_df.drop([0, 71455]) # remove row for Sun and AlphaCentauriA for convenience
    top_chosen = len(sorted_df)
    nearest_stars_table = sorted_df.head(top_chosen)
    
    return nearest_stars_table

def get_name(star_info):
    star_name = star_info['proper'] 
    if star_name is np.nan:
        star_name = star_info['gl']

    if star_name is np.nan:
        try:
            star_name = 'HIP_' + str(int(star_info['hip']))
        except:
            pass
            
    return star_name

# custom formatter:
def centuries_to_thousands(x, pos):
    return f'{x / 10:.0f}'

def find(star_table, upper_bound=7.0):
    results = []

    print(f'Stars with minimum achieved distance <= {upper_bound} light years:\n')

    for num, row in star_table.iterrows():
        star_info = row.to_dict()
        star_name = get_name(star_info)
        star_pos = [star_info['x'], star_info['y'], star_info['z']]
        star_velocity = [star_info['vx'], star_info['vy'], star_info['vz']]

        t_min, d_min = hyperbola_minimum(star_pos, star_velocity)

        if d_min <= upper_bound:
            t_thousands = t_min/10
            print(f'{star_name} : d = {d_min:.2f} for t = {t_thousands:.1f} * 10^3 years')

            results.append({'star_name': star_name, 
                            'star_pos': star_pos,
                            'star_velocity': star_velocity,
                            't_min': t_min,
                            'd_min': d_min})
    
    print(f'\n{len(results)} stars in total.')
    return results

def draw_interval(stars, t_start, t_end, label='', lw=1.0, fs=5):
    fig, ax = plt.subplots()

    for d in stars:        
        times = range(10*t_start, 10*t_end)
        tx, ty = -10, -0.25

        distances = []
        for t in times: # time in centuries for better resolution
            distances.append(dist([0, 0, 0], get_coord(d['star_pos'], d['star_velocity'], t)))
                
        plt.plot(times, distances, color='black', linewidth=lw, zorder=-1)
        if times[0] <= d['t_min'] and d['t_min'] <= times[-1]:
            plt.text(d['t_min'] + tx, d['d_min'] + ty, d['star_name'], fontsize=fs, color='black')

    formatter = tkr.FuncFormatter(centuries_to_thousands)
    ax.xaxis.set_major_formatter(formatter)

    plt.xlabel("Time [thousands of years]")
    plt.ylabel("Distance to Sun [in light years]")

    plt.ylim([0, 10])
    plt.grid(True)

    plt.savefig(f'nearest_stars_{label}.svg', bbox_inches='tight')

if __name__ == '__main__':
    df = pd.read_csv('hygdata_v41.csv')
    nearest_stars_table = prepare_data(df)

    stars = find(nearest_stars_table, 7.0) 

    draw_interval(stars, -2000, -20, 'Past', lw = 0.5, fs = 3)
    draw_interval(stars, -20, 100, 'Standard', lw = 1.0, fs = 3)
    draw_interval(stars, 100, 2000, 'Long', lw = 0.5, fs = 3)
