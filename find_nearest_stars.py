import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
import numpy as np
import pandas as pd 

from mini_star_lib import light_years_per_parsec

pd.set_option('display.max_rows', None)

def get_coord(x0, v, t):
    return [x0[0] + v[0]*t, x0[1] + v[1]*t, x0[2] + v[2]*t]

def dist(p, q):
    return np.sqrt((p[0]-q[0])**2 + (p[1]-q[1])**2 + (p[2]-q[2])**2)

def hyperbola_minimum(star_pos, star_velocity):
    P, Q, R = star_pos
    A, B, C = star_velocity
    Vn, t_min = (A*A + B*B + C*C), 0.0
    
    if Vn >= 1e-8: 
        t_min = -(A*P + B*Q + C*R)/Vn
    d_min = np.sqrt((P+A*t_min)**2 + (Q+B*t_min)**2 + (R+C*t_min)**2)
    return t_min, d_min

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

def centuries_to_thousands(x, pos):
    return f'{x / 10:.0f}'

def find(star_table, upper_bound=5.0):
    results = []
    print(f'Stars with minimum achieved distance <= {upper_bound} light years:\n')

    for num, row in star_table.iterrows():
        star_info = row.to_dict()
        star_name, star_pos = get_name(star_info), [star_info['x'], star_info['y'], star_info['z']]
        star_velocity = [star_info['vx'], star_info['vy'], star_info['vz']]
        
        t_min, d_min = hyperbola_minimum(star_pos, star_velocity)

        if d_min <= upper_bound:
            t_thousands = t_min / 10
            print(f'{star_name} : d = {d_min:.2f} for t = {t_thousands:.1f} * 10^3 years')
            results.append({'star_name': star_name, 'star_pos': star_pos,
                            'star_velocity': star_velocity, 't_min': t_min, 'd_min': d_min})
    
    print(f'\n{len(results)} stars in total.')
    return results

def text_overlap(text_pos, pt):
    for pos in text_pos:
        if np.abs(pos[0]-pt[0]) <= 50 and np.abs(pos[1]-pt[1]) <= 0.3:
            return True
    return False

def draw_interval(stars, t_start, t_end, label='', lw=1.0, fs=5, max_dist=8.0):
    fig, ax = plt.subplots()
    text_pos = []
    
    for i, d in enumerate(stars):        
        times, distances = range(10*t_start, 10*t_end), []
        tx, ty = -10, -0.15
        for t in times: # time in centuries for better resolution
            distances.append(dist([0, 0, 0], get_coord(d['star_pos'], d['star_velocity'], t)))
        distances, plot_visible = np.array(distances), False

        if np.any((distances >= 0) & (distances <= max_dist)):
            plt.plot(times, distances, color='black', linewidth=lw, zorder=2*i+3)
            plot_visible = True

        if times[0] <= d['t_min'] and d['t_min'] <= times[-1] and plot_visible:
            ttx, tty = d['t_min'] + tx, d['d_min'] + ty
            if text_overlap(text_pos, (ttx, tty)):
                tty -= 0.15
            plt.text(ttx, tty, d['star_name'], fontsize=fs, color='#1f1f1f', zorder=2*(i+len(stars)+1)+3)
            text_pos.append((ttx, tty))

    ax.xaxis.set_major_formatter(tkr.FuncFormatter(centuries_to_thousands))
    plt.xlabel("Time [thousands of years]")
    plt.ylabel("Distance to Sun [in light years]")
    plt.ylim([0, 8])
    plt.grid(True, zorder=0)
    plt.savefig(f'nearest_stars_{label}.svg', bbox_inches='tight')
