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
    # TODO?
    if star_name is np.nan:
        star_name = star_info['gl']
    if star_name is np.nan:
        try:
            star_name = 'HIP_' + str(int(star_info['hip']))
        except:
            pass            
    return star_name

def find(star_table, central_star_name='Sun', upper_bound=5.0):
    print(f'Stars with minimum achieved distance <= {upper_bound} light years to {central_star_name}:\n')
    results, t = [], star_table.loc[star_table['proper'] == central_star_name]

    px, py, pz = list(t['x'])[0], list(t['y'])[0], list(t['z'])[0]
    vx, vy, vz = list(t['vx'])[0], list(t['vy'])[0], list(t['vz'])[0]
    
    star_table['x'] = star_table['x'] - px
    star_table['y'] = star_table['y'] - py
    star_table['z'] = star_table['z'] - pz
    
    star_table['vx'] = star_table['vx'] - vx
    star_table['vy'] = star_table['vy'] - vy
    star_table['vz'] = star_table['vz'] - vz

    star_table['dot_product'] = pd.eval("star_table.x*star_table.vx + star_table.y*star_table.vy + star_table.z*star_table.vz")   
    star_table['v_norm']      = pd.eval("star_table.vx**2 + star_table.vy**2 + star_table.vz**2") 
    star_table['t_min']       = pd.eval("-star_table.dot_product / star_table.v_norm")
    star_table['d_min']       = pd.eval("(star_table.x**2 + star_table.y**2 + star_table.z**2 - star_table.dot_product**2 / star_table.v_norm)**0.5")
    
    star_table.dropna()
    star_table = star_table[star_table['d_min'] <= upper_bound]

    print(f'filtered star_table length: {len(star_table)}')
    print(star_table)

    for num, row in star_table.iterrows():
        star_info = row.to_dict()

        if star_info['proper'] == central_star_name:
            continue

        star_name = get_name(star_info)
        
        star_pos      = [star_info['x'], star_info['y'], star_info['z']]
        star_velocity = [star_info['vx'], star_info['vy'], star_info['vz']]
        
        #t_min, d_min = hyperbola_minimum(star_pos, star_velocity)
        results.append({'star_name': star_name, 'star_pos': star_pos,
                            'star_velocity': star_velocity, 't_min': star_info['t_min'], 'd_min': star_info['d_min']})
        d_min, t_thousands = star_info['d_min'], star_info['t_min'] / 10
        print(f'{star_name} : d = {d_min:.2f} for t = {t_thousands:.1f} * 10^3 years')

    print(f'\n{len(results)} stars in total.')    
    return results

def text_overlap(text_pos, pt):
    for pos in text_pos:
        if np.abs(pos[0]-pt[0]) <= 50 and np.abs(pos[1]-pt[1]) <= 0.3:
            return True
    return False

def draw_interval(stars, t_start, t_end, central_star_name='Sun', label='', lw=1.0, fs=5, max_dist=8.0):
    fig, ax = plt.subplots()
    plt.suptitle(f'Nearest stars for {central_star_name}')

    text_pos, tx, ty = [], -10, -0.15
    
    for i, d in enumerate(stars):        
        times, distances = range(10*t_start, 10*t_end), []
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

    ax.xaxis.set_major_formatter(tkr.FuncFormatter(lambda x, pos: f'{x / 10:.0f}'))
    plt.xlabel('Time [thousands of years]')
    plt.ylabel(f'Distance to star [in light years]')
    plt.ylim([0, max_dist])
    plt.grid(True, zorder=0)
    plt.savefig(f'nearest_stars_{label}.svg', bbox_inches='tight')
