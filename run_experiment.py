import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
import numpy as np
import pandas as pd 
from find_nearest_stars import find, draw_interval
from mini_star_lib import get_uvw, light_years_per_parsec, kms_to_lyc, Gaia_EDR3_ID_to_common

def prepare_HYG_data(df):
    sorted_df = df.sort_values(by='dist')[['id', 'hip', 'gl', 'proper', 'dist', 'x', 'y', 'z', 'vx', 'vy', 'vz']]
    for k in ['dist', 'x', 'y', 'z', 'vx', 'vy', 'vz']:
        sorted_df[k] = light_years_per_parsec * sorted_df[k]
    for k in ['vx', 'vy', 'vz']:
        sorted_df[k] = 100 * sorted_df[k] # velocity to ly / century
    sorted_df = sorted_df.drop([0, 71455]) # remove row for Sun and AlphaCentauriA for convenience
    return sorted_df 

def HYG_experiment(path_to_HYG):
    df = pd.read_csv(path_to_HYG)
    nearest_stars_table = prepare_HYG_data(df)
    stars = find(nearest_stars_table, 5.0) 
    draw_interval(stars, -20, 100, 'HYG_test', lw = 1.0, fs = 3)

def GCNS_experiment_with_patch(path_to_GCNS, patch=True):
    column_names = ['GaiaEDR3','RAdeg','e_RAdeg','DEdeg','e_DEdeg',
                'Plx','e_Plx','pmRA','e_pmRA','pmDE',
                'e_pmDE','Gmag','RFG','BPmag','RFBP',
                'RPmag','RFRP','E(BP/RP)','RUWE','IPDfmp',
                'RV','e_RV','r_RV','f_RV','GCNSprob',
                'WDprob','Dist1','Dist16','Dist50','Dist84',
                'xcoord50','xcoord16','xcoord84','ycoord50','ycoord16',
                'ycoord84','zcoord50','zcoord16','zcoord84','Uvel50',
                'Uvel16','Uvel84','Vvel50','Vvel16','Vvel84',
                'Wvel50','Wvel16','Wvel84','GUNN','r_GUNN',
                'gmag','e_gmag','rmag','e_rmag','imag',
                'e_imag','zmag','e_zmag','2MASS','Jmag',
                'e_Jmag','Hmag','e_Hmag','Ksmag','e_Ksmag',
                'WISE','W1mag','e_W1mag','W2mag','e_W2mag',
                'W3mag','e_W3mag','W4mag','e_W4mag']

    df_gcns = pd.read_csv(path_to_GCNS, header=None)
    df_gcns.columns = column_names
    
    if patch:
        cnt_nans, fixed = df_gcns['RV'].isna().sum(), 0
        part = (100*cnt_nans/len(df_gcns))
        print(f'Missed RV values: {cnt_nans} / {len(df_gcns)} ({part:.2f}% of data)')
        
        df_survey = pd.read_csv('sosdr1.csv', engine='pyarrow', dtype_backend='pyarrow')
        df_rvstdcat = pd.read_csv('rvstdcat.csv')

        survey_map   = dict(zip(df_survey['sosdr1_gaiaSourceId'], df_survey['sosdr1_RVcor_merged']))
        rvstdcat_map = dict(zip(df_rvstdcat['Source'], df_rvstdcat['RV']))

        for index, row in df_gcns.iterrows():
            star_info = row.to_dict()

            if np.isnan(star_info['RV']):
                success, Gaia_id = 0, star_info['GaiaEDR3'] # cross-match by GaiaEDR3 id
                if Gaia_id in rvstdcat_map:
                    success = 1
                    df_gcns.iloc[index, df_gcns.columns.get_loc('RV')] = rvstdcat_map[Gaia_id]
                if success == 0 and Gaia_id in survey_map:
                    success = 1
                    df_gcns.iloc[index, df_gcns.columns.get_loc('RV')] = survey_map[Gaia_id]
                fixed += success

        part = (100*(cnt_nans-fixed)/len(df_gcns))
        print(f'After patch missed RV values: {cnt_nans - fixed} / {len(df_gcns)} ({part:.2f}% of data)')
        
        star_names = []
        for index, row in df_gcns.iterrows():
            star_info = row.to_dict()
            star = star_info['GaiaEDR3']
            
            if star in Gaia_EDR3_ID_to_common:
                star = Gaia_EDR3_ID_to_common[star]
            if np.isnan(star_info['Uvel50']) and not np.isnan(star_info['RV']): # fix by restored RV
                V = get_uvw(star_info)
                df_gcns.iloc[index, df_gcns.columns.get_loc('Uvel50')] = V[0]
                df_gcns.iloc[index, df_gcns.columns.get_loc('Vvel50')] = V[1]
                df_gcns.iloc[index, df_gcns.columns.get_loc('Wvel50')] = V[2]
            star_names.append(star)

    df_new = pd.DataFrame({
                           'proper' : star_names,
                           'dist': light_years_per_parsec * 1000 * df_gcns['Dist50'].to_numpy(),
                           'x': light_years_per_parsec * df_gcns['xcoord50'].to_numpy(),
                           'y': light_years_per_parsec * df_gcns['ycoord50'].to_numpy(),
                           'z': light_years_per_parsec * df_gcns['zcoord50'].to_numpy(),
                           'vx': kms_to_lyc * df_gcns['Uvel50'].to_numpy(),
                           'vy': kms_to_lyc * df_gcns['Vvel50'].to_numpy(),
                           'vz': kms_to_lyc * df_gcns['Wvel50'].to_numpy()
                          })
    
    stars = find(df_new, 5.0)
    draw_interval(stars, -20, 100, 'GCNS_test', lw = 1.0, fs = 3)

if __name__ == '__main__':
    HYG_experiment('./HYG.csv')
    GCNS_experiment_with_patch('./GCNS.csv')
