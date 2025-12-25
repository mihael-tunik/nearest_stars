import numpy as np

Gaia_EDR3_ID_to_common = {
                          5853498713190525696: 'Alpha Centauri A/B', 
                          4472832130942575872: 'Barnards Star',
                          3864972938605115520: 'Wolf 359',
                          762815470562110464 : 'Lalande 21185',
                          4810594479418041856: 'Kapteyns Star',
                          1926461164913660160: 'Ross 248',
                          1129149723913123456: 'Gliese 445',
                          2552928187080872832: 'Van Maanens',
                          6131272458806020736: 'L 326-21',
                          6697578465310949376: 'Gliese 783',
                          2011565220332867584: 'HIP 117795',
                          5231593594752514304: 'HD 95123B',
                          1952802469918554368: '2MASS J2146+3813',
                          2452378776434477184: 'Tau Ceti'
                         }

light_years_per_parsec = 3.2615637769
kms_to_lyc     = 0.000333564095198152
AU_to_km       = 149597870.7
year_to_s      = 31557600
speed_to_kms   = AU_to_km / year_to_s # = 4.74047 milliarcseconds / year to km / s

# from article 'Reconsidering the Galactic coordinate system'
N_B = [[-0.066988739410, -0.872755765850, -0.483538914637],
       [ 0.492728466081, -0.450346958020,  0.744584633279],
       [-0.867600811149, -0.188374601732,  0.460199784785]]

X_0 = [[ 0.999925679496,  0.011181483239,  0.004859003772],
       [-0.011181483221,  0.999937484893, -0.000027170294],
       [-0.004859003815, -0.000027162595,  0.999988194602]]

# compute star velocity (U, V, W) by RA, DE, Plx, proper motion and radial velocity
# requires transformation between the equatorial and galactic systems A_G matrix inside
# see 'The Hipparcos and Tycho Catalogues' and 'Reconsidering the Galactic coordinate system' for details
def get_uvw(star_info):
    alpha = star_info['RAdeg'] * (np.pi / 180)
    delta, omega = star_info['DEdeg'] * (np.pi / 180), star_info['Plx']
    mu_a, mu_d, vr  = star_info['pmRA'], star_info['pmDE'], star_info['RV']
    
    v = [speed_to_kms * mu_a / omega, speed_to_kms * mu_d / omega, vr ]
    A_G_t = np.array(N_B) @ np.array(X_0) # more precise?
    A = [[-np.sin(alpha), -np.sin(delta)*np.cos(alpha)  , np.cos(delta)*np.cos(alpha) ],
         [np.cos(alpha) , -np.sin(delta)*np.sin(alpha)  , np.cos(delta)*np.sin(alpha) ],
         [0             ,  np.cos(delta)                , np.sin(delta)               ]]
    return A_G_t @ np.array(A) @ np.array(v)
