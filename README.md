## About project
This repository contains the code for reproducing classic result from stellar kinematics: computed distances from the nearest stars to the solar system.

<img src="nearest_stars_GCNS_test.svg" alt="drawing" width="100%"/>

Code uses two main data sources with rectangular coordinates and stellar velocity vectors.
- HYG star database (110K stars), contains data from Hipparcos, Yale Bright Star and Gliese catalogues, includes distant stars >1000 light years away;
- GCNS star database (330K stars), contains stars within 100 parsec radius from the solar system;

Surprisingly, GCNS star database do not contain some nearby stars which is not as good for purposes of this project.
In this code, it is shown how some velocity vectors can be recovered from auxiliary data sources for GCNS.

## A bit of mathematics
Lets choose rectangular coordinate system with center in the sun.
Stellar velocity vectors are considered constant.

Then distance to the sun from any star with position $(P, Q, R)$ and vector velocity $(A, B, C)$ 
in our linear model is given by

$$d = \sqrt{(P+A \cdot t)^2 + (Q+B \cdot t)^2 + (R+C \cdot t)^2}.$$

This equation is actually a upper branch of _hyperbola_ in $(t, d)$-coordinates, since equation
$$d^2 - (A^2 + B^2 + C^2) \cdot t^2 - 2 \cdot (A \cdot P + B \cdot Q + C \cdot R) \cdot t = P^2 + Q^2 + R^2$$

is equivalent to

$$d^2 - (\sqrt{A^2 + B^2 + C^2} \cdot t + \frac{A \cdot P + B \cdot Q + C \cdot R}{\sqrt{A^2 + B^2 + C^2}})^2 = (P^2 + Q^2 + R^2) - \frac{(A \cdot P + B \cdot Q + C \cdot R)^2}{A^2 + B^2 + C^2}$$

by "extracting the square" trick.

Last formula means that $d$ is minimized when 
$$t = -\frac{A \cdot P + B \cdot Q + C \cdot R}{A^2 + B^2 + C^2}.$$

Given star *positions* and *velocities*, let's use this expression for testing.

> It is easy to prove (by Cauchy-Schwarz inequality) that $(P^2 + Q^2 + R^2) - \frac{(A \cdot P + B \cdot Q + C \cdot R)^2}{A^2 + B^2 + C^2}$ is non-negative and verify that our formula for distance is always upper branch of hyperbola for non-zero velocity in general case when sun does not meet with other star on trajectory.

## Example
Find HYG .csv dataset and place it in the same folder. 
```python
if __name__ == '__main__':
    HYG_experiment('./HYG.csv')
    # GCNS_experiment_with_patch('./GCNS.csv')
```
Then run (details are in the code):
```
python3 run_experiment.py
```
One can also use this code for any dataset with computed rectangular coordinates and velocities.

## How to use in code
In function 
```python
def find(star_table, upper_bound=5.0):
    ...
```
input pandas dataframe should have columns: _'x', 'y', 'z', 'vx', 'vy', 'vz'_ for computations and _'proper'_ column with star names.
And units inside _star\_table_ dataframe should be translated to _light years_ for coordinates and _light years per century_ for velocities

Use it in your code as follows:
```python
stars = find(stars_table, upper_bound) 
draw_interval(stars, t_start, t_end, label, ...) # time should be in thousands of years
```
This will print stars with $t_{min}, d_{min}$ and draw a picture with distances.

## How to restore velocity vectors
Unfortunately, 75% of stars in GCNS data do not have vector velocity (including Barnard's star, one of the nearest to our solar system).
However, the data contains columns with measured right ascention, declination, proper motions and more. 

If $(\alpha, \delta, \pi)$ data was given and proper motions $\mu_{\delta}, \mu_{\alpha^*}$ with radial velocity $v_r$ was measured, then velocity vector in Galactic coordinate system can be computed as follows:
```python
def get_uvw(star_info):
    alpha = star_info['RAdeg'] * (np.pi / 180)
    delta, omega = star_info['DEdeg'] * (np.pi / 180), star_info['Plx']
    mu_a, mu_d, vr  = star_info['pmRA'], star_info['pmDE'], star_info['RV']
    
    v = [speed_to_kms * mu_a / omega, speed_to_kms * mu_d / omega, vr ]
    A_G_t = np.array(N_B) @ np.array(X_0)
    A = [[-np.sin(alpha), -np.sin(delta)*np.cos(alpha)  , np.cos(delta)*np.cos(alpha) ],
         [np.cos(alpha) , -np.sin(delta)*np.sin(alpha)  , np.cos(delta)*np.sin(alpha) ],
         [0             ,  np.cos(delta)                , np.sin(delta)               ]]
    return A_G_t @ np.array(A) @ np.array(v)
```
The main problem is RV is often missed together with vector velocity.
That means dataset should be patched with RV values in order to get more stars for our experiment.

In experiment, _sosdr1.csv_ data from "Survey of surveys I" (11 541 195 stars) and _'rvstdcat.csv'_ dataset (4 813 stars) with averaged $v_r$ from "Gaia DR2 radial velocity standard stars catalog" was used.

## Results (HYG data)
For HYG star database program had found 12 stars with minimum achieved distances less than 5 light years.

<img src="nearest_stars_HYG_test.svg" alt="drawing" width="100%"/>

As an interesting fact, program found a star GJ 2005 with achieved minimum distance 2.42 ly in just 30 000 years from now.
For comparison current minimum distance is about 4.25 light years for Proxima Centauri.

## Results (GCNS data)
Results for GCNS data are expected to be more precize. 
Difference can be seen, for example, by estimates for famous Gl 710 star.
HYG data usage lead to minimal distance estimate equal to 1.22 ly in 1.437 million years. 
However, GCNS data leads to estimate 0.17 ly for in 1.282 million years.

For augmented version of GCNS dataset this program had found 44 stars with minimum achieved distances less than 5 light years.
