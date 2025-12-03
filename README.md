## About project
<img src="nearest_stars_Standard.svg" alt="drawing" width="100%"/>

This repository contains the code for reproducing classic result from stellar kinematics.
Code uses HYG star database (110K stars!) in order to obtain rectangular coordinates and velocity vectors.
Also this repository is notable example of pretty mathematics and datascience.

## How to use
Find HYG .csv dataset and place it in the same folder. Then run (details are in the code):
```
python3 find_nearest_stars.py
```
One can also use this code for any dataset with computed rectangular coordinates and velocities.
In function 
```
def find(star_table, upper_bound=7.0):
    ...
```
Units inside _star\_table_ pandas dataframe should be translated to _light years_ for coordinates and _light years per century_ for velocities.

## A bit of mathematics
Lets choose rectangular coordinate system with center in the sun.
Stellar velocity vectors are considered constant.

Then distance to the sun from any star with position $(P,Q,R)$ and vector velocity $(A, B, C)$ 
in our linear model is given by

$$d = \sqrt{(P+A \cdot t)^2 + (Q+B \cdot t)^2 + (R+C \cdot t)^2}.$$

This equation is actually a upper branch of _hyperbola_, since equation
$$d^2 - (A^2 + B^2 + C^2) \cdot t^2 - 2 \cdot (A \cdot P + B \cdot Q + C \cdot R) \cdot t = P^2 + Q^2 + R^2$$

is equivalent to

$$d^2 - (\sqrt{A^2 + B^2 + C^2} \cdot t + \frac{A \cdot P + B \cdot Q + C \cdot R}{\sqrt{A^2 + B^2 + C^2}})^2 = (P^2 + Q^2 + R^2) - \frac{(A \cdot P + B \cdot Q + C \cdot R)^2}{A^2 + B^2 + C^2}$$

by "extracting the square" trick.

Last formula means that $d$ is minimized when 
$$t = -\frac{A \cdot P + B \cdot Q + C \cdot R}{A^2 + B^2 + C^2}.$$

Given star positions and velocities, let's use this expression for testing.

> It is easy to check (also by "extracting the square") $$(P^2 + Q^2 + R^2) - \frac{(A \cdot P + B \cdot Q + C \cdot R)^2}{A^2 + B^2 + C^2} \geq 0$$ and verify it is exactly upper branch of hyperbola for $(P, Q, R) \neq (0, 0, 0)$ and $(A, B, C) \neq (0, 0, 0)$.

## Results
My programm had found 29 stars with minimum less than 7 light years.

<img src="nearest_stars_Past.svg" alt="drawing" width="48%"/>
<img src="nearest_stars_Long.svg" alt="drawing" width="48%"/>

Star Gl 710 will approach the sun as close as 1.21 ly in 1.437 million years.
Star with HIP identifier 38965 had minimum distance 1.86 ly 1.087 million years ago.

Also program found a star GJ 2005 with achieved minimum distance 2.41 ly in 30 000 years from now.
For compare current minimum distance is about 4 light years for Proxima Centauri.
