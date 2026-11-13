"""Modern orbital elements from JPL, used only for the 'modern' ablation variants.
Source: E.M. Standish & J.G. Williams, 'Keplerian Elements for Approximate Positions of the Major Planets',
        Table 2a (valid 3000 BC - 3000 AD), https://ssd.jpl.nasa.gov/planets/approx_pos.html
Values are evaluated at 1500 CE (T = -5 centuries from J2000), the epoch of Tantrasangraha.
"""
T_1500 = -5.0

# name: (a0, a_rate, e0, e_rate, perihelion_longitude0, perihelion_rate)   [au, au/cy, -, /cy, deg, deg/cy]
TABLE_2A = {
    "earth":   (1.00000018, -0.00000003, 0.01673163, -0.00003661, 102.93005885, 0.31795260),
    "venus":   (0.72332102, -0.00000026, 0.00676399, -0.00005107, 131.76755713, 0.05679648),
    "jupiter": (5.20248019, -0.00002864, 0.04853590,  0.00018026,  14.27495244, 0.18199196),
    "saturn":  (9.54149883, -0.00003065, 0.05550825, -0.00032044,  92.86136063, 0.54179478),
}

# Mean longitude L (deg, deg/century) from Table 2a, and the extra Table 2b terms (b, c, s, f) for Jupiter and Saturn,
# which include the Jupiter-Saturn 'great inequality'.
MEAN_LONGITUDE = {"earth": (100.46691572, 35999.37306329), "venus": (181.97970850, 58517.81560260),
                  "jupiter": (34.33479152, 3034.90371757), "saturn": (50.07571329, 1222.11494724)}
TABLE_2B = {"jupiter": (-0.00012452, 0.06064060, -0.35635438, 38.35125000),
            "saturn":  ( 0.00025899, -0.13434469, 0.87320147, 38.35125000)}

def mean_longitude_J2000(name, T, with_2b=True):
    """Mean longitude in the J2000 ecliptic frame at T centuries from J2000 (array OK)."""
    import numpy as np
    L0, dL = MEAN_LONGITUDE[name]
    L = L0 + dL * T
    if with_2b and name in TABLE_2B:
        b, c, s, f = TABLE_2B[name]
        L = L + b * T**2 + c * np.cos(np.radians(f * T)) + s * np.sin(np.radians(f * T))
    return L % 360

def elements(name, T=T_1500):
    a0, da, e0, de, p0, dp = TABLE_2A[name]
    return {"a": a0 + da * T, "e": e0 + de * T, "perihelion_J2000": p0 + dp * T}

if __name__ == "__main__":
    for n in TABLE_2A:
        el = elements(n)
        print(f"{n:8s} a = {el['a']:.5f} au   e = {el['e']:.5f}   2e = {2*el['e']:.4f}   "
              f"perihelion (J2000 frame) = {el['perihelion_J2000']:.2f} deg")