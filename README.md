# About

*** IN TESTING ***

moontran _is intended to_ be a coupled lunar spectral irradiance and empirical atmospheric radiative transfer model for estimating lunar spectral irradiance at the surface of the earth for a specified time and position (latitude and longitude).

The model estimates top-of-atmosphere (TOA) lunar spectral irradiance using Miller and Turner's (2009) lunar spectral irradiance model but with astronomical terms (e.g. solar and lunar coordinates, moon phase angle, true moon zenith angle from the earth's surface) estimated using Skyfield (Rhodes 2019) with JPL Planetary Development Ephemeris 440 (Park et al. 2021). The RADTRAN model (Gregg and Carder 1990) with a cloud-cover modification (Gregg 2002, Gregg et al. 2009) is then used to estimate spectral transmission of TOA lunar irradiance to the earth's surface. The cloud-cover modification is Slingo's (1989) delta-Eddington approximation for two-stream spectral clound transmission.
