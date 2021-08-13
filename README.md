# Description

*** IN DEVELOPMENT -- August 12, 2021 ***

*moontran* is a coupled lunar spectral irradiance and empirical atmospheric radiative transfer model for estimating lunar spectral irradiance at the surface of the earth, and just below the sea surface, for a specified time and location (latitude and longitude), given environmental conditions.

The model estimates top-of-atmosphere (TOA) lunar spectral irradiance using Miller and Turner's (2009) lunar spectral irradiance model but with astronomical terms (e.g. solar and lunar coordinates, moon phase angle, true moon zenith angle from the earth's surface) calculated using [Skyfield](https://github.com/skyfielders/python-skyfield/) (Rhodes, 2019) with JPL Planetary Development Ephemeris 440 (Park et al., 2021). The RADTRAN model (Gregg and Carder, 1990) with a cloud-cover modification (Gregg, 2002; Gregg and Casey, 2009) is then used to estimate spectral transmission of TOA lunar irradiance to the earth's surface. The cloud-cover modification is Slingo's (1989) delta-Eddington approximation for two-stream spectral clound transmission.

## Data sets

The *moontran* repository includes files containing data sets from several sources that are necessary to estimate lunar spectral irradiances. The contents of each file are described here:

<b>lunar_irrad__1AU_MeanME_350_700.csv</b>
<br>Spectral top of atmosphere (TOA) lunar irradiances for 350-700 nm wavelengths (1 nm resolution); a subset of the 202-2800 nm TOA spectral irradiances from Miller and Turner (2009).
- w_v_um: Wavelength in micrometers.
- lunar_E_mwm2um: Top of atmosphere lunar irradiance in mW<sup>-2</sup>&micro;m<sup>-1</sup>.

<b>radtran_params.csv</b>
<br>Spectral absorption coefficients and model parameters for the atmosphere and clouds.
- w_v_nm: Wavelength in nanometers
- solar_E_wm2um: Top of atomsphere solar irradiance in Wm<sup>-2</sup>&micro;m<sup>-1</sup> (not used in moontran; Bird and Riordan, 1986).
- a_ozone: Atmospheric absorption coefficient for ozone (Bird and Riordan, 1986).
- a_umg: Atmospheric absorption coefficient for uniformly mixed gasses (Bird and Riordan, 1986).
- a_water: Atmospheric absorption absorption coefficient for water vapor (Bird and Riordan, 1986).
- cloud_a_i: _a<sub>i</sub>_ parameter for delta-Eddington approximation for two-stream approach for cloud transmission (Slingo, 1989).
- cloud_b_i: _b<sub>i</sub>_ parameter for delta-Eddington approximation for two-stream approach for cloud transmission (Slingo, 1989).
- cloud_c_i: _c<sub>i</sub>_ parameter for delta-Eddington approximation for two-stream approach for cloud transmission (Slingo, 1989).
- cloud_d_i: _d<sub>i</sub>_ parameter for delta-Eddington approximation for two-stream approach for cloud transmission (Slingo, 1989).
- cloud_e_i: _e<sub>i</sub>_ parameter for delta-Eddington approximation for two-stream approach for cloud transmission (Slingo, 1989).
- cloud_f_i: _f<sub>i</sub>_ parameter for delta-Eddington approximation for two-stream approach for cloud transmission (Slingo, 1989).

# Usage

Running *moontran* requires using a JSON file that specifies all of the arguments that are needed to run the model. For an example of a JSON input file, see the 'moontran_input.json' file in the repository.

In the JSON file, arguments that specify time (e.g., month, hour), space (latitude, longitude), and environmental conditions (e.g., water_vapor, wind_speed) must have length equal to one (1) or the number of observations (n_obs). Arguments that specify paths to data sources for constants, user-specified constants (e.g., cloud_droplet_radius), and mode for the model (e.g., adj_to_utc, cloud_modification) must have length equal to one (1). 

Tip for R Users: The 'rjson' package converts lists to JSON strings that can be written to a file and used by moontran().

# References

- Bird, R.E., Riordan, C., 1986. Simple solar spectral model for direct and diffuse irradiance on horizontal and tilted planes at the Earth’s surface for cloudless atmospheres. J. Clim. Appl. Meteorol. 25, 87–97. https://doi.org/10.1175/1520-0450(1986)025<0087:SSSMFD>2.0.CO;2
- Gregg, W.W., Carder, K.L., 1990. A simple spectral solar irradiance model for cloudless maritime atmospheres. Limnol. Oceanogr. 35, 1657–1675. https://doi.org/10.4319/lo.1990.35.8.1657
- Gregg, W.W., Casey, N.W., 2009. Skill assessment of a spectral ocean-atmosphere radiative model. J. Mar. Syst. 76, 49–63. https://doi.org/10.1016/j.jmarsys.2008.05.007<br>
- Miller, S.D., Turner, R.E., 2009. A dynamic lunar spectral irradiance data set for NPOESS/VIIRS day/night band nighttime environmental applications. IEEE Trans. Geosci. Remote Sens. 47, 2316–2329. https://doi.org/10.1109/TGRS.2009.2012696<br>
- Park, R.S., Folkner, W.M., Williams, J.G., Boggs, D.H., 2021. The JPL Planetary and Lunar Ephemerides DE440 and DE441. Astron. J. 161, 105. https://doi.org/10.3847/1538-3881/abd414<br>
- Rhodes, B., 2019. Skyfield: High precision research-grade positions for planets and Earth satellites generator. Version 1.39. ascl:1907.024. https://github.com/skyfielders/python-skyfield<br>
- Slingo, A., 1989. A GCM Parameterization for the Shortwave Radiative Properties of Water Clouds. J. Atmos. Sci. 46, 1419–1427. [https://doi.org/10.1175/1520-0469(1989)046<1419:AGPFTS>2.0.CO;2](https://doi.org/10.1175/1520-0469(1989)046<1419:AGPFTS>2.0.CO;2)<br>

# Legal Disclaimer

This repository is a software product and is not official communication of the National Oceanic and Atmospheric Administration (NOAA), or the United States Department of Commerce (DOC). All NOAA GitHub project code is provided on an 'as is' basis and the user assumes responsibility for its use. Any claims against the DOC or DOC bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation, or favoring by the DOC. The DOC seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by the DOC or the United States Government.
