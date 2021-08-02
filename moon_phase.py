""" Model spectral lunar irradiance at the Earth's surface for ocean ecosystems

This model calculates spectral lunar irradiance at a user specified time
and location at the surface of the Earth. Skyfield (Rhodes 2019) is used
to calculate time-varying geometry of the sun, Moon, and Earth relative
to latitude and longitude on the Earth. Skyfield positions are based on
JPL Plantary and Lunar Development Ephemeris 440 (Park et al. 2021)
Lunar irradiance is calculated using phase-corrected top of atmosphere
(TOA) incident lunar irradiances from the Miller and Turner (2009) model
and look-up tables. Transmission of TOA lunar irradiance through just below the sea surface
uses the RADTRAN spectral irradiance model for marine systems (Gregg and Carder 1990) with
cloud shortwave radiative transfer modifications described for OASIM
(Gregg 2002, Gregg and Casey 2009) and WASI (Gege 2012). Cloud shortwave radiative transfer
modifications use Slingo's (1989) delta-Eddington approximation for
two-stream spectral cloud transmission.

Gege, P., 2012. Analytic model for the direct and diffuse components of downwelling
    spectral irradiance in water. Appl. Opt. 51, 1407–1419.
    https://doi.org/10.1364/AO.51.001407
Gregg, W.W., Carder, K.L., 1990. A simple spectral solar irradiance model for
    cloudless maritime atmospheres. Limnol. Oceanogr. 35, 1657–1675.
    https://doi.org/10.4319/lo.1990.35.8.1657
Gregg, W.W., Casey, N.W., 2009. Skill assessment of a spectral ocean-atmosphere
    radiative model. J. Mar. Syst. 76, 49–63.
    https://doi.org/10.1016/j.jmarsys.2008.05.007
Miller, S.D., Turner, R.E., 2009. A dynamic lunar spectral irradiance data set for
    NPOESS/VIIRS day/night band nighttime environmental applications.
    IEEE Trans. Geosci. Remote Sens. 47, 2316–2329.
    https://doi.org/10.1109/TGRS.2009.2012696
Park, R.S., Folkner, W.M., Williams, J.G., Boggs, D.H., 2021. The JPL Planetary
    and Lunar Ephemerides DE440 and DE441. Astron. J. 161, 105.
    https://doi.org/10.3847/1538-3881/abd414
Rhodes, B., 2019. Skyfield: High precision research-grade positions for planets
    and Earth satellites generator. Version 1.39. ascl:1907.024.
    https://rhodesmill.org/skyfield/
    https://github.com/skyfielders/python-skyfield
Slingo, A., 1989. A GCM Parameterization for the Shortwave Radiative Properties of
    Water Clouds. J. Atmos. Sci. 46, 1419–1427.
    https://doi.org/10.1175/1520-0469(1989)046<1419:AGPFTS>2.0.CO;2

"""

# import sys
import datetime
import pytz
import numpy as np
import pandas as pd
from skyfield.api import load
from skyfield.framelib import ecliptic_frame
from skyfield.api import N, E, wgs84

# list_of_args = sys.argv
# print(list_of_args[0])
# print(list_of_args[1])

# Input *TEST*
longitude = np.array([-123, -123.0])
latitude = np.array([47.4, 47.4])
year = np.array([2021, 2021])
hour = np.array([22, 20])
minute = np.array([30, 00])
month = np.array([7, 7])
day = np.array([23, 23])
surface_pressure = np.array([1009, 1009])
water_vapor = np.array([1.42])
wind_speed = np.array([0])
alpha = np.array([1.14])
air_mass = np.array([10])
relative_humidity = np.array([50])
ozone = "vanheuklon"
adj_to_utc = True
only_toa = False
cloud_modification = True
visual_range = np.array([8])
aerosol_scale_height = np.array([1.0])
cloud_droplet_radius = np.array([11.8])
liquid_water_path = np.array([5.0])
spectral_path = 'lunar_irrad_1AU_MeanME_350_700.csv'
radtran_params_path = 'radtran_params.csv'

# Set constants
km_per_au = 149598073.0
mean_earthsun_dist = 149598022.6071 / km_per_au
mean_earthmoon_dist = 384400.0 / km_per_au
radius_earth = 6378.140 / km_per_au
water_refraction = 1.341
CONST = (1/(6.626176E-34*299792458.0))*1.0E-10
mWum_to_Wnm = 1.0E-6

# Load TOA mean lunar spectral irradiance table, distance variables (Miller and Turner 2009)
spectral_df = pd.read_csv(spectral_path, header = 0)
w_v = spectral_df['wavelength'].drop_duplicates()
n_wavelength = len(w_v)
spectral_df['lunar_phase'] = np.repeat(np.arange(0.0,181.0,1.0), n_wavelength)

# Load Gregg and Carder (1990) and Slingo (1989) parameters
radtran_df = pd.read_csv(radtran_params_path, header = 0)

# Load JPL Planetary and Lunar Development Ephemeris 440 (Park et al. 2021)
eph = load('de440.bsp')
sun, moon, earth = eph['sun'], eph['moon'], eph['earth']

# Set up time zones
utc = pytz.timezone("UTC")

# Load timescale
ts = load.timescale()

# Output
hhour = []
t = []
phase = []
altitude = []
azimuth = []
moon_zenith = []
earthmoon_dist = []
earthsun_dist = []
moonsun_dist = []
lunar_toa = np.empty(shape = [len(longitude), n_wavelength])

# Calculate positions, distances, phases, zenith, and altitude
for ii in range(len(longitude)):
    
    # Set UTC Time
    hhour.append(hour[ii] + minute[ii] / 60.0)
    dt = datetime.datetime(year[ii], month[ii], day[ii], hour[ii], minute[ii])
    
    if adj_to_utc:
        # Correction based on eastings (west longitude = negative)
        dt = dt + datetime.timedelta(hours = abs(longitude[ii]) / 15.0)
        dt  = dt.replace(tzinfo = utc)
        
    t.append(ts.utc(dt.year, dt.month, dt.day, dt.hour, dt.minute))

    # Moon Phase (with waxing and waning - 0-360 degrees)
    e = earth.at(t[ii])
    _, slon, _ = e.observe(sun).apparent().frame_latlon(ecliptic_frame)
    _, mlon, _ = e.observe(moon).apparent().frame_latlon(ecliptic_frame)
    phase.append((mlon.degrees - slon.degrees) % 360.0)

    # Earth (Observer)-Moon Distance
    loc = earth + wgs84.latlon(latitude[ii] * N, longitude[ii] * E)
    astrometric = loc.at(t[ii]).observe(moon)
    alt, az, d_moon = astrometric.apparent().altaz()
    altitude.append(alt)
    azimuth.append(az)
    earthmoon_dist.append(d_moon)

    # Moon zenith angle
    moon_zenith.append(90.0-alt.degrees)
    moon_cos_zenith = np.cos(np.deg2rad(moon_zenith[ii]))

    # Earth-Sun Distance
    alt_sun, az_sun, d_sun = loc.at(t[ii]).observe(sun).apparent().altaz()

    # Moon-Sun Distance
    moonsun_dist.append(np.sqrt(sum((sun.at(t[ii]).position.au - moon.at(t[ii]).position.au)**2.0)))
    earthsun_dist.append(d_sun)

    # MT2009 Moon Phase (no waxing/waning - 0-180 degrees)
    mt_phase = abs(phase[ii] - 180.0)
    cos_phase_angle = np.cos(np.deg2rad(mt_phase))

    # MT2009 scaling factor
    t_1 = mean_earthsun_dist**2.0 + mean_earthmoon_dist**2.0 + 2.0*mean_earthmoon_dist*mean_earthsun_dist*cos_phase_angle
    t_2 = d_sun.au**2.0 + d_moon.au**2.0 + 2.0*d_moon.au*d_sun.au*cos_phase_angle
    t_3 = ((mean_earthmoon_dist - radius_earth) / (d_moon.au - radius_earth))**2.0
    moon_etc_factor = (t_1 / t_2) * t_3

    # MT2009 Lunar irradiance at phase with fractional phase degree weighting
    frac = mt_phase % 1
    lunar_toa[ii,:] = moon_etc_factor * ((1.0-frac)*spectral_df['E'][np.floor(mt_phase) == spectral_df['lunar_phase']].to_numpy() +
     (frac)*spectral_df['E'][np.ceil(mt_phase) == spectral_df['lunar_phase']].to_numpy())

    # (GC 13) Geometric airmass path length for zenith <= 90 degrees(Karsten and Young 1989)
    if moon_zenith[ii] <= 90:
        geometric_airmass = 1.0 / (moon_cos_zenith + 0.50572 * (96.07995 - moon_zenith[ii])**-1.6364)
    else:
        geometric_airmass = np.sqrt((707.89 * moon_cos_zenith)**2.0 + 2.0*707.89 + 1.0) - 707.89 * moon_cos_zenith

    # (GC 14) Effective airmass of ozone (Paltridge and Platt 1976)
    amoz = (1+22/6370)/(moon_cos_zenith**2 + 2*(22/6370))**0.5

    # (GC 16) Pressure-corrected airmass
    pressure_corr_airmass = geometric_airmass * (surface_pressure[ii]/1013)

    # Total ozone scale height (Van Heuklon 1979)
    if ozone == "vanheuklon":
        if latitude[ii] > 0:
            oz_par = np.array([0.15, 0.04, -30.0, 3.0, 1.28])
        else:
            oz_par = np.array([0.1, 0.03, 152.625, 2.0, 1.25])

        if (longitude[ii] > 0) & (latitude[ii] > 0):
            oz_par = np.append(oz_par, [20.0])
        elif (longitude[ii] > 0) & (latitude[ii] <= 0):
            oz_par = np.append(oz_par, [-75.0])
        else:
            oz_par = np.append(oz_par, [0.0])

        total_ozone = 0.235+(oz_par[0] + oz_par[1]*np.sin(0.9865*(dt.timetuple().tm_yday+oz_par[2])*0.017453)+0.02*np.sin(0.017453*(oz_par[3])*(oz_par[5])))*(np.sin(oz_par[4]*0.017453*latitude[ii]))**2.0
    else:
        total_ozone = ozone

    # (17) Transmission: total ozone (absorption coefficients from Inn and Tanaka [1953])
    transmission_toz = np.exp(-1*radtran_df['a_ozone']*amoz*total_ozone)

    # (18) Transmission: uniform mixed gasses (Bird and Riordan 1986)
    transmission_umg = np.exp(-1.41*radtran_df['a_umg']*pressure_corr_airmass/(1+118.3*radtran_df['a_umg']*pressure_corr_airmass)**0.45) 
  
    # (19) Transmission: Water (Bird and Riordan 1986)
    transmission_water = np.exp(-0.2385*radtran_df['a_water']*water_vapor*geometric_airmass/(1+20.07*radtran_df['a_water']*water_vapor*geometric_airmass)**0.45)

    # (35) Asymmetry parameter for scattering phase function
    if alpha < 0.0:
        asym <- 0.82
    elif alpha > 1.2:
        asym = 1.2
    else:
        asym = -0.1417*alpha+0.82

    # Forward scattering probability
    b_3 = np.log(1-asym)
    b_1 = b_3*(1.459 + b_3*(0.1595 + b_3*0.4129))
    b_2 = b_3*(0.0783 + b_3*(-0.3824 - b_3*0.5874))
    fsp = 1 - 0.5*np.exp((b_1 + b_2*moon_cos_zenith)*moon_cos_zenith)
  
    # (36) Single-scattering albedo
    omega_a = (-0.0032*air_mass + 0.972) * np.exp(3.06*relative_humidity*1.0E-4)

    # (39-43) Sea-foam reflectance as a function of wind stress (Trenberth et al. 1989)
    if wind_speed <= 4.0:
        rho_f = 0.0
    elif (wind_speed > 4.0) & (wind_speed <= 7.0):
        drag = (0.62 + 1.56 * wind_speed**-1.0)*1E-3
        rho_f = (2.2E-5* 1.2E3 * drag * wind_speed**2)-4E-4
    else:
        drag = (0.49 + 0.065 * wind_speed) * 1E-3
        rho_f = (4.5E-5 * 1.2E3 * drag - 4E-5) * wind_speed**2.0

    # Future foam spectral correction 
    # foam_spectral_correction <- 1
    # rho_f <- rho_f * foam_spectral_correction

    # (44-47) Direct specular reflectance
    rho_dsp = 0

    if moon_zenith[ii] >= 40:
        if wind_speed > 2:
            rho_dsp = 0.0253 * np.exp((-7.14E-4 * wind_speed + 0.0618)*(moon_zenith[ii] - 40.0))
        else:
            moon_z_rad = np.deg2rad(moon_zenith[ii])
            refracted_angle = np.arcsin(np.sin(moon_z_rad) / water_refraction)
            rho_dsp = 0.5 * (0.5 * np.sin(moon_z_rad)**2.0 / np.sin(moon_z_rad + refracted_angle)**2+ np.tan(moon_z_rad - refracted_angle)**2.0 / np.tan(moon_z_rad + refracted_angle)**2.0)

    # Diffuse specular reflectance (Burt 1954)
    if wind_speed > 4:
        rho_ssp = 0.057
    else:
        rho_ssp = 0.066

    # (37) Surface reflectance (direct)
    rho_direct = rho_dsp + rho_f
    
    # (38) Surface reflectance (diffuse)
    rho_diffuse = rho_ssp + rho_f

    # Cloudless sky = 100% transmisson
    transmission_cloud_direct = 1.0
    transmission_cloud_diffuse = 1.0
    transmission_cloud_total_direct = 1.0

    if cloud_modification:
        # Cloud Optical Depth (Slingo 1989 - Eqn. 1) with a mean cloud droplet radius
        tau_clouds = liquid_water_path * (radtran_df['cloud_a_i'] + radtran_df['cloud_b_i']/cloud_droplet_radius)

        # Single Scatter Albedo (Slingo 1989 - Eqn. 2)
        omega_clouds = 1.0 - radtran_df['cloud_c_i'] - radtran_df['cloud_d_i'] * cloud_droplet_radius

        # Asymmetry parameter (g) (Slingo 1989 - Eqn. 3)
        asym_clouds = radtran_df['cloud_e_i'] + radtran_df['cloud_f_i'] * cloud_droplet_radius

        # Fraction of scattered diffuse radiation scattered backwards (Slingo 1989 - Eqn. 6)
        beta_diffuse = 3.0/7.0 * (1.0 - asym_clouds)

        # Fraction of direct radiation scattered backwards (Slingo 1989 - Eqn. 7)
        beta_direct = 0.5 - 3.0 * moon_cos_zenith * asym_clouds / (4.0 * (1.0 + asym_clouds))

        # Fraction of scattered direct flux emerging at zenith angles close to incident beam (Slingo 1989 - Eqn. 8)
        f_foreward = asym_clouds**2.0

        # Reciprocal of cosine of diffuse upward (Slingo 1989 - Eqn. 9)
        u_1 = 7.0/4.0

        # Reciprocal of cosine of diffuse downward (Slingo 1989 - Eqn. 10)
        u_2 = 7.0/4.0*(1.0 - (1.0 - omega_clouds)/(7.0 * omega_clouds * beta_diffuse))

        # (Slingo 1989 - Eqn. 11)
        a_1 = u_1 * (1.0 - omega_clouds * (1.0 - beta_diffuse))

        # (Slingo 1989 - Eqn. 12)
        a_2 = u_2 * omega_clouds * beta_diffuse

        # (Slingo 1989 - Eqn. 13)
        a_3 = (1.0 - f_foreward) * omega_clouds * beta_direct

        # (Slingo 1989 - Eqn. 14)
        a_4 = (1.0 - f_foreward) * omega_clouds * (1-beta_direct)

        # (Slingo 1989 - Eqn. 15)
        epsilon = np.sqrt(abs(a_1**2 - a_2**2))

        # (Slingo 1989 - Eqn. 16)
        M = a_2 / (a_1 + epsilon)

        # (Slingo 1989 - Eqn. 17)
        EE = np.exp(-1.0 * epsilon * tau_clouds)

        # (Slingo 1989 - Eqn. 18)
        gamma_1 = ((1.0 - omega_clouds * f_foreward) * a_3 - moon_cos_zenith *(a_1*a_3 + a_2*a_4))/((1.0-omega_clouds*f_foreward)**2.0 - epsilon**2.0 * moon_cos_zenith**2.0)

        # (Slingo 1989 - Eqn. 19)
        gamma_2 = -1.0 * ((1.0 - omega_clouds * f_foreward) * a_4 - moon_cos_zenith * (a_1*a_4 + a_2*a_3)) / ((1 - omega_clouds * f_foreward)**2.0 - epsilon**2.0 * moon_cos_zenith**2.0)

        # Direct cloud transmission (Slingo 1989 - Eqn. 20)
        transmission_cloud_direct = np.exp(-1.0 * (1.0 - omega_clouds * f_foreward) * tau_clouds / moon_cos_zenith)

        # Diffuse reflectivity for diffuse incident radiation (Slingo 1989 - Eqn. 21)
        ref_diffuse = M*(1.0 - EE**2.0)/(1 - EE**2.0 * M**2.0)

        # Diffuse transmissiviity for diffuse incident radiation (Slingo 1989 - Eqn. 22)
        transmission_cloud_diffuse = EE * (1.0 - M**2.0)/(1 - EE**2 * M**2)

        # Diffuse reflectivity for direct incident radiation (Slingo 1989 - Eqn. 23)
        ref_direct = -1.0 * gamma_2 * ref_diffuse - gamma_1 * transmission_cloud_diffuse * transmission_cloud_direct + gamma_1

        # Diffuse transmissivity for direct incident radiation
        transmission_cloud_direct_incident = -1.0 * gamma_2 * transmission_cloud_diffuse - gamma_1 * transmission_cloud_diffuse * ref_diffuse + gamma_2 * transmission_cloud_direct

        # Total transmission to direct
        transmission_cloud_total_direct = transmission_cloud_direct + transmission_cloud_direct_incident

        # Ignore aerosol and Rayleign transmission if clouds are used
        transmission_rayleigh = 1
        transmission_aerosol = 1
        transmission_as = 0
        transmission_aa = 1
    else:
        # (28) Kochsmieder formula for aerosol beam attenuation at 550 (Fitzgerald 1989)
        ca_550 = 3.91 / visual_range

        # (29) Aerosol optical thickness
        taua_550 = ca_550 * aerosol_scale_height

        # (Rearranging 27) - Calculate turbidity coefficient (beta)
        beta = taua_550 * 0.55**alpha
        
        # (27) Wavelength-dependent aerosol optical thickness
        tau_a = beta * (w_v**-1.0) 

        # (15) Transmittance: Rayleigh (Bird and Riordan 1986)
        transmission_rayleigh = np.exp(-1.0 * pressure_corr_airmass / (w_v**4.0 * (115.6406 - 1.335 / w_v**2.0)))

        # (30) Transmission: Aerosol 
        transmission_aerosol = np.exp(-1.0 * geometric_airmass * tau_a)

        # (7) Transmission: Aerosol scattering
        transmission_as = np.exp(-1.0 * omega_a * tau_a * geometric_airmass)

        # (5) Transmission: Aerosol absorption
        transmission_aa = np.exp(-1.0 * (1.0 - omega_a) * tau_a * geometric_airmass)

    # (9) Direct downwelling irradiance
    direct_surface_wm2 = lunar_toa[ii,:] * moon_etc_factor * moon_cos_zenith * transmission_aerosol * transmission_umg * transmission_water * transmission_toz * transmission_rayleigh * transmission_cloud_total_direct * mWum_to_Wnm
    direct_subsurface_wm2 = direct_surface_wm2 * (1.0 - rho_direct)

    # (4) Diffuse component from Rayleigh scattering
    diffuse_rayleigh = lunar_toa[ii,:] * moon_etc_factor * moon_cos_zenith *transmission_toz * transmission_umg * transmission_water * transmission_aa * (1.0 - transmission_rayleigh**0.95)*0.5 * mWum_to_Wnm
    diffuse_scattering = lunar_toa[ii,:] * moon_etc_factor * moon_cos_zenith * transmission_toz * transmission_umg * transmission_water * transmission_rayleigh * (1.0 - transmission_as) * fsp * transmission_cloud_diffuse * mWum_to_Wnm
  
    # (10) Diffuse downwelling irradiance
    diffuse_surface_wm2 = diffuse_rayleigh + diffuse_scattering
    diffuse_subsurface_wm2 = diffuse_surface_wm2 * (1.0 - rho_diffuse)

    print("Moon zenith (degrees): %s"  %moon_zenith[ii])
    print("Moon phase (deg): %s"  %mt_phase)
    print("Geometric air mass: %s"  %geometric_airmass)
    print("Press. corr. air mass: %s"  %pressure_corr_airmass)
    print("Effective airmass of ozone: %s"  %amoz)
    print("TOA Irradiance: %s"  %np.sum(lunar_toa[ii,:]*mWum_to_Wnm))
    print("Diffuse surface: %s"  %np.sum(diffuse_surface_wm2))
    print("Direct surface: %s" %np.sum(direct_surface_wm2))
    print("Moon factor: %s" %moon_etc_factor)
    print("Transmission (ozone) %s" %np.min(transmission_toz))
    print("Transmission (umg) %s" %np.min(transmission_umg))
    print("Transmission (water) %s" %np.min(transmission_water))
    print("Transmission (aa) %s" %np.min(transmission_aa))
    print("Transmission (Rayleigh) %s" %np.min(transmission_rayleigh))
    print("Transmission (cloud direct) %s" %np.min(transmission_cloud_total_direct))
    print("Transmission (cloud diffuse) %s" %np.min(transmission_cloud_diffuse))
    print("Moon cosine zenith %s" %moon_cos_zenith)
    print("----")

    # Convert watts to photon flux density
    diffuse_surface_pfd = w_v * diffuse_surface_wm2 * CONST
    diffuse_subsurface_pfd = w_v * diffuse_subsurface_wm2 * CONST
  
    direct_surface_pfd = w_v * direct_surface_wm2 * CONST
    direct_subsurface_pfd = w_v * direct_subsurface_wm2 * CONST

if only_toa:
    print("Moon zenith (degrees): %s"  %moon_zenith[ii])
    print("Moon phase (deg): %s"  %mt_phase)
    print("Geometric air mass: %s"  %geometric_airmass)
    print("Press. corr. air mass: %s"  %pressure_corr_airmass)
    print("Effective airmass of ozone: %s"  %amoz)
