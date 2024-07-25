# utils.py
import numpy as np
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import Latitude, Longitude, EarthLocation, AltAz, SkyCoord

def load_star_data(csv_file='sextantCalculate/NavigationalStars.csv'):
    NavStars = np.genfromtxt(csv_file, delimiter=',', dtype=None, encoding='utf-8')
    NavStars = np.delete(NavStars, 0, axis=0)  # delete the label row
    return NavStars

def calculate_position(star_observations, csv_file='sextantCalculate/NavigationalStars.csv'):
    NavStars = load_star_data(csv_file)
    
    numObs = len(star_observations)
    A = np.zeros([numObs, 3])  # initialize the A matrix
    for i in range(numObs):  # populate the A matrix
        
        # unpack observation
        starNum = star_observations[i][0]
        Ho = star_observations[i][1]
        obsTime = Time(star_observations[i][2], format='unix')
        
        # read star's declination
        dec = Latitude(NavStars[starNum][3], unit=u.deg)  # Declination is in the 4th column
        
        # calculate Greenwich Hour Angle (GHA) of star
        ERA = Time.earth_rotation_angle(obsTime, longitude=0.0)  # earth rotation angle at obsTime, also called GHA_star in some books
        SHA = Longitude(360 - NavStars[starNum][2], unit=u.deg)  # Sideral hour angle (RA is in the 3rd column)
        GHA = Longitude(ERA + SHA, unit=u.deg)  # greenwich hour angle 
        
        A[i] = [np.cos(dec)*np.cos(GHA), np.cos(dec)*np.sin(GHA), np.sin(dec)] / np.sin(Ho)
    
    B = np.ones((numObs, 1))  # create B vector
    
    # convert to least squares problem
    pseudoA = np.matmul(np.transpose(A), A)
    pseudoB = np.matmul(np.transpose(A), B)
    
    # solve
    X_star = np.linalg.solve(pseudoA, pseudoB) * u.km  # solve Ax=B

    x, y, z = X_star[0], X_star[1], X_star[2]  # unpack results
    estLat = Latitude(np.arctan2(z, np.hypot(x, y))[0], unit=u.deg)
    estLon = -Longitude(np.arctan2(y, x)[0], unit=u.deg, wrap_angle=180*u.deg)
    
    return estLat.degree, estLon.degree
