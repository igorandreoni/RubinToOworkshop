#!/usr/bin/env python

import rubin_nights.dayobs_utils as rn_dayobs
import numpy as np
import healpy as hp
import matplotlib.dates as mdates
import matplotlib.units as munits
from ligo.skymap.io.fits import read_sky_map
from ligo.skymap.postprocess import find_greedy_credible_levels
import ligo.skymap.plot
import getpass
import os
import sys
from astropy.time import Time, TimeDelta
import pandas as pd
import copy
import argparse
import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(description="Run LSST survey sim on one skymap")

parser.add_argument(
    "--save_dir",
    type=str,
    required=False,
    help="The path to save the resulting opsims.",
    default="/scratch/smacbr/testSims/",
)

parser.add_argument(
    "--force_strategy",
    type=str,
    required=False,
    help="The type of strategy to enforce on the resulting simulation. Options are 'gold', 'silver', 'custom1', and 'custom2'",
    default="silver",
)

args = parser.parse_args()

force_strategy = args.force_strategy
save_dir = args.save_dir
superevent = "S251112cm"

username = getpass.getuser()
location = os.getenv("EXTERNAL_INSTANCE_URL", "")

# This is the location of the RUBIN_SIM_DATA_DIR on my science cluster
# TODO: Update this path to the location of rubin sim data on your machine
os.environ["RUBIN_SIM_DATA_DIR"] = "/shares/soares-santos.physik.uzh/rubin_sim_data"

# If you are outside of an RSP, just use USDF and your own USDF-RSP token
# TODO: Update this path
# See https://rsp.lsst.io/guides/auth/creating-user-tokens.html
tokenfile = os.path.join(os.path.expanduser("~"), ".usdf_token")
site = "usdf"

# Configure some paths
ts_config_scheduler_root = None

# TODO: Clone ts_config_scheduler and update the path here
if ts_config_scheduler_root is None:
    # Just make a new clone for ts_config_scheduler
    ts_config_scheduler_root = os.path.join(
        os.path.expanduser("~"), "ts_config_scheduler"
    )
    sys.path.insert(0, os.path.join(os.path.expanduser("~"), "lsst_survey_sim"))

    do_git_stuff = True

assert isinstance(ts_config_scheduler_root, str), "Please set ts_config_scheduler_root"

logger.info("Defining functions to be used for handling the ToO object")


def getSimNightsFromSkymapMetadata():
    # I am leaving this at 10 days for now. You can modify this however you like.
    return 10


def get_90_percent_binary_map(file_path, target_nside=32):
    """
    Reads a FITS skymap and returns a binary map at target_nside
    representing the 90% localization area.
    """
    # 1. Read the sky map
    skymap = read_sky_map(file_path, moc=False)
    prob = skymap[0]

    # 2. Calculate the credible level map
    credible_levels = find_greedy_credible_levels(prob)

    # 3. Create the binary mask at the original resolution (90% area)
    # This is what we do on the summit
    binary_mask_high_res = (credible_levels <= 0.9).astype(float)

    # 4. Resample to target_nside (default 32)
    # This is what we do on the summit
    binary_map_low_res = hp.ud_grade(binary_mask_high_res, target_nside)

    # 5. Ensure it remains binary and inclusive.
    binary_map_final = (binary_map_low_res > 0).astype(int)

    return binary_map_final


def forcedType(strategyType):
    """
    Maps a human-readable strategy name to the specific scheduler strategy to be executed.

    Args:
        strategyType (str): The strategy level (e.g., "gold", "silver").

    Returns:
        str or int: Returns "GW_case_B" for gold, "GW_case_D" for silver,
            or -1 if the strategy is unrecognized.
    """
    if strategyType == "gold":
        return "GW_case_B"
    elif strategyType == "silver":
        return "GW_case_D"
    elif strategyType == "custom1":
        return "GW_case_custom1"
    elif strategyType == "custom2":
        return "GW_case_custom2"
    else:
        return -1


def selectRandomAng(myMap):
    """
    Selects a random sky coordinate (RA/Dec) from non-zero pixels in a HEALPix map.

    Args:
        myMap (np.ndarray): A HEALPix map array.

    Returns:
        tuple: A tuple containing (latitude, longitude) in radians,
            specifically (theta - pi/2, 2*pi - phi).
    """
    theta, phi = hp.pix2ang(32, np.random.choice(np.argwhere(myMap > 0)[0]))
    return theta - np.pi / 2, 2 * np.pi - phi


def makeTargetoO_object_LIGO(event_id, targetTime, stratType, _id):
    """
    Creates a TargetoO object from a LIGO skymap.

    Args:
        filePath (str): Path to the LIGO skymap file.
        targetTime (float): The MJD time for the target event.
        stratType (str): The observation strategy type (e.g., "gold").
        _id (int/str): Unique identifier for the TargetoO object.

    Returns:
        TargetoO: An instantiated TargetoO object populated with footprint,
            coordinates, and alert metadata.
    """
    mjdTime = targetTime
    duration = getSimNightsFromSkymapMetadata()
    footprint = get_90_percent_binary_map(
        f"https://gracedb.ligo.org/api/superevents/{event_id}/files/bayestar.fits.gz"
    )
    alertType = forcedType(stratType)
    ra_rad_c, dec_rad_c = selectRandomAng(footprint)
    return TargetoO(
        _id,
        footprint,
        mjdTime,
        duration,
        ra_rad_center=ra_rad_c,
        dec_rad_center=dec_rad_c,
        too_type=alertType,
    )


def getDayObsFromDateTime(dateTime, dateJust=0):
    """
    Converts a NumPy datetime64 object into an integer format (YYYYMMDD).

    Args:
        dateTime (np.datetime64): The base datetime object.
        dateJust (int, optional): Number of days to offset the date. Defaults to 0.

    Returns:
        int: The date represented as an integer (e.g., 20350101).
    """
    dt = dateTime + np.timedelta64(dateJust, "D")
    formatted = np.datetime_as_string(dt, unit="D").replace("-", "")
    return int(formatted)


def split_regions(x):
    """
    Parses and flattens a list of comma-separated region strings into a sorted list.

    Args:
        x (iterable): A collection of strings containing region names
            (e.g., ["Region A, Region B", "Region C"]).

    Returns:
        list: A sorted list of unique, whitespace-stripped region names.
    """
    regions = set()
    for k in x:  # .target_name:
        regions = regions.union(set([kk.replace(" ", "") for kk in k.split(",")]))
    regions = list(regions)
    regions.sort()
    return regions


eventTriggerTime = 60991  # The MJD that S251112cm went off
day_obs = getDayObsFromDateTime(
    np.datetime64(Time(eventTriggerTime, format="mjd").to_datetime()), dateJust=-2
)
sim_nights = getSimNightsFromSkymapMetadata()

import getpass
import os
import warnings
import copy
import logging

logging.getLogger("lsst_survey_sim").setLevel(logging.INFO)
logging.getLogger("rubin_nights").setLevel(logging.INFO)

import numpy as np
import sqlite3
import healpy as hp
import matplotlib.pyplot as plt

import colorcet as cc
import skyproj

import datetime
from astropy.coordinates import SkyCoord
import astropy.units as u

from rubin_scheduler.site_models import Almanac
from rubin_scheduler.scheduler import sim_runner
from rubin_scheduler.scheduler.schedulers import CoreScheduler
from rubin_scheduler.scheduler.features import Conditions
from rubin_scheduler.scheduler.utils import TargetoO, SimTargetooServer
from rubin_scheduler.utils import (
    ddf_locations,
    angular_separation,
    approx_ra_dec2_alt_az,
    Site,
    SURVEY_START_MJD,
)

import rubin_sim.maf as maf
from rubin_sim.data import get_baseline

from rubin_nights import connections
import rubin_nights.lfa_data as rn_lfa
import rubin_nights.dayobs_utils as rn_dayobs
import rubin_nights.plot_utils as rn_plots
import rubin_nights.augment_visits as augment_visits
import rubin_nights.rubin_scheduler_addons as rn_sch
import rubin_nights.rubin_sim_addons as rn_sim
import rubin_nights.observatory_status as observatory_status
import rubin_nights.scriptqueue as scriptqueue
import rubin_nights.scriptqueue_formatting as scriptqueue_formatting

from rubin_scheduler.scheduler.surveys import ToOScriptedSurvey

import importlib

from lsst_survey_sim import lsst_support, simulate_lsst, plot

band_colors = rn_plots.PlotStyles.band_colors

# today_dayobs is the day of today - useful for checking if we're running in the past or not
today_dayobs = rn_dayobs.day_obs_str_to_int(rn_dayobs.today_day_obs())

day_obs_time = rn_dayobs.day_obs_to_time(day_obs)
next_day_obs_time = rn_dayobs.day_obs_to_time(day_obs) + TimeDelta(1, format="jd")

next_day_obs = rn_dayobs.day_obs_str_to_int(
    rn_dayobs.time_to_day_obs(next_day_obs_time)
)

# Some parameters relating to downtime setup for model observatory
# In general, if only simulating one night, probably want no downtime or clouds
day_downtime = day_obs
if sim_nights == 1:
    add_downtime = False
    real_downtime = False
    add_clouds = False
else:
    add_downtime = True
    real_downtime = True
    add_clouds = True

logger.info("Setting up to with downtimes like:")
logger.info("add_downtime =", add_downtime)
logger.info("real_downtime = ", real_downtime)
logger.info("add_clouds = ", add_clouds)

survey_start = SURVEY_START_MJD
programs = ["BLOCK-365", "BLOCK-407", "BLOCK-408"]
logger.info("Potential science programs : ", programs)

sunset, sunrise = rn_dayobs.day_obs_sunset_sunrise(day_obs, sun_alt=-12)
logger.info(f"Simulation for {sim_nights} nights starting on :")
logger.info(
    f"DayObs {day_obs}, -12 deg sunset {sunset.iso}, -12 deg sunrise {sunrise.iso}"
)

sunset18, sunrise18 = rn_dayobs.day_obs_sunset_sunrise(day_obs, sun_alt=-18)
logger.info(
    f"DayObs {day_obs}, -18 deg sunset {sunset18.iso}, -18 deg sunrise {sunrise18.iso}"
)

# Git checkout ts_config_scheduler and set the configs to use.
# .. have to sort something out about what happens if the repo isn't clean ..
do_git_stuff = False
if do_git_stuff:
    ts_commit = "develop"
    # ts_commit = current_commit
    simulate_lsst.get_configuration(ts_commit, clone_path=ts_config_scheduler_root)
config_script_path = os.path.join(
    ts_config_scheduler_root,
    "Scheduler/feature_scheduler/maintel/",
    "fbs_config_lsst_survey.py",
)
config_ddf_script_path = os.path.join(
    ts_config_scheduler_root, "Scheduler/ddf_gen/lsst_ddf_gen_block_407.py"
)

import pickle

try_toos = True
manual_too = True

if try_toos and (not manual_too):
    lookback = TimeDelta(8, format="jd")
    toos = simulate_lsst.fetch_too_events(sunset - lookback, sunrise)
    if toos is not None:
        with open("test_too.p", "wb") as f:
            pickle.dump(toos, f)
elif try_toos and manual_too:
    # Make the too object here
    toos = []
    toos.append(
        makeTargetoO_object_LIGO(superevent, eventTriggerTime, force_strategy, 0)
    )

else:
    try:
        with open("test_too.p", "rb") as f:
            toos = pickle.load(f)
    except FileNotFoundError:
        toos = None

toos = toos[-1:]  # use only one ToO
logger.info(f"ToO type: {toos[0].too_type}")
logger.info(f"Number of ToOs: {len(toos)}")

# Wrap up ToOs for sim target server - only last currently relevant, since we want to observe one object
too_server = SimTargetooServer(toos)

# Configure band scheduler
band_scheduler = simulate_lsst.setup_band_scheduler()

# Check what bands it expects for day_obs?
conditions = Conditions(mjd=sunset.mjd)

# Configure scheduler and add initial visits
starting_scheduler, _, nside = simulate_lsst.setup_scheduler(
    config_script_path=config_script_path,
    config_ddf_script_path=config_ddf_script_path,
    day_obs=day_obs,
    band_scheduler=band_scheduler,
    too_server=too_server,
    initial_opsim=None,
)

# # Add the new surveys

additional_too_surveys = []

# Use no safety masks - this is fine for survey performance comparison
masks = []

# ## custom1

times = np.array([0, 2, 4, 24, 48, 72], float)
bands_at_times = ["ugri", "ugri", "ugri", "gri", "gri", "gri"]
nvis = [4, 4, 4, 6, 6, 6]
exptimes = [30.0, 30.0, 30.0, 30.0, 30.0, 30.0]
additional_too_surveys.append(
    ToOScriptedSurvey(
        masks,
        nside=32,
        followup_footprint=np.ones(hp.nside2npix(32)),
        times=times,
        bands_at_times=bands_at_times,
        nvis=nvis,
        exptimes=exptimes,
        detailers=None,
        too_types_to_follow=["GW_case_custom1"],
        survey_name="ToO, GW_case_custom1",
        target_name_base="GW_case_custom1",
        observation_reason="gw_case_custom1",
        science_program="BLOCK-407",
        flushtime=48,
        n_snaps=1,
    )
)

# # custom2

times = np.array([0, 2, 4, 24, 48, 72], float)
bands_at_times = ["gr", "gr", "gri", "gri", "gri", "gri"]
nvis = [4, 4, 4, 6, 6, 6]
exptimes = [30.0, 30.0, 30.0, 30.0, 30.0, 30.0]
additional_too_surveys.append(
    ToOScriptedSurvey(
        masks,
        nside=32,
        followup_footprint=np.ones(hp.nside2npix(32)),
        times=times,
        bands_at_times=bands_at_times,
        nvis=nvis,
        exptimes=exptimes,
        detailers=None,
        too_types_to_follow=["GW_case_custom2"],
        survey_name="ToO, GW_case_custom2",
        target_name_base="GW_case_custom2",
        observation_reason="gw_case_custom2",
        science_program="BLOCK-407",
        flushtime=48,
        n_snaps=1,
    )
)

for s in additional_too_surveys:
    starting_scheduler.survey_lists[0].append(s)

# Could consider placing a single seeing value
swap_seeing = False
if swap_seeing:
    # The fill value for missing DIMM values is 1
    seeing = 1.0
else:
    seeing = None

# Instantiate the observatory
starting_observatory, survey_info = simulate_lsst.setup_observatory(
    day_obs=day_downtime,
    nside=nside,
    add_downtime=add_downtime,
    add_clouds=add_clouds,
    seeing=seeing,
    real_downtime=False,  # Leaving this as false to avoid downtime, and needing an opsim file
    too_server=too_server,
)

mask_days = 30
logger.info(f"Modeled availability over next {mask_days} nights")
mask = np.where(
    survey_info["dayobsmjd"] < rn_dayobs.day_obs_to_time(day_obs).mjd + mask_days
)
logger.info("Total nighttime {}, ".format(survey_info["hours_in_night"][mask].sum()))
logger.info("total downtime {}, ".format(survey_info["downtime_per_night"][mask].sum()))
logger.info(
    "available time {}".format(
        survey_info["hours_in_night"][mask].sum()
        - survey_info["downtime_per_night"][mask].sum()
    )
)

reset = True
if reset:
    logger.info("Resetting the observatory and scheduler instance")
    observatory = copy.deepcopy(starting_observatory)
    scheduler = copy.deepcopy(starting_scheduler)

# Run the simulation for some nights
logger.info("Running simulation")
observations, scheduler, observatory, rewards, obs_rewards, survey_info = (
    simulate_lsst.run_sim(
        scheduler=scheduler,
        band_scheduler=band_scheduler,
        observatory=observatory,
        survey_info=survey_info,
        day_obs=day_obs,
        sim_nights=sim_nights,
        keep_rewards=False,
    )
)

# Save the observations
obsDF = pd.DataFrame(observations)
savePath = os.path.join(save_dir, f"{superevent}_obs.h5")
logger.info(f"Save dir: {savePath}")
obsDF.to_hdf(savePath, key="visits")
logger.info(f"Observations saved to {savePath}.")
logger.info(f"Script finished for ToO object {superevent}. Exiting.")
