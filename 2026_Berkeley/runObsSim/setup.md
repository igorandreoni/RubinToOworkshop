# How to run these sims

Contact author: [@seanmacb](https://github.com/seanmacb) on github, or [my email](mailto:sean.macbride@physik.uzh.ch)

1. Create your environment:
```
conda create --name <env> --file requirements.txt
```
2. With your environment active, run `scheduler_download_data`.
3. Download the necessary data for the `rubin_scheduler`. [This link](https://rubin-scheduler.lsst.io/data-download.html#data-download) has more information about what is needed to download. The above command,`scheduler_download_data`, should cover mostly everything, but you may need extra `sky_brightness` files, which can be found [here](https://s3df.slac.stanford.edu/data/rubin/sim-data/sims_skybrightness_pre/h5_2023_09_12/). You only will need the `sky_brightness` files for the date range (in MJD) that you care about.
4. If you haven't already using the requirements file, install `ts_fbs_utils`:
```
pip install --user git+https://github.com/lsst-ts/ts_fbs_utils
```
5. Create a [usdf token](https://usdf-rsp.slac.stanford.edu/settings/tokens), and store it somewhere on your machine (remember the path for later).
6. Clone `ts_config_scheduler`:
```
git clone git@github.com:lsst-ts/ts_config_scheduler.git
```
7. In `runSim_S251112cm.py`, update the path to your token file, the path to `ts_config_scheduler`, and the path to the location of rubin sim data. All of these items are marked in the code by `# TODO`.
8. Run the script using the following syntax:
```
runSim_S251112cm.py [-h] [--save_dir SAVE_DIR] [--force_strategy FORCE_STRATEGY]

options:
  --save_dir SAVE_DIR   The path to save the resulting opsims.
  --force_strategy FORCE_STRATEGY
                        The type of strategy to enforce on the resulting simulation. Options are 'gold', 'silver', 'custom1', and 'custom2'
```
9. As you desire, change the desired strategy using the `--force_strategy` flag, and lines 395:443.