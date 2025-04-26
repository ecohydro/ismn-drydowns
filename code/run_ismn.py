#!usr/bin/env python
# -*- coding: utf-8 -*-
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

__author__ = 'Bryn Morgan'
__contact__ = 'bryn.morgan@geog.ucsb.edu'
__copyright__ = '(c) Bryn Morgan 2024'

__license__ = 'MIT'
__date__ = 'Tue 09 Apr 24 14:23:58'
__version__ = '1.0'
__status__ = 'initial release'
__url__ = ''

"""

Name:           run_ismn.py
Compatibility:  Python 3.7.0
Description:    Description of what program does

URL:            https://

Requires:       list of libraries required

Dev ToDo:       None

AUTHOR:         Bryn Morgan
ORGANIZATION:   University of California, Santa Barbara
Contact:        bryn.morgan@geog.ucsb.edu
Copyright:      (c) Bryn Morgan 2024


"""


# IMPORTS

import multiprocessing as mp
from configparser import ConfigParser

import drydowns
from drydowns import Agent, TowerAgent
from drydowns.config import Config
from drydowns.mylogger import getLogger
from drydowns.utils import is_true

import time

from soil.ismnagent import ISMNAgent

# from .ismnagent import ISMNAgent
# Create a logger
log = getLogger(__name__)

agent_dict = {
    'ISMN' : ISMNAgent,
    'FLUXNET' : TowerAgent,
    'SMAP' : Agent
}

def main():
    """Main execution script ot run the drydown analysis"""
    start = time.perf_counter()
    log.info("--- Initializing the model ---")

    # _______________________________________________________________________________________________
    # Read config
    # cfg = ConfigParser()
    # cfg.read("config.ini")
    config = Config('soil/config.ini').config
    cfg = config[config.get('RUN','type')]

    # cfg_model = cfg["MODEL"]

    # Initiate agent
    # if cfg['DATA']['data_type'] != 'SMAP':
    # TODO: Get object type from config
    # if cfg.name != 'SMAP':
    #     agent = TowerAgent(cfg=cfg)
    # else:
    #     agent = Agent(cfg=cfg)
    agent = agent_dict[cfg.name](cfg=cfg)
    agent.initialize()

    # _______________________________________________________________________________________________
    # Define serial/parallel mode
    # run_mode = cfg["MODEL"]["run_mode"]
    run_mode = cfg["run_mode"]
    log.info(f"--- Analysis started with {run_mode} mode ---")

    # _______________________________________________________________________________________________
    # Verbose models to run
    log.info(f"Running the following models:")
    # if is_true(cfg["MODEL"]["exponential_model"]):
    #     log.info(f"Exponential model")
    # if is_true(cfg["MODEL"]["q_model"]):
    #     log.info(f"q model")
    # if is_true(cfg["MODEL"]["sigmoid_model"]):
    #     log.info(f"Sigmoid model")
    for mod_name in ['exponential_model', 'q_model', 'sigmoid_model']:
        if cfg.getboolean(mod_name):
            log.info(f"{mod_name}")
    # [m for m in ['exponential_model', 'q_model', 'sigmoid_model'] if cfg_model.getboolean(m)]

    # Run the model
    if run_mode == "serial":
        results = agent.run(agent.data_ids[500])
    
    elif run_mode == "parallel":
        # nprocess = int(cfg["MULTIPROCESSING"]["nprocess"])
        # nprocess = cfg.getint("MULTIPROCESSING", "nprocess")
        nprocess = cfg.getint("nprocess")
        with mp.Pool(nprocess) as pool:
            results = list(pool.imap(agent.run, agent.data_ids))
        pool.close()
        pool.join()
    else:
        log.info(
            "run_mode in config is invalid: should be either 'serial' or 'parallel'"
        )

    # _______________________________________________________________________________________________
    # Finalize the model
    log.info(f"--- Finished analysis ---")

    if run_mode == "serial":
        if results.empty:
            log.info("No results are returned")
        else:
            try:
                agent.finalize(results)
            except NameError:
                log.info("No results are returned")

    elif run_mode == "parallel":
        if not results:
            log.info("No results are returned")
        else:
            try:
                agent.finalize(results)
            except NameError:
                log.info("No results are returned")

    end = time.perf_counter()
    log.info(f"Run took : {(end - start):.6f} seconds")


if __name__ == "__main__":
    main()
