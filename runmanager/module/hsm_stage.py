#!/usr/bin/env python3
'''Shared helpers for staging files off HSM (tape) via ghils/head, and for
requesting bulk staging through the hstage daemon.'''

import logging
import os
import shlex
import subprocess

logger = logging.getLogger('__main__').getChild(__name__)

HSTAGE_PATH = '/ghi/fs02/hstage/requests'
HSTAGE_NMIN = 20
HSTAGE_NMAX = 10000


#______________________________________________________________________________
def ghils_staged(path):
    ''' Return True if path is already staged on disk (ghils reports G, B or BP). '''
    try:
        proc = subprocess.run(shlex.split(f'ghils {path}'),
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    except OSError as e:
        logger.warning(f'ghils failed for {path}: {e}')
        return False
    out = proc.stdout.decode().split()
    return bool(out) and out[0] in ('G', 'B', 'BP')


#______________________________________________________________________________
def prefetch_path(path):
    ''' Trigger HSM staging for a single path in the background
    ("head" recalls it). '''
    if path is None:
        return
    if not os.path.exists(path):
        logger.warning(f'Prefetch: file not found, skipping: {path}')
        return
    if ghils_staged(path):
        logger.debug(f'Prefetch: already staged: {path}')
        return
    logger.info(f'Prefetching (staging) {path}')
    try:
        subprocess.Popen(shlex.split(f'head {path}'),
                          stdout=subprocess.DEVNULL,
                          stderr=subprocess.DEVNULL)
    except OSError as e:
        logger.warning(f'Prefetch failed for {path}: {e}')


#______________________________________________________________________________
def register_bulk_stage(tag, paths):
    ''' Request staging of many files at once via the hstage daemon, by
    dropping a list file under HSTAGE_PATH (same convention as
    RunManager.register_staging). Raises ValueError if too many files. '''
    if len(paths) > HSTAGE_NMAX:
        raise ValueError(
            f'Too many files to stage at once ({len(paths)} > {HSTAGE_NMAX}).')
    stage_list_path = os.path.join(
        HSTAGE_PATH, tag + '.lst')
    with open(stage_list_path, 'w') as f:
        for p in paths:
            f.write(p + '\n')
    return stage_list_path
