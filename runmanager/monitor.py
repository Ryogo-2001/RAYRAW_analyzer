#!/usr/bin/env python3
#______________________________________________________________________________
import argparse
import fcntl
import json
import logging
import logging.config
import os
import sys
import time
import yaml
from datetime import datetime

top_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(top_dir, 'module'))

import singlerun

try:
  from rich.console import Console, Group
  from rich.table import Table
  from rich.live import Live
  from rich.panel import Panel
  from rich.text import Text
  from rich import box
  _USE_RICH = True
except ImportError:
  _USE_RICH = False
  import utility
  from utility import pycolor as cl

DISPLAY_PERIOD = 5  # second
logger = logging.getLogger(__name__)
parsed = None

# ── status → (rich style, symbol) ─────────────────────────────────────────────
_STATUS_MAP = [
  ('done',       'bold green',   '✓'),
  ('running',    'bold cyan',    '▶'),
  ('merging',    'bold yellow',  '⟳'),
  ('failed',     'bold red',     '✗'),
  ('terminated', 'red',          '■'),
  ('staged',     'white',        '○'),
  ('purged',     'dim',          '◌'),
  ('initialized','dim',          '·'),
  ('unknown',    'bold magenta', '?'),
]

def _status_style(status):
  for key, style, sym in _STATUS_MAP:
    if key in status:
      return style, sym
  return 'white', '·'

#______________________________________________________________________________
# rich version
#______________________________________________________________________________
def _build_renderable(info, skip_done):
  counts = dict(total=0, done=0, running=0, failed=0, other=0)

  table = Table(
    box=box.SIMPLE_HEAD,
    border_style='bright_black',
    header_style='bold white',
    show_edge=True,
    expand=False,
    padding=(0, 1),
  )
  table.add_column('KEY',          style='bold',       width=10,  no_wrap=True)
  table.add_column('STATUS',                           width=26,  no_wrap=True)
  table.add_column('PROG',         justify='center',   width=7,   no_wrap=True)
  table.add_column('BIN',          style='dim',        width=18,  no_wrap=True)
  table.add_column('CONF',         style='dim',        width=18,  no_wrap=True)
  table.add_column('DATA (#EVT)',                      width=26,  no_wrap=True)
  table.add_column('ROOT',         style='dim',        width=18,  no_wrap=True)
  table.add_column('TIME',         justify='right',    width=10,  no_wrap=True)

  for key, item in sorted(info.items(),
                          key=lambda x: (int(x[0]) if str(x[0]).isdigit()
                                         else x[0])):
    status = singlerun.SingleRun.decode_status(item)
    counts['total'] += 1
    if 'done'       in status: counts['done']    += 1
    elif 'running'  in status or 'merging' in status: counts['running'] += 1
    elif 'failed'   in status or 'terminated' in status: counts['failed'] += 1
    else:                      counts['other']   += 1

    if skip_done and 'done' in status:
      continue

    style, sym = _status_style(status)
    ptime = singlerun.SingleRun.decode_time(item)

    infile = None
    if 'data' in item and item['data'] is not None:
      infile = item['data']
    elif 'dstin' in item and item['dstin'] and len(item['dstin']) == 1:
      infile = item['dstin'][0]
    nev  = item.get('nev', '?')
    inbuf = f'{os.path.basename(str(infile))} ({nev})'

    nseg = item.get('nseg', 0)
    prog = item.get('prog', 0)
    prog_str = f'{prog}/{nseg}' if nseg else '-'

    status_text = Text()
    status_text.append(f'{sym} ', style=style)
    status_text.append(status,    style=style)

    table.add_row(
      key[:10],
      status_text,
      prog_str,
      os.path.basename(item['bin'])[-18:],
      os.path.basename(item['conf'])[-18:],
      inbuf[-26:],
      os.path.basename(item['root'])[-18:],
      ptime,
    )

  # ── header panel ────────────────────────────────────────────────────────────
  now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
  summary = Text()
  summary.append('  TOTAL ',  style='bold white')
  summary.append(f'{counts["total"]}   ', style='bold white')
  summary.append('DONE ',    style='bold green')
  summary.append(f'{counts["done"]}   ', style='bold green')
  summary.append('RUNNING ', style='bold cyan')
  summary.append(f'{counts["running"]}   ', style='bold cyan')
  summary.append('FAILED ',  style='bold red')
  summary.append(f'{counts["failed"]}   ', style='bold red')
  summary.append(f'  {now}',  style='dim')

  header = Panel(summary,
                 title='[bold bright_blue] Job Monitor [/]',
                 border_style='bright_blue',
                 padding=(0, 1))

  # ── footer ───────────────────────────────────────────────────────────────────
  footer = Text(f'  Ctrl-C to exit  ·  refresh every {DISPLAY_PERIOD}s',
                style='dim', justify='center')

  return Group(header, table, footer), counts

#______________________________________________________________________________
def _main_rich(path):
  console = Console()
  with Live(console=console, refresh_per_second=4, screen=True) as live:
    ptime = time.time()
    try:
      while True:
        try:
          with open(path, 'r') as f:
            try:    fcntl.flock(f.fileno(), fcntl.LOCK_SH)
            except IOError as e:
              logger.error(e); time.sleep(1); continue
            buff = f.read()
            fcntl.flock(f.fileno(), fcntl.LOCK_UN)
          info = json.loads(buff)
        except FileNotFoundError as e:
          logger.error(e); break
        except ValueError as e:
          logger.error(e); time.sleep(1); continue

        renderable, _ = _build_renderable(info, parsed.skip)
        live.update(renderable)

        dtime = DISPLAY_PERIOD - (time.time() - ptime)
        if dtime > 0:
          time.sleep(dtime)
        ptime = time.time()
    except KeyboardInterrupt:
      pass

#______________________________________________________________________________
# fallback: original plain-text version (no rich)
#______________________________________________________________________________
def _display_plain(filename):
  with open(filename, 'r') as f:
    try:    fcntl.flock(f.fileno(), fcntl.LOCK_SH)
    except IOError as e:
      logger.error(e); return False
    buff = f.read()
    fcntl.flock(f.fileno(), fcntl.LOCK_UN)
  try:
    info = json.loads(buff)
  except ValueError as e:
    logger.error(e); return False

  os.system('clear')
  header = (cl.reverce + cl.bold
            + 'KEY'.ljust(8)          + '  '
            + 'STAT'.ljust(18)        + '  '
            + 'BIN'.ljust(16)         + '  '
            + 'CONF'.ljust(16)        + '  '
            + 'DATA(#EVENT)'.ljust(24)+ '  '
            + 'ROOT'.ljust(16)        + '  '
            + 'TIME'.ljust(8)
            + cl.end)
  print(header)

  for key, item in sorted(info.items(),
                          key=lambda x: (int(x[0]) if str(x[0]).isdigit()
                                          else x[0])):
    status = singlerun.SingleRun.decode_status(item)
    if parsed.skip and 'done' in status:
      continue
    ptime = singlerun.SingleRun.decode_time(item)
    infile = None
    if 'data' in item and item['data'] is not None:
      infile = item['data']
    elif 'dstin' in item and len(item['dstin']) == 1:
      infile = item['dstin'][0]
    nev   = item.get('nev')
    inbuf = f'{os.path.basename(str(infile))} ({str(nev)})'
    row = (cl.bold + key[:8].ljust(8) + cl.end + '  '
           + f'{cl.reverce}{cl.bold}{cl.red}{status}{cl.end}'.ljust(16+20)
           + '  '
           + os.path.basename(item['bin'])[-16:].ljust(16)  + '  '
           + os.path.basename(item['conf'])[-16:].ljust(16) + '  '
           + inbuf[-24:].ljust(24)                          + '  '
           + os.path.basename(item['root'])[-16:].ljust(16) + '  '
           + ptime.rjust(8))
    print(row)
  return True

def _main_plain(path):
  ptime = time.time()
  try:
    while _display_plain(path):
      dtime = DISPLAY_PERIOD - (time.time() - ptime)
      if dtime > 0:
        time.sleep(dtime)
      ptime = time.time()
  except KeyboardInterrupt as e:
    logger.info(e)
  except FileNotFoundError as e:
    logger.error(e)

#______________________________________________________________________________
def main(path):
  if _USE_RICH:
    _main_rich(path)
  else:
    _main_plain(path)

#______________________________________________________________________________
if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-s', '--skip', help='skip runs with status of "done"',
                      action='store_true')
  parser.add_argument('stat_json', help='stat json file path')
  parsed, unparsed = parser.parse_known_args()
  log_conf = os.path.join(top_dir, 'logging_config.yml')
  with open(log_conf, 'r') as f:
    logging.config.dictConfig(yaml.safe_load(f))
  main(parsed.stat_json)
