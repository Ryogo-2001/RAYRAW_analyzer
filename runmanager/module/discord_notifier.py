#!/usr/bin/env python3

import json
import logging
import os
import urllib.request
import urllib.error

logger = logging.getLogger('__main__').getChild(__name__)

WEBHOOK_URL_ENV = 'DISCORD_WEBHOOK_URL'
WEBHOOK_URL_FILE = os.path.expanduser('~/.config/hyps/discord_webhook')

#______________________________________________________________________________
def _get_webhook_url():
    url = os.environ.get(WEBHOOK_URL_ENV)
    if url:
        return url.strip()
    if os.path.isfile(WEBHOOK_URL_FILE):
        with open(WEBHOOK_URL_FILE) as f:
            url = f.read().strip()
        if url:
            return url
    return None

#______________________________________________________________________________
def notify(message: str, username: str = 'HYPS RunManager'):
    url = _get_webhook_url()
    if not url:
        logger.debug('Discord webhook URL not set — skipping notification')
        return
    payload = json.dumps({'username': username, 'content': message}).encode('utf-8')
    req = urllib.request.Request(
        url,
        data=payload,
        headers={'Content-Type': 'application/json',
                 'User-Agent': 'Python/3'},
        method='POST')
    try:
        with urllib.request.urlopen(req, timeout=10) as resp:
            if resp.status not in (200, 204):
                logger.warning(f'Discord webhook returned status {resp.status}')
    except urllib.error.URLError as e:
        logger.warning(f'Discord webhook failed: {e}')
