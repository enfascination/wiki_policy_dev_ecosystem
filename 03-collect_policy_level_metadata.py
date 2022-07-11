#!/usr/bin/env python3

import requests
import pandas as pd
import numpy as np
import datetime
import re
import urllib.parse

## this is the only line you'll need to change to get this to work in 
sitelink_df=pd.read_csv('generated_data/policy_page_links.csv')
# sitelink_df=pd.read_csv('/content/drive/MyDrive/Wikipedia Policy Development BIGSSS/policy_page_links.csv')

# get the page name
sitelink_df['title'] = sitelink_df.url.apply(lambda x: urllib.parse.unquote(re.sub('^.*/wiki/', '', x)))

sitelink_df['base_api'] = 'https://' + sitelink_df.lang + '.wikipedia.org/w/api.php?titles=' + sitelink_df.title + '&format=json'

timestamp_list=[]
for base_api in sitelink_df['base_api']:
    api_link = base_api + '&action=query&prop=revisions&rvlimit=1&rvprop=timestamp&rvdir=newer'
    r = requests.get(api_link)
    data = r.json()
    if '-1' in data['query']['pages']: ## page is missing
        timestamp = timestamp_list.append('NA')
    else:
        timestamp=list(data['query']['pages'].values())[0]['revisions'][0]['timestamp']
    timestamp_list.append(timestamp)

