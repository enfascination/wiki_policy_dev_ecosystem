#!/usr/bin/env python3

import requests
import pandas as pd
import numpy as np
import datetime
import re
import urllib.parse
from datetime import datetime 

## this is the only line you'll need to change to get this to work in 
sitelink_df=pd.read_csv('generated_data/policy_page_links.csv')
# sitelink_df=pd.read_csv('/content/drive/MyDrive/Wikipedia Policy Development BIGSSS/policy_page_links.csv')

# get the page name
sitelink_df['title'] = sitelink_df.url.apply(lambda x: urllib.parse.unquote(re.sub('^.*/wiki/', '', x)))

timestamp_list = []
for row in sitelink_df[['lang', 'title']].iterrows():
    lang = row[1]['lang']
    title = row[1]['title']

    params = {'titles' : title,
              'format' : 'json',
              'action' : 'query',
              'prop' : 'revisions',
              'rvlimit' : 1,
              'rvprop' : 'timestamp',
              'rvdir' : 'newer' }

    base_url = 'https://' + lang + '.wikipedia.org/w/api.php'
   
    # gather the data
    r = requests.get(base_url, params=params)
    data = r.json()

    print(lang, title)
    if len(data['query']['pages']) != 1:
        print("ERROR: more than one page id returned")
        break

    # iterate through each pageid
    for pageid in data['query']['pages']:
        page_info = data['query']['pages'][pageid]
        if 'missing' in page_info:
            timestamp = 'NA'
        else:
            timestamp = page_info['revisions'][0]['timestamp']

        timestamp_list.append(timestamp)

## append the timestamp
sitelink_df['timestamp'] = timestamp_list

# drop the NA which is a missing page
sitelink_df = sitelink_df[sitelink_df['timestamp'] != "NA"]

## create order within wiki
df_new = pd.DataFrame()
for lang, df_subset in sitelink_df.groupby('lang'):
    df_subset = df_subset.sort_values(by="timestamp")
    df_subset["order_within_wiki"] = range(0, df_subset.shape[0])
    df_new = df_new.append(df_subset)
sitelink_df = df_new

## create order across wiki
df_new = pd.DataFrame()
for qid, df_subset in sitelink_df.groupby('QID'):
    df_subset = df_subset.sort_values(by="timestamp")
    df_subset["order_across_wikis"] = range(0, df_subset.shape[0])
    df_new = df_new.append(df_subset)
sitelink_df = df_new

# write out the post03 data processing
sitelink_df.to_csv('generated_data/policy_page_links-post03.csv', index=False)

# old code to split groups intos groups of 50 titles by wiki
# title_max = 1
# for lang, df_subset in sitelink_df.groupby('lang'):
#     title_batches = [df_subset['title'][i:i+title_max] for i in range(0, df_subset.shape[0], title_max)]
#     for batch in title_batches:
#         title_slug = "|".join([x for x in batch])
