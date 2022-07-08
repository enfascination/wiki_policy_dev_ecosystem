BIGSSS project on how wiki language editions evolve policies and influence each other

1.

Run SPARQL query to get the list of policy pages in the Wikidata category
"Wikimedia policy pages and guidelines" (Q4656150):
https://www.wikidata.org/wiki/Q4656150

  Script: 01-sparql_query.txt
  Query: https://w.wiki/5RCg
  Output: generated_data/policy_page_sitelinks-sparql.csv (symlink)

2.

  Get the full list of policy pages for each of the policies listed in Wikidata
  collected from the query in #1.

  Script: 02-collect_sitelinks_for_all_policy_articles.py 
  Inputs: generated_data/policy_page_sitelinks-sparql.csv (#1)
  Outputs: edge_list.csv policy_page_links.csv
