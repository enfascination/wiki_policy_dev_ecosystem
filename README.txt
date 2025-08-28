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

3. 

  Enrich the metadata

  Script: 03-collect_policy_level_metadata.py
  Inputs: generated_data/policy_page_links.csv
  Outputs: generated_data/policy_page_links-post03.csv

4.

  Test the hypothesis that there are more than one policy
  development path, representing the problem as a sequence clustering. Script
  is is a stochastic numerical optimizer performing many fits for an optimal
  allocation of sequences to k  groups. This can be run for k={1,2,3,4}, possibly higher, without real edits (perhaps turning some blocks of code on and off. 

  Script: `04-cluster_fit.r`
  Inputs: some version of `generated_data/policy_page_links-post03.csv`
  Outputs: `.Rdata` files of arrays of the population of fits rather,
  ac[k]_results.Rdata`

5.

  Generate analytics on sequence results, for robustness and presentation. to
  perform analyses to understand the characteristics of the population of
  clusters. important for descriptives and robustness checks. Latter for now produces a few simple descriptive plots.

  Scripts: `05-analyze_clusters.r`, `06-figures.r`
  Inputs: The output of `04-cluster_fit.r`, `.Rdata` files 
  Outputs: textual outputs to shell and some images


