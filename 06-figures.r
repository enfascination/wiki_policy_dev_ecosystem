#install.packages("forcats")
library(ggplot2)
library(tidyverse)
library(readr)
library(forcats)

policy_adoption <- read_delim("generated_data/policy_page_links.csv", col_names=T, delim=',')
dim(policy_adoption)
#str(policy_adoption)
#spec(policy_adoption)

policy_hist <- ggplot(policy_adoption, aes(fct_infreq(policy_name_en))) + geom_bar()+ theme(axis.title.x = element_text(""), axis.text.x = element_text(angle = 80, hjust=1))  + xlab('')

print(policy_hist)
ggsave(policy_hist, filename="0_policy_freqs.png")


pa_gini <- policy_adoption %>% group_by(policy_name_en) %>% summarize(a = n()) %>% arrange(a)
pa_gini$f <- pa_gini$a / sum(pa_gini$a)
pa_gini$c <- cumsum(pa_gini$f)
pa_gini$n <- 1:nrow(pa_gini)
policy_gini <- ggplot(pa_gini) + geom_line(aes(x=n, y=c)) + geom_line(color="grey", aes(x=n, y=seq(0,1,length.out=nrow(pa_gini))))
print(policy_gini)
ggsave(policy_gini, filename="0_policy_gini.png", width=3)



