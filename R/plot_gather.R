
# fraction of sample identified bar char ----------------------------------

plot_gather_classified <- function(gather_df) {
  # generate a dataframe that contains information about the percent of the sample that was unclassified.
  # then, bind this new dataframe to the sourmash_taxonomy dataframe within the plot
  # that way, this information is not carried on beyond this plot.
  gather_df <- read_gather(Sys.glob("tests/testthat/*gather*.csv"))
  sourmash_taxonomy_unclassified <- sourmash_taxonomy %>%
    group_by(query_name) %>%
    summarize(sum_f_unique_weighted = sum(f_unique_weighted)) %>%
    mutate(f_unique_weighted = 1 - sum_f_unique_weighted,
           database = "unclassified") %>%
    select(-sum_f_unique_weighted)

  ggplot(sourmash_taxonomy %>%
           bind_rows(sourmash_taxonomy_unclassified) %>%
           mutate(database = factor(database, levels = c("archaea", "bacteria", "fungi", "protozoa", "viral", "unclassified"))),
         aes(x = query_name, y = f_unique_weighted, fill = database)) +
    geom_col() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_brewer(palette = "Set2")
}

# upset plots -------------------------------------------------------------

from_gather_to_upset_df <- function(gather_df){}

plot_gather_upset <- function(upset_df){

}
