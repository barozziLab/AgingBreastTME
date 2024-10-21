#~~~~~~~~~~~~
# JI function
#############

jaccard_index <- function(v1, v2) {
  
  i <- intersect(v1, v2)
  u <- union(v1, v2)
  
  ji <- length(i) / length(u)
  
  return(ji)
  
}

#~~~~~~~~~~~~~~~~~~~~~
# LIANA custom dotplot
######################
# adapted from liana_dotplot function

liana_custom_dotplot <- function (liana_res, source_groups = NULL, target_groups = NULL, plot.title = "Source",
                                  specificity = "natmi.edge_specificity", magnitude = "sca.LRscore", 
                                  y.label = "Interactions (Ligand -> Receptor)", size.label = "Interaction\nSpecificity", 
                                  colour.label = "Expression\nMagnitude", size_range = c(1, 10))
{
  entities <- c("ligand.complex", "receptor.complex")
  liana_mod <- liana_res
  liana_mod %<>% rename(magnitude = !!magnitude) %>% rename(specificity = !!specificity) %>% 
    unite(entities, col = "interaction", sep = " -> ") %>% 
    unite(c("source", "target"), col = "source_target", remove = FALSE)
  interactions_order <- liana_mod %>% pull("interaction") %>% 
    unique()
  liana_mod %<>% mutate(interaction = factor(interaction, levels = rev(interactions_order))) %>% 
    mutate(across(where(is.character), as.factor))
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                 "#0072B2", "#D55E00", "#CC79A7", "#DF69A7")
  suppressWarnings(ggplot(liana_mod, aes(x = target, y = interaction, 
                                         colour = magnitude, size = specificity, group = target)) + 
                     geom_point() + 
                     scale_color_gradientn(colours = viridis::viridis(20), limits = c(0.8, 0.982)) + 
                     scale_size_continuous(range = size_range, limits = c(0.04, 1)) + 
                     facet_grid(. ~source, space = "free", scales = "free", switch = "y") + 
                     labs(y = y.label, colour = colour.label, size = size.label, 
                          x = "Target", title = plot.title) + 
                     theme_bw(base_size = 20) + 
                     facet_wrap(~age) +
                     theme(legend.text = element_text(size = 16), 
                           axis.text.x = element_text(colour = cbPalette[1:length(unique(liana_mod$source))], 
                                                      size = 20, angle = 90, hjust = 1, vjust = 0.5), 
                           axis.title.x = element_text(colour = "gray6"), 
                           axis.text.y = element_text(size = 18, vjust = 0.5), 
                           legend.title = element_text(size = 18), panel.spacing = unit(0.1, "lines"), 
                           strip.background = element_rect(fill = NA), 
                           plot.title = element_text(vjust = 0, hjust = 0.5, colour = "gray6"), 
                           strip.text = element_text(size = 24, colour = "gray6")))
}


#~~~~~~~~~~~~~~~~~~~~~~~
# LIANA dotplot function
########################

liana_overview_dotplot <- function(liana_res, age = "young", summary_stat = "Mean", color_palette, level = "subcompartment") {
  
  stats <- liana_res %>% 
    group_by(age, source, target) %>% 
    summarize(Mean = mean(interaction_score),
              Median = median(interaction_score),
              n_interactions = n())
  
  plt <- ggplot(stats %>% filter(age == !!age), aes(x = source, y = target)) +
    geom_point(aes(size = n_interactions, fill = !!sym(summary_stat)), shape = 21) +
    theme_classic() +
    xlab(paste0("Sender (", level, ")")) +
    ylab(paste0("Receiver (", level, ")")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_distiller(palette = color_palette, direction = 1) +
    scale_size(range = c(0.5, 5), name = "Number of interactions") +
    guides(fill = guide_colorbar(title = paste0(summary_stat, " interaction score"))) +
    ggtitle(paste0("Predicted interactions in ", age, " mice"))
  
  return(plt)
  
}