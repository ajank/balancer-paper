# copied from /g/korbel/shared/projects/drosophila_balancer/analyses/figures/scripts

# colors.R
# Using mostly colors from here:
# http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=7

my_colors = c(del       = "deepskyblue2",
              dup       = "#e41a1c",        # 228,26,28
              wildtype  = "#4daf4a",
              balancer  = "#377eb8",
              breakpoint = "#984ea3",
              common    = "#999999",
              hetero    = "darkorange",
              error     = "firebrick1")
my_theme =  theme_minimal() + theme(strip.background = element_rect(), plot.title = element_text(hjust = 0.5,face = "bold"))