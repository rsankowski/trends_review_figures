library(rhdf5)
library(tidyverse)
library(spatstat)

#load cartana data
#chen et al cell 2020
lst <- list()

chen_folder <- file.path("data", "Cartana", "cartana_chen_et_al") 
chen_files <- list.files(chen_folder)

chen <- map(chen_files, function(x) {
  a <- read_tsv(file.path(chen_folder, x)) 
  lst[[x]] <- a
})

names(chen) <- gsub("Cartana_genecalls_|\\..*", "",chen_files)

chen_sum <- chen %>%
  bind_rows(.id="table") %>% 
    mutate(species=case_when(
      grepl("HS", table) ~ "Human",
      TRUE ~ "Mouse"
    )) %>%
  group_by(species, table) %>%
  summarise(xrange=max(X)-min(X),
            yrange=max(Y)-min(Y),
            area=xrange*yrange*.32^2/10e6,
            probe_count=n())

nndist <- data.frame()

nndist <- map(names(chen), function(x) {
  a <- distinct(chen[[1]][,3:4])
  bind_rows(nndist, data.frame(table=x, nndist=min(nndist(a$X, a$Y)*0.32)))
}) %>%
  bind_rows()

#load microglia data
micr <- read.csv2(file.path("data", "cell_counts_final_nn.csv")) %>%
  filter(Antigen=="Iba1" & Region %in% c("Cortex", "Marklager")) %>% 
  distinct(Image, .keep_all=T) %>% 
  na.omit() %>% 
  mutate(Mean_cell_dist=as.numeric(Mean_cell_dist))

#plots
df <- data.frame(Protocol=c("ST", "Visium", "Slide-seq", "HDST","Cartana"),
                 Dotsize=c(0.1, 0.055, 0.01, 0.002, 0.00032))

df <- df %>% mutate(Dots_per_1mmsq=1/Dotsize^2,
                    Protocol=factor(Protocol, levels = c("ST", "Visium", "Slide-seq", "HDST","Cartana")),
                    Resolution=round(sqrt(Dotsize^2+Dotsize^2)*1000,2))

df %>% ggplot(aes(Protocol,Dots_per_1mmsq, group=1, fill=Protocol, label=round(Dots_per_1mmsq,2))) +
  geom_point(size=20, pch=21) +
  expand_limits(y=0) +
  scale_y_log10() +
  theme_linedraw() +
  theme(text=element_text(size=25),
        panel.grid.minor = element_blank()) +
  labs(title="Maximal number of spatially distinct data points", y="Max. dots per 1mm^2") +
  coord_flip() +
  geom_text(nudge_x = .2, size=5) +
  scale_fill_brewer(guide=F, palette = "Set3")

ggsave(file.path("plots", "spatial_protocols_number_of_dots_per_mmsq.pdf"), useDingbats=F)

#expected resolution
df %>%
  ggplot(aes(Protocol, Resolution, fill=Protocol, label=Resolution)) +
  geom_point(size=20, pch=21) +
  expand_limits(y=0) +
  scale_y_log10() +
  theme_linedraw() +
  theme(text=element_text(size=25),
        panel.grid.minor = element_blank()) +
  labs(title="Spatial resolution", y="Resolution (µm)") +
  coord_flip() +
  geom_text(nudge_x = .2, size=5) +
  scale_fill_brewer(guide=F, palette = "Set3")
ggsave(file.path("plots", "spatial_protocols_expected_resolution.pdf"), useDingbats=F)

