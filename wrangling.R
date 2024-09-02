library(baad.data)
library(geodata)
library(sp)
library(tidyverse)

FR_occurrences = read_csv("https://raw.githubusercontent.com/ashleylang/BAM/v.01/FR_occurrences.csv")
FR_measurements = read_csv("https://raw.githubusercontent.com/ashleylang/BAM/v.01/FR_measurements.csv") %>%
  pivot_wider(names_from = measurementType, 
              values_from = measurementValue) %>% 
  rename(Myc_type = 'Mycorrhiza type')


#summarise FungalRoot database  
fungal_root <- left_join(FR_occurrences, FR_measurements, by = "CoreID") %>%
  dplyr::select(order, family, scientificName, Myc_type) %>%
  separate(scientificName, c("genus", "species"), extra = "drop", fill = "right") %>%
  group_by(order, family, genus, species, Myc_type) %>%
  summarise(n = n())

#make Myc type groupings: 
unique(fungal_root$Myc_type)

ECMS = c("EcM, AM undetermined", "EcM, no AM colonization")
ERCS = c("ErM, AM", "ErM, EcM", "ErM")

fungal_root$myc_group <- ifelse(fungal_root$Myc_type == "AM", "AM", ifelse(fungal_root$Myc_type == "EcM,AM", "ECM/AM", ifelse(fungal_root$Myc_type %in% ECMS, "ECM", ifelse(fungal_root$Myc_type %in% ERCS, "ERC", "Other"))))

#Choosing the most common mycorrhizal type associated with each species of plant
fungal_root_sp <- fungal_root %>%
  group_by(order, family, genus, species, myc_group) %>%
  summarise(number= sum(n)) %>%
  slice_max(order_by = number, n = 1) %>% 
  arrange(genus, species) %>% 
  filter(myc_group != "ECM/AM" & myc_group != "Other" & myc_group != "ERC") %>% 
  group_by(order, family, genus, species) %>% 
  mutate(dupe = n()>1) %>% 
  filter(dupe==F)

## Read in BAAD database, remove groups with poor coverage: myc types other than AM and ECM, and decidious gymnosperms
baad <- baad.data::baad_data()
dict <- as.data.frame(baad$dictionary)
baad_df <- as.data.frame(baad$data) %>% 
  dplyr::select(studyName, latitude, longitude, species, vegetation, map, mat, pft, a.lf, h.t, d.bh, m.lf, m.st, m.so, m.rt, m.to, 	ma.ilf) %>%
  separate(species, c("genus", "species"), extra = "drop", fill = "right") %>%
  left_join(fungal_root_sp, by = c("genus", "species")) %>% 
  mutate(LmTm = m.lf/ m.to,
         RmTm = m.rt/m.to,
         SmTm= m.st/m.to,
         pft = as.factor(pft)) %>% 
  filter( myc_group != "ERC" & myc_group != "Other" & pft != "DG")

###Extract MAT/MAP from WorldClim----
r = worldclim_global(var = "bio", res = 10, path = "wc", version = "2.1")
r <- r[[c(1,12)]]

coords <- data.frame(x=baad_df$longitude, y=baad_df$latitude) %>%
  na.omit() %>%
  distinct()
points <- vect(coords,
               geom=c("x", "y"),
               crs = "EPSG:4326")
values <- terra::extract(r,points)
df <- cbind.data.frame(coords,values) %>%
  rename(latitude = y, longitude = x,
         Temp = wc2.1_10m_bio_1, Prec = wc2.1_10m_bio_12) %>% 
  select(-ID)


###filtering----
#Now the full version of the data with baad_df, myc types, and climate:
full_df = baad_df %>% 
  left_join(df, by=c("latitude", "longitude")) %>%
  filter(RmTm>0, h.t >.5, myc_group != "ECM/AM", family!= "Ericaceae") %>% 
  unite(study_species, studyName, genus, species, sep="_", remove=F) %>% 
  mutate(log_ht = log(h.t),
         leaf_habit= case_when(pft=="EA" | pft== "EG" ~ "evergreen",
                               pft=="DA" ~ "deciduous"),
         evo_group= case_when(pft=="EA" ~ "angiosperm",
                              pft== "EG" ~ "gymnosperm",
                              pft=="DA" ~ "angiosperm") ) 

###climate and geography of dataset----

AM_ECM=c("#C49F50", "#91BBA8")

sub <- full_df %>%
  group_by(Temp, Prec, leaf_habit, myc_group, longitude, latitude) %>%
  summarise(n = n())

clim_space <- ggplot(sub, aes(x= Temp, y=Prec))+
  geom_point(aes(colour=myc_group,shape=leaf_habit, size=n), alpha=0.65)+
  scale_colour_manual(labels = c("AM", "ECM"),values=AM_ECM)+
  scale_shape_manual(labels = c("Deciduous", "Evergreen"), values=c(16,17))+
  labs(x=expression("Mean Annual Temperature ("*degree*C*")"), y="Mean Annual\nPrecipitation (mm)")+
  theme(legend.title=element_blank(),legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, "mm"), legend.spacing.y = unit(0, "mm"), 
        legend.text=element_text(size=8), legend.position = c(.04, .85), 
        axis.title=element_text(size=11), axis.text=element_text(size=10))+
  guides(color = guide_legend(override.aes = list(size=3), order=1), shape = guide_legend(override.aes = list(size=4), order=2))+
  scale_size(range = c(3,8), guide="none")
clim_space

#map
world <- map_data("world")

map=ggplot(data=world)+
  geom_polygon(aes(x=long, y=lat, group=group), fill="white", color="black", size=0.2)+
  coord_fixed(1.3)+
  theme(axis.text=element_blank(), axis.line = element_blank(), 
        axis.ticks=element_blank(), axis.title=element_blank(), 
        panel.background = element_rect(fill = "white"), 
        panel.border = element_blank(), legend.position="none", 
        plot.margin = unit(c(.01,.01,.01,.01), "lines"))+
  geom_point(aes(x = longitude, y = latitude,color=myc_group), data = sub, size = 1)+
  scale_colour_manual(values=AM_ECM)


####model prep-----
#make clean dataset for models 
full_df_mod <- full_df %>%
  dplyr::select(RmTm, LmTm, SmTm, m.to, log_ht, leaf_habit, myc_group, Temp, Prec, study_species, family, order) %>%
  drop_na() %>% 
  separate(study_species, into=c("Study", "Genus", "Species"), sep="_", remove=F) %>% 
  unite(SppName, c(Genus, Species), sep="_") %>% 
  mutate(.before = m.to,
         tot = rowSums(.[,1:3]),
         roots = RmTm/tot,
         shoots = SmTm/tot,
         leaves = LmTm/tot)
#1429 observations

full_df_mod %>% group_by(myc_group) %>% summarise(fam_num = n_distinct(SppName))
#species by myc type : 34 AM, 21 ECM

unique(full_df_mod$order)
#16 orders

full_df_mod %>% group_by(myc_group) %>% summarise(order_num = n_distinct(order))
#orders by myc type : 15 AM, 4 ECM

full_df_mod %>% group_by(leaf_habit) %>% summarise(order_num = n_distinct(order))
#orders by leaf habit: 18 decidous, 16 evergreen

full_df_mod %>% group_by(leaf_habit, myc_group) %>% summarise(n = n())
#group numbers

#check variable distributions
hist(full_df_mod$RmTm)
hist(full_df_mod[full_df_mod$myc_group=="AM",]$RmTm)
hist(full_df_mod$LmTm)
hist(full_df_mod$SmTm)
hist(full_df_mod$m.to) #needs transformation
full_df_mod$log_totbio = log(full_df_mod$m.to)
hist(full_df_mod$log_totbio)
