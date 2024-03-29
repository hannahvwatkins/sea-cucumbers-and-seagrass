#This script is for the data processing for Watkins et al., 2023 in Ecosphere
#Here we load the various csvs containing raw data for this project and 
#compile them into cleaned datasets - one for the 3 month long experiment 
#and its associated data (including the blade heights, grunt counts, sand 
#depths, and elemental analysis), another for this experiment but with the blade 
#heights condensed into a single mean and sd for each quadrat (so each 
#observation in the df corresponds to a single quadrat at a single time point), 
#and the third is for the 10 day growth rate experiment
#note that these datasets also include data from a separate experiment - these
#data are removed prior to the analyses conducted for this paper
#Please forgive the convoluted and inefficient ways I transformed these data, 
#much of this code comes from the early stages of my MSc and it turns out we
#picked extremely inconvenient ways to organize our inputted data!

#load in necessary packages
library(tidyverse) #for general data processing
library(janitor) #for cleaning column names quickly

#PART 1: TREATMENTS, PATCHES, AND CUKE DENSITIES--------------------------------
#locations of patches
locations <- read_csv("Data/cuke_patch_full_list.csv") %>% 
  #fix a couple naming discrepancies
  mutate(name = case_when(name == "76" ~ "76S",
                          name == "96S" ~ "96",
                          TRUE ~ name)) %>% 
  #select the variables of interest
  dplyr::select(patch = name, lat, lon)

#measurements of patch sizes
patches <- read_csv("Data/patch_sizes_raw.csv") %>% 
  #select variables of interest
  dplyr::select(patch, seagrass_perimeter, 
                patch_perimeter1:patch_perim_remaining) %>% 
  #separate the patch perim remaining into the correct number of cols since we
  #listed anything more than 4 individual reef sections in a single column
  separate(patch_perim_remaining, sep = ", ", 
           into = c("patch_perimeter5",
                    "patch_perimeter6",
                    "patch_perimeter7",	
                    "patch_perimeter8",
                    "patch_perimeter9",	
                    "patch_perimeter10",
                    "patch_perimeter11",	
                    "patch_perimeter12",
                    "patch_perimeter13",	
                    "patch_perimeter14",
                    "patch_perimeter15",	
                    "patch_perimeter16",
                    "patch_perimeter17",
                    "patch_perimeter18",	
                    "patch_perimeter19")) %>% 
  #convert values to numeric instead of characters
  mutate(across(c(patch_perimeter5:patch_perimeter19), as.numeric),
         #convert all perimeters to radii
         across(c(seagrass_perimeter:patch_perimeter19), 
                ~ .x/(2*pi), 
                .names = "radius_{.col}"),
         #then areas
         across(c(radius_seagrass_perimeter:radius_patch_perimeter19),
                ~ (.x ^2) * pi,
                .names = "area_{.col}")) %>% 
         #then add up all of the areas for each row
  rowwise() %>% 
  mutate(total_reef_area = 
           #couldn't figure out how to streamline this part so writing out all
           #the names because it's faster - womp womp
           sum(c(area_radius_patch_perimeter1,
                 area_radius_patch_perimeter2,
                 area_radius_patch_perimeter3,
                 area_radius_patch_perimeter4,
                 area_radius_patch_perimeter5,
                 area_radius_patch_perimeter6,
                 area_radius_patch_perimeter7,
                 area_radius_patch_perimeter8,
                 area_radius_patch_perimeter9,
                 area_radius_patch_perimeter10,
                 area_radius_patch_perimeter11,
                 area_radius_patch_perimeter12,
                 area_radius_patch_perimeter13,
                 area_radius_patch_perimeter14,
                 area_radius_patch_perimeter15,
                 area_radius_patch_perimeter16,
                 area_radius_patch_perimeter17,
                 area_radius_patch_perimeter18,
                 area_radius_patch_perimeter19), 
               na.rm = TRUE),
         true_seagrass_area = 
           area_radius_seagrass_perimeter - total_reef_area,
         #just rename the total area encompassed by the halo including the reef
         #to something more intuitive 
         patch_size = area_radius_seagrass_perimeter) %>% 
  #then just select the parts we care about
  dplyr::select(patch, total_reef_area, true_seagrass_area, patch_size)

#treatment assignments and original + manipulated cuke densities
treatments <- read_csv("Data/treatment_changes_raw.csv") %>% 
  #connect to the patch information and location info
  left_join(patches, by = "patch") %>% 
  left_join(locations, by = "patch") %>% 
  #convert spp-specific counts into total cukes per m2
  mutate(original_cuke_density = (ft_original + dd_original)/true_seagrass_area,
         final_cuke_density = (ft_final + dd_final)/true_seagrass_area,
         log_original_cuke_density = log(original_cuke_density),
         #deal with zeros - 50% of lowest non-zero value - evidence for this 
         #being the appropriate choice is in the analysis script (TLDR, a 
         #sensitivity analysis shows that the results are robust to the choice
         #of constant here)
         log_final_cuke_density = log(final_cuke_density + 0.5*0.0150220),
         total_final = dd_final + ft_final,
         #excretion rates based on Munger et al 2023 - this is in umol NH4+/h
         total_ft_excretion_initial = ft_original * 12.0,
         total_dd_excretion_initial = dd_original * 15.6,
         total_ft_excretion_final = ft_final * 12.0,
         total_dd_excretion_final = dd_final * 15.6,
         total_cuke_excretion_initial = total_ft_excretion_initial + 
           total_dd_excretion_initial,
         total_cuke_excretion_final = total_ft_excretion_final + 
           total_dd_excretion_final) %>% 
  #and then keep only the columns we need to streamline the final csvs
  dplyr::select(patch, treatment, original_cuke_density, final_cuke_density,
                true_seagrass_area, total_reef_area, log_original_cuke_density, 
                log_final_cuke_density, lat, lon, total_cuke_excretion_initial,
                total_cuke_excretion_final, patch_size)


#PART 2: MAIN SAMPLING------------------------------------------------------
#blade heights from main sampling
heights <- read_csv("Data/blade_height_raw.csv",
                    col_types = cols(blade_num = col_number())) %>% 
  #create new cols for transect number and quadrat number
  separate(quadrat, sep = c(1, 2, 3), 
           into = c(NA, "transect", NA, "quadrat")) %>% 
  #fill in missing blade numbers to create unique IDs for each row - the actual
  #numbers used don't matter, we just need a unique identifier for each row,
  #which in this case is the combo of patch, period, transect, quadrat, and 
  #blade_num
  group_by(patch, sampling_period, transect, quadrat) %>% 
  mutate(blade_num = row_number(),
         transect = as.numeric(transect),
         quadrat = as.numeric(quadrat)) %>% 
  ungroup()

#calculate a mean height for each quadrat for analyses at a quadrat level
grouped_heights <- heights %>% 
  group_by(patch, sampling_period, transect, quadrat) %>% 
  summarize(mean_blade_height = mean(height_cm, na.rm = TRUE),
            sd_blade_height = sd(height_cm, na.rm = TRUE),
            n_blade_height = n(),
            se_blade_height = sd_blade_height/sqrt(n_blade_height))

#grunt surveys
grunts <- read_csv("Data/grunt_surveys_raw.csv",
                   col_types = cols(percent_on_patch = col_character())) %>% 
  #we need to convert length to mass which we can do with the equation from
  #Fiona Francis' thesis using Rock Sound grunt data. The relationship is 
  #w = a*L^b, where a is 0.012 and b is 3.072. Weight will be in g and length 
  #in cm
  mutate(mass = 0.012 * (total_length_cm^3.072),
         log10_mass = log10(mass)) %>% 
  #and we can extract the weight-excretion relationship estimated by Allgeier 
  #et al 2015 (PNAS)- note that I put NH4 here, but it's actually ugrams of N/h.
  #to avoid confusion with the percent N stuff later, I convert to umol
  #of NH4 below (to make these values comparable to the cukes)
  mutate(log10_nh4 = 1.65 + (0.76 * log10_mass),
         #calculate value for all fish in that size class
         nh4 = (10^log10_nh4)*total_number) %>% 
  #then we can group by patch and add up the total amount of nitrogen (in 
  #ug/hour)
  group_by(patch) %>% 
  summarize(total_grunt_nh4_rate_ug = sum(nh4),
            #then convert ugrams of nitrogen to umol - since there is one 
            #molecule of N per molecule of NH4, this will be the same value as 
            #umol of NH4 to allow for easy comparison with the cukes
            #molar mass of N is 14.01 (so 14.01ug per umol)
            total_grunt_nh4_rate = total_grunt_nh4_rate_ug/14.01) %>% 
  ungroup()


#info on where each sample came from for the elemental analysis
elem_sample <- read_csv("Data/elemental_samples_raw.csv") %>% 
  #remove empty cols
  dplyr::select(!(X11:X28))
#nitrogen estimates for each sample
elem <- read_csv("Data/elemental_composition_raw.csv") %>% 
  dplyr::select(!(X4:X20)) %>% 
  #remove empty rows at top
  slice(-(1:3)) %>% 
  #convert top row to column names
  row_to_names(1) %>% 
  clean_names() %>% 
  #link elemental composition with the info about where sample came from
  left_join(elem_sample, by = c("sample" = "sample_id")) %>% 
  rename(sampling_period = time_point,
         patch = site) %>% 
  #remove "t"s from each sampling period
  separate(sampling_period, sep = 1, into = c(NA, "sampling_period")) %>% 
  #remove "q"s from each quadrat
  separate(quadrat, sep = 1, into = c(NA, "quadrat")) %>% 
  dplyr::select(patch, sampling_period, quadrat, 
                percent_total_n, percent_total_c) %>% 
  #put in correct format
  mutate(sampling_period = as.numeric(sampling_period),
         quadrat = as.numeric(quadrat),
         percent_total_n = as.numeric(percent_total_n)/100,
         percent_total_c = as.numeric(percent_total_c)/100)

#sand depths from RG's project
sand <- read_csv("Data/sand_depths_raw.csv") %>% 
  group_by(patch) %>% 
  #RG took six measures radiating from each patch, but only the first three
  #from each actually correspond to our quadrats, so we'll take those
  slice_min(distance_patch, n = 3) %>% 
  mutate(quadrat = row_number()) %>% 
  ungroup() %>% 
  rowwise() %>% 
  #take the mean value from the four corners of each quadrat
  mutate(mean_sand_depth = mean(c(depth1, depth2, depth3, depth4)),
         sd_sand_depth = sd(c(depth1, depth2, depth3, depth4))) %>% 
  ungroup() %>% 
  #and now we can remove the individual measures since they aren't particularly
  #helpful and will just make the df more crowded. We'll also remove distance
  #and counts since they are already in the sampling df from the same quadrats
  dplyr::select(!depth1:depth4 & !distance_patch & !shoot_count)
  

#full destructive sampling data (from main experiment)
#destructive sounds dramatic but it's just the data we had to actually disturb
#(i.e. cut/core) the seagrass for
sampling <- read_csv("Data/destructive_sampling_raw.csv") %>% 
  #lets remove the columns that aren't needed for analysis
  dplyr::select(!c(syringe_1, syringe_2, core_bag_num)) %>% 
  #now add in each of the other dfs
  left_join(elem, by = c("patch", "sampling_period", "quadrat")) %>% 
  #we'll use the one time point of grunts as an approximation of the consistent
  #fish supply (since Fi's work suggests that the abundances stay stable
  #over the summer)
  left_join(grunts, by = "patch") %>% 
  #even though the sand depths were only measured on T3, we'll use them as an
  #approximation for the depth at all three transects since there's a consistent
  #relationship between distance and depth
  left_join(sand, by = c("patch", "quadrat")) %>% 
  #and finally add in the info on cuke densities
  left_join(treatments, by = c("patch", "treatment")) %>% 
  filter(treatment!='ft' & treatment != 'dd') %>% 
  #make excretion per m2, and scale all continuous variables, as well as
  #deal with zeros for logging - 50% of lowest non-zero value (the test for this
  #being an appropriate choice is in the analysis script)
  mutate(grunt_pee_m2 = (total_grunt_nh4_rate/true_seagrass_area),
         cuke_pee_final_m2 = (total_cuke_excretion_final/true_seagrass_area),
         log_grunt_pee_m2 = log(grunt_pee_m2 + (0.5 * 0.4946114)),
         log_cuke_pee_m2 = log(cuke_pee_final_m2 + (0.5 * 0.1802642)),
         #scale all predictors
         scale_log_og_cukes = scale(log_original_cuke_density),
         scale_log_cukes_final = scale(log_final_cuke_density),
         scale_log_grunt_pee_m2 = scale(log_grunt_pee_m2),
         scale_log_cuke_pee_m2 = scale(log_cuke_pee_m2),
         scale_distance = scale(distance_patch),
         scale_sand = scale(mean_sand_depth),
         scale_seagrass_area = scale(true_seagrass_area),
         scale_reef_area = scale(total_reef_area),
         scale_patch_size = scale(patch_size),
         #we want to code treatment as continuous here, because there is a
         #direct linear relationship between the three levels since they are
         #proportions of original biomass removed
         cont_treatment = case_when(treatment == "natural" ~ 0,
                                    treatment == "half" ~ -0.5,
                                    treatment == "zero" ~ -1),
         scale_cont_treatment = scale(cont_treatment))

#one version where each row is a quadrat
full_sampling <- sampling %>% 
  left_join(grouped_heights, 
            by = c("patch", "sampling_period", "transect", "quadrat"))
#write_csv(full_sampling, "Data/main_sampling_clean.csv")

#another version where each row is an individual blade from the blade height
#data
full_sampling_all_heights <- sampling %>% 
  #need the right_join here so it doesn't get rid of the NAs
  right_join(heights, by = c("patch", "sampling_period", "transect", "quadrat"))
#write_csv(full_sampling_all_heights, 
#          "Data/main_sampling_clean_with_heights.csv")

#PART 3: GROWTH EXPERIENT------------------------------------------------------
#this is the least tidy format ever (damn you, past Hannah), so we'll pivot
#the three sections separately and join them
treatments2 <- treatments %>% 
  dplyr::select(patch, treatment, original_cuke_density, final_cuke_density, 
                true_seagrass_area, total_reef_area, log_original_cuke_density,
                log_final_cuke_density, lat, lon, total_cuke_excretion_final,
                patch_size)

#info on the shoots cut at each quadrat
growth_exp_2 <- read_csv("Data/seagrass_shaving_raw.csv") %>% 
  dplyr::select(patch, shoots_q1_cut:shoots_q3_cut) %>% 
  pivot_longer(shoots_q1_cut:shoots_q3_cut, names_to = "quadrat") %>% 
  rename(shoots_cut = value) %>% 
  mutate(quadrat = as.numeric(case_when(quadrat == "shoots_q1_cut" ~ "1",
                                        quadrat == "shoots_q2_cut" ~ "2",
                                        quadrat == "shoots_q3_cut" ~ "3")))

#info on the shoots found at each quadrat after 9 days
growth_exp_3 <- read_csv("Data/seagrass_shaving_raw.csv") %>% 
  dplyr::select(patch, date_final, shoots_q1_found:shoots_q3_found) %>% 
  pivot_longer(shoots_q1_found:shoots_q3_found, names_to = "quadrat") %>% 
  rename(shoots_found = value) %>% 
  mutate(quadrat = as.numeric(case_when(quadrat == "shoots_q1_found" ~ "1",
                                        quadrat == "shoots_q2_found" ~ "2",
                                        quadrat == "shoots_q3_found" ~ "3")))

#info on the quadrats
growth_exp <- read_csv("Data/seagrass_shaving_raw.csv") %>% 
  dplyr::select(patch:distance_patch3) %>% 
  pivot_longer(distance_patch1:distance_patch3, names_to = "quadrat") %>% 
  rename(distance = value) %>% 
  mutate(quadrat = 
           as.numeric(case_when(quadrat == "distance_patch1" ~ "1",
                                quadrat == "distance_patch2" ~ "2",
                                quadrat == "distance_patch3" ~ "3"))) %>% 
  left_join(growth_exp_2, by = c("patch", "quadrat")) %>% 
  left_join(growth_exp_3, by = c("patch", "quadrat"))

#and now we'll combine that with the individual blade measurements
growths <- read_csv("Data/seagrass_growth_raw.csv") %>% 
  left_join(growth_exp, by = c("patch", "quadrat")) %>% 
  left_join(treatments2, by = "patch") %>% 
  #and we can also add in the approximate sand depths along the transect
  left_join(sand, by = c("patch", "quadrat")) %>% 
  #and the grunt data
  left_join(grunts, by = "patch") %>% 
  filter(treatment!='ft' & treatment != 'dd') %>% 
  #deal with zeros for logging - 50% of lowest non-zero value
  mutate(grunt_pee_m2 = (total_grunt_nh4_rate/true_seagrass_area),
         cuke_pee_final_m2 = (total_cuke_excretion_final/true_seagrass_area),
         log_grunt_pee_m2 = log(grunt_pee_m2 + (0.5 * 0.4946114)),
         log_cuke_pee_m2 = log(cuke_pee_final_m2 + (0.5 * 0.1802642)),
         scale_log_cuke_density_final = scale(log_final_cuke_density),
         scale_log_og_cukes = scale(log_original_cuke_density),
         scale_log_grunt_pee_m2 = scale(log_grunt_pee_m2),
         scale_log_cuke_pee_m2 = scale(log_cuke_pee_m2),
         scale_distance = scale(distance),
         scale_sand = scale(mean_sand_depth),
         scale_seagrass_area = scale(true_seagrass_area),
         scale_reef_area = scale(total_reef_area),
         cont_treatment = case_when(treatment == "natural" ~ 0,
                                    treatment == "half" ~ -0.5,
                                    treatment == "zero" ~ -1),
         scale_cont_treatment = scale(cont_treatment),
         scale_patch_size = scale(patch_size))

#and save
#write_csv(growths, "Data/seagrass_growth_clean.csv")

#PART 4: REFORMATTING-----------------------------------------------------------
#to make the analysis easier, we can reformat the full sampling dataset
#we've done that in four ways here - creating a new column that shows the 
#initial values for each metric in each quadrat that we can use as an offset,
#the same data but with a single mean value for the three quadrats at each
#distance from the patch, and then two more dfs where we convert the response
#to a proportional change over time since that's a lot easier and more 
#intuitive to model

#Offset by quadrat
#we'll take the initial values at each site and treat them as an offset (which
#we can also use to calculate a proportional change)
#since our questions are about how seagrass changed after manipulation, this is
#the easiest way to look at that.
t0 <- full_sampling %>% 
  filter(sampling_period == 0) %>% 
  #create a unique ID for each tq
  unite(c(transect, quadrat), col = "quad_id", remove = FALSE) %>% 
  #only keep what we need
  transmute(patch = patch, 
            quad_id = quad_id, 
            shoot_count_offset = shoot_count,
            above_ground_offset = above_ground,
            below_ground_offset = below_ground,
            mean_height_offset = mean_blade_height,
            sd_height_offset = sd_blade_height,
            percent_total_n_offset = percent_total_n,
            observer_offset = observer)
offset_df <- full_sampling %>% 
  filter(sampling_period != 0) %>% 
  unite(c(transect, quadrat), col = "quad_id", remove = FALSE) %>%
  #this pairs each observation with it's pre-manipulation value
  left_join(t0, by = c("patch", "quad_id"))
#write_csv(offset_df, "Data/sampling_data_offset.csv")

#offset averaged across Q1, Q2, and Q3
offset_df_avg <- offset_df %>% 
  group_by(sampling_period, patch, quadrat) %>% 
  summarize(treatment = first(treatment),
            heading = first(heading),
            mean_distance = mean(distance_patch),
            mean_shoot_count = mean(shoot_count, na.rm = TRUE),
            mean_height = mean(mean_blade_height, na.rm = TRUE),
            mean_above = mean(above_ground, na.rm = TRUE),
            mean_below = mean(below_ground, na.rm = TRUE),
            mean_total_n = mean(percent_total_n, na.rm = TRUE),
            mean_shoot_count_offset = mean(shoot_count_offset, na.rm = TRUE),
            mean_height_offset = mean(mean_height_offset, na.rm = TRUE),
            mean_above_offset = mean(above_ground_offset, na.rm = TRUE),
            mean_below_offset = mean(below_ground_offset, na.rm = TRUE),
            mean_total_n_offset = mean(percent_total_n_offset, na.rm = TRUE),
            total_grunt_nh4_rate = first(total_grunt_nh4_rate),
            grunt_pee_m2 = first(grunt_pee_m2),
            cuke_pee_final_m2 = first(cuke_pee_final_m2),
            log_grunt_pee_m2 = first(log_grunt_pee_m2),
            log_cuke_pee_m2 = first(log_cuke_pee_m2),
            mean_sand_depth = first(mean_sand_depth),
            original_cuke_density = first(original_cuke_density),
            final_cuke_density = first(final_cuke_density),
            log_original_cuke_density = first(log_original_cuke_density),
            log_final_cuke_density = first(log_final_cuke_density),
            true_seagrass_area = first(true_seagrass_area),
            scale_log_og_cukes = first(scale_log_og_cukes),
            scale_log_cukes_final = first(scale_log_cukes_final),
            scale_log_grunt_pee_m2 = first(scale_log_grunt_pee_m2),
            scale_log_cuke_pee_m2 = first(scale_log_cuke_pee_m2),
            scale_distance = first(scale_distance),
            scale_sand = first(scale_sand),
            scale_seagrass_area = first(scale_seagrass_area),
            scale_reef_area = first(scale_reef_area),
            scale_patch_size = first(scale_patch_size),
            cont_treatment = first(cont_treatment),
            scale_cont_treatment = first(scale_cont_treatment),
            lat = first(lat),
            long = first(lon)) %>% 
  ungroup()
#write_csv(offset_df_avg, "Data/sampling_data_averaged_offset.csv")

#proportions by quadrat
proportions_df <- offset_df %>% 
  mutate(shoot_prop = (shoot_count)/shoot_count_offset,
         height_prop = (mean_blade_height)/mean_height_offset,
         above_prop = (above_ground)/above_ground_offset,
         below_prop = (below_ground)/below_ground_offset,
         percent_total_n_prop = (percent_total_n)/percent_total_n_offset)
#write_csv(proportions_df, "Data/sampling_data_proportions.csv")

#proportions averaged
proportions_df_avg <- offset_df_avg %>% 
  mutate(mean_shoot_prop = (mean_shoot_count)/mean_shoot_count_offset,
         mean_height_prop = (mean_height)/mean_height_offset,
         mean_above_prop = (mean_above)/mean_above_offset,
         mean_below_prop = (mean_below)/mean_below_offset,
         mean_percent_total_n_prop = (mean_total_n)/mean_total_n_offset)
#write_csv(proportions_df_avg, "Data/sampling_data_averaged_proportions.csv")
