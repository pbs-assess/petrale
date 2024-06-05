## Data used in this analysis was retrieved using a development branch of the gfdata package
## these functions pull all data, including sets that have been deemed unusable due to survey domain changes
# remotes::install_github("pbs-assess/gfdata", ref = "trials")
# library(gfdata)
#
# dd <- get_survey_sets2("Petrale Sole", ssid = NULL,
#                        remove_false_zeros = TRUE, usability = NULL)
#
# ds <- get_survey_samples2("Petrale Sole",
#                           include_event_info = TRUE,
#                           unsorted_only = FALSE)
#
# saveRDS(dd, "dataIn/all-survey-sets-petrale.rds")
# saveRDS(ds, "dataIn/all-survey-samples-petrale.rds")

library(tidyverse)

# these survey_abbrev will be used
surveys_included <- c("MSSM QCS", "MSSM WCVI",
                      "HS MSA", "SYN HS", "SYN QCS", "SYN WCHG", "SYN WCVI")


dset <- readRDS("dataIn/all-survey-sets-petrale.rds") %>%
  # this removes duplications in MSSM surveys
  filter(
    survey_abbrev %in% surveys_included,
    ## some MSSM sets are in both as QCS and WCVI
    !(survey_series_id == 6 & latitude < 50),
    !(survey_series_id == 7 & latitude > 50)
  ) %>%
  mutate(
    survey_abbrev = ifelse(survey_abbrev == "MSSM" & latitude < 50, "MSSM WCVI",
                           ifelse(survey_abbrev == "MSSM" & latitude > 50, "MSSM QCS", survey_abbrev)
    ),
    survey_type = as.factor(
      case_when(
        survey_abbrev == "HS MSA"~"MSA",
        survey_abbrev %in% c("MSSM WCVI", "MSSM QCS") & year>2002 & year<=2005~"MSSM<=05",
        survey_abbrev %in% c("MSSM WCVI", "MSSM QCS") & year>2005~"MSSM>05",
        survey_abbrev %in% c("MSSM WCVI", "MSSM QCS") & year <= 2002~"MSSM <03",
        survey_abbrev %in% c("SYN HS", "SYN QCS", "SYN WCHG", "SYN WCVI")~"SYN",
        TRUE~survey_abbrev
      ))
  ) %>% distinct()

dsamp <- readRDS("dataIn/all-survey-samples-petrale.rds") %>%
  filter(survey_abbrev %in% surveys_included) %>%
  mutate(
    survey_abbrev = ifelse(survey_abbrev == "MSSM" & latitude < 50, "MSSM WCVI",
                           ifelse(survey_abbrev == "MSSM" & latitude > 50, "MSSM QCS", survey_abbrev)
    ),
    survey_type = as.factor(
      case_when(
        survey_abbrev == "HS MSA"~"MSA",
        survey_abbrev %in% c("MSSM WCVI", "MSSM QCS") & year>2002 & year<=2005~"MSSM<=05",
        survey_abbrev %in% c("MSSM WCVI", "MSSM QCS") & year>2005~"MSSM>05",
        survey_abbrev %in% c("MSSM WCVI", "MSSM QCS") & year <= 2002~"MSSM <03",
        survey_abbrev %in% c("SYN HS", "SYN QCS", "SYN WCHG", "SYN WCVI")~"SYN",
        TRUE~survey_abbrev
      ))) %>%
  distinct()

# correct error in year due to false bottom contact time
dset[dset$fishing_event_id == 886005,]$year <- 2005

# save combined processed data
saveRDS(dset, "data-generated/all-sets-used.rds")
saveRDS(dsamp, "data-generated/all-samples-used.rds")

# data checks
unique(dsamp$survey_abbrev)
check_for_duplicates <- dsamp[duplicated(dsamp$specimen_id), ]

