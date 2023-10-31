# Energetic trade-offs in migration decision-making, reproductive effort, and subsequent parental care in a long-distance migratory bird. *Proceedings of the Royal Society B*.
## Schindler, A. R., A. D. Fox, C. K. Wikle, B. M. Ballard, A. J. Walsh, S. B. A. Kelly, L. Cao, L. R. Griffin, M. D. Weegman. 

#### Model code and associated data files to quanitfy the effects of reproductive preparation during spring migration on subsequent breeding outcomes, breeding outcomes on autumn migration characteristics, and autumn migration characteristics on subsequent parental survival in Greenland white-fronted geese (*Anser albifrons flavirostris*). 
___
### Authors
Alexander R. Schindler  
Department of Biology, University of Saskatchewan, Saskatoon, SK, Canada

Anthony D. Fox  
Department of Ecoscience, Aarhus University, Aarhus, Denmark

Christopher K. Wikle  
Department of Statistics, University of Missouri, Columbia, MO, USA

Bart M. Ballard  
Caesar Kleberg Wildlife Research Institute, Texas A&M University-Kingsville, Kingsville, TX, USA

Alyn J. Walsh  
National Parks and Wildlife Service, Dublin, Ireland

Seán B. A. Kelly  
National Parks and Wildlife Service, Dublin, Ireland

Lei Cao  
State Key Laboratory of Urban and Regional Ecology, Research Center for Eco-Environmental Sciences, Chinese Academy of Sciences, Beijing, China
University of Chinese Academy of Sciences, Beijing, China

Larry R. Griffin  
Wildfowl & Wetlands Trust, Gloucester, UK
ECO-LG Limited, Dumfries, UK

Mitch D. Weegman  
Department of Biology, University of Saskatchewan, Saskatoon, SK, Canada
___
### Files
- `full_annual_cycle_model_categorical.R` - Includes all code to read in data files and run the full annual cycle model with breeding outcome as a categorical variable (i.e., successful attempt, failed attempt, deferral).  
- `full_annual_cycle_model_bernoulli.R` - Includes all code to read in data files and run the full annual cycle model with breeding outcome as a Bernoulli variable (i.e., success, fail).  
- `spring_dat.csv` - All data for the spring sub-model (summarised to each bird/year/sub-season).  
- `autumn_dat.csv` - All data for the autumn sub-model  (summarised to each bird/year/sub-season).  

### Data file column names
- `id` - numerical identifier for each bird marked in the study
- `year` - numerical identifier for each year in the study (i.e., 1-5 corresponds with 2018-2022)
- `sub_season` - numerical identifier for the sub-season in the corresponding season
- `breeding_outcome` - 1 = successful breeding attempt, 2 = failed breeding attempt, 3 = breeding deferral
- `breeding_success` - 1 = breeding success, 0 = breeding failure
- `first_day` - day of year (i.e., days since 31 Dec the previous year) for the first day of the corresponding bird/year/sub-season
- `log_ODBA` - log-transformed overall dynamic body acceleration
- `num_feed_fixes` - number of ACC fixes classified as feeding
- `num_ACC_fixes` - total number of ACC fixes
- `mean_precip` - mean daily cumulative precipitation
- `prop_days_below_freezing` - proportion of days below freezing
- `prop_storm_days` - proportion of days with a severe storm
- `prop_grass` - proportion of time spent in grasslands
- `prop_ag` - proportion of time spent in (non-grassland) agricultural land
- `prop_bog` - proportion of time spent in peat bogs
- `survival` (autumn data only) - autumn survival

