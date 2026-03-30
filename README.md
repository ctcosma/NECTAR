# Purpose
This repository hosts an example implementation of the **Network-Enhanced Conservation Tool for Analysis and Recommendation (NECTAR)** framework, as descibed in "Local interaction networks reconstructed from global biodiversity data improve pollinator restoration decision making" by Baiotto and Cosma et al. (2026). 

NECTAR is a reproducible framework for predicting large-scale plant-pollinator interaction metawebs using existing interaction data and phylogenetic information, and downscaling these large-scale metawebs to regional metawebs using spatial co-occurrence and the overlap of key life-stage events (e.g., for flower visitation interactions, regional flowering phenology and pollinator adult flight period). While the code below is our application of NECTAR to California plant-pollinator communities, it can be extended to other regions and pollinator taxa with relative ease by re-configuring the extent of input occurrence, interaction, and environmental data, the spatial boundaries used to downscale to regional metawebs, and the "last say" authority for name harmonization in the region of interest (e.g., for plants in California, we use Jepson as the ultimate taxonomic authority).

# Organization
NECTAR is organized into 12 steps that are intended to be completed in order. The results presented in Baiotto and Cosma generally rely on the code in steps 01-07 and 12 and skip steps 08-11, which are various data summary/aggregation steps to support the implementation of NECTAR on the California Native Plant Societies' Calscape Pollinator Companion tool (beta version available at: https://pollinators.calscape.org/pollinator-companion/welcome).

📁 01_Data_Download - Automated data download (e.g., for plant occurrence and supporting environmental/ecological data from Calflora) and descriptions for manual data download and very basic cleaning/formatting.

📁 02_Occurrence_Cleaning - Harmonization of occurrence records from all sources and spatial, temporal, taxonomic, and column cleaning generally following established occurrence cleaning workflows (e.g., biodiversity data cleaning workflow - `bdc`)

📁 03_Checklist_Generation - Generation of taxon-specific and regionally-explicit checklists of plant and pollinator species using external sources/authorities (e.g., using Jepson to constrain to only California native plants) or occurrence data (when there is not a good external authority, e.g., hoverflies).

📁 04_SDMs - Ensemble species distribution models for all species in our checklists using a three-tiered approach based on data availability.

📁 05_Phenology - Estimation of region-level phenometrics for key life stage events using a three-tiered approach based on data availability.

📁 06_Interaction_Data - Downloading, scraping, and aggregation/cleaning of plant-pollinator interaction data from multiple sources.

📁 07_Interaction_Predictions - Generation of large-scale metawebs using simple phylogenetic constraints and graph embedding with transfer learning, and downscaling of large-scale metawebs to regional metawebs using spatial and phenological co-occurrence for potentially interacting pairs identified in the large-scale metaweb.

📁 08_Alpha_Diversity - Summarization of pollinator support for all native plants included, broken up by ecoregion (this is the foundation of the ranked draw selection of pollinator-supporting native plants)

📁 09_Beta_Diversity - Identification of complimentary pollinator diversity companion plants (plants that support the most additional pollinator species) for each plant in each ecoregion using the predicted intercation networks and constraining companion plants to only those found in the same ecosystem(s) as the focal plant.

📁 10_Phenological_Diversity - Identification of complimentary phenological companion plants (plants that add the most additional contiguous flowering time) for each plant in each ecoregion using the estimated phenometrics and constraining companion plants to only those found in the same ecosystem(s) as the focal plant.

📁 11_Crop_Support - Interaction predictions for crops grown in California to aid in the identification of native plants that would best serve as companion plants for each crop (in each region), based on the diversity of polliators that we expect to visit the crops that also would be likely to be supported by each native plant.

📁 12_Analyses - Analyses conducted specifically for Baiotto and Cosma et al. (2026), including comparison of pollinator support native plants identified with NECTAR to random selection of plants and external lists (e.g., Xerces Society regional plant lists and Pollinator Partnership ecoregional planting guides), ecosystem restoration scenarios, and network/data summarisation. 

📁 99_Supporting - Supporting script files contain functions re-used many times (or for methods that could be easily used in other applications - e.g., phenometric estimated), including the single objective genetic algorithm we use for identifying complimentary sets of native plants in the optimized planting scenarios, taxonomic name harmonization, and phenometrics.
