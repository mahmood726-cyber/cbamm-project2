CBAMM: Comprehensive Bayesian and Advanced Meta-Analysis Methods
<!-- badges: start -->
Show Image
<!-- badges: end -->
CBAMM (Comprehensive Bayesian and Advanced Meta-Analysis Methods) is an R package that implements a unified framework for advanced meta-analysis methods, including transportability weighting, robust variance estimation, Bayesian analysis with model stacking, and novel conflict detection.
Key Features

🎯 Transportability Weighting: Entropy balancing and covariate adjustment for target populations
🔧 Robust Methods: HKSJ adjustment and robust variance estimation (RVE)
📊 Bayesian Analysis: brms and JAGS integration with model stacking
🔍 Sensitivity Analysis: PET-PEESE, publication bias tests, missing studies analysis
⚠️ Conflict Detection: K-means clustering to identify discordant studies
🤖 Smart Advisor: Automated methodological recommendations
📈 Rich Visualizations: Forest plots, funnel plots, diagnostic plots

Installation
Install the development version from GitHub:
r# Install devtools if you haven't already
install.packages("devtools")

# Install CBAMM
devtools::install_github("yourusername/cbamm")
Optional Dependencies
For full functionality, install these optional packages:
r# For transportability weighting
install.packages("WeightIt")

# For Bayesian analysis  
install.packages(c("brms", "rjags", "loo", "posterior"))

# For robust variance estimation
install.packages("clubSandwich")

# For enhanced publication bias testing
install.packages("weightr")

# For clustering diagnostics
install.packages("cluster")
Quick Start
rlibrary(cbamm)

# Simulate example data
data <- simulate_cbamm_data(n_rct = 15, n_obs = 10, n_mr = 5)

# Run comprehensive analysis
results <- cbamm(data)

# View results
print(results)
summary(results)
plot(results)
Advanced Usage
Custom Configuration
r# Configure specific methods
config <- cbamm_config(
  methods = list(
    transport = TRUE,
    bayesian = TRUE,
    conflict_detection = TRUE,
    pet_peese = TRUE
  ),
  estimators = c("REML", "PM"),
  bayesian = list(
    method = "brms",
    chains = 4,
    iter = 2000
  )
)

results <- cbamm(data, config = config)
Transportability to Target Population
r# Define target population
target_pop <- list(
  age_mean = 65,
  female_pct = 0.6,
  bmi_mean = 28,
  charlson = 1.5
)

# Run analysis with transportability weighting
results <- cbamm(data, target_population = target_pop)
Extract Specific Results
r# Meta-analysis results
meta_results <- results$meta_results

# Bayesian results
bayesian_results <- results$bayesian_results

# Advisor recommendations  
advisor <- results$advisor_recommendations

# Conflict detection
conflicts <- results$conflict_detection
Package Structure
The package is organized into modular components:

Core Analysis: robust_rma(), meta-analysis with HKSJ
Transportability: compute_transport_weights(), entropy balancing
Bayesian Methods: run_bayesian_analysis(), brms/JAGS integration
Sensitivity: pet_peese(), publication bias detection
Conflict Detection: detect_study_conflicts(), K-means clustering
Visualization: generate_cbamm_plots(), comprehensive plotting

Data Requirements
Your data should include these columns:
Required:

study_id: Unique study identifier
yi: Effect sizes (log hazard ratios, log odds ratios, etc.)
se: Standard errors

Recommended:

study_type: Study design ("RCT", "OBS", "MR")
grade: GRADE quality assessment ("High", "Moderate", "Low", "Very low")

For Transportability:

age_mean: Mean age
female_pct: Proportion female
bmi_mean: Mean BMI
charlson: Charlson comorbidity index

Methods Implemented
Meta-Analysis Methods

Random-effects meta-analysis (REML, DL, PM, ML estimators)
Hartung-Knapp-Sidik-Jonkman (HKSJ) adjustment
Robust variance estimation for dependent studies

Transportability

Entropy balancing with WeightIt
Robust entropy optimization fallback
GRADE quality score weighting

Bayesian Analysis

brms integration with model stacking
JAGS fallback implementation
LOO-based model comparison

Sensitivity Analysis

PET-PEESE for publication bias
Egger's test, Begg's test
Trim-and-fill, selection models
Missing studies sensitivity

Novel Features

K-means conflict detection with silhouette optimization
Automated advisor recommendations
Multiverse analysis framework

Citation
If you use CBAMM in your research, please cite:
[Your Citation Information]
Contributing
Contributions are welcome! Please see our contributing guidelines.
License
This project is licensed under the MIT License - see the LICENSE file for details.
Support

📫 Issues: GitHub Issues
📖 Documentation: Package Website
💬 Discussions: GitHub Discussions
