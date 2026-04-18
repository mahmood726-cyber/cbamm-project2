# CBAMM: Comprehensive Bayesian and Advanced Meta-Analysis Methods

[![R-CMD-check](https://github.com/cbamm-dev/cbamm/workflows/R-CMD-check/badge.svg)](https://github.com/cbamm-dev/cbamm/actions)

## Overview

CBAMM provides sophisticated meta-analysis capabilities including network meta-analysis, Bayesian methods, and advanced sensitivity analyses within a unified R framework.

## Key Features

- **Comprehensive Meta-Analysis**: Advanced meta-analytical methods with intelligent method selection
- **Network Meta-Analysis**: Multi-treatment comparisons with ranking and inconsistency assessment
- **Bayesian Integration**: Hierarchical modeling with model stacking capabilities
- **Publication Bias Assessment**: PET-PEESE correction and selection model approaches
- **Transport Analysis**: External validity assessment through transportability weighting
- **Sensitivity Analysis**: Comprehensive robustness evaluation across methodological choices
- **Utility Functions**: Effect size conversion, formatting, and data simulation tools

## Installation

```r
# From CRAN (when available)
install.packages("cbamm")

# Development version
devtools::install_github("cbamm-dev/cbamm")
```

## Quick Start

### Basic Meta-Analysis

```r
library(cbamm)

# Generate example data
data <- simulate_cbamm_data()

# Perform comprehensive meta-analysis
result <- cbamm(data)
summary(result)
plot(result)
```

### Network Meta-Analysis

```r
# Load network data
data(network_example_data)

# Perform network meta-analysis
nma <- network_meta_analysis(
  data = network_example_data,
  studies = "study_id", 
  treatments = "treatment",
  effect_size = "effect", 
  std_error = "se"
)

# View results
print(nma$results$treatment_effects)
plot(nma, type = "forest")
```

### Utility Functions

```r
# Effect size conversion
convert_effect_size(c(0.2, 0.5, 0.8))

# Statistical formatting
format_p(c(0.001, 0.023, 0.156))
format_ci(0.25, 0.10, 0.40)

# Heterogeneity assessment
calculate_i_squared(effect_sizes, standard_errors)
```

## Documentation

- **Main Vignette**: `vignette("cbamm-introduction")` - Comprehensive guide with examples
- **Function Help**: `help(package = "cbamm")` - Individual function documentation
- **GitHub Pages**: [Full Documentation](https://cbamm-dev.github.io/cbamm/) (when available)

## Data Requirements

### Standard Meta-Analysis
- `study_id`: Study identifiers
- `yi`: Effect sizes  
- `se`: Standard errors

### Network Meta-Analysis
- User-specified column names for studies, treatments, effects, and standard errors
- Connected network structure required
- Minimum three treatments

## Advanced Features

### Configuration System
```r
# Customize analysis approaches
config <- cbamm_config()
config <- update_config(config, bayesian = TRUE, transport_weighting = TRUE)
result <- cbamm(data, config = config)
```

### Publication Bias Assessment
```r
# PET-PEESE analysis
bias_result <- pet_peese(data)

# Comprehensive bias sensitivity
bias_sensitivity <- run_publication_bias_sensitivity(data)
```

### Transport Analysis
```r
# Assess external validity
transport_weights <- compute_transport_weights(data, target_population)
transport_result <- cbamm(data, transport_weights = transport_weights)
```

## Performance

CBAMM provides multiple performance modes:
- **Fast**: `cbamm_fast()` for large datasets
- **Standard**: `cbamm()` for comprehensive analysis
- **Custom**: Configure specific methods through `cbamm_config()`

## Citation

When using CBAMM in research, please cite:

```
CBAMM Development Team (2025). CBAMM: Comprehensive Bayesian and Advanced 
Meta-Analysis Methods. R package version 0.1.0.
```

## Getting Help

- **Issues**: [GitHub Issues](https://github.com/cbamm-dev/cbamm/issues)
- **Discussions**: [GitHub Discussions](https://github.com/cbamm-dev/cbamm/discussions)
- **Email**: cbamm.dev@example.com

## License

GPL-3 License. See [LICENSE](LICENSE) for details.

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

