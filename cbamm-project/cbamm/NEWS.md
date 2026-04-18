# CBAMM 0.1.0

## Major Features

* **Comprehensive Meta-Analysis**: Core meta-analysis functionality with intelligent method selection
* **Network Meta-Analysis**: Multi-treatment comparisons with automatic reference selection and connectivity assessment
* **Bayesian Integration**: Hierarchical modeling with model stacking capabilities
* **Publication Bias Assessment**: PET-PEESE correction and comprehensive sensitivity analyses
* **Transport Analysis**: External validity assessment through transportability weighting
* **Utility Functions**: Comprehensive effect size conversion, statistical formatting, and data simulation

## Key Functions

### Core Analysis
* `cbamm()` - Comprehensive meta-analysis with automated method selection
* `cbamm_fast()` - Optimized analysis for large datasets
* `cbamm_config()` - Flexible configuration system
* `network_meta_analysis()` - Multi-treatment network meta-analysis

### Utility Functions
* `simulate_cbamm_data()` - Generate realistic meta-analysis datasets
* `convert_effect_size()` - Effect size metric conversion
* `format_p()`, `format_ci()` - Statistical result formatting
* `calculate_i_squared()` - Heterogeneity assessment
* `detect_outliers_iqr()` - Outlier identification

### Specialized Methods
* `pet_peese()` - Publication bias correction
* `robust_rma()` - Robust variance estimation
* `compute_transport_weights()` - External validity weighting
* `run_publication_bias_sensitivity()` - Comprehensive bias assessment

## Documentation

* Comprehensive vignette with practical examples
* Complete function documentation with usage examples
* Professional README with quick start guide
* Extensive utility functions for research workflows

## Performance

* Intelligent method selection based on data characteristics
* Multiple performance modes for different dataset sizes
* Robust error handling and input validation
* Integration with established R meta-analysis ecosystem

## Data Support

* Flexible data input formats
* Automatic data validation and structure checking
* Support for both traditional and network meta-analysis designs
* Comprehensive example datasets included

This initial release provides a solid foundation for advanced meta-analysis research with plans for continued methodological expansion in future versions.

