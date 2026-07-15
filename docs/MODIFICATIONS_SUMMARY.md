# SnPM Analysis Tool - Modifications Summary

## Overview
This document summarizes the modifications made to the SnPM Analysis Tool to add:
1. **Pearson Correlation** option (in addition to existing Spearman)
2. **Covariate Control** functionality for partial correlation analysis

## Modified Files

### 1. SnPMAnalysisGui.m
**Purpose**: Main GUI interface for the analysis tool

**Changes Made**:
- **Line 22-29**: Added covariate panel UI components
  - `CovariatePanel` - Panel for covariate selection
  - `CovariateFileButton` - Button to browse for covariate CSV file
  - `CovariateFileLabel` - Label showing selected covariate file
  - `UseCovariatesCheckBox` - Checkbox to enable/disable covariate use

- **Line 57**: Added `covariate_file` property to store covariate file path

- **Line 162-207**: Added `CovariateFileButtonPushed` callback function
  - Validates CSV file format (minimum 2 columns)
  - Displays preview of covariate columns
  - Auto-enables covariate checkbox when file is loaded
  - Provides user feedback on number of subjects and covariates

- **Line 188-189**: Modified `RunAnalysisButtonPushed` to pass covariate parameters
  - `params.covariate_file` - Path to covariate file
  - `params.use_covariates` - Boolean flag for covariate usage

- **Line 287-289**: Modified `ResetButtonPushed` to clear covariate settings

- **Line 395**: Updated GridLayout row heights to accommodate covariate panel

- **Line 488-515**: Added covariate panel UI creation in `createComponents`

- **Line 560**: **Added 'correlationP' to comparison dropdown** - KEY FEATURE
  - Users can now select Pearson correlation in addition to Spearman

### 2. core_snpm_analysis.m
**Purpose**: Core analysis engine

**Changes Made**:
- **Line 1-18**: Updated function documentation to include covariate parameters

- **Line 38-48**: Added covariate parameter extraction with defaults
  - `covariate_file` - defaults to empty string
  - `use_covariates` - defaults to false

- **Line 50-103**: **Added covariate loading and processing section**
  - Reads covariate CSV file
  - Validates format (Subject ID + covariates)
  - Converts subject IDs to standardized format
  - Handles errors gracefully with warnings

- **Line 140-169**: **Added covariate-subject matching logic**
  - Matches covariate subjects with data subjects
  - Excludes subjects without covariate data
  - Updates subject match info structure
  - Reports matching statistics

- **Line 265-272**: **Modified permutation analysis to pass covariates**
  - Passes covariate matrix to TFCE function for correlation analyses
  - Maintains backward compatibility (works without covariates)

- **Line 428-447**: **Added covariate information to results**
  - Reports number and names of covariates used
  - Shows subject matching statistics
  - Stores covariate data in results structure

### 3. dependencies/partial_correlation.m
**Purpose**: New helper function for partial correlation

**Features**:
- Implements residualization method for partial correlation
- Removes linear effects of covariates from both neural data and behavioral variables
- Handles NaN values appropriately
- Works channel-by-channel for efficiency

**Algorithm**:
```matlab
For each variable (behavioral and each channel):
  1. Create design matrix: X = [ones(n,1), covariates]
  2. Fit linear model: Y = X*beta + residuals
  3. Extract residuals: residuals = Y - X*beta
  4. Use residuals for correlation analysis
```

## New Test Data

### test_data/ Directory
Created comprehensive test dataset with:

1. **data1_condition1.csv** (15 subjects × 10 channels)
   - Simulated EEG power values for condition 1
   - Subject IDs: sub001-sub015

2. **data2_condition2.csv** (15 subjects × 10 channels)
   - Simulated EEG power values for condition 2
   - Matched subjects with condition 1
   - Systematically higher values to simulate effect

3. **covariates.csv** (15 subjects × 4 covariates)
   - Age (25-36 years)
   - Gender (0=Female, 1=Male)
   - Education (14-20 years)
   - Handedness (1=Right, -1=Left)

4. **README.md**
   - Comprehensive testing guide
   - 5 test scenarios
   - Troubleshooting tips
   - Data format specifications

## Usage Guide

### Basic Correlation Analysis (No Covariates)

1. Launch GUI: `SnPMAnalysisGui`
2. Select Data 1 file (e.g., `test_data/data1_condition1.csv`)
3. Select Data 2 file (e.g., `test_data/data2_condition2.csv`)
4. Choose output directory
5. Select comparison type:
   - **correlationS** - Spearman correlation
   - **correlationP** - Pearson correlation (NEW!)
6. Set parameters (channels, permutations, etc.)
7. Click "Run Analysis"

### Partial Correlation with Covariates (NEW!)

1. Follow steps 1-4 above
2. Click "Select Covariate CSV File..." button
3. Choose covariate file (e.g., `test_data/covariates.csv`)
4. Verify covariate preview dialog
5. Ensure "Include Covariates in Analysis" is checked
6. Select correlation type (correlationS or correlationP)
7. Click "Run Analysis"

**Result**: Analysis will compute partial correlations controlling for all covariates in the CSV file

## Covariate CSV File Format

```csv
Subject,Covariate1,Covariate2,Covariate3,...
sub001,value1,value2,value3,...
sub002,value1,value2,value3,...
...
```

**Requirements**:
- First column MUST be Subject IDs (matching data files)
- Remaining columns are covariates (numeric values)
- Column names will be displayed in results
- All subjects in data files should have covariate data

## Technical Implementation Details

### Partial Correlation Method

The implementation uses the **residualization approach**:

1. **For behavioral variable**:
   ```matlab
   beta = [1, covariates] \ behavioral_data
   residual_behavioral = behavioral_data - [1, covariates] * beta
   ```

2. **For each EEG channel**:
   ```matlab
   beta = [1, covariates] \ channel_data
   residual_channel = channel_data - [1, covariates] * beta
   ```

3. **Compute correlation on residuals**:
   ```matlab
   correlation = corr(residual_behavioral, residual_channel)
   ```

This is mathematically equivalent to partial correlation but more computationally efficient for permutation testing.

### Permutation Testing with Covariates

The permutation procedure:
1. Residualize data and behavioral variables
2. Permute subject labels on residualized data
3. Compute correlation on each permutation
4. Build null distribution
5. Calculate p-values

This maintains the covariate relationships while testing the partial correlation.

## Backward Compatibility

All modifications maintain **full backward compatibility**:
- Existing analyses without covariates work unchanged
- Covariate functionality is optional (checkbox must be enabled)
- Default behavior (no covariates) is identical to original version
- All existing comparison types (pairedT, unpairedT, etc.) unaffected

## Benefits of These Modifications

### 1. Pearson Correlation Option
- **Flexibility**: Choose appropriate correlation for data characteristics
- **Linear relationships**: Pearson more powerful for linear associations
- **Comparison**: Can compare Pearson vs Spearman results

### 2. Covariate Control
- **Confound control**: Remove effects of age, gender, education, etc.
- **Increased validity**: More accurate correlation estimates
- **Publication quality**: Meets standards for controlling confounds
- **Flexibility**: User-defined covariates via CSV import

### 3. User-Friendly Implementation
- **Simple CSV format**: Easy to prepare covariate files
- **Automatic validation**: File format checking with helpful error messages
- **Visual feedback**: Preview of loaded covariates
- **Clear reporting**: Covariate information in results

## Testing Recommendations

### Quick Test (Development)
- Use test_data files
- 1000 permutations
- Runtime: <1 minute

### Full Test (Validation)
- Use test_data files
- 10000 permutations
- Compare results with/without covariates
- Verify both Pearson and Spearman

### Production Use
- Real EEG data (164 or 178 channels)
- 10000+ permutations
- Appropriate covariates for your study
- Document covariate choices in methods

## Known Limitations

1. **Covariate file format**: Must be CSV with specific structure
2. **Subject matching**: All subjects must have covariate data (or will be excluded)
3. **Numeric covariates**: Categorical variables must be dummy-coded
4. **Linear relationships**: Assumes linear relationships between covariates and data
5. **Correlation only**: Covariate control currently only for correlation analyses

## Future Enhancements

Potential improvements for future versions:
1. Covariate control for group comparisons (t-tests)
2. Automatic dummy coding for categorical variables
3. Covariate selection (choose which covariates to include)
4. Interaction terms between covariates
5. Non-linear covariate relationships
6. Multiple covariate files for different subject groups

## References

### Theoretical Background
- **Maris & Oostenveld (2007)**: Nonparametric statistical testing of EEG- and MEG-data. *Journal of Neuroscience Methods*, 164, 177-190.
- **Nichols & Holmes (2001)**: Nonparametric permutation tests for functional neuroimaging. *Human Brain Mapping*, 15, 1-25.

### Partial Correlation
- **Cohen et al. (2003)**: Applied Multiple Regression/Correlation Analysis for the Behavioral Sciences. Chapter on partial correlation.

## Support

For questions or issues:
1. Check `test_data/README.md` for testing guidance
2. Review this document for implementation details
3. Examine example CSV files in `test_data/` folder
4. Verify covariate file format matches specifications

## Version History

**Version 2.0** (2025-01-14)
- Added Pearson correlation option
- Implemented covariate control functionality
- Created comprehensive test dataset
- Enhanced documentation

**Version 1.0** (Original)
- Basic SnPM analysis with Spearman correlation
- Group comparisons (paired/unpaired t-tests)
- TFCE and cluster-based corrections