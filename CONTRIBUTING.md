# Contributing to MastingTrends

Thank you for considering contributing to MastingTrends! This document provides guidelines for contributing to the project.

## Ways to Contribute

1. **Report bugs**: Open an issue describing the bug and how to reproduce it
2. **Suggest enhancements**: Open an issue describing your idea
3. **Submit pull requests**: Fix bugs or add features
4. **Improve documentation**: Help make the documentation clearer
5. **Share your analyses**: Contribute example datasets or case studies

## Development Setup

1. Fork the repository
2. Clone your fork:
   ```bash
   git clone https://github.com/YOUR-USERNAME/MastingTrends.git
   cd MastingTrends
   ```

3. Install R and recommended packages:
   ```r
   install.packages(c("ggplot2", "readxl", "Kendall", "trend"))
   ```

4. Run the test script to ensure everything works:
   ```r
   source("scripts/test_functions.R")
   ```

## Code Style Guidelines

### R Code Style

- Follow the [tidyverse style guide](https://style.tidyverse.org/)
- Use meaningful variable and function names
- Add comments for complex logic
- Include roxygen2-style documentation for functions

Example function documentation:
```r
#' Brief description of function
#'
#' More detailed description
#'
#' @param param1 Description of parameter 1
#' @param param2 Description of parameter 2
#' @return Description of return value
#' @export
function_name <- function(param1, param2) {
  # Function body
}
```

### General Guidelines

- **Keep functions focused**: Each function should do one thing well
- **Handle errors gracefully**: Check inputs and provide helpful error messages
- **Maintain backward compatibility**: Don't break existing APIs unless necessary
- **Write tests**: Add tests for new functionality
- **Update documentation**: Keep README and other docs in sync with code changes

## Submitting Changes

1. Create a new branch for your changes:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. Make your changes
3. Test your changes thoroughly
4. Commit with clear, descriptive messages:
   ```bash
   git commit -m "Add feature: brief description"
   ```

5. Push to your fork:
   ```bash
   git push origin feature/your-feature-name
   ```

6. Open a pull request with:
   - Clear description of changes
   - Rationale for the changes
   - Any relevant issue numbers

## Testing

Before submitting a pull request:

1. Run the test script:
   ```r
   source("scripts/test_functions.R")
   ```

2. Test with example data:
   ```r
   source("scripts/example_analysis.R")
   ```

3. Verify that existing functionality still works
4. Test edge cases (missing data, single species, etc.)

## Adding New Features

When adding new features:

1. **New analysis methods**: Add to `scripts/trend_analysis.R`
2. **New visualizations**: Add to `scripts/visualization.R`
3. **New data formats**: Add to `scripts/data_utils.R`
4. **New workflows**: Create a new script in `scripts/`

Make sure to:
- Document the new feature
- Add examples to `example_analysis.R`
- Update the README if needed

## Reporting Bugs

When reporting bugs, include:

1. **Description**: Clear description of the bug
2. **Steps to reproduce**: Minimal code to reproduce the issue
3. **Expected behavior**: What you expected to happen
4. **Actual behavior**: What actually happened
5. **Environment**: R version, OS, package versions
6. **Data**: If possible, a minimal example dataset

Example bug report:
```
**Description**: calculate_cv() fails with error when all values are zero

**Steps to reproduce**:
data <- data.frame(Year = 2000:2010, Harvest = 0)
calculate_cv(data$Harvest)

**Expected**: Should return NA or handle gracefully
**Actual**: Error: division by zero

**Environment**: R 4.2.0, Ubuntu 22.04
```

## Feature Requests

When requesting features, include:

1. **Use case**: Why you need this feature
2. **Proposed solution**: How you think it should work
3. **Alternatives**: Other approaches you've considered
4. **Examples**: Example code or workflows

## Documentation

Good documentation is crucial. When contributing:

- Update README.md for major changes
- Update QUICKSTART.md for usage changes
- Add inline comments for complex code
- Include function documentation
- Provide examples

## Community Guidelines

- Be respectful and inclusive
- Provide constructive feedback
- Help others when you can
- Stay on topic
- Follow the code of conduct

## Questions?

If you have questions about contributing:

1. Check existing issues and documentation
2. Open a new issue with your question
3. Tag it with "question"

## License

By contributing, you agree that your contributions will be licensed under the same license as the project (see LICENSE file).

## Acknowledgments

Contributors will be acknowledged in:
- README.md (for significant contributions)
- Git commit history (all contributions)
- Release notes (for features in releases)

Thank you for contributing to MastingTrends!
