library(covr)

# Calculate coverage
cov <- package_coverage()

# Print coverage report
report(cov)

# Create coverage report in HTML
report(cov, file = "coverage-report.html")

# Get coverage percentage
percent_coverage <- percent_coverage(cov)

cat(sprintf("\nTotal coverage: %.1f%%\n", percent_coverage))
