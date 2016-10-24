# percentile-bootstrap
Percentile bootstrap function for paired-samples tests

This repository contains a function for performing a percentile bootstrap based on trimmed means, to be used on paired-samples data.
The percentile bootstrap generates B bootstrap samples from the condition 1/condition 2 difference scores to derive an X% confidence interval of the difference.
For single sample tests it creates a bootstrap sample using the single sample data.

A function is also provide to evaluate the Type 1 error rate of the test using different sample sizes.
These Type 1 error rates can be compared to those from Student's paired-samples t test and Yuen's t.
