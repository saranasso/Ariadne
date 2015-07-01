# Ariadne
=========

This is the Matlab code related to the research project published in the paper:

"Ariadne’s thread: a robust software solution leading to automated absolute and relative quantification of SRM data"

http://pubs.acs.org/doi/abs/10.1021/pr500996s

Summary
<br/>
Selected reaction monitoring (SRM) MS is a highly selective and sensitive technique to quantify protein abundances in complex biological samples. To enhance the pace of SRM large studies, a validated, robust method to fully automate absolute quantification, and to substitute for interactive evaluation, would be valuable. To address this demand, we present Ariadne, a Matlab® software.

To quantify monitored targets, Ariadne exploits metadata imported from the transition lists, and targets can be filtered according to mProphet output. Signal processing and statistical learning approaches are combined to compute peptide quantifications. To robustly estimate absolute abundances, the external calibration curve method is applied, ensuring linearity over the measured dynamic range.

Ariadne was benchmarked against mProphet and Skyline by comparing its quantification performance on three different dilution series, featuring either noisy/smooth traces without background or smooth traces with complex background. Results, evaluated as efficiency, linearity, accuracy, and precision of quantification, showed that Ariadne’s performance is independent of data smoothness and complex background presence, and that Ariadne outperforms mProphet on the noisier dataset and improved twofold Skyline’s accuracy and precision for the lowest abundant dilution with complex background. Remarkably, Ariadne could statistically distinguish from each other all different abundances, discriminating dilutions as low as 0.1 and 0.2 fmol. These results suggest that Ariadne offers reliable and automated analysis of large-scale SRM differential expression studies.

