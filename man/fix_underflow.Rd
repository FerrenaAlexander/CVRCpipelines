% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gsea.R
\name{fix_underflow}
\alias{fix_underflow}
\title{Fix underflow in a sorted list of -log pvalues}
\usage{
fix_underflow(scores, logFC_vec, verbose = F)
}
\arguments{
\item{scores}{named vector of -log or -log10(pvalues) times sign(LFC), sorted hi to low. For the +inf and -inf values, one suggestion is to sort them by logFoldChange. names are gene sybols or IDs, values are -logP values * sign of LFC}

\item{logFC_vec}{named vector of some effect size, such as log2fc, by which the underflow P value DEGs may be sorted, or else the sorting may be arbitrary. names are gene symbols or IDs, values are L2FC (or some other effect size such as L2FC * pct.diff)}
}
\value{
a vector of P value scores, where underflow -logP values have been ranked appropriated
}
\description{
Given a list of P values from Seurat Wilcox or EdgeR test, one should multiply P values times -log10 times sign of lfc. With this input, sometimes underflow occurs. This function will fix this underflow (ie for input to GSEA) by adding +1 to all Inf values and -1 to all -Inf values.
}
