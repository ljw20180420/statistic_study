#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import (
    ttest_ind,
    mannwhitneyu,
    ttest_rel,
    wilcoxon,
    binomtest,
    bws_test,
    ranksums,
    brunnermunzel,
)
from scipy.stats import PermutationMethod
from process import get_data

df = get_data("files_for_templated")
a = df["spycas9"]
result_dict = {"b": [], "method": [], "alternative": [], "p-value": []}
for cas in ["spymac", "ispymac"]:
    b = df[cas]

    #####################################
    # The t-test on independent samples (or the so-called unpaired t-test).
    #####################################

    # equal_var=True performs normal unpaired t-test.
    # nan_policy="omit" omit missing data.
    # alternative="greater" means a is smaller than b, which decrease p-value.
    result_dict["b"].append(cas)
    result_dict["method"].append("t-test-ind")
    alternative = "greater"
    result_dict["alternative"].append(alternative)
    result_dict["p-value"].append(
        ttest_ind(
            a, b, equal_var=True, nan_policy="omit", alternative=alternative
        ).pvalue.item()
    )

    # equal_var=False performs Welch’s t-test.
    result_dict["b"].append(cas)
    result_dict["method"].append("Welch-test")
    alternative = "greater"
    result_dict["alternative"].append(alternative)
    result_dict["p-value"].append(
        ttest_ind(
            a, b, equal_var=False, nan_policy="omit", alternative=alternative
        ).pvalue.item()
    )

    # trim=0.05 performs a trimmed (Yuen’s) t-test, which is an extension of Welch’s t-test and help to filter outliers.
    result_dict["b"].append(cas)
    result_dict["method"].append("Yuen-test")
    alternative = "greater"
    result_dict["alternative"].append(alternative)
    result_dict["p-value"].append(
        ttest_ind(
            a, b, nan_policy="omit", alternative=alternative, trim=0.05
        ).pvalue.item()
    )

    # # method=PermutationMethod(n_resamples=9999) performs permutation test over t-statistics.
    # result_dict["b"].append(cas)
    # result_dict["method"].append("permute-t-test")
    # alternative = "greater"
    # result_dict["alternative"].append(alternative)
    # result_dict["p-value"].append(
    #     ttest_ind(
    #         a,
    #         b,
    #         nan_policy="omit",
    #         alternative=alternative,
    #         method=PermutationMethod(n_resamples=9999),
    #     ).pvalue.item()
    # )

    #####################################
    # The U-test on independent samples (or non-parametric version of the t-test for independent samples.)
    #####################################

    # The Mann-Whitney U test requires that the underlying distributions of two samples are the same. It does not care what the underlying distribution is.
    result_dict["b"].append(cas)
    result_dict["method"].append("U-test")
    alternative = "greater"
    result_dict["alternative"].append(alternative)
    result_dict["p-value"].append(
        mannwhitneyu(a, b, alternative=alternative, nan_policy="omit").pvalue.item()
    )

    # # method=PermutationMethod(n_resamples=9999) performs permutation test over U-statistics.
    # result_dict["b"].append(cas)
    # result_dict["method"].append("permute-U-test")
    # alternative = "greater"
    # result_dict["alternative"].append(alternative)
    # result_dict["p-value"].append(
    #     mannwhitneyu(
    #         a,
    #         b,
    #         alternative=alternative,
    #         nan_policy="omit",
    #         method=PermutationMethod(n_resamples=9999),
    #     ).pvalue.item()
    # )

    ######################################
    # The binomial test.
    ######################################
    result_dict["b"].append(cas)
    result_dict["method"].append("binomial-test")
    alternative = "greater"
    result_dict["alternative"].append(alternative)
    result_dict["p-value"].append(
        binomtest(
            (a > b).sum().item(),
            ((a > b).sum().item() + (a < b).sum().item()),
            p=0.5,
            alternative=alternative,
        ).pvalue.item()
    )

    # plt.close()
    # (a - b).hist(bins=30).get_figure().savefig(f"{cas}_diff_hist.png")

    #######################################
    # The t-test on two related samples (or the so-called paired t-test).
    #######################################
    result_dict["b"].append(cas)
    result_dict["method"].append("t-test-rel")
    alternative = "greater"
    result_dict["alternative"].append(alternative)
    result_dict["p-value"].append(
        ttest_rel(a, b, nan_policy="omit", alternative=alternative).pvalue.item()
    )

    ###########################################
    # The Wilcoxon signed-rank test (non-parametric version of the paired t-test).
    ###########################################
    result_dict["b"].append(cas)
    result_dict["method"].append("Wilcoxon-signed-rank-test")
    alternative = "greater"
    result_dict["alternative"].append(alternative)
    result_dict["p-value"].append(
        wilcoxon(
            a,
            b,
            alternative=alternative,
            method="auto",
            nan_policy="omit",
        ).pvalue.item()
    )

    # ############################################
    # # The Baumgartner-Weiss-Schindler test on two independent samples (cannot handle missing data).
    # ############################################
    # result_dict["b"].append(cas)
    # result_dict["method"].append("Baumgartner-Weiss-Schindler-test")
    # alternative = "greater"
    # result_dict["alternative"].append(alternative)
    # bws_test(a, b, alternative=alternative)

    ################################################
    # The Wilcoxon rank-sum test for two samples. The null hypothesis that two sets of measurements are drawn from the same distribution.
    ################################################
    result_dict["b"].append(cas)
    result_dict["method"].append("Wilcoxon-rank-sum-test")
    alternative = "greater"
    result_dict["alternative"].append(alternative)
    result_dict["p-value"].append(
        ranksums(a, b, alternative=alternative, nan_policy="omit").pvalue.item()
    )

    #################################################
    # The Brunner-Munzel test on two independent samples. It does not require the equal variance (like the Wilcoxon-Mann-Whitney’s U-test,) and same distribution.
    #################################################
    result_dict["b"].append(cas)
    result_dict["method"].append("Brunner-Munzel-test")
    alternative = "greater"
    result_dict["alternative"].append(alternative)
    result_dict["p-value"].append(
        brunnermunzel(a, b, alternative=alternative, nan_policy="omit").pvalue.item()
    )


result_df = pd.DataFrame(result_dict)
result_df.pivot(
    columns="b", index=["method", "alternative"], values="p-value"
).reset_index().to_csv("result.csv", index=False)
