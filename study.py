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
from plotnine import (
    ggplot,
    aes,
    geom_violin,
    geom_sina,
    geom_boxplot,
    geom_point,
    scale_x_discrete,
    scale_color_manual,
    scale_fill_manual,
)

methods = ["t-test-ind", "Welch-test", "Yuen-test", "t-test-rel"]
# methods = [
#     "t-test-ind",
#     "Welch-test",
#     "Yuen-test",
#     "t-test-rel",
#     "U-test",
#     "binomial-test",
#     "Wilcoxon-signed-rank-test",
#     "Wilcoxon-rank-sum-test",
#     "Brunner-Munzel-test",
# ]
for target in ["tem_1", "tem_2", "tem_3", "tem_4"]:
    df = get_data(
        data_dir="tem_nofilter",
        target=target,
        sum_x_thres=100,
        mut_fre_x_thres={"spycas9": 0.08, "spymac": 0.02, "ispymac": 0.02},
    )
    df.to_csv(f"{target}.csv", index=False)
    (
        ggplot(
            df.melt(
                id_vars="sgrna",
                value_vars=["spycas9", "spymac", "ispymac"],
                var_name="cas",
                value_name=target,
            ),
            mapping=aes(x="cas", y=target),
        )
        + geom_violin()
        + geom_sina(mapping=aes(fill="cas", color="cas"), alpha=0.01)
        + geom_boxplot(width=0.05, outlier_alpha=0.0)
        + scale_x_discrete(limits=["spycas9", "spymac", "ispymac"])
        + scale_color_manual(
            values={"spycas9": "#FF0000", "spymac": "#00FF00", "ispymac": "#0000FF"}
        )
        + scale_fill_manual(
            values={"spycas9": "#FF0000", "spymac": "#00FF00", "ispymac": "#0000FF"}
            # values=["#FF0000", "#00FF00", "#0000FF"],
        )
    ).save(f"{target}.png")
    for cas in ["spycas9", "spymac", "ispymac"]:
        (
            ggplot(
                df[[cas]].fillna(0.0).assign(sgRNA_idx=lambda df: np.arange(len(df))),
                mapping=aes(x=cas, y="sgRNA_idx"),
            )
            + geom_point()
        ).save(f"{cas}_{target}.png")

    a = df["spycas9"]
    result_dict = {"b": [], "method": [], "alternative": [], "p-value": []}
    for cas in ["spymac", "ispymac"]:
        for alternative in ["less", "greater"]:
            b = df[cas]

            #####################################
            # The t-test on independent samples (or the so-called unpaired t-test).
            #####################################

            # equal_var=True performs normal unpaired t-test.
            # nan_policy="omit" omit missing data.
            result_dict["b"].append(cas)
            result_dict["method"].append("t-test-ind")
            result_dict["alternative"].append(alternative)
            result_dict["p-value"].append(
                ttest_ind(
                    a, b, equal_var=True, nan_policy="omit", alternative=alternative
                ).pvalue.item()
            )

            # equal_var=False performs Welch’s t-test.
            result_dict["b"].append(cas)
            result_dict["method"].append("Welch-test")
            result_dict["alternative"].append(alternative)
            result_dict["p-value"].append(
                ttest_ind(
                    a, b, equal_var=False, nan_policy="omit", alternative=alternative
                ).pvalue.item()
            )

            # trim=0.05 performs a trimmed (Yuen’s) t-test, which is an extension of Welch’s t-test and help to filter outliers.
            result_dict["b"].append(cas)
            result_dict["method"].append("Yuen-test")
            result_dict["alternative"].append(alternative)
            result_dict["p-value"].append(
                ttest_ind(
                    a, b, nan_policy="omit", alternative=alternative, trim=0.05
                ).pvalue.item()
            )

            # # method=PermutationMethod(n_resamples=9999) performs permutation test over t-statistics.
            # result_dict["b"].append(cas)
            # result_dict["method"].append("permute-t-test")
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
            result_dict["alternative"].append(alternative)
            result_dict["p-value"].append(
                mannwhitneyu(
                    a, b, alternative=alternative, nan_policy="omit"
                ).pvalue.item()
            )

            # # method=PermutationMethod(n_resamples=9999) performs permutation test over U-statistics.
            # result_dict["b"].append(cas)
            # result_dict["method"].append("permute-U-test")
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
            result_dict["alternative"].append(alternative)
            result_dict["p-value"].append(
                ttest_rel(
                    a, b, nan_policy="omit", alternative=alternative
                ).pvalue.item()
            )

            ###########################################
            # The Wilcoxon signed-rank test (non-parametric version of the paired t-test).
            ###########################################
            result_dict["b"].append(cas)
            result_dict["method"].append("Wilcoxon-signed-rank-test")
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
            # result_dict["alternative"].append(alternative)
            # bws_test(a, b, alternative=alternative)

            ################################################
            # The Wilcoxon rank-sum test for two samples. The null hypothesis that two sets of measurements are drawn from the same distribution.
            ################################################
            result_dict["b"].append(cas)
            result_dict["method"].append("Wilcoxon-rank-sum-test")
            result_dict["alternative"].append(alternative)
            result_dict["p-value"].append(
                ranksums(a, b, alternative=alternative, nan_policy="omit").pvalue.item()
            )

            #################################################
            # The Brunner-Munzel test on two independent samples. It does not require the equal variance (like the Wilcoxon-Mann-Whitney’s U-test,) and same distribution.
            #################################################
            result_dict["b"].append(cas)
            result_dict["method"].append("Brunner-Munzel-test")
            result_dict["alternative"].append(alternative)
            result_dict["p-value"].append(
                brunnermunzel(
                    a, b, alternative=alternative, nan_policy="omit"
                ).pvalue.item()
            )

    result_df = pd.DataFrame(result_dict).query("method in @methods")
    result_df.pivot(
        columns="b", index=["method", "alternative"], values="p-value"
    ).reset_index().to_csv(f"result_{target}.csv", index=False)
