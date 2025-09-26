#!/usr/bin/env python

import os
import re
import pandas as pd
import pathlib


def get_data(data_dir: os.PathLike) -> pd.DataFrame:
    data_dir = pathlib.Path(os.fspath(data_dir))
    dfs = []
    for file in os.listdir(data_dir):
        if (
            not re.search(r"^36t-", file)
            and not re.search(r"^B2-", file)
            and not re.search(r"^A2-", file)
            and not re.search(r"^D2-", file)
            and not re.search(r"^i10t-", file)
        ):
            continue

        df = pd.read_csv(data_dir / file, header=0)
        if re.search(r"^36t-", file) or re.search(r"^B2-", file):
            df["cas"] = "spymac"
        elif re.search(r"^A2-", file) or re.search(r"^D2-", file):
            df["cas"] = "spycas9"
        elif re.search(r"^i10t-", file):
            df["cas"] = "ispymac"
        else:
            raise ValueError("Unexpected sample")

        mat = re.search(r"([^-]+)-[^-\d](\d)n?+-(\d)-wt(\d)", file)
        df["cellline"] = mat.group(1)
        df["chip"] = int(mat.group(2))
        df["bio_rep"] = int(mat.group(3))
        df["tech_rep"] = int(mat.group(4))

        dfs.append(df)

    # df = (
    #     pd.concat(dfs)
    #     .value_counts(["sgrna", "cas", "cellline", "chip", "bio_rep", "tech_rep"])
    #     .reset_index()
    #     .query("count > 1")[["cas", "count"]]
    #     .value_counts()
    # )

    # df = (
    #     pd.concat(dfs)
    #     .groupby(["sgrna", "cas", "cellline", "chip", "bio_rep", "tech_rep"])["tem_1"]
    #     .min()
    #     .reset_index(drop=False)
    #     .groupby(["sgrna", "cas"])["tem_1"]
    #     .mean()
    #     .reset_index(drop=False)
    #     .pivot(
    #         columns=["cas"],
    #         values="tem_1",
    #         index="sgrna",
    #     )
    # ).reset_index(drop=False)

    df = (
        (pd.concat(dfs).groupby(["sgrna", "cas"])["tem_1"].mean())
        .reset_index(drop=False)
        .pivot(
            columns=["cas"],
            values="tem_1",
            index="sgrna",
        )
    )

    return df
