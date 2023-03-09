"""
Created April 2021
@author: Harikrishnan Ramadasan
- HR Essentialome analysis
"""

import pandas as pd

def find_range(x):
    if x["Start"] <= x["Site"] <= x["Stop"]:
        return x["Site"]


def find_leng(x):
    return len(x["Count"])


def get_count(x) -> dict:
    """

    :param x: rows of dataframe
    :return: dictionary with PTM, counts per position
    """
    counts = dict()
    for i in x["Real"]:
        counts[int(i)] = counts.get(int(i), 0) + 1  # get all ptms for each site

    sorted_counts = dict(sorted(counts.items()))

    return sorted_counts


def split_count(x):
    """
    :param x: rows of dataframe
    :return: Total PTMs
    """
    new = []
    for i in x["Count"]:
        j = i.split(":")[1]
        new.append(int(j))
    return sum(new)


def filter_matches(zeros: bool = False):
    df = pd.read_csv("dbPTM2021_AllPTMs.txt", sep="\t")
    print(df)
    df1 = pd.read_csv("AllHRPs.txt", sep="\t")
    print(df1)
    df2 = pd.merge(df, df1, left_on="UniprotID", right_on="UniprotID")

    df2["Real"] = df2.apply(find_range, axis=1)
    df2 = df2[~df2["Real"].isna()]  # removes all rows not falling in range

    df2 = df2.groupby(["UniprotID", "AminoAcid", "Start", "Length", "Stop"])[
        'Real'].agg(list).reset_index()
    df2["HR_Id"] = df2["UniprotID"] + "_" + df2["AminoAcid"] + "_" + df2["Start"].astype(
        str)  # HR_id is a unique identifier for each HR

    if zeros:
        df2 = df2[['UniprotID', 'HR_Id', 'Length']]
        df2.to_csv("HR_PTM_Aggregate_2021_zeros.txt",
                   sep="\t")  # gives all rows that doesn't fall in range,ie positions where PTM dont occur

        return df2

    df2["Count"] = df2.apply(get_count, axis=1)
    df2["Number_of_Positions"] = df2.apply(find_leng, axis=1)
    df2['Count'] = df2['Count'].astype(str).str.strip(
        '{}').str.split(', ')
    df2["No_of_Subs"] = df2.apply(split_count, axis=1)
    df2 = df2[["UniprotID", "Length",
               "Number_of_Positions", "No_of_Subs", "HR_Id"]]

    df2.loc[len(df2), ['UniprotID', "HR_Id", 'Length', "Number_of_Positions", "No_of_Subs"]] = ['Total', (
            df2["Number_of_Positions"].sum(
            ) + df2["No_of_Subs"].sum()) / df2["Length"].sum(),
                                                                                                df2["Length"].sum(),
                                                                                                df2[
                                                                                                    "Number_of_Positions"].sum(),
                                                                                                df2[
                                                                                                    "No_of_Subs"].sum()]  # calculate final rows for density and number of PTMS and positions

    df2.to_csv("HR_PTM_Aggregate_2021_mod.txt", sep="\t")

    return df2
