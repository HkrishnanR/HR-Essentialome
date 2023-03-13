"""
Created December 2021
@author: Harikrishnan Ramadasan <harikrishnan@students.iisertirupati.ac.in>
- HR Essentialome analysis
"""

import numpy as np
import pandas as pd
from math import isnan


def calculate_age():
    df2 = pd.read_csv(".../Taxons.txt", sep="\t")  # Path to file with all taxa
    d_l = df2["OMACode"].tolist()
    dfr = df2[df2["Taxon"].isin(tuple(["Fungi", "Viridiplanteae"]))]  # Separate out plants and fungus
    dfrr = dfr[dfr["Taxon"] == "Fungi"]
    dr = dfr[dfr["Taxon"] == "Viridiplanteae"]

    FG = dfrr["OMACode"].tolist()
    PL = dr["OMACode"].tolist()

    df = pd.read_csv("../oma_pairs_5.txt", sep="\t")  # file with orthology information

    print(df)
    df1 = pd.read_csv("../HRP_Typeadded", sep="\t")  # File with HR information
    df1 = df1[["Identifier", "AminoAcid", "TypeofHRP", "Code"]]
    d = dict(zip(df1["Identifier"], df1["TypeofHRP"]))
    print(df1)
    df["TypeA"] = df["ProteinA"].map(d)  # Map the type of HRs
    df["TypeB"] = df["ProteinB"].map(d)  # Map the type of HRs
    df["Code"] = [x[:5] for x in df["ProteinB"]]  # Extract OMA code of each protein
    df["TypeB"] = df["TypeB"].replace(np.nan, "NoHRP")
    df["TypeA"] = df["TypeA"].replace(np.nan, "NoHRP")
    k = dict(zip(df["ProteinA"], df["TypeA"]))
    print(df)

    x_l = df["Code"].unique().tolist()
    x_l = x_l + ["HUMAN"]

    dd = df.groupby('ProteinA').apply(lambda x: dict(zip(x["Code"], x["TypeB"]))).reset_index().rename(
        columns={0: "new"})  # Groupby proteins and get oma code and type of HR
    dff = pd.concat([dd.drop(['new'], axis=1), dd["new"].apply(pd.Series)], axis=1) # Split list to rows and append back to the original dataframe
    dff["HUMAN"] = dff["ProteinA"].map(k) # Map type of HR
    kl = [i for i in d_l if i in x_l]
    kl = ["ProteinA"] + kl
    dff = dff[kl] # Extracting columns of interest
    dffd = dff
    dff = dff.replace('NoHRP', np.nan)
    dff = dff.assign(last=dff.iloc[:, 1:].apply(lambda x: x.last_valid_index(), axis=1))  # Get last index where value not zero
    final = dff[["ProteinA", "last"]]
    print(final)

    ll = dict(zip(df2["OMACode"], df2["Taxon"]))
    final["Sp"] = final["last"].map(ll)  # map OMAcode to Taxon
    final["Sp"] = final["Sp"] + "_Specific"
    FF = final
    final = final[final["Sp"].isin(tuple(["Fungi_Specific", "Viridiplanteae_Specific"]))]
    dffd = dffd[dffd["ProteinA"].isin(tuple(final["ProteinA"].tolist()))]
    drd = dffd.set_index('ProteinA').T.to_dict('dict')

    for i, j in drd.items():
        for k in j.copy():
            if type(j[k]) == float and isnan(j[k]):
                del j[k]
    drr = {}
    for i in drd:
        drr[i] = []

    for i, j in drd.items():
        for l in j:
            for a in drr:
                if a == i:
                    drr[i].append(l)

    dfx = pd.DataFrame(drr.items(), columns=["Protein", "Species"])
    for i, rows in dfx.iterrows():  # matching each condition to Fungi_plant, Plant_specific and Fungi_specific
        if (len(list(set(dfx["Species"][i]) & set(PL))) >= 1) and (len(list(set(dfx["Species"][i]) & set(FG)))) >= 1:
            dfx.at[i, "Sp"] = "Fungi_plant"
        if (len(list(set(dfx["Species"][i]) & set(PL))) >= 1) and (len(list(set(dfx["Species"][i]) & set(FG)))) == 0:
            dfx.at[i, "Sp"] = "Plant_Specific"
        if (len(list(set(dfx["Species"][i]) & set(PL))) == 0) and (len(list(set(dfx["Species"][i]) & set(FG)))) >= 1:
            dfx.at[i, "Sp"] = "Fungi_Specific"

    ddx = dict(zip(dfx["Protein"], dfx["Sp"]))

    FF.loc[FF['ProteinA'].isin(ddx.keys()), 'Sp'] = FF['ProteinA'].map(ddx) # Selecting proteins falling to the three categories mentioned above

    FF.to_csv("HR_Age.txt", sep="\t")



