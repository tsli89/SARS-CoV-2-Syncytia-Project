import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
plt.close('all')
import openpyxl

def main():

    print("\nBeginning Program")
    args = sys.argv[1:]
    if len(args) != 2:
        raise IOError("Script must be called with three arguments specifying the metadata, WHO-COVID-19-global-data and  vaccination-data csv files.")
    metafile = args[0]
    cases = args[1]
    print("Reading metadata from %s" % metafile)

    if not os.path.isfile(metafile):
        raise IOERROR("Cannot find metadata file of %s in the current directory." % metafile)

    dataraw = pd.read_csv(metafile, sep="\t", usecols=["Virus name", "Collection date", "Host", "Clade", "Pango lineage", "Variant", "AA Substitutions", "Is complete?"])

    dataraw = dataraw[~dataraw['Collection date'].isin(["2020", "2021", "2022", "2023"])] # remove inaccurate/ambiguous dates ("2020", "2021",etc.)
    dataraw2 = dataraw
    dataraw3 = dataraw2[dataraw2["Host"] == 'Human'] # limit to human strains only, no bat, pangolin, mink, environmental, etc.
    dataraw4 = dataraw3[dataraw3["AA Substitutions"].str.find(")") > 0] # remove if unclosed parentheses
    data = dataraw4[dataraw4["Is complete?"] == True]
    data = data.set_index(data["Collection date"])

    data['Collection date'] = pd.to_datetime(data['Collection date'], format='ISO8601')

    data_time = pd.to_datetime(data.index, format='ISO8601')
    datetime_index = pd.DatetimeIndex(data_time.values)
    data = data.set_index(datetime_index)
    del data['Host']


    print("Reading metadata from %s" % cases)
    if not os.path.isfile(cases):
        raise IOERROR("Cannot find data file of %s in the current directory." % cases)
    casesraw = pd.read_csv(cases, usecols=["Date_reported", "New_cases"])
    cases = casesraw.set_index(casesraw["Date_reported"])
    cases['Date_reported'] = pd.to_datetime(cases['Date_reported'], format='ISO8601')
    cases_time = pd.to_datetime(cases.index, format='ISO8601')
    cases_time_index = pd.DatetimeIndex(cases_time.values)
    cases = cases.set_index(cases_time_index)
    cases = cases.sort_index()

    stats = pd.DataFrame(
        {
            "Variant": [
                "Original (No Spike Mutations)",
                "D614G Only",
                "Alpha",
                "Beta",
                "Gamma",
                "Delta",
                "BA.1",
                "BA.2",
                "BA.4",
                "BA.5",
                "BA.2.12.1",
                "BA.2.75",
                "BQ.1",
                "XBB",
            ],
            "First Appearance": ['', '', '', '', '', '', '', '', '', '', '', '', '', ''],
            "Last Appearance":  ['', '', '', '', '', '', '', '', '', '', '', '', '', ''],
            "Delta":            ['', '', '', '', '', '', '', '', '', '', '', '', '', ''],
            "Mean":             ['', '', '', '', '', '', '', '', '', '', '', '', '', ''],
            "Median":           ['', '', '', '', '', '', '', '', '', '', '', '', '', ''],
            "Max Incidence":    ['', '', '', '', '', '', '', '', '', '', '', '', '', ''],
            "Date of Max":      ['', '', '', '', '', '', '', '', '', '', '', '', '', ''],
        }
    )
    stats = stats.set_index(stats["Variant"])

    noSpike = data[data["AA Substitutions"].str.find("Spike") <= 0] # showing some omicrons and delta with no spike mutations. Exclude these.
    noSpike = noSpike[noSpike["Variant"].isnull()] # remove Omicron, etc, with no spike mutations
    stats.loc["Original (No Spike Mutations)", "First Appearance"] = noSpike["Collection date"].min()
    stats.loc["Original (No Spike Mutations)","Last Appearance"]  = noSpike["Collection date"].max()
    stats.loc["Original (No Spike Mutations)","Mean"]  = noSpike["Collection date"].mean()
    stats.loc["Original (No Spike Mutations)","Median"]  = noSpike["Collection date"].median()
    stats.loc["Original (No Spike Mutations)","Delta"]   = stats.loc["Original (No Spike Mutations)","Last Appearance"] - stats.loc["Original (No Spike Mutations)","First Appearance"]

    noSpikeTable = noSpike[["Collection date"]]
    noSpikeTable["Original"] = 1

    originala = data[data["Clade"] == 'O']
    original = originala[originala["Pango lineage"] == "Unassigned"]


    D614Gonlya = data[data["AA Substitutions"].str.count("Spike") == 1]
    D614Gonly = D614Gonlya[D614Gonlya["AA Substitutions"].str.count("D614G") == 1]
    stats.loc["D614G Only", "First Appearance"] = D614Gonly["Collection date"].min()
    stats.loc["D614G Only","Last Appearance"]  = D614Gonly["Collection date"].max()
    stats.loc["D614G Only","Mean"]  = D614Gonly["Collection date"].mean()
    stats.loc["D614G Only","Median"]  = D614Gonly["Collection date"].median()
    stats.loc["D614G Only","Delta"]   = stats.loc["D614G Only","Last Appearance"] - stats.loc["D614G Only","First Appearance"]

    D614GonlyTable = D614Gonly[["Collection date"]]
    D614GonlyTable["D614G Only"] = 1

    alpha = data[data["Variant"] == 'Former VOC Alpha GRY (B.1.1.7+Q.*) first detected in the UK']
    stats.loc["Alpha", "First Appearance"] = alpha["Collection date"].min()
    stats.loc["Alpha","Last Appearance"]  = alpha["Collection date"].max()
    stats.loc["Alpha","Mean"]  = alpha["Collection date"].mean()
    stats.loc["Alpha","Median"]  = alpha["Collection date"].median()
    stats.loc["Alpha","Delta"]   = stats.loc["Alpha","Last Appearance"] - stats.loc["Alpha","First Appearance"]
    alphaTable = alpha[["Collection date"]]
    alphaTable["Alpha"] = 1

    beta = data[data["Variant"] == 'Former VOC Beta GH/501Y.V2 (B.1.351+B.1.351.2+B.1.351.3) first detected in South Africa']
    stats.loc["Beta", "First Appearance"] = beta["Collection date"].min()
    stats.loc["Beta","Last Appearance"]  = beta["Collection date"].max()
    stats.loc["Beta","Mean"]  = beta["Collection date"].mean()
    stats.loc["Beta","Median"]  = beta["Collection date"].median()
    stats.loc["Beta","Delta"]   = stats.loc["Beta","Last Appearance"] - stats.loc["Beta","First Appearance"]
    betaTable = beta[["Collection date"]]
    betaTable["Beta"] = 1


    gamma = data[data["Variant"] == 'Former VOC Gamma GR/501Y.V3 (P.1+P.1.*) first detected in Brazil/Japan']
    stats.loc["Gamma", "First Appearance"] = gamma["Collection date"].min()
    stats.loc["Gamma","Last Appearance"]  = gamma["Collection date"].max()
    stats.loc["Gamma","Mean"]  = gamma["Collection date"].mean()
    stats.loc["Gamma","Median"]  = gamma["Collection date"].median()
    stats.loc["Gamma","Delta"]   = stats.loc["Gamma","Last Appearance"] - stats.loc["Gamma","First Appearance"]
    gammaTable = gamma[["Collection date"]]
    gammaTable["Gamma"] = 1

    delta = data[data["Variant"] == 'Former VOC Delta GK (B.1.617.2+AY.*) first detected in India']
    stats.loc["Delta", "First Appearance"] = delta["Collection date"].min()
    stats.loc["Delta","Last Appearance"]  = delta["Collection date"].max()
    stats.loc["Delta","Mean"]  = delta["Collection date"].mean()
    stats.loc["Delta","Median"]  = delta["Collection date"].median()
    stats.loc["Delta","Delta"]   = stats.loc["Delta","Last Appearance"] - stats.loc["Delta","First Appearance"]
    deltaTable = delta[["Collection date"]]
    deltaTable["Delta"] = 1

    ba1 = data[data["Pango lineage"] == 'BA.1']
    stats.loc["BA.1", "First Appearance"] = ba1["Collection date"].min()
    stats.loc["BA.1","Last Appearance"]  = ba1["Collection date"].max()
    stats.loc["BA.1","Mean"]  = ba1["Collection date"].mean()
    stats.loc["BA.1","Median"]  = ba1["Collection date"].median()
    stats.loc["BA.1","Delta"]   = stats.loc["BA.1","Last Appearance"] - stats.loc["BA.1","First Appearance"]
    ba1Table = ba1[["Collection date"]]
    ba1Table["BA.1"] = 1

    ba2 = data[data["Pango lineage"] == 'BA.2']
    stats.loc["BA.2", "First Appearance"] = ba2["Collection date"].min()
    stats.loc["BA.2","Last Appearance"]  = ba2["Collection date"].max()
    stats.loc["BA.2","Mean"]  = ba2["Collection date"].mean()
    stats.loc["BA.2","Median"]  = ba2["Collection date"].median()
    stats.loc["BA.2","Delta"]   = stats.loc["BA.2","Last Appearance"] - stats.loc["BA.2","First Appearance"]
    ba2Table = ba2[["Collection date"]]
    ba2Table["BA.2"] = 1

    ba4 = data[data["Pango lineage"] == 'BA.4']
    stats.loc["BA.4", "First Appearance"] = ba4["Collection date"].min()
    stats.loc["BA.4","Last Appearance"]  = ba4["Collection date"].max()
    stats.loc["BA.4","Mean"]  = ba4["Collection date"].mean()
    stats.loc["BA.4","Median"]  = ba4["Collection date"].median()
    stats.loc["BA.4","Delta"]   = stats.loc["BA.4","Last Appearance"] - stats.loc["BA.4","First Appearance"]
    ba4Table = ba4[["Collection date"]]
    ba4Table["BA.4"] = 1

    ba5 = data[data["Pango lineage"] == 'BA.5']
    stats.loc["BA.5", "First Appearance"] = ba5["Collection date"].min()
    stats.loc["BA.5","Last Appearance"]  = ba5["Collection date"].max()
    stats.loc["BA.5","Mean"]  = ba5["Collection date"].mean()
    stats.loc["BA.5","Median"]  = ba5["Collection date"].median()
    stats.loc["BA.5","Delta"]   = stats.loc["BA.5","Last Appearance"] - stats.loc["BA.5","First Appearance"]
    ba5Table = ba5[["Collection date"]]
    ba5Table["BA.5"] = 1

    ba2121 = data[data["Pango lineage"] == 'BA.2.12.1']
    stats.loc["BA.2.12.1", "First Appearance"] = ba2121["Collection date"].min()
    stats.loc["BA.2.12.1","Last Appearance"]  = ba2121["Collection date"].max()
    stats.loc["BA.2.12.1","Mean"]  = ba2121["Collection date"].mean()
    stats.loc["BA.2.12.1","Median"]  = ba2121["Collection date"].median()
    stats.loc["BA.2.12.1","Delta"]   = stats.loc["BA.2.12.1","Last Appearance"] - stats.loc["BA.2.12.1","First Appearance"]
    ba2121Table = ba2121[["Collection date"]]
    ba2121Table["BA.2.12.1"] = 1

    ba275 = data[data["Pango lineage"] == 'BA.2.75']
    stats.loc["BA.2.75", "First Appearance"] = ba275["Collection date"].min()
    stats.loc["BA.2.75","Last Appearance"]  = ba275["Collection date"].max()
    stats.loc["BA.2.75","Mean"]  = ba275["Collection date"].mean()
    stats.loc["BA.2.75","Median"]  = ba275["Collection date"].median()
    stats.loc["BA.2.75","Delta"]   = stats.loc["BA.2.75","Last Appearance"] - stats.loc["BA.2.75","First Appearance"]
    ba275Table = ba275[["Collection date"]]
    ba275Table["BA.2.75"] = 1

    bq1 = data[data["Pango lineage"] == 'BQ.1']
    stats.loc["BQ.1", "First Appearance"] = bq1["Collection date"].min()
    stats.loc["BQ.1","Last Appearance"]  = bq1["Collection date"].max()
    stats.loc["BQ.1","Mean"]  = bq1["Collection date"].mean()
    stats.loc["BQ.1","Median"]  = bq1["Collection date"].median()
    stats.loc["BQ.1","Delta"]   = stats.loc["BQ.1","Last Appearance"] - stats.loc["BQ.1","First Appearance"]
    bq1Table = bq1[["Collection date"]]
    bq1Table["BQ.1"] = 1

    xbb = data[data["Pango lineage"].str.startswith("XBB") == True]
    stats.loc["XBB", "First Appearance"] = xbb["Collection date"].min()
    stats.loc["XBB","Last Appearance"]  = xbb["Collection date"].max()
    stats.loc["XBB","Mean"]  = xbb["Collection date"].mean()
    stats.loc["XBB","Median"]  = xbb["Collection date"].median()
    stats.loc["XBB","Delta"]   = stats.loc["XBB","Last Appearance"] - stats.loc["XBB","First Appearance"]
    xbbTable = xbb[["Collection date"]]
    xbbTable["XBB"] = 1

    result = pd.concat([noSpikeTable, D614GonlyTable, alphaTable, betaTable, gammaTable, deltaTable, ba1Table, ba2Table, ba4Table, ba5Table, ba2121Table, ba275Table, bq1Table, xbbTable], axis=0)

    casesresult = pd.concat([result, cases], axis=0)
    del casesresult['Collection date']
    del casesresult['Date_reported']
    casesresult = casesresult.sort_index()
    casesresult = casesresult.resample("2W").sum()

    del result['Collection date']
    result2 = result.sort_index()
    result3 = result2.resample("2W").sum()

    stats.loc["Original (No Spike Mutations)","Max Incidence"]  = result3["Original"].max()
    stats.loc["Original (No Spike Mutations)","Date of Max"]  = result3.index[result3["Original"] == result3["Original"].max()].tolist()
    stats.loc["D614G Only","Max Incidence"]  = result3["D614G Only"].max()
    stats.loc["D614G Only", "Date of Max"] = result3.index[result3["D614G Only"] == result3["D614G Only"].max()].tolist()
    stats.loc["Alpha","Max Incidence"]  = result3["Alpha"].max()
    stats.loc["Alpha","Date of Max"]  = result3.index[result3["Alpha"] == result3["Alpha"].max()].tolist()
    stats.loc["Beta","Max Incidence"]  = result3["Beta"].max()
    stats.loc["Beta","Date of Max"]  = result3.index[result3["Beta"] == result3["Beta"].max()].tolist()
    stats.loc["Gamma","Max Incidence"]  = result3["Gamma"].max()
    stats.loc["Gamma","Date of Max"]  = result3.index[result3["Gamma"] == result3["Gamma"].max()].tolist()
    stats.loc["Delta","Max Incidence"]  = result3["Delta"].max()
    stats.loc["Delta","Date of Max"]  = result3.index[result3["Delta"] == result3["Delta"].max()].tolist()
    stats.loc["BA.1","Max Incidence"]  = result3["BA.1"].max()
    stats.loc["BA.1","Date of Max"]  = result3.index[result3["BA.1"] == result3["BA.1"].max()].tolist()
    stats.loc["BA.2","Max Incidence"]  = result3["BA.2"].max()
    stats.loc["BA.2","Date of Max"]  = result3.index[result3["BA.2"] == result3["BA.2"].max()].tolist()
    stats.loc["BA.4","Max Incidence"]  = result3["BA.4"].max()
    stats.loc["BA.4","Date of Max"]  = result3.index[result3["BA.4"] == result3["BA.4"].max()].tolist()
    stats.loc["BA.5","Max Incidence"]  = result3["BA.5"].max()
    stats.loc["BA.5","Date of Max"]  = result3.index[result3["BA.5"] == result3["BA.5"].max()].tolist()
    stats.loc["BA.2.12.1","Max Incidence"]  = result3["BA.2.12.1"].max()
    stats.loc["BA.2.12.1","Date of Max"]  = result3.index[result3["BA.2.12.1"] == result3["BA.2.12.1"].max()].tolist()
    stats.loc["BA.2.75","Max Incidence"]  = result3["BA.2.75"].max()
    stats.loc["BA.2.75","Date of Max"]  = result3.index[result3["BA.2.75"] == result3["BA.2.75"].max()].tolist()
    stats.loc["BQ.1","Max Incidence"]  = result3["BQ.1"].max()
    stats.loc["BQ.1","Date of Max"]  = result3.index[result3["BQ.1"] == result3["BQ.1"].max()].tolist()
    stats.loc["XBB","Max Incidence"]  = result3["XBB"].max()
    stats.loc["XBB","Date of Max"]  = result3.index[result3["XBB"] == result3["XBB"].max()].tolist()

    fractionTable = result3[["Original", "D614G Only", "Alpha", "Beta", "Gamma", "Delta", "BA.1", "BA.2", "BA.4", "BA.5", "BA.2.12.1", "BA.2.75", "BQ.1", "XBB"]]
    fractionTable["Total"] = fractionTable["Original"] + fractionTable["D614G Only"] + fractionTable["Alpha"] + fractionTable["Beta"] + fractionTable["Gamma"] + fractionTable["Delta"] + fractionTable["BA.1"] + fractionTable["BA.2"] + fractionTable["BA.4"] + fractionTable["BA.5"] + fractionTable["BA.2.12.1"] + fractionTable["BA.2.75"] + fractionTable["BQ.1"] + fractionTable["XBB"]
    fractionTable["Original"] = (fractionTable["Original"] / fractionTable["Total"])
    fractionTable["D614G Only"] = (fractionTable["D614G Only"] / fractionTable["Total"])
    fractionTable["Alpha"] = (fractionTable["Alpha"] / fractionTable["Total"])
    fractionTable["Beta"] = (fractionTable["Beta"] / fractionTable["Total"])
    fractionTable["Gamma"] = (fractionTable["Gamma"] / fractionTable["Total"])
    fractionTable["Delta"] = (fractionTable["Delta"] / fractionTable["Total"])
    fractionTable["BA.1"] = (fractionTable["BA.1"] / fractionTable["Total"])
    fractionTable["BA.2"] = (fractionTable["BA.2"] / fractionTable["Total"])
    fractionTable["BA.4"] = (fractionTable["BA.4"] / fractionTable["Total"])
    fractionTable["BA.5"] = (fractionTable["BA.5"] / fractionTable["Total"])
    fractionTable["BA.2.12.1"] = (fractionTable["BA.2.12.1"] / fractionTable["Total"])
    fractionTable["BA.2.75"] = (fractionTable["BA.2.75"] / fractionTable["Total"])
    fractionTable["BQ.1"] = (fractionTable["BQ.1"] / fractionTable["Total"])
    fractionTable["XBB"] = (fractionTable["XBB"] / fractionTable["Total"])
    del fractionTable["Total"]

    ax = fractionTable.plot.area(legend=False)

    casesresult["New_cases"].plot(secondary_y=True, color="White");

    ax.set_ylabel("Variant fraction of Total");
    ax.right_ax.set_ylabel("New cases");

    result3.to_excel("Variant_counts.xlsx", sheet_name="Variant counts", index=True)
    stats.to_excel("stats.xlsx", sheet_name="stats", index=False)
    fractionTable.to_excel("Variant_fractions.xlsx", sheet_name="Variant fractions", index=True)
    casesresult.to_excel("New_cases.xlsx", sheet_name="New_cases", index=True)
#    ax = result3.plot.area()

    test = 12

main()


