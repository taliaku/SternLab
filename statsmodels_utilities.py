#! /powerapps/share/python-anaconda-3.6/bin/python
import statsmodels.api as sm


def parse_summary_MH_test(summary):
    #Mantel-Haenszel test
    res = {}
    sum_csv = summary.as_csv()
    sum_csv = sum_csv.split("\n")
    pooled_odds = [i.strip() for i in sum_csv[1].split(",")]
    res["pooled_odds_estimate"] = pooled_odds[1]
    res["pooled_odds_LCB"] = pooled_odds[2]
    res["pooled_odds_UCB"] = pooled_odds[3]
    pooled_log_odds = [i.strip() for i in sum_csv[2].split(",")]
    res["pooled_log_odds_estimate"] = pooled_log_odds[1]
    res["pooled_log_odds_LCB"] = pooled_log_odds[2]
    res["pooled_log_odds_UCB"] = pooled_log_odds[3]
    pooled_risk_ratio = [i.strip() for i in sum_csv[3].split(",")]
    res["pooled_risk_ratio"] =  pooled_risk_ratio[1]
    OR_test_1 =  [i.strip() for i in sum_csv[6].split(",")]
    res["OR_1_statistic"] = OR_test_1[1]
    res["OR_1_pvalue"] = float(OR_test_1[2])
    OR_constant =  [i.strip() for i in sum_csv[7].split(",")]
    res["OR_constant_statistic"] = OR_constant[1]
    res["OR_constant_pvalue"] = OR_constant[2]
    res["table_num"] = [i.strip() for i in sum_csv[9].split(",")][1]
    res["min_n"] = [i.strip() for i in sum_csv[10].split(",")][1]
    res["max_n"] = [i.strip() for i in sum_csv[11].split(",")][1]
    res["avg_n"] = [i.strip() for i in sum_csv[12].split(",")][1]
    res["total_n"] = [i.strip() for i in sum_csv[13].split(",")][1]

    return res