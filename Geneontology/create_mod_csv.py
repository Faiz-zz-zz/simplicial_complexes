import json
import csv
simp_data = json.loads(open("all_pred.json").read())
ppi_data = json.loads(open("all_pred_ppi.json").read())


method_map = {}
global_comb = None
for method in ["random_forest_classifier", "ranodm_forest_regressor", "linear_regression", "logistic_regression"]:
    print(method)
    go_ids = list(set(simp_data[method].keys()) & set(ppi_data[method].keys()))
    average_effec = {} # (comb): {simp: (cal, cnt), ppi: (val, cnt)}
    for go_id in go_ids:
        preds_simp, pred_ppi = simp_data[method][go_id], ppi_data[method][go_id]
        combs = preds_simp.keys()
        global_comb = combs
        for comb in combs:
            if comb not in average_effec:
                average_effec[comb] = {"simp": [0, 0], "ppi": [0, 0]}
            average_effec[comb]["simp"][0] += preds_simp[comb]["rms"]
            average_effec[comb]["simp"][1] += 1
            average_effec[comb]["ppi"][0] += pred_ppi[comb]["rms"]
            average_effec[comb]["ppi"][1] += 1
    method_map[method] = average_effec

with open("metrics.csv", "w") as outfile:
    fieldnames = ["Method"]
    for comb in global_comb:
        fieldnames.append(comb + " PPI")
        fieldnames.append(comb + " SIMP")
        fieldnames.append(comb + " Gain")        

    writer = csv.DictWriter(outfile, fieldnames=fieldnames)
    writer.writeheader()
    for method in method_map:
        result = {name: None for name in fieldnames}
        result["Method"] = method
        for comb in global_comb:
            result[comb + " PPI"] = method_map[method][comb]["ppi"][0] / method_map[method][comb]["ppi"][1]
            result[comb + " SIMP"] = method_map[method][comb]["simp"][0] / method_map[method][comb]["simp"][1]
            result[comb + " Gain"] = "{0:.2f}%".format(((result[comb + " SIMP"] - result[comb + " PPI"]) / abs(result[comb + " PPI"])) * 100)
        writer.writerow(result)




    

