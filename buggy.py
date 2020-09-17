patients = [[70, 1.8], [80, 1.9], [150, 1.7]]

def calculate_bmis(wheight, hweight):
    return wheight / (hweight ** 2)

for patient in patients:
    wheight, hweight = patient[1]
    bims = calculate_bmis(hwight, wheight)
    print("Patient's BMI is: %f" % bims)

