patients = [[70, 1.8], [80, 1.9], [150, 1.7], [90, 2.8]]

def calculate_bmi(weight, height):
    return weight / (height ** 2)

for patient in patients:

    weight, height = patient
    bmi = calculate_bmi(weight, height)
    print("Patient's BMI is: %f" % bmi)

