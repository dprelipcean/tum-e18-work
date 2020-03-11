import csv


with open('./CP_MC_data_SPD.CP_MC') as f:
    data = list()
    for row in f:
        mother, bachelor, isobar = row.split()

        try:
            data.append([float(mother), float(bachelor), float(isobar)])
        except ValueError:
            data.pop(-1)
            continue


with open('./CP_MC_data_SPD_sorted.CP_MC', mode='w') as f:
    writer = csv.writer(f, delimiter=' ')
    data.sort(key=lambda x: x[1])  # Sort by bachelor mass
    for row in data:
        writer.writerow(row)
