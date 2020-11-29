
# clean values in Engineer's estimate column 
fp = open('data.csv', 'r')
fp_r = open('cleanedData.csv', 'w')

ESTIMATE_COL = 3
fp.readline()
for line in fp:
    values = line.replace('\n', '').split(',')

    estimates = values[ESTIMATE_COL].split()
    estimates_numeric = [] 

    for i in estimates:
        estimates_numeric.append(float(i.replace('%', '').replace('\"', '')))

    values.append(min(estimates_numeric))
    output = '' 
    for i, v in enumerate(values):
        output += str(v)

        if i < len(values)-1:
            output += ', '
    output += '\n' 
    fp_r.write(output)

fp_r.close()




