# from surfleetDataRead import dataRead

def dataRead(path):

    # Input: A full file path
    # *** MAKE SURE FILE PATH IS A RAW STRING
    # File can only be either comma seperated or space seperated

    # Output: A dictionary with each data column in file stored as lists of float values.

    f = open(path)
    DATA = f.read()
    f.close()
    
    data = {
        'master' : DATA.splitlines()
    }

    sep = ' '
    if ',' in data['master'][1]:
        sep = ','

    for j in range(0,len(data['master'][1].split(sep))):
        data['{}'.format(j)] = [float(data['master'][i].split(sep)[j]) for i in range(1,len(data['master']))]

    data.pop('master')

    return data