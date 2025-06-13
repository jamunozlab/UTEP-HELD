def read_dos(path_dos):
    column1 = []
    column2 = []
    with open(path_dos, 'r') as f:
        content = f.readlines()
        for line in content:
            # print(line)
            parts = line.strip().split()
            if len(parts) == 2: 
                column1.append(parts[0])
                column2.append(parts[1])
    return column1, column2