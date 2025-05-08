import os

def write_poscar(data, lat_par, path=".", filename="POSCAR"):
    full_path = os.path.join(path, filename)
    with open(full_path, 'w') as file:
        file.write(f"{data['System']}\n")
        file.write(f"{lat_par:.2f}\n")  # Use passed-in lat_par here
        for vector in data['Lattice Vectors']:
            file.write(" ".join(f"{x:.16f}" for x in vector) + "\n")
        file.write(" ".join(data['Elements'][0].split()) + "\n")  # ensure space between "Ni" and "Ti"
        file.write(" ".join(str(count) for count in data['Element Counts']) + "\n")
        file.write(f"{data['Coordinate Type']}\n")
        for coord in data['Coordinates']:
            file.write(" ".join(f"{x:.16f}" for x in coord) + "\n")
