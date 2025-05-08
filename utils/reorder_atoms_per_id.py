def reorder_and_renumber_atoms(input_file_path, output_file_path):
    """
    Ensures equal number of Ni (type 1) and Ti (type 2) atoms.
    Converts excess atoms from the majority type to the minority type and renumbers IDs.
    """
    with open(input_file_path, 'r') as input_file:
        lines = input_file.readlines()

    # Find start of atoms
    atom_section_start = next(i for i, line in enumerate(lines) if line.strip() == "Atoms") + 2
    atom_lines = [line.strip() for line in lines[atom_section_start:] if line.strip()]

    id_1_lines = [line for line in atom_lines if line.split()[1] == "1"]
    id_2_lines = [line for line in atom_lines if line.split()[1] == "2"]

    n1, n2 = len(id_1_lines), len(id_2_lines)
    print(f'Before correction: nickel: {n1}, titanium: {n2}')

    # Balance atom types
    if n1 > n2:
        diff = (n1 - n2) // 2
        # convert atom_type 2 to type 1
        to_convert = id_1_lines[-diff:]
        id_1_lines = id_1_lines[:-diff]
        converted = []
        for line in to_convert:
            parts = line.split()
            parts[1] = "2"
            converted.append(" ".join(parts))
        id_2_lines.extend(converted)
        # convert atom_type 1 to type 2
    elif n2 > n1:
        diff = (n2 - n1) // 2
        to_convert = id_2_lines[-diff:]
        id_2_lines = id_2_lines[:-diff]
        converted = []
        for line in to_convert:
            parts = line.split()
            parts[1] = "1"
            converted.append(" ".join(parts))
        id_1_lines.extend(converted)

    # Ensure final match
    n1_final, n2_final = len(id_1_lines), len(id_2_lines)
    if n1_final != n2_final:
        raise ValueError(f"Mismatch still exists after correction: Ni = {n1_final}, Ti = {n2_final}")

    print(f'After correction: nickel: {n1_final}, titanium: {n2_final}')

    organized_lines = id_1_lines + id_2_lines
    renumbered_lines = []
    for new_id, line in enumerate(organized_lines, start=1):
        parts = line.split()
        parts[0] = str(new_id)
        renumbered_lines.append(" ".join(parts))

    with open(output_file_path, 'w') as output_file:
        output_file.write("\n".join(renumbered_lines) + "\n")
    return n1_final,n2_final