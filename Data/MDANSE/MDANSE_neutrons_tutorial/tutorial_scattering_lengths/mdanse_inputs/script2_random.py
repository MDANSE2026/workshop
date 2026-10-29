import random

N_ATOMS = 2048  # this should match the number of atoms in your trajectory.

isotopes = ["Atm1"] * int(N_ATOMS/2) + ["Atm2"] * int(N_ATOMS/2)

mapping = {}
for i in range(N_ATOMS):
    random_element = random.choice(isotopes)
    isotopes.remove(random_element)
    mapping[str(i)] = random_element

print(str(mapping).replace("'", '"'))
