from rdkit import Chem


def substructure_search(mols: list, mol: str) -> list:
    """
    Returns list of molecules with specified substructure or messages that there are no such molecules.
    ARGUMENTS:
        - mols: list of molecules in SMILES notation
        - mol: sought substructure in SMILES notation
    RETURNS:
        list of molecules with specified substructure in SMILES notation
    """
    molecule = Chem.MolFromSmiles(mol)
    match = [smiles_str for smiles_str in mols
             if Chem.MolFromSmiles(smiles_str).HasSubstructMatch(molecule)]

    if not match:
        print(f"The substructure '{mol}' is not found in any of the proposed molecules.")
    return match


substructure = 'c1ccccc1'
molecules = ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]
result = substructure_search(molecules, substructure)
print(f"The substructure '{substructure}' is found in the following molecules: {result}")
