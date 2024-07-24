from fastapi import FastAPI
from rdkit import Chem


app = FastAPI()

mol_db = [{"molecule_id": 1, "molecule_structure": "CCO"},
          {"molecule_id": 2, "molecule_structure": "c1ccccc1"},
          {"molecule_id": 3, "molecule_structure": "CC(=O)O"},
          {"molecule_id": 4, "molecule_structure": "CC(=O)Oc1ccccc1C(=O)O"}
          ]


@app.post("/molecules/add")
def add_molecule():
    pass


# 5. List all molecules
@app.get("/molecules/all")
def get_all_molecules():
    return mol_db


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
    else:
        return match


substructure = 'c1ccccc1'
molecules = ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]
result = substructure_search(molecules, substructure)

if result is not None:
    print(f"The substructure '{substructure}' is found in the following molecules: {result}")
