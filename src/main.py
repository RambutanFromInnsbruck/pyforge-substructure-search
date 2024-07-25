from fastapi import FastAPI, HTTPException, status
from rdkit import Chem
from models import Molecule


app = FastAPI()

mols_db = [{"molecule_id": 1, "molecule_structure": "CCO"},
          {"molecule_id": 2, "molecule_structure": "c1ccccc1"},
          {"molecule_id": 3, "molecule_structure": "CC(=O)O"},
          {"molecule_id": 4, "molecule_structure": "CC(=O)Oc1ccccc1C(=O)O"}
          ]


# 0. Home page
@app.get("/", tags=["Start"], status_code=status.HTTP_200_OK)
def read_root():
    return {"message": "Welcome to Substructure Search Application"}

# 1. Add molecule (smiles) and its identifier
@app.post("/molecules/add", tags=["Molecules"], status_code=status.HTTP_201_CREATED)
def add_molecule(mol: Molecule):
    mols_db.append(mol)
    return mol

# 2. Get molecule by identifier
@app.get("/molecules/{molecule_id}", tags=["Molecules"], status_code=status.HTTP_200_OK)
def get_mol_by_id(molecule_id: int):
    for mol in mols_db:
        if mol["molecule_id"] == molecule_id:
            return mol
    raise HTTPException(status_code=404, detail="Molecule is not found")

# 3. Updating a molecule by identifier
@app.put("/molecules/{molecule_id}", tags=["Molecules"], status_code=status.HTTP_200_OK)
def update_molecule(molecule_id: int, updated_molecule: Molecule):
    for index, mol in enumerate(mols_db):
        if mol["molecule_id"] == molecule_id:
            mols_db[index] = updated_molecule.dict()
            return updated_molecule.dict()
    raise HTTPException(status_code=404, detail="Molecule is not found")

# 4. Delete a molecule by identifier
@app.delete("/molecules/{molecule_id}", tags=["Molecules"], status_code=status.HTTP_200_OK)
def delete_molecule(molecule_id: int):
    for index, mol in enumerate(mols_db):
        if mol["molecule_id"] == molecule_id:
            deleted_molecule = mols_db.pop()
            return deleted_molecule
    raise HTTPException(status_code=404, detail="Molecule is not found")

# 5. List all molecules
@app.get("/molecules", tags=["Molecules"], status_code=status.HTTP_200_OK)
def get_all_molecules():
    """
    Return list of all molecules
    """
    return mols_db

# 6. Substructure search for all added molecules
@app.get("/molecules/search/{search_molecule}", tags=["Molecules"], status_code=status.HTTP_200_OK)
def substructure_search(search_molecule: str):
    """
    Returns list of molecules with specified substructure or messages that there are no such molecules.
    ARGUMENTS:
        - **search_molecule**: sought substructure in SMILES notation
    RETURNS:
        list of molecules with specified substructure in SMILES notation
    """
    molecule = Chem.MolFromSmiles(search_molecule)
    match = [mol['molecule_structure'] for mol in mols_db
             if Chem.MolFromSmiles(mol['molecule_structure']).HasSubstructMatch(molecule)]
    return match
